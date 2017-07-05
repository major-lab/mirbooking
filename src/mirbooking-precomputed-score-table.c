#include "mirbooking-precomputed-score-table.h"

/*
 * The score table act as a pre-computed set of tables that hold the
 * hybridation probabilities for pair of sequences of various lengths.
 */
typedef gfloat MirbookingPrecomputedTable[16384][16384];

typedef struct
{
    GBytes                           *precomputed_table_bytes;
    const MirbookingPrecomputedTable *precomputed_table;
} MirbookingPrecomputedScoreTablePrivate;

struct _MirbookingPrecomputedScoreTable
{
    MirbookingScoreTable parent_instance;
    MirbookingPrecomputedScoreTablePrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingPrecomputedScoreTable, mirbooking_precomputed_score_table, MIRBOOKING_TYPE_SCORE_TABLE)

static void
mirbooking_precomputed_score_table_init (MirbookingPrecomputedScoreTable *self)
{
    self->priv = g_new0 (MirbookingPrecomputedScoreTablePrivate, 1);
}

static void
mirbooking_precomputed_score_table_finalize (GObject *object)
{
    MirbookingPrecomputedScoreTable *self = MIRBOOKING_PRECOMPUTED_SCORE_TABLE (object);

    if (self->priv->precomputed_table_bytes != NULL)
    {
        g_bytes_unref (self->priv->precomputed_table_bytes);
    }

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_precomputed_score_table_parent_class)->finalize (object);
}

static gsize
sequence_index (const gchar *seq, gsize seq_len)
{
    gint i;
    gsize index = 0;

    // note: 4**i = 2**2i = 2 << 2i - 1
    // for i = 0, we do it manually as it would result in a negative shift

    for (i = 0; i < seq_len; i++)
    {
        gint base = i == 0 ? 1 : (2 << (2 * i - 1));
        switch (seq[seq_len - i - 1])
        {
            case 'A':
                index += 0 * base;
                break;
            case 'C':
                index += 1 * base;
                break;
            case 'G':
                index += 2 * base;
                break;
            case 'T':
            case 'U':
                index += 3 * base;
                break;
            default:
                g_return_val_if_reached (0);
        }
    }

    return index;
}

static gfloat
compute_score (MirbookingScoreTable *self,
               MirbookingSequence   *a,
               gsize                 a_offset,
               MirbookingSequence   *b,
               gsize                 b_offset,
               gsize                 len,
               GError              **error)
{
    union
    {
        gint32 i;
        gfloat f;
    } ret;

    gsize i = sequence_index (mirbooking_sequence_get_subsequence (a, a_offset, len), len);
    gsize j = sequence_index (mirbooking_sequence_get_subsequence (b, b_offset, len), len);

    switch (len)
    {
        case 7:
            ret.f = (*MIRBOOKING_PRECOMPUTED_SCORE_TABLE (self)->priv->precomputed_table)[i][j];
            break;
        default:
            g_return_val_if_reached (0.0);
    }

    ret.i = GINT32_FROM_BE (ret.i);

    return ret.f;
}

enum
{
    PROP_SCORE_TABLE = 1
};

static void
mirbooking_precomputed_score_table_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    switch (property_id)
    {
        case PROP_SCORE_TABLE:
            g_value_set_boxed (value, MIRBOOKING_PRECOMPUTED_SCORE_TABLE (object)->priv->precomputed_table);
            break;
        default:
                g_assert_not_reached ();
    }
}

static void
mirbooking_precomputed_score_table_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingPrecomputedScoreTable *self = MIRBOOKING_PRECOMPUTED_SCORE_TABLE (object);

    GBytes *precomputed_table;

    switch (property_id)
    {
        case PROP_SCORE_TABLE:
            precomputed_table = g_value_get_boxed (value);
            g_return_if_fail (g_bytes_get_size (precomputed_table) == 16384 * 16384 * sizeof (gfloat));
            self->priv->precomputed_table_bytes = g_bytes_ref (precomputed_table);
            self->priv->precomputed_table = g_bytes_get_data (precomputed_table, NULL);
            break;
        default:
                g_assert_not_reached ();
    }
}

static void
mirbooking_precomputed_score_table_class_init (MirbookingPrecomputedScoreTableClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->get_property = mirbooking_precomputed_score_table_get_property;
    object_class->set_property = mirbooking_precomputed_score_table_set_property;
    object_class->finalize     = mirbooking_precomputed_score_table_finalize;

    klass->parent_class.compute_score = compute_score;

    g_object_class_install_property (object_class,
                                     PROP_SCORE_TABLE,
                                     g_param_spec_boxed ("score-table", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

/**
 * mirbooking_precomputed_score_table_new:
 * Returns: (transfer full)
 */
MirbookingPrecomputedScoreTable *
mirbooking_precomputed_score_table_new (gfloat *data)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_TABLE,
                         "score-table", g_bytes_new_take (data, sizeof (gfloat)),
                         NULL);
}

/**
 * mirbooking_precomputed_score_table_new_from_bytes:
 *
 * Returns: (transfer full)
 */
MirbookingPrecomputedScoreTable *
mirbooking_precomputed_score_table_new_from_bytes (GBytes *precomputed_table)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_TABLE,
                         "score-table", precomputed_table,
                         NULL);
}

/**
 * mirbooking_precomputed_score_table_new_from_stream:
 *
 * Returns: (transfer full): A #MirbookingPrecomputedScoreTable or %NULL if any
 * stream-related operation failed and @error will be set.
 */
MirbookingPrecomputedScoreTable *
mirbooking_precomputed_score_table_new_from_stream (GInputStream *stream, GError **error)
{
    g_autoptr (GMemoryOutputStream) out = G_MEMORY_OUTPUT_STREAM (g_memory_output_stream_new_resizable ());

    if (g_output_stream_splice (G_OUTPUT_STREAM (out),
                                stream,
                                G_OUTPUT_STREAM_SPLICE_CLOSE_SOURCE | G_OUTPUT_STREAM_SPLICE_CLOSE_TARGET,
                                NULL,
                                error))
    {
        return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_TABLE,
                             "score-table", g_memory_output_stream_steal_as_bytes (out),
                             NULL);
    }
    else
    {
        return NULL;
    }
}
