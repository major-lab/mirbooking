#include "mirbooking-precomputed-score-table.h"

typedef struct
{
    GBytes *precomputed_table_bytes;
    gsize   seed_offset;
    gsize   seed_length;
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

static gfloat
compute_score (MirbookingScoreTable *score_table,
               MirbookingMirna      *mirna,
               MirbookingTarget     *target,
               gsize                 position,
               GError              **error)
{
    union
    {
        gint32 i;
        gfloat f;
    } ret;

    MirbookingPrecomputedScoreTable *self = MIRBOOKING_PRECOMPUTED_SCORE_TABLE (score_table);

    gsize seed_offset = self->priv->seed_offset;
    gsize seed_len    = self->priv->seed_length;

    // the seed does not fit in the target
    if (position + seed_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)))
    {
        return 0.0f;
    }

    // the seed does not fit in the mirna
    if (seed_offset + seed_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (mirna)))
    {
        return 0.0f;
    }

    gsize  data_len;
    const gfloat *data = g_bytes_get_data (self->priv->precomputed_table_bytes, &data_len);

    gssize i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), seed_offset, seed_len);
    gssize j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position, seed_len);

    // either sequence index is undefined
    if (i == -1 || j == -1)
    {
        return 0.0f;
    }

    gsize k = i * (1l << 2 * seed_len) + j;

    g_return_val_if_fail (k * sizeof (gfloat) < data_len, 0.0f);

    ret.f = data[k];
    ret.i = GINT32_FROM_BE (ret.i);

    return ret.f;
}

enum
{
    PROP_SCORE_TABLE = 1,
    PROP_SEED_OFFSET,
    PROP_SEED_LENGTH
};

static void
mirbooking_precomputed_score_table_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    switch (property_id)
    {
        case PROP_SCORE_TABLE:
            g_value_set_boxed (value, MIRBOOKING_PRECOMPUTED_SCORE_TABLE (object)->priv->precomputed_table_bytes);
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
            self->priv->precomputed_table_bytes = g_bytes_ref (precomputed_table);
            break;
        case PROP_SEED_OFFSET:
            self->priv->seed_offset = g_value_get_uint (value);
            break;
        case PROP_SEED_LENGTH:
            self->priv->seed_length = g_value_get_uint (value);
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
    g_object_class_install_property (object_class,
                                     PROP_SEED_OFFSET,
                                     g_param_spec_uint ("seed-offset", "", "", 0, G_MAXUINT, MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_OFFSET, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SEED_LENGTH,
                                     g_param_spec_uint ("seed-length", "", "", 1, G_MAXUINT, MIRBOOKING_PRECOMPUTED_SCORE_TABLE_DEFAULT_SEED_LENGTH, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

/**
 * mirbooking_precomputed_score_table_new:
 * Returns: (transfer full)
 */
MirbookingPrecomputedScoreTable *
mirbooking_precomputed_score_table_new (gfloat *data, gsize seed_offset, gsize seed_len)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_TABLE,
                         "score-table", g_bytes_new_take (data, (1l << 2 * seed_len) * (1l << 2 * seed_len) * sizeof (gfloat)),
                         "seed-offset", seed_offset,
                         "seed-length", seed_len,
                         NULL);
}

/**
 * mirbooking_precomputed_score_table_new_from_bytes:
 *
 * Returns: (transfer full)
 */
MirbookingPrecomputedScoreTable *
mirbooking_precomputed_score_table_new_from_bytes (GBytes *precomputed_table, gsize seed_offset, gsize seed_len)
{
    return g_object_new (MIRBOOKING_TYPE_PRECOMPUTED_SCORE_TABLE,
                         "score-table", precomputed_table,
                         "seed-offset", seed_offset,
                         "seed-length", seed_len,
                         NULL);
}

/**
 * mirbooking_precomputed_score_table_new_from_stream:
 *
 * Returns: (transfer full): A #MirbookingPrecomputedScoreTable or %NULL if any
 * stream-related operation failed and @error will be set.
 */
MirbookingPrecomputedScoreTable *
mirbooking_precomputed_score_table_new_from_stream (GInputStream *stream, gsize seed_offset, gsize seed_len, GError **error)
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
                             "seed-offset", seed_offset,
                             "seed-length", seed_len,
                             NULL);
    }
    else
    {
        return NULL;
    }
}
