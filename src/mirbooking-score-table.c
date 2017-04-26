#include "mirbooking-score-table.h"

/*
 * The score table act as a pre-computed set of tables that hold the
 * hybridation probabilities for pair of sequences of various lengths.
 */

typedef gfloat MirbookingScoreTable3x3[64][64];
typedef gfloat MirbookingScoreTable4x4[256][256];
typedef gfloat MirbookingScoreTable7x7[16384][16384];

typedef struct
{
    const MirbookingScoreTable3x3 *table3;
    const MirbookingScoreTable4x4 *table4;
    GBytes                        *table7_bytes;
    const MirbookingScoreTable7x7 *table7;
} MirbookingScoreTablePrivate;

struct _MirbookingScoreTable
{
    GObject parent_instance;
    MirbookingScoreTablePrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingScoreTable, mirbooking_score_table, G_TYPE_OBJECT)

static void
mirbooking_score_table_init (MirbookingScoreTable *self)
{
    self->priv = g_new (MirbookingScoreTablePrivate, 1);
}

static void
mirbooking_score_table_finalize (GObject *object)
{
    MirbookingScoreTable *self = (MIRBOOKING_SCORE_TABLE (object));

    g_bytes_unref (self->priv->table7_bytes);
    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_score_table_parent_class)->finalize (object);
}

static void
mirbooking_score_table_class_init (MirbookingScoreTableClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->finalize = mirbooking_score_table_finalize;
}

MirbookingScoreTable *
mirbooking_score_table_new ()
{
    return g_object_new (MIRBOOKING_TYPE_SCORE_TABLE, NULL);
}

MirbookingScoreTable *
mirbooking_score_table_new_precomputed (GBytes *precomputed_table)
{
    MirbookingScoreTable *ret;

    ret = g_object_new (MIRBOOKING_TYPE_SCORE_TABLE, NULL);

    switch (g_bytes_get_size (precomputed_table))
    {
        case 16384 * 16384 * sizeof (gfloat):
            ret->priv->table7_bytes = g_bytes_ref (precomputed_table);
            ret->priv->table7       = g_bytes_get_data (precomputed_table, NULL);
            break;
        default:
            g_object_unref (ret);
            g_return_val_if_reached (NULL);
    }

    return ret;
}

static guint
sequence_index (const gchar *seq, gsize seq_len)
{
    gint i;
    guint index = 0;

    for (i = 0; i < seq_len; i++)
    {
        switch (seq[i])
        {
            case 'A':
                index += 4 << i;
                break;
            case 'C':
                index += 4 << i;
                break;
            case 'G':
                index += 4 << i;
                break;
            case 'T':
                index += 4 << i;
                break;
            default:
                g_return_val_if_reached (0);
        }
    }

    return index;
}

/**
 * mirbooking_score_table_compute_score:
 *
 * Compute the score for a pair of #MirbookingSequence objects, extracting a
 * substring of a fixed length and evaluating the score table at its
 * corresponding index.
 */
gfloat
mirbooking_score_table_compute_score (MirbookingScoreTable *self,
                                      MirbookingSequence   *a,
                                      gsize                 a_offset,
                                      MirbookingSequence   *b,
                                      gsize                 b_offset,
                                      gsize                 len)
{
    gfloat ret;

    gint i = sequence_index (mirbooking_sequence_get_subsequence (a, a_offset, len), len);
    gint j = sequence_index (mirbooking_sequence_get_subsequence (b, b_offset, len), len);

    switch (len)
    {
        case 7:
            ret = (gfloat) GUINT32_FROM_BE ((guint32) *self->priv->table7[i][j]);
            break;
        default:
            g_return_val_if_reached (0.0);
    }

    return ret;
}
