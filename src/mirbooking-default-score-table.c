#include "mirbooking-default-score-table.h"

#include <math.h>
#include <sparse.h>

#define R 1.987203611e-3
#define T 310.15

/*
 * For the duplex: CUCCAUC&GAUGGAG, ViennaRNA reports a free energy of -9.30
 * kcal/mol.
 *
 * Wee et al. measured a dissociation constant for a mouse Ago2 protein
 * carrying a guide miRNA with only the seed pairing of 26±2 pM, which
 * correspond to a free energy of -15.02 kcal/mol.
 *
 * On the other hand, Salomon et al. instead measured 15±2 pm, which correspond
 * to -15.36 kcal/mol.
 *
 * We thus impute the -6.06 kcal/mol to AGO2 entropic contribution.
 *
 * Reference: Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into
 * Domains with Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5
 * (November 21, 2012): 1055–67, https://doi.org/10.1016/j.cell.2012.10.036.
 * */
#define AGO2_SCORE (-6.06f)

typedef struct
{
    GBytes       *seed_scores_bytes;
    SparseMatrix  seed_scores; /* view over @seed_scores_bytes */
    gsize         seed_offset;
    gsize         seed_length;
    GBytes       *supplementary_scores_bytes;
    SparseMatrix  supplementary_scores;
} MirbookingDefaultScoreTablePrivate;

struct _MirbookingDefaultScoreTable
{
    MirbookingScoreTable parent_instance;
    MirbookingDefaultScoreTablePrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingDefaultScoreTable, mirbooking_default_score_table, MIRBOOKING_TYPE_SCORE_TABLE)

static void
mirbooking_default_score_table_init (MirbookingDefaultScoreTable *self)
{
    self->priv = g_new0 (MirbookingDefaultScoreTablePrivate, 1);
}

static void
_init_sparse_matrix_from_bytes (SparseMatrix *sm, GBytes *bytes)
{
    gsize *d = (gsize*) g_bytes_get_data (bytes, NULL);

    g_return_if_fail (d != NULL);

    gsize n       = *d;
    gsize nnz     = *(d + 1);
    gsize* rowptr = d + 2;
    gsize* colind = rowptr + n + 1;
    gfloat* data  = (gfloat*) (colind + nnz);

    /* initialize the sparse matrix */
    sm->storage        = SPARSE_MATRIX_STORAGE_CSR;
    sm->type           = SPARSE_MATRIX_TYPE_FLOAT;
    sm->shape[0]       = n;
    sm->shape[1]       = n;
    sm->s.csr.nnz      = nnz;
    sm->s.csr.colind   = colind;
    sm->s.csr.rowptr   = rowptr;
    sm->default_data.f = INFINITY;
    sm->data           = data;

}

static void
mirbooking_default_score_table_finalize (GObject *object)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (object);

    if (self->priv->seed_scores_bytes != NULL)
    {
        g_bytes_unref (self->priv->seed_scores_bytes);
    }

    if (self->priv->supplementary_scores_bytes != NULL)
    {
        g_bytes_unref (self->priv->supplementary_scores_bytes);
    }

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_default_score_table_parent_class)->finalize (object);
}

static gdouble
_compute_Kd (gdouble deltaG)
{
    return 1e12 * exp (deltaG / (R * T));
}

static gdouble
compute_score (MirbookingScoreTable *score_table,
               MirbookingMirna      *mirna,
               MirbookingTarget     *target,
               gsize                 position,
               GError              **error)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (score_table);

    gsize seed_offset = self->priv->seed_offset;
    gsize seed_len    = self->priv->seed_length;

    // the seed does not fit in the target
    if (position + seed_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)))
    {
        return INFINITY;
    }

    // the seed does not fit in the mirna
    if (seed_offset + seed_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (mirna)))
    {
        return INFINITY;
    }

    gssize seed_i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), seed_offset, seed_len);
    gssize seed_j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position, seed_len);

    // either sequence index is undefined
    if (seed_i == -1 || seed_j == -1)
    {
        return INFINITY;
    }

    gfloat seed_score = sparse_matrix_get_float (&self->priv->seed_scores, seed_i, seed_j);

    gfloat supplementary_score = 0.0f;
    gsize supplementary_offset = 12;
    gsize supplementary_len = 4;
    if (self->priv->supplementary_scores_bytes != NULL &&
        position >= (supplementary_offset - supplementary_len) &&
        mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (mirna)) >= (supplementary_offset + supplementary_len))
    {
        gsize supplementary_i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), supplementary_offset, supplementary_len);
        gsize supplementary_j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position - supplementary_offset + supplementary_len, supplementary_len);
        if (supplementary_i != -1 && supplementary_j != -1)
        {
            supplementary_score = sparse_matrix_get_float (&self->priv->supplementary_scores,
                                                           supplementary_i,
                                                           supplementary_j);
        }
    }

    return _compute_Kd (seed_score + supplementary_score + mirbooking_target_get_accessibility_score (target, position) + AGO2_SCORE);
}

static gboolean
compute_positions (MirbookingScoreTable  *score_table,
                   MirbookingMirna       *mirna,
                   MirbookingTarget      *target,
                   gsize                **positions,
                   gsize                 *positions_len,
                   GError               **error)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (score_table);

    gsize k = 0;

    gsize seed_offset = self->priv->seed_offset;
    gsize seed_len    = self->priv->seed_length;

    gssize i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), seed_offset, seed_len);

    // miRNA seed is undefined (thus no suitable targets)
    g_return_val_if_fail (i != -1, FALSE);

    gsize j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), 0, seed_len);

    gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target));

    gsize total_positions_len = seq_len - seed_len + 1;

    gsize  *_positions = NULL;

    gsize p;
    for (p = 0; p < total_positions_len; p++)
    {
        if (sparse_matrix_get_float (&self->priv->seed_scores, i, j) < INFINITY)
        {
            _positions = g_realloc (_positions, (k + 1) * sizeof (gsize));
            _positions[k++] = p;
        }

        if (p < seq_len - seed_len)
        {
            gsize out = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), p, 1);
            gsize in  = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), p+seed_len, 1);

            if (in == -1)
            {
                break; // FIXME
            }

            g_assert_cmpint (in, !=, -1);
            g_assert_cmpint (out, !=, -1);

            j -= out * (2l << (2 * (seed_len - 1) - 1));
            j *= 4;
            j += in;

            g_assert_cmpint (j, <, 16384);
            g_assert_cmpint (j, >=, 0);
        }
    }

    *positions     = _positions;
    *positions_len = k;

    return TRUE;
}

enum
{
    PROP_SEED_SCORES = 1,
    PROP_SEED_OFFSET,
    PROP_SEED_LENGTH,
    PROP_SUPPLEMENTARY_SCORES
};

static void
mirbooking_default_score_table_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    switch (property_id)
    {
        case PROP_SEED_SCORES:
            g_value_set_boxed (value, MIRBOOKING_DEFAULT_SCORE_TABLE (object)->priv->seed_scores_bytes);
            break;
        case PROP_SUPPLEMENTARY_SCORES:
            g_value_set_boxed (value, MIRBOOKING_DEFAULT_SCORE_TABLE (object)->priv->supplementary_scores_bytes);
        default:
                g_assert_not_reached ();
    }
}

static void
mirbooking_default_score_table_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (object);

    GBytes *score_table;

    switch (property_id)
    {
        case PROP_SEED_SCORES:
            score_table = g_value_get_boxed (value);
            self->priv->seed_scores_bytes = g_bytes_ref (score_table);
            _init_sparse_matrix_from_bytes (&self->priv->seed_scores,
                                            self->priv->seed_scores_bytes);
            break;
        case PROP_SEED_OFFSET:
            self->priv->seed_offset = g_value_get_uint (value);
            break;
        case PROP_SEED_LENGTH:
            self->priv->seed_length = g_value_get_uint (value);
            break;
        case PROP_SUPPLEMENTARY_SCORES:
            score_table = g_value_get_boxed (value);
            if (score_table != NULL)
            {
                self->priv->supplementary_scores_bytes = g_bytes_ref (score_table);
                _init_sparse_matrix_from_bytes (&self->priv->supplementary_scores,
                                                self->priv->supplementary_scores_bytes);
            }
            break;
        default:
            g_assert_not_reached ();
    }
}

static void
mirbooking_default_score_table_class_init (MirbookingDefaultScoreTableClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->get_property = mirbooking_default_score_table_get_property;
    object_class->set_property = mirbooking_default_score_table_set_property;
    object_class->finalize     = mirbooking_default_score_table_finalize;

    klass->parent_class.compute_score     = compute_score;
    klass->parent_class.compute_positions = compute_positions;

    g_object_class_install_property (object_class,
                                     PROP_SEED_SCORES,
                                     g_param_spec_boxed ("seed-scores", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SEED_OFFSET,
                                     g_param_spec_uint ("seed-offset", "", "", 0, G_MAXUINT, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SEED_OFFSET, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SEED_LENGTH,
                                     g_param_spec_uint ("seed-length", "", "", 1, G_MAXUINT, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SEED_LENGTH, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SUPPLEMENTARY_SCORES,
                                     g_param_spec_boxed ("supplementary-scores", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

/**
 * mirbooking_default_score_table_new:
 * Returns: (transfer full)
 */
MirbookingDefaultScoreTable *
mirbooking_default_score_table_new (gfloat *data, gsize seed_offset, gsize seed_len)
{
    g_return_val_if_fail (data != NULL, NULL);

    return g_object_new (MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE,
                         "seed-scores", g_bytes_new_take (data, 0),
                         "seed-offset", seed_offset,
                         "seed-length", seed_len,
                         NULL);
}

/**
 * mirbooking_default_score_table_new_from_bytes:
 *
 * Returns: (transfer full)
 */
MirbookingDefaultScoreTable *
mirbooking_default_score_table_new_from_bytes (GBytes *seed_scores, gsize seed_offset, gsize seed_len, GBytes *supp_scores)
{
    g_return_val_if_fail (g_bytes_get_data (seed_scores, NULL) != NULL, NULL);

    return g_object_new (MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE,
                         "seed-scores", seed_scores,
                         "seed-offset", seed_offset,
                         "seed-length", seed_len,
                         "supplementary-scores", supp_scores,
                         NULL);
}

/**
 * mirbooking_default_score_table_new_from_stream:
 *
 * Returns: (transfer full): A #MirbookingDefaultScoreTable or %NULL if any
 * stream-related operation failed and @error will be set.
 */
MirbookingDefaultScoreTable *
mirbooking_default_score_table_new_from_stream (GInputStream *stream, gsize seed_offset, gsize seed_len, GError **error)
{
    g_autoptr (GMemoryOutputStream) out = G_MEMORY_OUTPUT_STREAM (g_memory_output_stream_new_resizable ());

    if (g_output_stream_splice (G_OUTPUT_STREAM (out),
                                stream,
                                G_OUTPUT_STREAM_SPLICE_CLOSE_SOURCE | G_OUTPUT_STREAM_SPLICE_CLOSE_TARGET,
                                NULL,
                                error))
    {
        return g_object_new (MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE,
                             "seed-scores", g_memory_output_stream_steal_as_bytes (out),
                             "seed-offset", seed_offset,
                             "seed-length", seed_len,
                             NULL);
    }
    else
    {
        return NULL;
    }
}
