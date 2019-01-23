#include "mirbooking-default-score-table.h"

#include <math.h>
#include <sparse.h>

#define R 1.987203611e-3
#define T 310.15

#define SEED_OFFSET 1
#define SEED_LENGTH 7

/*
 * For the duplex: 'CUACCUC&GAGGUAG', ViennaRNA reports a binding energy of
 * -9.37 kcal/mol.
 *
 * Wee et al. measured a dissociation constant for a mouse Ago2 protein
 * carrying a guide miRNA with only the seed pairing of 26±2 pM, which
 * correspond to a free energy of -15.02 kcal/mol.
 *
 * On the other hand, Salomon et al. instead measured 15±2 pm, which correspond
 * to -15.36 kcal/mol.
 *
 * In addition, the seed setup experiment in both experiments had 'A', which
 * shoudld account for a -0.56 kcal/mol additional contribution (Schirle et al. 2015).
 *
 * We thus the latest value and impute the -5.43 kcal/mol gap to AGO2 entropic
 * contribution.
 *
 * Reference:
 * Liang Meng Wee et al., “Argonaute Divides Its RNA Guide into
 * Domains with Distinct Functions and RNA-Binding Properties,” Cell 151, no. 5
 * (November 21, 2012): 1055–67, https://doi.org/10.1016/j.cell.2012.10.036.
 *
 * Nicole T Schirle et al., “Water-Mediated Recognition of T1-Adenosine Anchors
 * Argonaute2 to MicroRNA Targets,” ed. Phillip D Zamore, ELife 4 (September
 * 11, 2015): e07646, https://doi.org/10.7554/eLife.07646.
 */

#define AGO2_SCORE (-5.43f)

typedef struct
{
    GBytes       *seed_scores_bytes;
    SparseMatrix  seed_scores; /* view over @seed_scores_bytes */
    MirbookingDefaultScoreTableSupplementaryModel supplementary_model;
    GBytes       *supplementary_scores_bytes;
    SparseMatrix  supplementary_scores;
    /* hints for filtering interactions */
    MirbookingDefaultScoreTableFilter filter;
    gpointer                          filter_user_data;
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

    gsize n             = *d;
    gsize nnz           = *(d + 1);
    const gsize* rowptr = d + 2;
    const gsize* colind = rowptr + n + 1;
    const gfloat* data  = (gfloat*) (colind + nnz);

    /* initialize the sparse matrix */
    sm->storage        = SPARSE_MATRIX_STORAGE_CSR;
    sm->type           = SPARSE_MATRIX_TYPE_FLOAT;
    sm->shape[0]       = n;
    sm->shape[1]       = n;
    sm->s.csr.nnz      = nnz;
    sm->s.csr.colind   = (gsize*) colind;
    sm->s.csr.rowptr   = (gsize*) rowptr;
    sm->default_data.f = INFINITY;
    sm->data           = (void*) data;
}

static void
mirbooking_default_score_table_constructed (GObject *object)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (object);

    _init_sparse_matrix_from_bytes (&self->priv->seed_scores,
                                    self->priv->seed_scores_bytes);

    if (self->priv->supplementary_scores_bytes != NULL)
    {
        _init_sparse_matrix_from_bytes (&self->priv->supplementary_scores,
                                        self->priv->supplementary_scores_bytes);
    }

    G_OBJECT_CLASS (mirbooking_default_score_table_parent_class)->constructed (object);
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

static gfloat
_get_subsequence_score (SparseMatrix     *scores,
                        MirbookingMirna  *mirna,
                        MirbookingTarget *target,
                        gsize             position,
                        gsize             mirna_position_offset,
                        gsize             subsequence_offset,
                        gsize             subsequence_len)
{
    g_return_val_if_fail ((1 << (2 * subsequence_len)) == scores->shape[0], INFINITY);

    // the subsequence does not fit in the target
    if (position + mirna_position_offset - subsequence_offset + subsequence_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)))
    {
        return INFINITY;
    }

    // the subsequence does not fit in the mirna
    if (subsequence_offset + subsequence_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (mirna)))
    {
        return INFINITY;
    }

    gsize subsequence_i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), subsequence_offset, subsequence_len);
    gsize subsequence_j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + mirna_position_offset - subsequence_offset, subsequence_len);

    if (subsequence_i == -1 || subsequence_j == -1)
    {
        return INFINITY;
    }

    gfloat score = sparse_matrix_get_float (scores,
                                            subsequence_i,
                                            subsequence_j);

    return score;
}

static gdouble
_compute_Kd (gdouble deltaG)
{
    return 1e12 * exp (deltaG / (R * T));
}

static gboolean
is_g_bulge (MirbookingTarget *target, gsize position)
{
    return position + 4 + 4 <= mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)) &&
        *mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), position + 3, 1) == 'G';
}

static gdouble
compute_score (MirbookingScoreTable *score_table,
               MirbookingMirna      *mirna,
               MirbookingTarget     *target,
               gsize                 position,
               GError              **error)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (score_table);

    /*
     * AGO2 has a slight preference for sites starting with 'A' at position t1.
     *
     * Reference:
     * Nicole T Schirle et al., “Water-Mediated Recognition of T1-Adenosine
     * Anchors Argonaute2 to MicroRNA Targets,” ed. Phillip D Zamore, ELife 4
     * (September 11, 2015): e07646, https://doi.org/10.7554/eLife.07646.
     */
    gfloat A_score = 0.0f;
    if (position + SEED_LENGTH + 1 <= mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)) &&
        *mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), position + SEED_LENGTH, 1) == 'A')
    {
        A_score = -0.56f;
    }

    gfloat seed_score = _get_subsequence_score (&self->priv->seed_scores,
                                                mirna,
                                                target,
                                                position,
                                                SEED_OFFSET,
                                                SEED_OFFSET,
                                                SEED_LENGTH);

    /*
     * We allow a 'G' nucleation bulge at position t5.
     *
     * There's a penalty of 1.2 kcal/mol for this non-canonical motif that was
     * taken from the authors free energy estimates.
     *
     * Reference:
     * Sung Wook Chi, Gregory J. Hannon, and Robert B. Darnell, “An Alternative
     * Mode of MicroRNA Target Recognition,” Nature Structural & Molecular
     * Biology 19, no. 3 (March 2012): 321–27, https://doi.org/10.1038/nsmb.2230.
     */
    if (is_g_bulge (target, position))
    {
        gsize i, j;
        i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), SEED_OFFSET, SEED_LENGTH);
        j = (1l << (2 * 4)) * mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position, 3) +
            mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + 4, 4);
        seed_score = MIN (seed_score, sparse_matrix_get_float (&self->priv->seed_scores, i, j) + 1.2f);
    }

    gfloat supplementary_score = 0;
    if (self->priv->supplementary_model == MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_3PRIME)
    {
        if (self->priv->supplementary_scores_bytes != NULL)
        {
            gint bulge;
            for (bulge = 0; bulge < 5; bulge++)
            {
                gfloat supplementary_score_ = _get_subsequence_score (&self->priv->supplementary_scores,
                                                                      mirna,
                                                                      target,
                                                                      position,
                                                                      4 - bulge,
                                                                      SEED_OFFSET + SEED_LENGTH + 4,
                                                                      4);

                supplementary_score = MIN (supplementary_score, supplementary_score_);
            }
        }
    }
    else if (self->priv->supplementary_model == MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_YAN_ET_AL_2018)
    {
        /*
         * Reference:
         * Yifei Yan et al., “The Sequence Features That Define Efficient and
         * Specific HAGO2-Dependent MiRNA Silencing Guides,” Nucleic Acids
         * Research, June 22, 2018, https://doi.org/10.1093/nar/gky546.
         */
        gfloat A_box = 0, B_box = 0, C_box = 0, D_box = 0;
        if (self->priv->supplementary_scores_bytes != NULL)
        {
            // A box
            A_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                            mirna,
                                            target,
                                            position,
                                            5,
                                            SEED_OFFSET + SEED_LENGTH,
                                            3);

            // B box
            B_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                            mirna,
                                            target,
                                            position,
                                            5,
                                            SEED_OFFSET + SEED_LENGTH + 3,
                                            3);

            // C box
            C_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                            mirna,
                                            target,
                                            position,
                                            5,
                                            SEED_OFFSET + SEED_LENGTH + 6,
                                            3);

            // TODO: D box (remaining nucleotides)

            // AAA::UUU
            gfloat mismatch_threshold = -1e-8; // sparse_matrix_get_float (&self->priv->supplementary_scores, 0, 63);

            if (B_box <= mismatch_threshold)
            {
                supplementary_score += B_box;
            }

            if (B_box <= mismatch_threshold && C_box <= mismatch_threshold)
            {
                supplementary_score += C_box;
            }

            if (A_box <= mismatch_threshold && B_box <= mismatch_threshold && C_box <= mismatch_threshold)
            {
                supplementary_score += A_box;
            }

            if (A_box <= mismatch_threshold && B_box <= mismatch_threshold && C_box <= mismatch_threshold && D_box <= mismatch_threshold)
            {
                supplementary_score += D_box;
            }
        }
    }

    return _compute_Kd (A_score + seed_score + supplementary_score + mirbooking_target_get_accessibility_score (target, position) + AGO2_SCORE);
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

    gssize i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), SEED_OFFSET, SEED_LENGTH);

    // miRNA seed is undefined (thus no suitable targets)
    if (i == -1)
    {
        *positions = NULL;
        *positions_len = 0;
        return TRUE;
    }

    gsize j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), 0, SEED_LENGTH);

    gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target));

    gsize total_positions_len = seq_len - SEED_LENGTH + 1;

    g_autofree gsize  *_positions = NULL;

    if (self->priv->filter != NULL && !self->priv->filter (self, mirna, target, -1, self->priv->filter_user_data))
    {
        *positions = NULL;
        *positions_len = 0;
        return TRUE;
    }

    gsize p;
    for (p = 0; p < total_positions_len; p++)
    {
        if ((sparse_matrix_get_float (&self->priv->seed_scores, i, j) < INFINITY || is_g_bulge (target, p)) &&
            mirbooking_target_get_accessibility_score (target, p) < INFINITY                                &&
            (self->priv->filter == NULL || self->priv->filter (self, mirna, target, p, self->priv->filter_user_data)))
        {
            _positions = g_realloc (_positions, (k + 1) * sizeof (gsize));
            _positions[k++] = p;
        }

        if (p < seq_len - SEED_LENGTH)
        {
            gsize out = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), p, 1);
            gsize in  = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), p+SEED_LENGTH, 1);

            if (in == -1)
            {
                break; // FIXME
            }

            g_assert_cmpint (in, !=, -1);
            g_assert_cmpint (out, !=, -1);

            j -= out * (2l << (2 * (SEED_LENGTH - 1) - 1));
            j *= 4;
            j += in;

            g_assert_cmpint (j, <, 16384);
            g_assert_cmpint (j, >=, 0);
        }
    }

    *positions     = g_steal_pointer (&_positions);
    *positions_len = k;

    return TRUE;
}

enum
{
    PROP_SEED_SCORES = 1,
    PROP_SUPPLEMENTARY_MODEL,
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
            break;
        case PROP_SUPPLEMENTARY_MODEL:
            self->priv->supplementary_model = g_value_get_uint (value);
            break;
        case PROP_SUPPLEMENTARY_SCORES:
            score_table = g_value_get_boxed (value);
            self->priv->supplementary_scores_bytes = score_table == NULL ? NULL : g_bytes_ref (score_table);
            break;
        default:
            g_assert_not_reached ();
    }
}

static void
mirbooking_default_score_table_class_init (MirbookingDefaultScoreTableClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->constructed  = mirbooking_default_score_table_constructed;
    object_class->get_property = mirbooking_default_score_table_get_property;
    object_class->set_property = mirbooking_default_score_table_set_property;
    object_class->finalize     = mirbooking_default_score_table_finalize;

    klass->parent_class.compute_score     = compute_score;
    klass->parent_class.compute_positions = compute_positions;

    g_object_class_install_property (object_class,
                                     PROP_SEED_SCORES,
                                     g_param_spec_boxed ("seed-scores", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SUPPLEMENTARY_MODEL,
                                     g_param_spec_uint ("supplementary-model", "", "", 0, 2, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SUPPLEMENTARY_SCORES,
                                     g_param_spec_boxed ("supplementary-scores", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

/**
 * mirbooking_default_score_table_new:
 *
 * Returns: (transfer full)
 */
MirbookingDefaultScoreTable *
mirbooking_default_score_table_new (GBytes *seed_scores, MirbookingDefaultScoreTableSupplementaryModel supp_model, GBytes *supp_scores)
{
    g_return_val_if_fail (g_bytes_get_data (seed_scores, NULL) != NULL, NULL);

    return g_object_new (MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE,
                         "seed-scores", seed_scores,
                         "supplementary-model", supp_model,
                         "supplementary-scores", supp_scores,
                         NULL);
}

void
mirbooking_default_score_table_set_filter (MirbookingDefaultScoreTable *self,
                                           MirbookingDefaultScoreTableFilter filter,
                                           gpointer user_data)
{
    self->priv->filter           = filter;
    self->priv->filter_user_data = user_data;
}
