#include "mirbooking-default-score-table.h"
#include "mirbooking-score-table-private.h"

#include <math.h>
#include <sparse.h>
#include <ctype.h>

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
    GDestroyNotify                    filter_user_data_destroy;
} MirbookingDefaultScoreTablePrivate;

struct _MirbookingDefaultScoreTable
{
    MirbookingScoreTable parent_instance;
    MirbookingDefaultScoreTablePrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingDefaultScoreTable, mirbooking_default_score_table, MIRBOOKING_TYPE_SCORE_TABLE)

/**
 * mirbooking_default_score_table_cutoff_filter:
 * @user_data: (type MirbookingDefaultScoreTableCutoffFilterUserData)
 */
gboolean
mirbooking_default_score_table_cutoff_filter  (MirbookingDefaultScoreTable *score_table,
                                               MirbookingMirna             *mirna,
                                               MirbookingTarget            *target,
                                               gssize                       position,
                                               gpointer                     user_data)
{
    MirbookingDefaultScoreTableCutoffFilterUserData *cutoff_filter_user_data = user_data;

    gdouble E0 = mirbooking_broker_get_sequence_quantity (cutoff_filter_user_data->broker, MIRBOOKING_SEQUENCE (mirna));
    gdouble S0 = mirbooking_broker_get_sequence_quantity (cutoff_filter_user_data->broker, MIRBOOKING_SEQUENCE (target));

    gdouble Km;
    if (position == -1)
    {
        Km = 0;
    }
    else
    {
        MirbookingScore score;
        mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (score_table),
                                              mirna,
                                              target,
                                              position,
                                              &score,
                                              NULL);
         Km = MIRBOOKING_SCORE_KM (score);
    }

    gdouble Z = E0 + S0 + Km;

    gdouble ES = ((Z - sqrt (pow (Z, 2) - 4 * E0 * S0)) / 2.0);

    return ES >= cutoff_filter_user_data->cutoff &&
        ((ES / S0) >= cutoff_filter_user_data->bound_fraction_cutoff) &&
        ((ES / E0) >= cutoff_filter_user_data->sponged_fraction_cutoff);
}

/**
 * mirbooking_default_score_table_blacklist_filter:
 * @user_data: (type MirbookingDefaultScoreTableBlacklistFilterUserdata)
 * The user_data pointer contains a GHashTable that acts as a hash set. It is
 * sufficient for an blacklisted interaction to appear as a key to this map
 * regardless of its associated value.
 */
gboolean
mirbooking_default_score_table_blacklist_filter (MirbookingDefaultScoreTable *score_table,
                                                 MirbookingMirna             *mirna,
                                                 MirbookingTarget            *target,
                                                 gssize                       position,
                                                 gpointer                     user_data)
{
    MirbookingDefaultScoreTableBlacklistFilterUserdata *ud = user_data;

    if (position == -1)
        return TRUE;

    MirbookingOccupant occupant = { .target = target, .position = position, .mirna = mirna };

    return !g_hash_table_contains (ud->blacklist, &occupant);
}

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

    if (self->priv->filter_user_data_destroy)
    {
        self->priv->filter_user_data_destroy (self->priv->filter_user_data);
    }

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_default_score_table_parent_class)->finalize (object);
}

static gsize
_pow4 (gsize n)
{
    return 1 << (2 * n);
}

static gfloat
_get_subsequence_score (SparseMatrix     *scores,
                        MirbookingMirna  *mirna,
                        MirbookingTarget *target,
                        gsize             position,
                        gsize             mirna_position,
                        gsize             subsequence_offset,
                        gsize             subsequence_len)
{
    g_return_val_if_fail (_pow4 (subsequence_len) == scores->shape[0], INFINITY);

    // the subsequence does not fit in the target
    if (position + mirna_position - subsequence_offset + 1 > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)))
    {
        return INFINITY;
    }

    // the subsequence does not fit in the mirna
    if (subsequence_offset + subsequence_len > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (mirna)))
    {
        return INFINITY;
    }

    gssize subsequence_i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), subsequence_offset, subsequence_len);
    gssize subsequence_j = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + mirna_position - subsequence_offset - subsequence_len + 1, subsequence_len);

    if (subsequence_i == -1 || subsequence_j == -1)
    {
        return INFINITY;
    }

    gfloat score = sparse_matrix_get_float (scores,
                                            subsequence_i,
                                            subsequence_j);

    return score;
}

static gboolean
is_g_bulge (MirbookingTarget *target, gsize position)
{
    return position + 4 + 4 <= mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)) &&
        toupper (*mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), position + 3, 1)) == 'G';
}

static gdouble
binding_energy (gdouble *G, gsize n)
{
    guint i;
    gdouble Gtot = 0, Z = 0;
    for (i = 0; i < n; i++)
    {
        if (isfinite (G[i]))
        {
            Gtot += exp (-G[i] / (R * T)) * G[i];
            Z += exp (-G[i] / (R * T));
        }
    }

    if (Z == 0)
    {
        return INFINITY;
    }

    return Gtot / Z;
}

static gdouble
binding_probability (gdouble *G, gsize n, gsize j)
{
    guint i;
    gdouble Z = 0;
    for (i = 0; i < n; i++)
    {
        if (isfinite (G[i]))
        {
            Z += exp (-G[i] / (R * T));
        }
    }

    if (Z == 0)
    {
        return INFINITY;
    }

    return exp (-G[j] / (R * T)) / Z;
}

static gboolean
compute_score (MirbookingScoreTable *score_table,
               MirbookingMirna      *mirna,
               MirbookingTarget     *target,
               gsize                 position,
               MirbookingScore      *score,
               GError              **error)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (score_table);

    MirbookingScore ret = {0};

    gdouble A_score = 0.0f;
    if (position + SEED_LENGTH + 1 <= mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)) &&
        toupper (*mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (target), position + SEED_LENGTH, 1)) == 'A')
    {
        A_score = T1_ADENOSINE_SCORE;
    }

    gfloat seed_score = _get_subsequence_score (&self->priv->seed_scores,
                                                mirna,
                                                target,
                                                position,
                                                SEED_OFFSET + SEED_LENGTH - 1,
                                                SEED_OFFSET,
                                                SEED_LENGTH);

    // Seed mismatches penalty to the arrival rate (Salomon et al. 2015)
    //
    // The two weights are obtained by maximum likelihood fit on a broad set of
    // experimental data. The expected values were around 0.4 and 0.9 for the 4
    // first and 3 last nucleotides of the seed respectively, which are quite
    // consistent with the obtained ones.
    gdouble w[7] = {0.38, 0.38, 0.38, 0.34, 0.71, 0.71, 0.71};

    const guint8 *_seed = mirbooking_sequence_get_subsequence (MIRBOOKING_SEQUENCE (mirna),
                                                               SEED_OFFSET,
                                                               SEED_LENGTH);

    guint8 seed[7];
    memcpy (seed, _seed, 7);

    const guint8 *target_seed_rc = mirbooking_sequence_get_subsequence_rc (MIRBOOKING_SEQUENCE (target),
                                                                           position ,
                                                                           SEED_LENGTH);

    // base forward rate which is modulated by seed mismatches
    ret.kf = KF;

    guint i;
    for (i = 0; i < 7; i++)
    {
        if ((toupper (seed[i]) == 'U' ? 'T' : toupper (seed[i])) != target_seed_rc[i])
        {
            ret.kf *= w[i];
        }
    }

    if (is_g_bulge (target, position))
    {
        gsize i, j;
        i = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), SEED_OFFSET, SEED_LENGTH);
        j = (1l << (2 * 4)) * mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position, 3) +
            mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + 4, 4);

        gdouble Z[2] = {seed_score, sparse_matrix_get_float (&self->priv->seed_scores, i, j) + G_BULGED_SEED_SCORE};
        seed_score = binding_energy (Z, 2);
    }

    gdouble supplementary_score = 0;
    if (self->priv->supplementary_model == MIRBOOKING_DEFAULT_SCORE_TABLE_SUPPLEMENTARY_MODEL_WEE_ET_AL_2012)
    {
        gfloat supplementary_score_ = _get_subsequence_score (&self->priv->supplementary_scores,
                                                              mirna,
                                                              target,
                                                              position,
                                                              SEED_OFFSET + SEED_LENGTH - 1,
                                                              SEED_OFFSET + SEED_LENGTH + 4,
                                                              4);
        gdouble z[2] = {0, supplementary_score_};
        supplementary_score = binding_energy (z, 2);

        // require at least the 3' supplementary bindings for cleavage
        ret.kcleave = binding_probability (z, 5, 1) * KCLEAVE;
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

        // A box
        A_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                        mirna,
                                        target,
                                        position,
                                        SEED_OFFSET + SEED_LENGTH - 1,
                                        SEED_OFFSET + SEED_LENGTH,
                                        3);

        // B box
        B_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                        mirna,
                                        target,
                                        position,
                                        SEED_OFFSET + SEED_LENGTH - 1,
                                        SEED_OFFSET + SEED_LENGTH + 3,
                                        3);

        // C box
        C_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                        mirna,
                                        target,
                                        position,
                                        SEED_OFFSET + SEED_LENGTH - 1,
                                        SEED_OFFSET + SEED_LENGTH + 6,
                                        3);

        // D box
        // TODO: consider all remaining nucleotides
        D_box = _get_subsequence_score (&self->priv->supplementary_scores,
                                        mirna,
                                        target,
                                        position,
                                        SEED_OFFSET + SEED_LENGTH - 1,
                                        SEED_OFFSET + SEED_LENGTH + 9,
                                        3);

        gdouble z[5] = {0,
            B_box,
            B_box + C_box,
            B_box + C_box + A_box,
            B_box + C_box + A_box + D_box};

        supplementary_score = binding_energy (z, 5);

        // at least A-box for cleavage
        ret.kcleave = (binding_probability (z, 5, 3) + binding_probability (z, 5, 4)) * KCLEAVE;
    }

    gdouble Kd = 1e12 * exp ((A_score + seed_score + supplementary_score + mirbooking_target_get_accessibility_score (target, position) + AGO2_SCORE) / (R * T));

    ret.kr = ret.kf * Kd;

    ret.krelease = ret.kr;

    // overall turnover rate
    ret.kcat = 1. / (1. / ret.kcleave + 1. / ret.krelease);

    *score = ret;

    return TRUE;
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

    g_autofree gsize *_positions = NULL;

    if (self->priv->filter != NULL &&
        !self->priv->filter (self, mirna, target, -1, self->priv->filter_user_data))
    {
        *positions = NULL;
        *positions_len = 0;
        return TRUE;
    }

    gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)) - SEED_LENGTH + 1;

    gsize p;
    gsize k = 0;
    for (p = 0; p < seq_len; p++)
    {
        if (self->priv->filter == NULL ||
            self->priv->filter (self, mirna, target, p, self->priv->filter_user_data))
        {
            MirbookingScore score;
            mirbooking_score_table_compute_score (MIRBOOKING_SCORE_TABLE (self), mirna, target, p, &score, NULL);
            if (MIRBOOKING_SCORE_IS_FINITE (score))
            {
                _positions = g_realloc (_positions, (k + 1) * sizeof (gsize));
                _positions[k++] = p;
            }
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
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
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
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
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
                                     g_param_spec_boxed ("seed-scores", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SUPPLEMENTARY_MODEL,
                                     g_param_spec_uint ("supplementary-model", "", "", 0, 2, MIRBOOKING_DEFAULT_SCORE_TABLE_DEFAULT_SUPPLEMENTARY_MODEL, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
    g_object_class_install_property (object_class,
                                     PROP_SUPPLEMENTARY_SCORES,
                                     g_param_spec_boxed ("supplementary-scores", "", "", G_TYPE_BYTES, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
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

/**
 * mirbooking_default_score_table_set_filter:
 */
void
mirbooking_default_score_table_set_filter (MirbookingDefaultScoreTable *self,
                                           MirbookingDefaultScoreTableFilter filter,
                                           gpointer user_data,
                                           GDestroyNotify destroy_func)
{
    self->priv->filter           = filter;
    self->priv->filter_user_data = user_data;
    self->priv->filter_user_data_destroy = destroy_func;
}
