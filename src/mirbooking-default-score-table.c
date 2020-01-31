#include "mirbooking-default-score-table.h"
#include "mirbooking-score-table-private.h"

#include <math.h>
#include <sparse.h>
#include <ctype.h>

typedef enum _StepType
{
    STEP_TYPE_MISMATCH,
    STEP_TYPE_MATCH,
    STEP_TYPE_GUIDE_BULGE,
    STEP_TYPE_TARGET_BULGE
} StepType;

/**
 * TODO: find a better name than step
 * @step_type: Type of step for the basepairing
 * @other_index: Index of the facing nucleotide in the other sequence
 */
typedef struct _Step
{
    StepType step_type;
    gsize    other_index;
} Step;

/**
 * Size of the dynamic programming window.
 */
#define L 21

/**
 * Shorthands for the work matrix.
 */
#define N  W[STEP_TYPE_MISMATCH]
#define M  W[STEP_TYPE_MATCH]
#define G  W[STEP_TYPE_GUIDE_BULGE]
#define _T W[STEP_TYPE_TARGET_BULGE]

/**
 * Shorthands for computing the minimum
 */
#define MIN3(a,b,c) MIN(a,MIN(b,c))
#define MIN4(a,b,c,d) MIN(a,MIN3(b,c,d))
#define MIN5(a,b,c,d,e) MIN(a,MIN4(b,c,d,e))

// intercept of the model (kcat/mol)
static gdouble kd_intercept = -9.721654427;

// cost for base pairing (kcal/mol)
// first index correspond to A, C, G and U
static gdouble E[4][L] =
{
    {
        -0.262183966,  0.105068373, -0.194341493, -0.081299257, -0.056447922,
        -0.239280099,  0.289289226,  0.186643095, -0.467251671, -0.229392879,
        -0.000300916, -0.149312855, -0.238999514, -0.099953415,  0.250045566,
        -0.30537667,   0.020955728, -0.070139814,  0.283916873, -1.143964436,
        -0.084198464
    },
    {
         0.320439139, -0.291260098, -0.263605051, -0.503914064, -0.262558287,
         0.429292601,  0.291326539, -0.09482782,  -0.034409021,  0.177011832,
        -0.047748145, -0.246900948,  0.181355517,  0.435540947, -0.161533569,
        -0.056905923,  0.101554701,  0.177882236,  0.249809152, -0.15524378,
         0
    },
    {
        -0.02541896,  0.074863363, -0.228212485,  0.382582093, 0.387863719,
         0.193156646, 0.189490511,  0.312623476, -0.238306128, 0.007632637,
         0.037705876, 0.202597916,  0.022289246,  0.180763037, 0.393816973,
        -0.10773831,  0.030143402,  0.018411937,  0.123897725, 0,
         0.16412786
    },
    {
        -0.032836213,  0.111328362, -1.417463397,  0.113766664,  0.04880133,
         0.551036696, -0.433057349,  0.354249053, -0.319907463, -0.288150253,
        -0.049861125, -0.10795833,   0.116736914,  0.223431981,  0.005152924,
        -0.067627737, -0.37731765,   0.031798643, -0.199313252, -0.015552099,
         0
    }
};

// bulge initiation and extension costs (kcal/mol)
static gdouble GBopening[L] =
{
    [0 ... 1]   = INFINITY,
    [2 ... 7]   = 0.120531131,
    [8 ... 12]  = 0.120531131,
    [13 ... 17] = 0.120531131,
    [18 ... 20] = INFINITY
};

static gdouble GBextension[L] =
{
    [0 ... 1]   = INFINITY,
    [2 ... 7]   = 0.122161608,
    [8 ... 12]  = 0.183027923,
    [13 ... 17] = 0,
    [18 ... 20] = INFINITY
};

static gdouble TBopening[L] =
{
    [0 ... 1]   = INFINITY,
    [2 ... 7]   = 1.004161049,
    [8 ... 12]  = 0.280616988,
    [13 ... 17] = 0.423176342,
    [18 ... 20] = INFINITY
};

static gdouble TBextension[L] =
{
    [0 ... 1]   = INFINITY,
    [2 ... 7]   = 0.029888898,
    [8 ... 12]  = 0.041007346,
    [13 ... 17] = 0.016248539,
    [18 ... 20] = INFINITY
};

static gdouble Pinit = 0.25766107;

// energy contribution weight
static gdouble kd_rnafold_relative_structure_multiplier = -0.204106906;

// ln(kcat)
// A <-> G or C <-> U
static gdouble kcat_transition[L] =
{
    -0.436733244, -0.319418493, -1.04651559,  -0.19822662,  -0.337674491,
    -0.529848392, -0.009803101, -0.582230978, -2.471311869, -2.343283254,
    -3.678198774, -2.331934398, -2.503647487, -0.852270673, -1.171758466,
    -1.181498281, -0.6811052,    0.03960261,   0.541476655,  0.720549257,
     0.251319632
};

// all other substitution
static gdouble kcat_transversion[L] =
{
    -0.010764119, -0.37042212,  -1.205981456, -0.404518801, 0.211674842,
    -1.347747417, -1.323046111, -1.201088586, -3.133902013, -5.408542199,
    -3.900980244, -1.021786583, -2.801551852, -0.803252792, -0.472283749,
    -0.604761391, -0.077204424,  0.795072416,  0.833168247,  0,
     0
};

static gdouble kcat_guide_bulge[L] =
{
      INFINITY,     0,           -1.631620059, -3.207753066, -3.007446355,
     -3.484551356, -5.105734757, -1.551700878, -3.205373271, -4.305629982,
     -3.15345819,  -3.38643224,  -3.890167923, -4.221593717, -1.875042509,
     -1.571151863, -0.693338747,  0,            0,            0,
      INFINITY
};

static gdouble kcat_target_bulge[L] =
{
     0,           -0.239746034, -1.246740243, -0.099418119, -0.167929447,
    -0.442086238, -1.280269319, -2.73455944,  -6.360570778, -6.24289057,
    -6.076824342, -3.31987613,  -1.37852684,  -1.677347536, -0.661427609,
    -0.069380976, 0.352987063,   0,            0.276728221,  0,
     INFINITY
};

typedef struct
{
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

    return ES >= cutoff_filter_user_data->cutoff && ((ES / S0) >= cutoff_filter_user_data->relative_cutoff);
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

static gboolean
is_complementary (gsize i, gsize j)
{
    // since A = 0, C = 1, G = 2 and U = 3
    return i + j == 3;
}

static gboolean
is_transition (gsize i, gsize j)
{
    // since A = 0, C = 1, G = 2 and U = 3
    return (i != j) && ((i + j == 2) || (i + j == 4));
}

static gboolean
is_transversion (gsize i, gsize j)
{
    return (i != j) && !is_transition (i, j);
}

static void
backtrack (gdouble (*W)[4][L][L],
           Step    (*solution_target)[L],
           gsize    i,
           Step    (*solution_mirna)[L],
           gsize    j)
{
    // find the minimum
    gdouble S = INFINITY;
    StepType k;
    for (k = STEP_TYPE_MISMATCH; k <= STEP_TYPE_TARGET_BULGE; k++)
    {
        if ((*W)[k][i][j] < S)
        {
            // min found
            S = (*W)[k][i][j];
            (*solution_target)[i] = (Step) {.step_type = k, .other_index = j};
            (*solution_mirna)[j]  = (Step) {.step_type = k, .other_index = i};
        }
    }

    // there should be at least one non-infinity solution
    g_assert_cmpfloat (S, !=, INFINITY);

    // backtrack based on the origin
    switch ((*solution_target)[i].step_type)
    {
        case STEP_TYPE_MATCH:
            if (i > 0 && j > 0)
                return backtrack (W, solution_target, i - 1, solution_mirna, j - 1);
            break;
        case STEP_TYPE_MISMATCH:
            if (i > 0 && j > 0)
                return backtrack (W, solution_target, i - 1, solution_mirna, j - 1);
            break;
        case STEP_TYPE_GUIDE_BULGE:
            if (j > 0)
                return backtrack (W, solution_target, i, solution_mirna, j - 1);
            break;
        case STEP_TYPE_TARGET_BULGE:
            if (i > 0)
                return backtrack (W, solution_target, i - 1, solution_mirna, j);
            break;
    }
}

static gboolean
compute_score (MirbookingScoreTable *score_table,
               MirbookingMirna      *mirna,
               MirbookingTarget     *target,
               gsize                 position,
               MirbookingScore      *score,
               GError              **error)
{
    MirbookingScore ret = {0};

    g_return_val_if_fail (L <= mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (mirna)), FALSE);

    if (position + L > mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target)))
    {
        ret.kr = INFINITY;
        *score = ret;
        return TRUE;
    }

    gint i, j;

    gdouble W[4][L][L];

    for (j = 0; j < L; j++)
    {
        gint t1  = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + L - 1, 1);
        M[0][j]  = E[t1][j];
        N[0][j]  = INFINITY;
        G[0][j]  = INFINITY;
        _T[0][j] = INFINITY;
    }

    for (i = 0; i < L; i++)
    {
        gint tj  = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + (L - i - 1), 1);
        M[i][0]  = E[tj][i];
        N[i][0]  = INFINITY;
        G[i][0]  = INFINITY;
        _T[i][0] = INFINITY;
    }

    for (i = 1; i < L; i++)
    {
        gssize ti = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + (L - i - 1), 1);
        g_assert_cmpint (ti, !=, -1);

        for (j = 1; j < L; j++)
        {
            gssize gj = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), j, 1);
            g_assert_cmpint (gj, !=, -1);

            if (is_complementary (ti, gj))
            {
                N[i][j]  = INFINITY;
                M[i][j]  = MIN5 (M[i - 1][j - 1], N[i - 1][j - 1] + Pinit, _T[i - 1][j - 1] + Pinit, G[i - 1][j - 1] + Pinit, Pinit);
                G[i][j]  = INFINITY;
                _T[i][j] = INFINITY;
                g_assert_cmpfloat (MIN4 (N[i][j], M[i][j], G[i][j], _T[i][j]), !=, INFINITY);
            }
            else
            {
                M[i][j]  = INFINITY;
                N[i][j]  = E[ti][j] + MIN (M[i - 1][j - 1], N[i - 1][j - 1]);
                G[i][j]  = MIN3 (M[i][j - 1] + GBopening[j], _T[i][j - 1] + GBextension[j], N[i][j - 1] + GBopening[j]);
                _T[i][j] = MIN3 (M[i - 1][j] + TBopening[j], _T[i - 1][j] + TBextension[j], N[i - 1][j] + TBopening[j]);
                g_assert_cmpfloat (MIN4 (N[i][j], M[i][j], G[i][j], _T[i][j]), !=, INFINITY);
            }
        }
    }

    gdouble deltaG = kd_intercept + MIN4 (M[L - 1][L - 1], N[L - 1][L - 1], G[L - 1][L - 1], _T[L - 1][L - 1]) + kd_rnafold_relative_structure_multiplier * 0.0 + mirbooking_target_get_accessibility_score (target, position);

    Step solution_target[L];
    Step solution_mirna[L];
    backtrack (&W, &solution_target, L - 1, &solution_mirna, L - 1);

    gdouble ln_kcat = 0;
    for (i = 0; i < L; i++)
    {
        gssize ti = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (target), position + (L - i - 1), 1);
        g_assert_cmpint (ti, !=, -1);

        gssize gj;

        switch (solution_target[i].step_type)
        {
            case STEP_TYPE_MISMATCH:
            case STEP_TYPE_MATCH:
                gj = mirbooking_sequence_get_subsequence_index (MIRBOOKING_SEQUENCE (mirna), solution_target[i].other_index, 1);
                g_assert_cmpint (gj, !=, -1);

                if (is_transition (ti, gj))
                {
                    ln_kcat += kcat_transition[i];
                }
                else if (is_transversion (ti, gj))
                {
                    ln_kcat += kcat_transversion[i];
                }
                break;
            case STEP_TYPE_TARGET_BULGE:
                // FIXME: check if offset is ok
                ln_kcat += kcat_target_bulge[i];
                break;
            case STEP_TYPE_GUIDE_BULGE:
                // FIXME: check if offset is ok
                ln_kcat += kcat_guide_bulge[i];
                break;
        }
    }

    g_print ("%f kcal/mol %.2e s^-1\n", deltaG, exp (ln_kcat));

    ret.kf = KF;
    ret.kr = ret.kf * 1e12 * exp (deltaG / (R * T));
    ret.kcleave = 0;
    ret.krelease = 0;
    ret.kcat = exp (ln_kcat);

    *score = ret;

    return TRUE;
}

static void
mirbooking_default_score_table_init (MirbookingDefaultScoreTable *self)
{
    self->priv = g_new0 (MirbookingDefaultScoreTablePrivate, 1);
}

static void
mirbooking_default_score_table_finalize (GObject *object)
{
    MirbookingDefaultScoreTable *self = MIRBOOKING_DEFAULT_SCORE_TABLE (object);
    g_free (self->priv);
    G_OBJECT_CLASS (mirbooking_default_score_table_parent_class)->finalize (object);
}

static void
mirbooking_default_score_table_class_init (MirbookingDefaultScoreTableClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->finalize = mirbooking_default_score_table_finalize;

    klass->parent_class.compute_score = compute_score;
}

/**
 * mirbooking_default_score_table_new:
 *
 * Returns: (transfer full)
 */
MirbookingDefaultScoreTable *
mirbooking_default_score_table_new ()
{
    return g_object_new (MIRBOOKING_TYPE_DEFAULT_SCORE_TABLE,
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
