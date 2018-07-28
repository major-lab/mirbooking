#include "mirbooking-broker.h"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <pb.h>
#include <string.h>
#if HAVE_OPENMP
#include <omp.h>
#endif
#include <sparse.h>
#include <stdio.h>

#define R 1.987203611e-3
#define T 310.15

typedef struct _IntegratorMeta
{
    const gchar *name;
    gint         steps;
    gdouble      c[4];
    gdouble      w[4];
} IntegratorMeta;

typedef gdouble IntegratorState[4];

#define PREDICT(x) (step_size*(integrator_meta.c[step] * x[step]))
#define CORRECT(x) (step_size*(integrator_meta.w[0] * x[0] + integrator_meta.w[1] * x[1] + integrator_meta.w[2] * x[2] + integrator_meta.w[3] * x[3]))

const IntegratorMeta INTEGRATOR_META[] =
{
     {"euler",       1, {0, 0,   0,   0}, {1,       0,       0,       0}},
     {"heuns",       2, {0, 1,   0,   0}, {0.5,     0.5,     0,       0}},
     {"runge-kutta", 4, {0, 0.5, 0.5, 1}, {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}}
};

typedef struct _MirbookingTargetSitesScores
{
    gfloat    *seed_scores;
    gsize     *positions;
    gsize      positions_len;
    GPtrArray *occupants;
} MirbookingTargetSitesScores;

void
mirbooking_target_sites_scores_clear (MirbookingTargetSitesScores *tss)
{
    g_free (tss->seed_scores);
    g_free (tss->positions);
    g_ptr_array_unref (tss->occupants);
}

typedef struct
{
    gdouble                     kappa;
    gdouble                     lambda;
    MirbookingBrokerIntegrator  integrator;
    gsize                       prime5_footprint;
    gsize                       prime3_footprint;
    MirbookingScoreTable       *score_table;

    GPtrArray                  *targets;
    GPtrArray                  *mirnas;

    GHashTable                 *quantification; // #MirbookingSequence -> #gfloat (initial quantity)

    /* all the target sites, stored contiguously */
    GArray     *target_sites;
    GHashTable *target_sites_by_target;
    GPtrArray  *target_sites_by_target_index;
    GArray     *target_sites_scores; // (target, mirna) -> positions, scores and occupants in the target
    MirbookingOccupant *occupants;

    /* state of the system */
    gdouble *E0; // number of mirnas
    gdouble *E;  // len(E0)
    gdouble *S0; // number of targets (sites are identical w.r.t. to initial concentration)
    gdouble *S;  // number of sites among targets (not aligned) thus much bigger than len(S0)
    gdouble *ES; // len(E) * len(S)
    gdouble *P;  // len(ES)

    /* state of the integrator */
    IntegratorState *dEdt;  // like E
    IntegratorState *dSdt;  // like S
    IntegratorState *dESdt; // len(E) * len(S)
    IntegratorState *dPdt;

    /* steady-state solver */
    SparseSolver *solver;
    SparseMatrix *J;         // len(ES) * len(ES)
    gdouble      *ES_delta;  // len(ES)
    gdouble      *dESdt_rhs; // len(ES)

    // footprint breaks convexity, which is necessary for Newton method, so we
    // solve a relaxed problem and iteratively refine the solution
    gsize effective_prime5_footprint;
    gsize effective_prime3_footprint;

} MirbookingBrokerPrivate;

struct _MirbookingBroker
{
    GObject                  parent_instance;
    MirbookingBrokerPrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingBroker, mirbooking_broker, G_TYPE_OBJECT)

static void
mirbooking_broker_init (MirbookingBroker *self)
{
    self->priv = g_new0 (MirbookingBrokerPrivate, 1);
    self->priv->targets = g_ptr_array_new_with_free_func (g_object_unref);
    self->priv->mirnas = g_ptr_array_new_with_free_func (g_object_unref);
    self->priv->quantification = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                   (GEqualFunc) mirbooking_sequence_equal);
}

enum
{
    PROP_KAPPA = 1,
    PROP_LAMBDA,
    PROP_5PRIME_FOOTPRINT,
    PROP_3PRIME_FOOTPRINT,
    PROP_SCORE_TABLE
};

static void
mirbooking_broker_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingBroker *self = MIRBOOKING_BROKER (object);

    switch (property_id)
    {
        case PROP_KAPPA:
            self->priv->kappa = g_value_get_double (value);
            break;
        case PROP_LAMBDA:
            self->priv->lambda = g_value_get_double (value);
            break;
        case PROP_5PRIME_FOOTPRINT:
            self->priv->prime5_footprint = g_value_get_uint (value);
            break;
        case PROP_3PRIME_FOOTPRINT:
            self->priv->prime3_footprint = g_value_get_uint (value);
            break;
        case PROP_SCORE_TABLE:
            self->priv->score_table = g_value_dup_object (value);
            break;
        default:
            g_assert_not_reached ();
    }
}

static void
mirbooking_broker_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    MirbookingBroker *self = MIRBOOKING_BROKER (object);

    switch (property_id)
    {
        case PROP_KAPPA:
            g_value_set_double (value, self->priv->kappa);
            break;
        case PROP_LAMBDA:
            g_value_set_double (value, self->priv->lambda);
            break;
        case PROP_5PRIME_FOOTPRINT:
            g_value_set_uint (value, self->priv->prime5_footprint);
            break;
        case PROP_3PRIME_FOOTPRINT:
            g_value_set_uint (value, self->priv->prime3_footprint);
            break;
        case PROP_SCORE_TABLE:
            g_value_set_object (value, self->priv->score_table);
        default:
            g_assert_not_reached ();
    }
}

static void
mirbooking_occupant_init (MirbookingOccupant* self, MirbookingMirna *mirna, gsize k, gfloat score)
{
    self->mirna    = g_object_ref (mirna);
    self->quantity = 0;
    self->k        = k;
    self->score    = score;
}


static void
mirbooking_occupant_clear (MirbookingOccupant *self)
{
    g_object_unref (self->mirna);
}

static void
mirbooking_target_site_clear (MirbookingTargetSite *self)
{
    g_object_unref (self->target);
    g_slist_free_full (self->occupants, (GDestroyNotify) mirbooking_occupant_clear);
}

static void
mirbooking_broker_finalize (GObject *object)
{
    MirbookingBroker *self = MIRBOOKING_BROKER (object);

    if (self->priv->score_table)
    {
        g_object_unref (self->priv->score_table);
    }

    g_hash_table_unref (self->priv->quantification);

    g_ptr_array_unref (self->priv->targets);
    g_ptr_array_unref (self->priv->mirnas);

    // iteration-specific stuff
    if (self->priv->target_sites)
    {
        g_array_free (self->priv->target_sites, TRUE);
        g_hash_table_unref (self->priv->target_sites_by_target);
        g_ptr_array_unref (self->priv->target_sites_by_target_index);
        g_array_free (self->priv->target_sites_scores, TRUE);
        g_free (self->priv->occupants);

        g_free (self->priv->E0);
        g_free (self->priv->E);
        g_free (self->priv->S0);
        g_free (self->priv->S);
        g_free (self->priv->ES);
        g_free (self->priv->P);

        g_free (self->priv->dEdt);
        g_free (self->priv->dSdt);
        g_free (self->priv->dESdt);
        g_free (self->priv->dPdt);

        sparse_solver_free (self->priv->solver);
        sparse_matrix_clear (self->priv->J);
        g_free (self->priv->ES_delta);
        g_free (self->priv->dESdt_rhs);
    }

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_broker_parent_class)->finalize (object);
}

static void
mirbooking_broker_class_init (MirbookingBrokerClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->set_property = mirbooking_broker_set_property;
    object_class->get_property = mirbooking_broker_get_property;
    object_class->finalize     = mirbooking_broker_finalize;

    g_object_class_install_property (object_class, PROP_KAPPA,
                                     g_param_spec_double ("kappa", "", "", 0, G_MAXDOUBLE, MIRBOOKING_BROKER_DEFAULT_KAPPA, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_LAMBDA,
                                     g_param_spec_double ("lambda", "", "", 0, G_MAXDOUBLE, MIRBOOKING_BROKER_DEFAULT_LAMBDA, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_5PRIME_FOOTPRINT,
                                     g_param_spec_uint ("prime5-footprint", "", "", 0, G_MAXUINT, MIRBOOKING_BROKER_DEFAULT_5PRIME_FOOTPRINT, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_3PRIME_FOOTPRINT,
                                     g_param_spec_uint ("prime3-footprint", "", "", 0, G_MAXUINT, MIRBOOKING_BROKER_DEFAULT_3PRIME_FOOTPRINT, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_SCORE_TABLE,
                                     g_param_spec_object ("score-table", "", "", MIRBOOKING_TYPE_SCORE_TABLE, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

/**
 * mirbooking_broker_new:
 *
 * Returns: (transfer full): A plain #Mirbooking instance
 */
MirbookingBroker *
mirbooking_broker_new (void)
{
    return g_object_new (MIRBOOKING_BROKER_TYPE, NULL);
}

void
mirbooking_broker_set_kappa (MirbookingBroker *self,
                             gdouble           kappa)
{
    self->priv->kappa = kappa;
}

void
mirbooking_broker_set_integrator (MirbookingBroker           *self,
                                  MirbookingBrokerIntegrator  integrator)
{
    self->priv->integrator = integrator;
}

void
mirbooking_broker_set_5prime_footprint (MirbookingBroker *self,
                                        gsize       footprint)
{
    self->priv->prime5_footprint = footprint;
}

void
mirbooking_broker_set_3prime_footprint (MirbookingBroker *self,
                                        gsize       footprint)
{
    self->priv->prime3_footprint = footprint;
}

/**
 * mirbooking_broker_get_score_table:
 * Obtain the #MirbookingScoreTable used by this for computing duplex scores.
 *
 * Returns: (transfer none)
 */
MirbookingScoreTable *
mirbooking_broker_get_score_table (MirbookingBroker *self)
{
    return self->priv->score_table;
}

void
mirbooking_broker_set_score_table (MirbookingBroker *self, MirbookingScoreTable *score_table)
{
    self->priv->score_table = score_table;
}

union gfloatptr
{
    gfloat   f;
    gpointer p;
};

static gfloat
gfloat_from_gpointer (gpointer ptr)
{
    union gfloatptr flt = { .p = ptr };
    return flt.f;
}

static gpointer
gpointer_from_gfloat (gfloat flt)
{
    union gfloatptr ptr = { .f = flt };
    return ptr.p;
}

/**
 * mirbooking_get_sequence_quantity:
 * @sequence: A #MirbookingSequence to retrieve quantity
 */
gfloat
mirbooking_broker_get_sequence_quantity (MirbookingBroker *self, MirbookingSequence *sequence)
{
    return gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, sequence));
}

/**
 * mirbooking_set_sequence_quantity:
 * @sequence: A #MirbookingSequence being quantified for the
 * upcoming execution
 *
 * Note that no new sequence can be added once #mirbooking_broker_evaluate has
 * been called.
 */
void
mirbooking_broker_set_sequence_quantity (MirbookingBroker *self, MirbookingSequence *sequence, gfloat quantity)
{
    g_return_if_fail (MIRBOOKING_IS_MIRNA (sequence) || MIRBOOKING_IS_TARGET (sequence));
    g_return_if_fail (self->priv->target_sites == NULL || g_hash_table_contains (self->priv->quantification, sequence));

    if (g_hash_table_insert (self->priv->quantification,
                             sequence,
                             gpointer_from_gfloat (quantity)));
    {
        if (MIRBOOKING_IS_MIRNA (sequence))
        {
            g_ptr_array_add (self->priv->mirnas,
                             g_object_ref (sequence));
        }
        else
        {
            g_ptr_array_add (self->priv->targets,
                             g_object_ref (sequence));
        }
    }
}

static int
sequence_cmp_desc (MirbookingSequence **a, MirbookingSequence **b, MirbookingBroker *mirbooking)
{
    gfloat a_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, *a));
    gfloat b_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, *b));
    return (a_quantity < b_quantity) - (a_quantity > b_quantity);
}

/**
 * Compute the footprint window in which two microRNA can have overlapping
 * footprint at this position.
 */
static void
_mirbooking_broker_get_footprint_window (MirbookingBroker            *self,
                                         const MirbookingTargetSite  *target_site,
                                         gsize                        prime5_footprint,
                                         gsize                        prime3_footprint,
                                         const MirbookingTargetSite **from_target_site,
                                         const MirbookingTargetSite **to_target_site)
{
    gsize window = prime5_footprint + prime3_footprint;
    const MirbookingTargetSite *_to_target_site;

    // find the lower target site
    *from_target_site = target_site - MIN (window, target_site->position);

    // find the upper target site
    _to_target_site = MIN (target_site + window,
                           &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len - 1));

    // we might overlap preceeding or following target sites
    while (_to_target_site->target != target_site->target)
    {
        --_to_target_site;
    }

    *to_target_site = _to_target_site;
}

gdouble
_mirbooking_broker_get_target_site_vacancy (MirbookingBroker           *self,
                                            const MirbookingTargetSite *target_site,
                                            gsize                       prime5_footprint,
                                            gsize                       prime3_footprint,
                                            gdouble                     S0)
{
    gdouble vacancy = 1;

    const MirbookingTargetSite *from_target_site, *to_target_site;
    _mirbooking_broker_get_footprint_window (self,
                                             target_site,
                                             prime5_footprint,
                                             prime3_footprint,
                                             &from_target_site,
                                             &to_target_site);

    // minimize vacancy around the footprint
    const MirbookingTargetSite *nearby_target_site;
    for (nearby_target_site = from_target_site; nearby_target_site <= to_target_site; nearby_target_site++)
    {
        vacancy *= 1 - (nearby_target_site->quantity / S0);
    }

    return vacancy;
}

/**
 * mirbooking_broker_get_target_site_vacancy:
 * Returns: The number of unoccupied copies of the target site
 */
gdouble
mirbooking_broker_get_target_site_vacancy (MirbookingBroker *self, const MirbookingTargetSite *target_site)
{
    MirbookingTarget *target = target_site->target;
    gfloat available_target_quantity = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, target));
    return _mirbooking_broker_get_target_site_vacancy (self,
                                                       target_site,
                                                       self->priv->prime5_footprint,
                                                       self->priv->prime3_footprint,
                                                       available_target_quantity);
}

typedef struct _MirbookingScoredTargetSite
{
    MirbookingTargetSite *target_site;
    gfloat                score;
} MirbookingScoredTargetSite;

/**
 * 'EVALUATE' will evaluate the partial derivatives and store them for the
 * given 'step'.
 * 'STEP' will perform the step described by 'step'.
 * 'UPDATE' will combine the previous, individual step computations into a
 * final update
 */
typedef enum _MirbookingBrokerIterMode
{
    MIRBOOKING_BROKER_ITER_MODE_EVALUATE,

    /* numerical integrator */
    MIRBOOKING_BROKER_ITER_MODE_STEP,

    /* just update */
    MIRBOOKING_BROKER_ITER_MODE_UPDATE,
} MirbookingBrokerIterMode;

static gboolean
_mirbooking_broker_prepare_step (MirbookingBroker *self)
{
    guint64 prepare_begin = g_get_monotonic_time ();

    gsize target_sites_len = 0;

    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (self->priv->score_table != NULL, FALSE);

    // sort internal targets and mirnas in descending quantity
    g_ptr_array_sort_with_data (self->priv->targets, (GCompareDataFunc) sequence_cmp_desc, self);
    g_ptr_array_sort_with_data (self->priv->mirnas, (GCompareDataFunc) sequence_cmp_desc, self);

    gint i;
    for (i = 0; i < self->priv->targets->len; i++)
    {
        target_sites_len += mirbooking_sequence_get_sequence_length (g_ptr_array_index (self->priv->targets, i));
    }

    // prepare an contiguous array
    self->priv->target_sites = g_array_sized_new (FALSE,
                                                  FALSE,
                                                  sizeof (MirbookingTargetSite),
                                                  target_sites_len);

    // automatically clear the target sites
    g_array_set_clear_func (self->priv->target_sites,
                            (GDestroyNotify) mirbooking_target_site_clear);

    // bookkeep each target site
    self->priv->target_sites_by_target = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                           (GEqualFunc) mirbooking_sequence_equal);

    // intialize sites
    for (i = 0; i < self->priv->targets->len; i++)
    {
        MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);

        g_hash_table_insert (self->priv->target_sites_by_target,
                             target,
                             &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len));

        gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target));

        gsize position;
        for (position = 0; position < seq_len; position++)
        {
            MirbookingTargetSite target_site;
            target_site.target    = g_object_ref (target);
            target_site.position  = position;
            target_site.occupants = NULL;
            target_site.quantity  = 0;
            g_array_append_val (self->priv->target_sites, target_site);
        }
    }

    // memoize in a array-friendly way the target sites
    self->priv->target_sites_by_target_index = g_ptr_array_new ();
    self->priv->S0 = g_new (gdouble, self->priv->targets->len);
    for (i = 0; i < self->priv->targets->len; i++)
    {
        gdouble q = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification,
                                                               g_ptr_array_index (self->priv->targets, i)));
        g_ptr_array_add (self->priv->target_sites_by_target_index,
                         g_hash_table_lookup (self->priv->target_sites_by_target, g_ptr_array_index (self->priv->targets, i)));
        self->priv->S0[i] = q;
    }

    // memoize unassigned quantities
    self->priv->E0 = g_new (gdouble, self->priv->mirnas->len);
    self->priv->E = g_new (gdouble, self->priv->mirnas->len);
    gint j;
    for (j = 0; j < self->priv->mirnas->len; j++)
    {
        gdouble q = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification,
                                                               g_ptr_array_index (self->priv->mirnas, j)));
        self->priv->E0[j] = q;
        self->priv->E[j]  = q;
    }

    // memoize score vectors
    self->priv->target_sites_scores = g_array_sized_new (FALSE,
                                                         FALSE,
                                                         sizeof (MirbookingTargetSitesScores),
                                                         self->priv->targets->len * self->priv->mirnas->len);
    g_array_set_clear_func (self->priv->target_sites_scores,
                            (GDestroyNotify) mirbooking_target_sites_scores_clear);

    // compute scores
    gsize occupants_len = 0;
    #pragma omp parallel for collapse(2) reduction(+:occupants_len)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);
            MirbookingMirna *mirna   = g_ptr_array_index (self->priv->mirnas, j);
            MirbookingTargetSitesScores *seed_scores = &g_array_index (self->priv->target_sites_scores, MirbookingTargetSitesScores, i * self->priv->mirnas->len + j);
            seed_scores->seed_scores = mirbooking_score_table_compute_scores (self->priv->score_table,
                                                                              mirna,
                                                                              target,
                                                                              &seed_scores->positions,
                                                                              &seed_scores->positions_len,
                                                                              NULL);
            occupants_len += seed_scores->positions_len;
        }
    }

    // pre-allocate occupants in contiguous memory
    self->priv->occupants = g_new (MirbookingOccupant, occupants_len);

    // create occupants
    gsize k = 0;
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingMirna *mirna = g_ptr_array_index (self->priv->mirnas, j);

            MirbookingTargetSite *target_sites = g_ptr_array_index (self->priv->target_sites_by_target_index,
                                                                    i);

            MirbookingTargetSitesScores *seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                       MirbookingTargetSitesScores,
                                                                       i * self->priv->mirnas->len + j);

            seed_scores->occupants = g_ptr_array_sized_new (seed_scores->positions_len);
            gint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];
                MirbookingOccupant *occupant = &self->priv->occupants[k];
                mirbooking_occupant_init (occupant, mirna, k++, seed_scores->seed_scores[p]);
                target_site->occupants = g_slist_prepend (target_site->occupants, occupant);
                g_ptr_array_add (seed_scores->occupants, occupant);
            }
        }
    }

    g_debug ("Number of duplexes: %lu", k);

    // count nnz entries in the Jacobian
    gsize nnz = 0;
    #pragma omp parallel for collapse(2) reduction(+:nnz)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTargetSite *target_sites = g_ptr_array_index (self->priv->target_sites_by_target_index,
                                                                    i);
            MirbookingTargetSitesScores *seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                       MirbookingTargetSitesScores,
                                                                       i * self->priv->mirnas->len + j);

            // footprint interactions
            gint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                // substitute targets
                gint z;
                for (z = 0; z < self->priv->targets->len; z++)
                {
                    MirbookingTargetSitesScores *alternative_seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                                           MirbookingTargetSitesScores,
                                                                                           z * self->priv->mirnas->len + j);
                    nnz += alternative_seed_scores->positions_len;
                }

                // substitute miRNAs (in the footprint)
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];

                const MirbookingTargetSite *from_target_site, *to_target_site;
                _mirbooking_broker_get_footprint_window (self,
                                                         target_site,
                                                         self->priv->prime5_footprint,
                                                         self->priv->prime3_footprint,
                                                         &from_target_site,
                                                         &to_target_site);

                const MirbookingTargetSite *ts;
                for (ts = from_target_site; ts <= to_target_site; ts++)
                {
                    GSList *occupant_list;
                    for (occupant_list = ts->occupants; occupant_list != NULL; occupant_list = occupant_list->next)
                    {
                        MirbookingOccupant *occupant = occupant_list->data;
                        if (occupant->mirna != g_ptr_array_index (self->priv->mirnas, j))
                        {
                            nnz += 1;
                        }
                    }
                }
            }
        }
    }

    g_debug ("nnz: %lu, sparsity: %f%%", nnz,  100.0f - 100.0f * (gdouble) nnz / (gdouble) (k*k));

    // allocate memory for the integrator and the solver
    self->priv->S  = g_new0 (gdouble, k);
    self->priv->ES = g_new0 (gdouble, k);
    self->priv->P  = g_new0 (gdouble, k);

    self->priv->dEdt  = g_new0 (IntegratorState, k);
    self->priv->dSdt  = g_new0 (IntegratorState, k);
    self->priv->dESdt = g_new0 (IntegratorState, k);
    self->priv->dPdt  = g_new0 (IntegratorState, k);

    self->priv->solver = sparse_solver_new (SPARSE_SOLVER_METHOD_SUPERLU);

    self->priv->J = g_new0 (SparseMatrix, 1);

    size_t shape[2] = {k, k};
    sparse_matrix_init (self->priv->J,
                        SPARSE_MATRIX_STORAGE_CSR,
                        shape,
                        nnz);
    self->priv->ES_delta = g_new0 (gdouble, k);
    self->priv->dESdt_rhs = g_new0 (gdouble, k);

    self->priv->effective_prime5_footprint = 0;
    self->priv->effective_prime3_footprint = 0;

    // initialize the sparse slots beforehand because it is not thread-safe and
    // we want to keep the in order for fast access
    // FIXME: this is slow..?
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTargetSite *target_sites = g_ptr_array_index (self->priv->target_sites_by_target_index,
                                                                    i);
            MirbookingTargetSitesScores *seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                       MirbookingTargetSitesScores,
                                                                       i * self->priv->mirnas->len + j);

            gint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                // footprint interactions
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];
                MirbookingOccupant *occupant = g_ptr_array_index (seed_scores->occupants, p);

                // substitute target
                gint z;
                for (z = 0; z < self->priv->targets->len; z++)
                {
                    MirbookingTargetSitesScores *alternative_seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                                           MirbookingTargetSitesScores,
                                                                                           z * self->priv->mirnas->len + j);
                    gint w;
                    for (w = 0; w < alternative_seed_scores->occupants->len; w++)
                    {
                        MirbookingOccupant *other_occupant = g_ptr_array_index (alternative_seed_scores->occupants,
                                                                                w);

                        sparse_matrix_set_value (self->priv->J,
                                                 occupant->k,
                                                 other_occupant->k,
                                                 0);
                    }
                }

                // substitute miRNAs (in the footprint)
                const MirbookingTargetSite *from_target_site, *to_target_site;
                _mirbooking_broker_get_footprint_window (self,
                                                         target_site,
                                                         self->priv->prime5_footprint,
                                                         self->priv->prime3_footprint,
                                                         &from_target_site,
                                                         &to_target_site);

                const MirbookingTargetSite *ts;
                for (ts = from_target_site; ts <= to_target_site; ts++)
                {
                    GSList *occupant_list;
                    for (occupant_list = ts->occupants; occupant_list != NULL; occupant_list = occupant_list->next)
                    {
                        MirbookingOccupant *other_occupant = occupant_list->data;
                        if (other_occupant->mirna != g_ptr_array_index (self->priv->mirnas, j))
                        {
                            sparse_matrix_set_value (self->priv->J,
                                                     occupant->k,
                                                     other_occupant->k,
                                                     0);
                        }
                    }
                }
            }
        }
    }

    g_debug ("Prepared the first step in %lums", 1000 * (g_get_monotonic_time () - prepare_begin) / G_USEC_PER_SEC);

    return TRUE;
}

static gboolean
_mirbooking_broker_step (MirbookingBroker             *self,
                         MirbookingBrokerStepMode      step_mode,
                         MirbookingBrokerIterMode      iter_mode,
                         gdouble                       step_size,
                         gint                          step,
                         gdouble                      *norm,
                         GError                      **error)
{
    gint i, j;
    gdouble _norm = 0;

    IntegratorMeta integrator_meta = INTEGRATOR_META[self->priv->integrator];

    gsize prime5_footprint = self->priv->effective_prime5_footprint;
    gsize prime3_footprint = self->priv->effective_prime3_footprint;

    guint64 step_begin = g_get_monotonic_time ();

    #pragma omp parallel for collapse(2) reduction(+:_norm)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);

            MirbookingTargetSite *target_sites = g_ptr_array_index (self->priv->target_sites_by_target_index,
                                                                    i);

            g_assert (target_sites->target == target);
            g_assert_cmpint (target_sites->position, ==, 0);

            MirbookingMirna *mirna = g_ptr_array_index (self->priv->mirnas, j);

            // fetch free energies for candidate MREs
            MirbookingTargetSitesScores *seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                       MirbookingTargetSitesScores,
                                                                       self->priv->mirnas->len * i + j);

            gint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                MirbookingOccupant *occupant = g_ptr_array_index (seed_scores->occupants, p);
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];

                g_assert (target_site->target == target);
                g_assert_cmpint (target_site->position, ==, seed_scores->positions[p]);
                g_assert (occupant->mirna == mirna);

                if (iter_mode == MIRBOOKING_BROKER_ITER_MODE_EVALUATE)
                {
                    gfloat score = occupant->score;

                    // The dissociation constant is derived from the duplex's Gibbs
                    // free energy
                    gdouble duplex_Kd = 1e9 * exp (score / (R*T)); // nM

                    // compute the dissociation constant
                    gdouble Kd = duplex_Kd / self->priv->kappa; // nM -> FPKM

                    // compute the vacancy
                    self->priv->S[occupant->k] = self->priv->S0[i] * _mirbooking_broker_get_target_site_vacancy (self,
                                                                                                                 target_site,
                                                                                                                 prime5_footprint,
                                                                                                                 prime3_footprint,
                                                                                                                 self->priv->S0[i]);

                    self->priv->ES[occupant->k] = occupant->quantity;
                    self->priv->P[occupant->k]  = occupant->cleaved_quantity;

                    // Here we apply a Michaelis-Menten kinetics
                    gdouble E  = self->priv->E[j];
                    gdouble S  = self->priv->S[occupant->k];
                    gdouble ES = self->priv->ES[occupant->k];
                    gdouble P  = self->priv->P[occupant->k];

                    const gdouble kf = self->priv->lambda / Kd;
                    const gdouble kr = self->priv->lambda * Kd;

                    // TODO: we need to implement a catalytic model
                    const gdouble kcat = self->priv->lambda * 0;

                    // dE/dt
                    gdouble dEdt  = -kf * E * S + kr * ES + kcat * ES;
                    gdouble dSdt  = -kf * E * S + kr * ES;
                    gdouble dESdt =  kf * E * S - kr * ES - kcat * ES;
                    gdouble dPdt  =                         kcat * ES;

                    _norm += pow (dESdt, 2) + pow (dEdt, 2) + pow (dSdt, 2) + pow (dPdt, 2);

                    if (step_mode == MIRBOOKING_BROKER_STEP_MODE_INTEGRATE)
                    {
                        self->priv->dEdt[occupant->k][step]  = dEdt;
                        self->priv->dSdt[occupant->k][step]  = dSdt;
                        self->priv->dESdt[occupant->k][step] = dESdt;
                        self->priv->dPdt[occupant->k][step]  = dPdt;
                    }
                    else if (step_mode == MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE)
                    {
                        // TODO: build a fast index for retrieving the Jacobian entries

                        // substitute target
                        gint z;
                        for (z = 0; z < self->priv->targets->len; z++)
                        {
                            MirbookingTargetSitesScores *alternative_seed_scores = &g_array_index (self->priv->target_sites_scores,
                                                                                                   MirbookingTargetSitesScores,
                                                                                                   z * self->priv->mirnas->len + j);
                            gint w;
                            for (w = 0; w < alternative_seed_scores->occupants->len; w++)
                            {
                                MirbookingOccupant *other_occupant = g_ptr_array_index (alternative_seed_scores->occupants, w);

                                g_assert (other_occupant->mirna == g_ptr_array_index (self->priv->mirnas, j));

                                gdouble dEdES = -1; // always
                                gboolean ofp = ABS (seed_scores->positions[p] - alternative_seed_scores->positions[w]) <= (prime3_footprint + prime5_footprint);
                                gdouble dSdES = (z == i && ofp) ? -S / (self->priv->S0[i] - target_site->quantity) : 0;

                                gdouble dESdES = kf * (E * dSdES + S * dEdES) - kr * (occupant->k == other_occupant->k ? 1 : 0);
                                sparse_matrix_set_value (self->priv->J,
                                                         occupant->k,
                                                         other_occupant->k,
                                                         dESdES);
                            }
                        }

                        // substitute miRNA (in the footprint)
                        const MirbookingTargetSite *from_target_site, *to_target_site, *ts;
                        _mirbooking_broker_get_footprint_window (self,
                                                                 target_site,
                                                                 prime5_footprint,
                                                                 prime3_footprint,
                                                                 &from_target_site,
                                                                 &to_target_site);
                        for (ts = from_target_site; ts <= to_target_site; ts++)
                        {
                            GSList *other_occupant_list;
                            for (other_occupant_list = ts->occupants; other_occupant_list != NULL; other_occupant_list = other_occupant_list->next)
                            {
                                MirbookingOccupant *other_occupant = other_occupant_list->data;

                                gdouble dEdES = occupant->mirna == other_occupant->mirna ? -1 : 0;

                                gdouble dSdES = -S / (self->priv->S0[i] - ts->quantity);

                                gdouble dESdES = kf * (E * dSdES + S * dEdES) - kr * (occupant->k == other_occupant->k ? 1 : 0);

                                sparse_matrix_set_value (self->priv->J,
                                                         occupant->k,
                                                         other_occupant->k,
                                                         dESdES);
                            }
                        }

                        self->priv->dESdt_rhs[occupant->k] = -dESdt;
                    }
                    else
                    {
                        g_assert_not_reached ();
                    }
                }
                else if (iter_mode == MIRBOOKING_BROKER_ITER_MODE_STEP)
                {
                    // step mode
                    // d[E]/dt
                    #pragma omp atomic update
                    self->priv->E[j] += PREDICT (self->priv->dEdt[occupant->k]);

                    // d[S]/dt
                    #pragma omp atomic update
                    target_site->quantity -= PREDICT (self->priv->dSdt[occupant->k]);

                    // d[ES]/dt
                    occupant->quantity += PREDICT (self->priv->dESdt[occupant->k]);

                    // d[P]/dt
                    #pragma omp atomic update
                    target_site->cleaved_quantity += PREDICT (self->priv->dPdt[occupant->k]);
                    occupant->cleaved_quantity    += PREDICT (self->priv->dPdt[occupant->k]);
                }
                else if (iter_mode == MIRBOOKING_BROKER_ITER_MODE_UPDATE)
                {
                    if (step_mode == MIRBOOKING_BROKER_STEP_MODE_INTEGRATE)
                    {
                        // d[E]/dt
                        #pragma omp atomic update
                        self->priv->E[j] += CORRECT (self->priv->dEdt[occupant->k]);

                        // d[S]/dt
                        #pragma omp atomic update
                        target_site->quantity -= CORRECT (self->priv->dSdt[occupant->k]);

                        // d[ES]/dt
                        occupant->quantity += CORRECT (self->priv->dESdt[occupant->k]);

                        // d[P]/dt
                        #pragma omp atomic update
                        target_site->cleaved_quantity += CORRECT (self->priv->dPdt[occupant->k]);
                        occupant->cleaved_quantity    += CORRECT (self->priv->dPdt[occupant->k]);
                    }
                    else if (step_mode == MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE)
                    {
                        // d[E]/dt
                        #pragma omp atomic update
                        self->priv->E[j] -=  step_size * self->priv->ES_delta[occupant->k];

                        // d[S]/dt
                        #pragma omp atomic update
                        target_site->quantity +=  step_size * self->priv->ES_delta[occupant->k];

                        // d[ES]/dt
                        occupant->quantity += step_size * self->priv->ES_delta[occupant->k];
                    }
                    else
                    {
                        g_assert_not_reached ();
                    }
                }
                else
                {
                    g_assert_not_reached ();
                }
            }
        }
    }

    if (iter_mode == MIRBOOKING_BROKER_ITER_MODE_EVALUATE)
    {
        g_debug ("Evaluated in %lums", 1000 * (g_get_monotonic_time () - step_begin) / G_USEC_PER_SEC);
    }
    else
    {
        g_debug ("Stepped in %lums", 1000 * (g_get_monotonic_time () - step_begin) / G_USEC_PER_SEC);
    }

    if (norm)
    {
        *norm = sqrt (_norm);
    }

    // refinement
    // if we're sufficiently close, we increase the footprint
    // FIXME: provide API for tuning this parameter
    if (norm                                                                  &&
        *norm <= 1e4                                                          &&
        self->priv->effective_prime5_footprint < self->priv->prime5_footprint &&
        self->priv->effective_prime3_footprint < self->priv->prime3_footprint)
    {
        self->priv->effective_prime5_footprint = self->priv->prime5_footprint;
        self->priv->effective_prime3_footprint = self->priv->prime3_footprint;
        g_debug ("Effective footprint is now [%lu, %lu]",
                 self->priv->effective_prime5_footprint,
                 self->priv->effective_prime3_footprint);

        // reevaluate to get the correct norm
        return _mirbooking_broker_step (self,
                                        step_mode,
                                        iter_mode,
                                        step_size,
                                        step,
                                        norm,
                                        error);
    }

    if (step_mode == MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE &&
        iter_mode == MIRBOOKING_BROKER_ITER_MODE_EVALUATE)
    {
        guint64 solve_begin = g_get_monotonic_time ();
        gboolean ret;
        ret = sparse_solver_solve (self->priv->solver,
                                   self->priv->J,
                                   self->priv->ES_delta,
                                   self->priv->dESdt_rhs);

        if (!ret)
        {
            g_set_error_literal (error,
                                 MIRBOOKING_ERROR,
                                 MIRBOOKING_ERROR_FAILED,
                                 "The solve step has failed.");
            return FALSE;
        }

        g_debug ("Solved in %lums", 1000 * (g_get_monotonic_time () - solve_begin) / G_USEC_PER_SEC);
    }

    return TRUE;
}

/**
 * mirbooking_broker_evaluate:
 * @norm (out): The norm of the system.
 *
 * Evaluate the current state of the system.
 *
 * Returns: %TRUE if the evaluation is successful, otherwise %FALSE and @error
 * is set accordingly
 */
gboolean
mirbooking_broker_evaluate (MirbookingBroker          *self,
                            MirbookingBrokerStepMode   step_mode,
                            gdouble                   *norm,
                            GError                   **error)
{
    static gsize init;
    if (g_once_init_enter (&init))
    {
        g_return_val_if_fail (_mirbooking_broker_prepare_step (self),
                              FALSE);
        g_once_init_leave (&init, 1);
    }

    // initial step with criterion evalation
    return _mirbooking_broker_step (self,
                                    step_mode,
                                    MIRBOOKING_BROKER_ITER_MODE_EVALUATE,
                                    0,
                                    0,
                                    norm,
                                    error);
}

/**
 * mirbooking_broker_step:
 *
 * Perform a step based on the last evaluation from #mirbooking_broker_evaluate.
 *
 * Returns: %TRUE on success, otherwise @error is set
 */
gboolean
mirbooking_broker_step (MirbookingBroker         *self,
                        MirbookingBrokerStepMode  step_mode,
                        gdouble                   step_size,
                        GError                   **error)
{
    IntegratorMeta integrator_meta = INTEGRATOR_META[self->priv->integrator];

    /* multi-step integrator methods */
    if (step_mode == MIRBOOKING_BROKER_STEP_MODE_INTEGRATE)
    {
        gint step;
        for (step = 1; step < integrator_meta.steps; step++)
        {
            // take the last predicted step
            g_return_val_if_fail (_mirbooking_broker_step (self,
                                                           step_mode,
                                                           MIRBOOKING_BROKER_ITER_MODE_STEP,
                                                           step_size,
                                                           step - 1,
                                                           NULL,
                                                           error), FALSE);

            // evaluate this position
            g_return_val_if_fail (_mirbooking_broker_step (self,
                                                           step_mode,
                                                           MIRBOOKING_BROKER_ITER_MODE_EVALUATE,
                                                           step_size,
                                                           step,
                                                           NULL,
                                                           error), FALSE);

            // restore to the original state with a negative step
            g_return_val_if_fail (_mirbooking_broker_step (self,
                                                           step_mode,
                                                           MIRBOOKING_BROKER_ITER_MODE_STEP,
                                                           -step_size,
                                                           step - 1,
                                                           NULL,
                                                           error), FALSE);
        }
    }

    // final update
    return _mirbooking_broker_step (self,
                                    step_mode,
                                    MIRBOOKING_BROKER_ITER_MODE_UPDATE,
                                    step_size,
                                    0,
                                    NULL,
                                    error);
}

/**
 * mirbooking_broker_get_target_sites:
 *
 * Obtain the computed #MirbookingTargetSite array by this #MirbookingBroker.
 *
 * Returns: (element-type MirbookingTargetSite) (transfer none): A view of the
 * computed #MirbookingTargetSite
 */
GArray *
mirbooking_broker_get_target_sites (MirbookingBroker *self)
{
    g_return_val_if_fail (self != NULL, NULL);
    g_return_val_if_fail (self->priv->target_sites != NULL, NULL);

    return self->priv->target_sites;
}

static gdouble
mean_silencing_by_number_of_sites (guint k)
{
    return 1 - pow (2, -0.0392 * k + 0.0054);
}

/**
 * mirbooking_broker_get_target_silencing:
 * @target: The #MirbookingTarget for which we are retrieving the silencing
 *
 * Compute the silencing of the target according to the Poisson-Binomial
 * distribution.
 *
 * First we estimate the number of occupant per target using the
 * Poisson-Binomial PMF and then we use an empirical regression for computing
 * the expected miRNA-induced silencing by weighting each discrete outcomes.
 *
 * For individual miRNA-induced silencing, this value should be distributed
 * proportionally to the individual occupancy.
 *
 * Returns: The silencing computed across all the occupied sites
 */
gdouble
mirbooking_broker_get_target_silencing (MirbookingBroker *self, MirbookingTarget *target)
{
    gdouble target_silencing = 0;
    MirbookingTargetSite *target_site = g_hash_table_lookup (self->priv->target_sites_by_target, target);
    gfloat target_quantity = mirbooking_broker_get_sequence_quantity (self, MIRBOOKING_SEQUENCE (target));

    g_return_val_if_fail (target_site != NULL, 0.0f);

    g_autoptr (GArray) probability_by_position = g_array_new (FALSE, FALSE, sizeof (gdouble));
    while (target_site < &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len) &&
           target_site->target == target)
    {
        if (target_site->quantity > 0)
        {
            gdouble proba = target_site->quantity / target_quantity;
            g_array_append_val (probability_by_position, proba);
        }

        ++target_site;
    }

    // TODO: compute the joint probability with footprint
    PoissonBinomial pb;
    pb_init (&pb,
             (gdouble*) probability_by_position->data,
             probability_by_position->len);
    guint k;
    #pragma omp parallel for reduction(+:target_silencing)
    for (k = 0; k <= probability_by_position->len; k++)
    {
        target_silencing += pb_pmf (&pb, k) * mean_silencing_by_number_of_sites (k);
    }

    pb_destroy (&pb);

    return target_silencing;
}
