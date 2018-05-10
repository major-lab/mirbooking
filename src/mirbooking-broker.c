#include "mirbooking-broker.h"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <pb.h>
#include <string.h>
#if HAVE_OPENMP
#include <omp.h>
#else
static gdouble
omp_get_wtime ()
{
    return (gdouble) g_get_monotonic_time () / (gdouble) G_USEC_PER_SEC;
}
#endif

#define R 1.987203611e-3
#define T 310.15

typedef struct _IntegratorMeta {
    gchar   *name;
    gint     steps;
    gdouble  c[4];
    gdouble  w[4];
} IntegratorMeta;

IntegratorMeta INTEGRATOR_META[] =
{
     {"euler",       1, {0, 0,   0,   0}, {1,       0,       0,       0}},
     {"heuns",       2, {0, 1,   0,   0}, {0.5,     0.5,     0,       0}},
     {"range-kutta", 4, {0, 0.5, 0.5, 1}, {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}}
};

typedef struct
{
    gdouble                    kappa;
    MirbookingBrokerIntegrator integrator;
    gdouble                     step_size;
    gdouble                     tolerance;
    guint64                     max_iterations;
    gsize                       prime5_footprint;
    gsize                       prime3_footprint;
    MirbookingScoreTable       *score_table;

    GPtrArray                  *targets;
    GPtrArray                  *mirnas;

    GHashTable                 *quantification; // #MirbookingSequence -> #gfloat (initial quantity)

    /* all the target sites, stored contiguously */
    GArray                     *target_sites;
    GHashTable                 *target_sites_by_target;
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
    PROP_CUTOFF,
    PROP_STEP_SIZE,
    PROP_MOMENTUM,
    PROP_TOLERANCE,
    PROP_MAX_ITERATIONS,
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
        case PROP_STEP_SIZE:
            self->priv->step_size = g_value_get_double (value);
            break;
        case PROP_TOLERANCE:
            self->priv->tolerance = g_value_get_double (value);
            break;
        case PROP_MAX_ITERATIONS:
            self->priv->max_iterations = g_value_get_uint64 (value);
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
        case PROP_STEP_SIZE:
            g_value_set_double (value, self->priv->step_size);
            break;
        case PROP_TOLERANCE:
            g_value_set_double (value, self->priv->tolerance);
            break;
        case PROP_MAX_ITERATIONS:
            g_value_set_uint64 (value, self->priv->max_iterations);
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

static MirbookingOccupant *
mirbooking_occupant_new (MirbookingMirna *mirna, gdouble quantity)
{
    MirbookingOccupant *ret = g_slice_new0 (MirbookingOccupant);
    ret->mirna              = g_object_ref (mirna);
    ret->quantity           = quantity;
    return ret;
}

static void
mirbooking_occupant_free (MirbookingOccupant *self)
{
    g_object_unref (self->mirna);
    g_slice_free (MirbookingOccupant, self);
}

static void
mirbooking_target_site_clear (MirbookingTargetSite *self)
{
    g_object_unref (self->target);
    g_slist_free_full (self->occupants, (GDestroyNotify) mirbooking_occupant_free);
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

    if (self->priv->target_sites)
    {
        g_array_free (self->priv->target_sites, TRUE);
        g_hash_table_unref (self->priv->target_sites_by_target);
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
    g_object_class_install_property (object_class, PROP_STEP_SIZE,
                                     g_param_spec_double ("step-size", "", "", 0, 1, MIRBOOKING_BROKER_DEFAULT_STEP_SIZE, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_TOLERANCE,
                                     g_param_spec_double ("tolerance", "", "", 0, G_MAXDOUBLE, MIRBOOKING_BROKER_DEFAULT_TOLERANCE, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_MAX_ITERATIONS,
                                     g_param_spec_uint64 ("max-iterations", "", "", 0, G_MAXUINT64, MIRBOOKING_BROKER_DEFAULT_MAX_ITERATIONS, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
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
mirbooking_broker_set_step_size (MirbookingBroker *self,
                                 gdouble           step_size)
{
    self->priv->step_size = step_size;
}

void
mirbooking_broker_set_tolerance (MirbookingBroker *self,
                                 gdouble           tolerance)
{
    self->priv->tolerance = tolerance;
}

void
mirbooking_broker_set_max_iterations (MirbookingBroker *self,
                                      guint64           max_iterations)
{
    self->priv->max_iterations = max_iterations;
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
 */
void
mirbooking_broker_set_sequence_quantity (MirbookingBroker *self, MirbookingSequence *sequence, gfloat quantity)
{
    g_return_if_fail (MIRBOOKING_IS_MIRNA (sequence) || MIRBOOKING_IS_TARGET (sequence));

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

gdouble
_mirbooking_broker_get_target_site_vacancy (MirbookingBroker           *self,
                                            const MirbookingTargetSite *target_site,
                                            gdouble                     available_target_quantity)
{
    gdouble vacancy = 1;

    gsize window = self->priv->prime5_footprint + self->priv->prime3_footprint;

    // find the lower target site
    const MirbookingTargetSite *from_target_site = target_site - MIN (window,
                                                                      target_site->position);

    // find the upper target site
    const MirbookingTargetSite *to_target_site = MIN (target_site + window + 1,
                                                      &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len));

    // minimize vacancy around the footprint
    const MirbookingTargetSite *nearby_target_site;
    for (nearby_target_site = from_target_site; nearby_target_site < to_target_site; nearby_target_site++)
    {
        // we might overlap preceeding or following target sites
        if (nearby_target_site->target == target_site->target)
        {
            vacancy *= 1 - (nearby_target_site->quantity / available_target_quantity);
        }
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
                                                       available_target_quantity);
}

typedef struct _MirbookingScoredTargetSite
{
    MirbookingTargetSite *target_site;
    gfloat                score;
} MirbookingScoredTargetSite;

typedef struct _MirbookingTargetSitesScores
{
    gfloat    *seed_scores;
    gsize     *positions;
    gsize      positions_len;
    GPtrArray *occupants;
} MirbookingTargetSitesScores;

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
    MIRBOOKING_BROKER_ITER_MODE_STEP,
    MIRBOOKING_BROKER_ITER_MODE_UPDATE
} MirbookingBrokerIterMode;

static gboolean
_mirbooking_broker_iter (MirbookingBroker             *self,
                         gdouble                      *available_target_quantities,
                         gdouble                      *available_mirna_quantities,
                         gdouble                      *unassigned_mirna_quantities,
                         GPtrArray                    *target_sites_by_target_index,
                         MirbookingTargetSitesScores  *target_sites_scores,
                         IntegratorMeta                integrator_meta,
                         gdouble                       step_size,
                         MirbookingBrokerIterMode      iter_mode,
                         gint                          step,
                         gdouble                      *jac_norm,
                         GError                      **error)
{
    gint i, j;
    gdouble _jac_norm = 0;

    #pragma omp parallel for collapse(2) reduction(+:_jac_norm)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);

            MirbookingTargetSite *target_sites = g_ptr_array_index (target_sites_by_target_index,
                                                                    i);

            g_assert (target_sites->target == target);
            g_assert_cmpint (target_sites->position, ==, 0);

            MirbookingMirna *mirna = g_ptr_array_index (self->priv->mirnas, j);

            gdouble available_target_quantity = available_target_quantities[i];
            gdouble available_mirna_quantity = available_mirna_quantities[j];

            // fetch free energies for candidate MREs
            MirbookingTargetSitesScores *seed_scores = &target_sites_scores[self->priv->mirnas->len * i + j];

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
                    gfloat seed_score = seed_scores->seed_scores[p];

                    // The dissociation constant is derived from the duplex's Gibbs
                    // free energy
                    gdouble duplex_Kd = 1e9 * exp (seed_score / (R*T)); // nM

                    // compute the dissociation constant
                    gdouble Kd = duplex_Kd / self->priv->kappa; // nM -> FPKM

                    // Here we apply a Michaelis-Menten kinetics
                    gdouble E = unassigned_mirna_quantities[j];
                    gdouble S = available_target_quantity * _mirbooking_broker_get_target_site_vacancy (self,
                                                                                                        target_site,
                                                                                                        available_target_quantity);

                    gdouble ES = occupant->quantity;

                    gdouble P = occupant->cleaved_quantity;

                    const gdouble kf = (1 / Kd);
                    const gdouble kr = Kd;

                    const gdouble kcat = 0;

                    // dE/dt
                    gdouble E_jac  = -kf * E * S + kr * ES + kcat * ES;
                    gdouble S_jac  = -kf * E * S + kr * ES;
                    gdouble ES_jac =  kf * E * S - kr * ES - kcat * ES;
                    gdouble P_jac  =                         kcat * ES;

                    _jac_norm += pow (ES_jac, 2) + pow (E_jac, 2) + pow (S_jac, 2) + pow (P_jac, 2);

                    occupant->E_jac[step]  = E_jac;
                    occupant->S_jac[step]  = S_jac;
                    occupant->ES_jac[step] = ES_jac;
                    occupant->P_jac[step]  = P_jac;
                }
                else if (iter_mode == MIRBOOKING_BROKER_ITER_MODE_STEP)
                {
                    // step mode
                    #define PREDICT(x) (step_size*(integrator_meta.c[step] * x[step]))

                    // d[E]/dt
                    #pragma omp atomic update
                    unassigned_mirna_quantities[j] += PREDICT (occupant->E_jac);

                    // d[S]/dt
                    #pragma omp atomic update
                    target_site->quantity -= PREDICT (occupant->S_jac);

                    // d[ES]/dt
                    occupant->quantity += PREDICT (occupant->ES_jac);

                    // d[P]/dt
                    #pragma omp atomic update
                    target_site->cleaved_quantity += PREDICT (occupant->P_jac);
                    occupant->cleaved_quantity    += PREDICT (occupant->P_jac);
                }
                else if (iter_mode == MIRBOOKING_BROKER_ITER_MODE_UPDATE)
                {
                    // final update
                    #define CORRECT(x) (step_size*(integrator_meta.w[0] * x[0] + integrator_meta.w[1] * x[1] + integrator_meta.w[2] * x[2] + integrator_meta.w[3] * x[3]))

                    // d[E]/dt
                    #pragma omp atomic update
                    unassigned_mirna_quantities[j] += CORRECT (occupant->E_jac);

                    // d[S]/dt
                    #pragma omp atomic update
                    target_site->quantity -= CORRECT (occupant->S_jac);

                    // d[ES]/dt
                    occupant->quantity += CORRECT (occupant->ES_jac);

                    // d[P]/dt
                    #pragma omp atomic update
                    target_site->cleaved_quantity += CORRECT (occupant->P_jac);
                    occupant->cleaved_quantity    += CORRECT (occupant->P_jac);
                }
                else
                {
                    g_assert_not_reached ();
                }
            }
        }
    }

    if (jac_norm)
    {
        *jac_norm = sqrt (_jac_norm);
    }

    return TRUE;
}

/**
 * mirbooking_broker_iter:
 * @step_size
 *
 * Perform a single iteration of the miRBooking algorithm
 *
 * Returns: %TRUE on success, otherwise @error is set
 */
gboolean
mirbooking_broker_iter (MirbookingBroker  *self,
                        gdouble            step_size,
                        GError           **error)
{
    return TRUE;
    /*
    return _mirbooking_broker_iter (self,
                                    step_size,
                                    0,
                                    NULL,
                                    error);
                                    */
}

gboolean
mirbooking_broker_run (MirbookingBroker *self, GError **error)
{
    gsize target_sites_len = 0;

    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (self->priv->score_table != NULL, FALSE);

    if (self->priv->target_sites != NULL)
    {
        g_set_error (error,
                     MIRBOOKING_ERROR,
                     MIRBOOKING_ERROR_FAILED,
                     "This broker has already run.");
        return FALSE;
    }

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

    self->priv->targets->len = self->priv->targets->len;
    self->priv->mirnas->len = self->priv->mirnas->len;

    // memoize in a array-friendly way the target sites
    g_autoptr (GPtrArray) target_sites_by_target_index = g_ptr_array_new ();
    gdouble *available_target_quantities = g_malloc (self->priv->targets->len * sizeof (gdouble));
    for (i = 0; i < self->priv->targets->len; i++)
    {
        gdouble q = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification,
                                                               g_ptr_array_index (self->priv->targets, i)));
        g_ptr_array_add (target_sites_by_target_index,
                         g_hash_table_lookup (self->priv->target_sites_by_target, g_ptr_array_index (self->priv->targets, i)));
        available_target_quantities[i] = q;
    }

    // memoize unassigned quantities
    gdouble *available_mirna_quantities = g_malloc (self->priv->mirnas->len * sizeof (gdouble));
    gdouble *unassigned_mirna_quantities = g_malloc (self->priv->mirnas->len * sizeof (gdouble));
    gint j;
    for (j = 0; j < self->priv->mirnas->len; j++)
    {
        gdouble q = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification,
                                                               g_ptr_array_index (self->priv->mirnas, j)));
        available_mirna_quantities[j]  = q;
        unassigned_mirna_quantities[j] = q;
    }

    // memoize score vectors
    MirbookingTargetSitesScores *target_sites_scores = g_new0 (MirbookingTargetSitesScores, self->priv->targets->len * self->priv->mirnas->len);

    #pragma omp parallel for collapse(2)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);
            MirbookingMirna *mirna   = g_ptr_array_index (self->priv->mirnas, j);
            MirbookingTargetSite *target_sites = g_hash_table_lookup (self->priv->target_sites_by_target,
                                                                      target);
            MirbookingTargetSitesScores *seed_scores = &target_sites_scores[i * self->priv->mirnas->len + j];
            seed_scores->seed_scores = mirbooking_score_table_compute_scores (self->priv->score_table,
                                                                              mirna,
                                                                              target,
                                                                              &seed_scores->positions,
                                                                              &seed_scores->positions_len,
                                                                              NULL);
            seed_scores->occupants = g_ptr_array_new ();
            gint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];

                MirbookingOccupant *occupant = mirbooking_occupant_new (mirna, 0);
                #pragma omp critical
                target_site->occupants = g_slist_prepend (target_site->occupants, occupant);
                g_ptr_array_add (seed_scores->occupants, occupant);
            }
        }
    }

    guint64 iteration = 0;
    gdouble jac_norm = G_MAXDOUBLE;

    IntegratorMeta integrator_meta = INTEGRATOR_META[self->priv->integrator];

    while (iteration < self->priv->max_iterations)
    {
        gdouble iteration_begin = omp_get_wtime ();
        iteration++;

        // initial step with criterion evalation
        _mirbooking_broker_iter (self,
                                 available_target_quantities,
                                 available_mirna_quantities,
                                 unassigned_mirna_quantities,
                                 target_sites_by_target_index,
                                 target_sites_scores,
                                 integrator_meta,
                                 self->priv->step_size,
                                 MIRBOOKING_BROKER_ITER_MODE_EVALUATE,
                                 0,
                                 &jac_norm,
                                 error);

        if (jac_norm <= MIRBOOKING_BROKER_DEFAULT_TOLERANCE)
        {
            break;
        }

        gint step;
        for (step = 1; step < integrator_meta.steps; step++)
        {
            // take the last predicted step
            _mirbooking_broker_iter (self,
                                     available_target_quantities,
                                     available_mirna_quantities,
                                     unassigned_mirna_quantities,
                                     target_sites_by_target_index,
                                     target_sites_scores,
                                     integrator_meta,
                                     self->priv->step_size,
                                     MIRBOOKING_BROKER_ITER_MODE_STEP,
                                     step - 1,
                                     NULL,
                                     error);

            // evaluate this position
            _mirbooking_broker_iter (self,
                                     available_target_quantities,
                                     available_mirna_quantities,
                                     unassigned_mirna_quantities,
                                     target_sites_by_target_index,
                                     target_sites_scores,
                                     integrator_meta,
                                     self->priv->step_size,
                                     MIRBOOKING_BROKER_ITER_MODE_EVALUATE,
                                     step,
                                     NULL,
                                     error);

            // restore to the original state with a negative step
            _mirbooking_broker_iter (self,
                                     available_target_quantities,
                                     available_mirna_quantities,
                                     unassigned_mirna_quantities,
                                     target_sites_by_target_index,
                                     target_sites_scores,
                                     integrator_meta,
                                     -self->priv->step_size,
                                     MIRBOOKING_BROKER_ITER_MODE_STEP,
                                     step - 1,
                                     NULL,
                                     error);
        }

        // final update
        _mirbooking_broker_iter (self,
                                 available_target_quantities,
                                 available_mirna_quantities,
                                 unassigned_mirna_quantities,
                                 target_sites_by_target_index,
                                 target_sites_scores,
                                 integrator_meta,
                                 self->priv->step_size,
                                 MIRBOOKING_BROKER_ITER_MODE_UPDATE,
                                 0,
                                 NULL,
                                 error);

        if (iteration % 1000 == 0)
        {
            g_debug ("iteration: %lu jac-norm: %.2e throughput: %.2f iter/sec",
                     iteration,
                     jac_norm,
                     1 / (omp_get_wtime () - iteration_begin));
        }
    }

    if (jac_norm > MIRBOOKING_BROKER_DEFAULT_TOLERANCE)
    {
        g_warning ("The maximum number of iterations were reached before complete convergence.");
    }
    else
    {
        g_debug ("Reached convergence in %lu iterations.", iteration);
    }

    g_free (target_sites_scores);
    g_free (available_mirna_quantities);
    g_free (unassigned_mirna_quantities);
    g_free (available_target_quantities);

    return TRUE;
}

static void
mirbooking_broker_run_in_thread (GTask        *task,
                                 gpointer      source_object,
                                 gpointer      task_data,
                                 GCancellable *cancellable)
{
    MirbookingBroker *broker = source_object;
    GError *err = NULL;

    if (mirbooking_broker_run (broker, &err))
    {
        g_task_return_boolean (task, TRUE);
    }
    else if (err != NULL)
    {
        g_task_return_error (task, err);
    }
    else
    {
        g_task_return_boolean (task, FALSE);
    }
}

/**
 * mirbooking_broker_run_async:
 *
 * Run the algorithm in a background thread via #GTask. The result can be
 * retrieved later on.
 */
void
mirbooking_broker_run_async (MirbookingBroker    *self,
                             GAsyncReadyCallback  callback,
                             gpointer             callback_data)
{
    GTask *task = g_task_new (self,
                              NULL,
                              callback,
                              callback_data);

    g_task_run_in_thread (task, mirbooking_broker_run_in_thread);
}

/**
 * mirbooking_broker_run_finish:
 */
gboolean
mirbooking_broker_run_finish (MirbookingBroker  *self,
                              GAsyncResult      *result,
                              GError           **error)
{
    g_return_val_if_fail (g_task_is_valid (result, self), FALSE);

    return g_task_propagate_boolean (G_TASK (result), error);
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
