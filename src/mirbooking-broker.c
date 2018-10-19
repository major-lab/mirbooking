#include "mirbooking-broker.h"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <pb.h>
#include <string.h>
#include <odeint.h>
#if HAVE_OPENMP
#include <omp.h>
#endif
#include <sparse.h>
#include <stdio.h>

typedef struct _MirbookingTargetSitesScores
{
    gsize     *positions;
    gsize      positions_len;
    GPtrArray *occupants;
} MirbookingTargetSitesScores;

void
mirbooking_target_sites_scores_clear (MirbookingTargetSitesScores *tss)
{
    g_free (tss->positions);
    g_ptr_array_unref (tss->occupants);
}

typedef struct
{
    gdouble                     kappa;
    gdouble                     lambda;
    gsize                       prime5_footprint;
    gsize                       prime3_footprint;
    MirbookingScoreTable       *score_table;

    GPtrArray                  *targets;
    GPtrArray                  *mirnas;

    GHashTable                 *quantification; // #MirbookingSequence -> #gfloat (initial quantity)

    /* whether or not the system has been initialized */
    gsize init;

    /* all the target sites, stored contiguously */
    GArray     *target_sites;
    GHashTable *target_sites_by_target;
    GPtrArray  *target_sites_by_target_index;
    GArray     *target_sites_scores; // (target, mirna) -> positions, scores and occupants in the target
    GArray     *occupants;

    /* state of the system */

    /* state of the system, which corresponds to the concatenation of the
     * various concentration vectors
     *
     * [E] [S] [ES] [P]
     */
    gdouble  t;
    gdouble *y;
    gsize    y_len;

    /* shortcuts over 'y' */
    gdouble *E;  // len(mirnas)
    gdouble *S;  // len(targets)
    gdouble *ES; // len(occupants)
    gdouble *P;  // len(S)

    /* odeint */
    OdeIntIntegrator *integrator;

    gdouble *F;

    /* shortcuts over 'F' */
    gdouble *dEdt;
    gdouble *dSdt;
    gdouble *dESdt;
    gdouble *dPdt;

    /* steady-state solver */
    SparseSolver *solver;
    /* for compactness and efficiency, the Jacobian is only defined over [ES] */
    SparseMatrix *J;        // len(ES) * len(ES)
    gdouble      *ES_delta; // len(ES)

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
    PROP_SCORE_TABLE,
    PROP_SPARSE_SOLVER
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
        case PROP_SPARSE_SOLVER:
            self->priv->solver = sparse_solver_new (g_value_get_enum (value));
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
mirbooking_occupant_init (MirbookingOccupant* self, MirbookingMirna *mirna, gdouble score, gdouble enzymatic_score)
{
    self->mirna = g_object_ref (mirna);
    self->score = score;
    self->enzymatic_score = enzymatic_score;
}


static void
mirbooking_occupant_clear (MirbookingOccupant *self)
{
    g_object_unref (self->mirna);
}

static size_t
_mirbooking_broker_get_occupant_index (MirbookingBroker *self, const MirbookingOccupant *occupant)
{
    return occupant - &g_array_index (self->priv->occupants, MirbookingOccupant, 0);
}

gdouble
_mirbooking_broker_get_occupant_quantity (MirbookingBroker *self, const MirbookingOccupant *occupant, const gdouble *ES)
{
    return ES[_mirbooking_broker_get_occupant_index (self, occupant)];
}

gdouble
mirbooking_broker_get_occupant_quantity (MirbookingBroker *self, const MirbookingOccupant *occupant)
{
    return _mirbooking_broker_get_occupant_quantity (self, occupant, self->priv->ES);
}

void
mirbooking_broker_set_occupant_quantity (MirbookingBroker *self, const MirbookingOccupant *occupant, gdouble quantity)
{
    g_return_if_fail (!self->priv->init);
    g_return_if_fail (quantity >= 0);

    self->priv->ES[_mirbooking_broker_get_occupant_index (self, occupant)] = quantity;
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
        g_array_free (self->priv->occupants, TRUE);

        g_free (self->priv->y);

        odeint_integrator_free (self->priv->integrator);

        sparse_solver_free (self->priv->solver);
        sparse_matrix_clear (self->priv->J);
        g_free (self->priv->ES_delta);
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
    g_object_class_install_property (object_class, PROP_SPARSE_SOLVER,
                                     g_param_spec_enum ("sparse-solver", "", "", MIRBOOKING_BROKER_SPARSE_SOLVER_ENUM, MIRBOOKING_BROKER_DEFAULT_SPARSE_SOLVER, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
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
mirbooking_broker_set_lambda (MirbookingBroker *self,
                              gdouble           lambda)
{
    self->priv->lambda = lambda;
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

void
mirbooking_broker_set_sparse_solver (MirbookingBroker *self,
                                     MirbookingBrokerSparseSolver sparse_solver)
{
    self->priv->solver = sparse_solver_new (sparse_solver);
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
    g_return_if_fail (self->priv->init == 0 || g_hash_table_contains (self->priv->quantification, sequence));
    g_return_if_fail (quantity > 0);

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

gdouble
mirbooking_broker_get_time (MirbookingBroker *self)
{
    return self->priv->t;
}

/**
 * mirbooking_broker_set_time:
 *
 * Set the initial time for the numerical integration.
 */
void
mirbooking_broker_set_time (MirbookingBroker *self, gdouble time)
{
    g_return_if_fail (!self->priv->init);
    self->priv->t = time;
}


static int
sequence_cmp_desc (MirbookingSequence **a, MirbookingSequence **b, MirbookingBroker *mirbooking)
{
    gfloat a_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, *a));
    gfloat b_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, *b));
    return (a_quantity < b_quantity) - (a_quantity > b_quantity);
}

/*
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
_mirbooking_broker_get_target_site_occupants_quantity (MirbookingBroker           *self,
                                                       const MirbookingTargetSite *target_site,
                                                       const gdouble              *ES)
{
    gdouble total_ES = 0;
    GSList *occupants_list;
    for (occupants_list = target_site->occupants; occupants_list != NULL; occupants_list = occupants_list->next)
    {
        MirbookingOccupant *occupant = occupants_list->data;
        total_ES += _mirbooking_broker_get_occupant_quantity (self, occupant, ES);
    }
    return total_ES;
}

gdouble
_mirbooking_broker_get_target_site_vacancy (MirbookingBroker           *self,
                                            const MirbookingTargetSite *target_site,
                                            gsize                       prime5_footprint,
                                            gsize                       prime3_footprint,
                                            gdouble                     St,
                                            const gdouble              *ES)
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
        vacancy *= 1 - (_mirbooking_broker_get_target_site_occupants_quantity (self, nearby_target_site, ES) / St);
    }

    g_assert_cmpfloat (vacancy, >=, 0);
    g_assert_cmpfloat (vacancy, <=, 1);

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
                                                       available_target_quantity,
                                                       self->priv->ES);
}

typedef struct _MirbookingScoredTargetSite
{
    MirbookingTargetSite *target_site;
    gfloat                score;
} MirbookingScoredTargetSite;

static gint
cmp_gsize (gconstpointer a, gconstpointer b)
{
    return *(const gsize*)a - *(const gsize*)b;
}

static gboolean
_mirbooking_broker_prepare_step (MirbookingBroker *self)
{
    guint64 prepare_begin = g_get_monotonic_time ();

    gsize target_sites_len = 0;

    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (self->priv->score_table != NULL, FALSE);
    g_return_val_if_fail (self->priv->solver != NULL, FALSE);

    // sort internal targets and mirnas in descending quantity
    g_ptr_array_sort_with_data (self->priv->targets, (GCompareDataFunc) sequence_cmp_desc, self);
    g_ptr_array_sort_with_data (self->priv->mirnas, (GCompareDataFunc) sequence_cmp_desc, self);

    gint i;
    #pragma omp parallel for reduction(+:target_sites_len)
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
            g_array_append_val (self->priv->target_sites, target_site);
        }
    }

    // memoize in a array-friendly way the target sites
    self->priv->target_sites_by_target_index = g_ptr_array_new ();
    for (i = 0; i < self->priv->targets->len; i++)
    {
        MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);
        MirbookingTargetSite *target_site = g_hash_table_lookup (self->priv->target_sites_by_target, target);
        g_ptr_array_add (self->priv->target_sites_by_target_index,
                         target_site);
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
    gint j;
    #pragma omp parallel for collapse(2) reduction(+:occupants_len)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);
            MirbookingMirna *mirna   = g_ptr_array_index (self->priv->mirnas, j);
            MirbookingTargetSitesScores *seed_scores = &g_array_index (self->priv->target_sites_scores, MirbookingTargetSitesScores, i * self->priv->mirnas->len + j);
            GError *err = NULL;
            if (mirbooking_score_table_compute_positions (self->priv->score_table,
                                                          mirna,
                                                          target,
                                                          &seed_scores->positions,
                                                          &seed_scores->positions_len,
                                                          &err))
            {
                occupants_len += seed_scores->positions_len;
            }
            else
            {
                g_critical ("%s", err->message);
            }
        }
    }

    // pre-allocate occupants in contiguous memory
    self->priv->occupants = g_array_sized_new (FALSE, FALSE, sizeof (MirbookingOccupant), occupants_len);

    // create occupants
    // FIXME: #pragma omp parallel for collapse(2)
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

            /* pre-allocate occupants */
            gsize occupants_offset;
            MirbookingOccupant _occupants[seed_scores->positions_len];
            #pragma omp critical
            {
                occupants_offset = self->priv->occupants->len;
                g_array_append_vals (self->priv->occupants,
                                     &_occupants,
                                     seed_scores->positions_len);
            }

            seed_scores->occupants = g_ptr_array_sized_new (seed_scores->positions_len);
            gint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];

                gdouble score = mirbooking_score_table_compute_score (self->priv->score_table,
                                                                      mirna,
                                                                      target_site->target,
                                                                      seed_scores->positions[p],
                                                                      NULL);

                gdouble enzymatic_score = mirbooking_score_table_compute_enzymatic_score (self->priv->score_table,
                                                                                          mirna,
                                                                                          target_site->target,
                                                                                          seed_scores->positions[p],
                                                                                          NULL);

                gsize k = occupants_offset + p;
                MirbookingOccupant *occupant = &g_array_index (self->priv->occupants, MirbookingOccupant, k);

                mirbooking_occupant_init (occupant,
                                          mirna,
                                          score / self->priv->kappa,
                                          enzymatic_score / self->priv->kappa);

                target_site->occupants = g_slist_prepend (target_site->occupants, occupant);
                g_ptr_array_add (seed_scores->occupants, occupant);
            }
        }
    }

    g_assert_cmpint (self->priv->occupants->len, ==, occupants_len);

    g_debug ("Number of duplexes: %u", self->priv->occupants->len);

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

                // substitute miRNAs (excluding this one as we consider it as a substitute target)
                nnz += g_slist_length (target_sites[seed_scores->positions[p]].occupants) - 1;
            }
        }
    }

    g_debug ("nnz: %lu, sparsity: %f%%", nnz,  100.0f - 100.0f * (gdouble) nnz / (gdouble) (self->priv->occupants->len * self->priv->occupants->len));

    self->priv->y_len = self->priv->mirnas->len + self->priv->targets->len + self->priv->occupants->len + self->priv->targets->len;

    // state of the system
    self->priv->y = g_new0 (gdouble, self->priv->y_len);

    // add shortcuts
    self->priv->E  = self->priv->y;
    self->priv->S  = self->priv->E  + self->priv->mirnas->len;
    self->priv->ES = self->priv->S  + self->priv->targets->len;
    self->priv->P  = self->priv->ES + self->priv->occupants->len;

    // setup initial conditions

    #pragma omp parallel for
    for (j = 0; j < self->priv->mirnas->len; j++)
    {
        gdouble q = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification,
                                                               g_ptr_array_index (self->priv->mirnas, j)));
        self->priv->E[j]  = q;
    }

    #pragma omp parallel for
    for (i = 0; i < self->priv->targets->len; i++)
    {
        MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);
        gdouble q = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification,
                                                               target));
        self->priv->S[i] = q;
    }

    // allocate memory for the integrator and the solver
    self->priv->F = g_new0 (gdouble, self->priv->y_len);

    // add shortcuts
    self->priv->dEdt  = self->priv->F;
    self->priv->dSdt  = self->priv->dEdt  + self->priv->mirnas->len;
    self->priv->dESdt = self->priv->dSdt  + self->priv->targets->len;
    self->priv->dPdt  = self->priv->dESdt + self->priv->occupants->len;

    self->priv->J = g_new0 (SparseMatrix, 1);
    self->priv->ES_delta = g_new0 (gdouble, self->priv->occupants->len);

    // integrator
    self->priv->integrator = odeint_integrator_new (ODEINT_METHOD_DORMAND_PRINCE,
                                                    &self->priv->t,
                                                    self->priv->y,
                                                    self->priv->y_len,
                                                    ODEINT_INTEGRATOR_DEFAULT_RTOL,
                                                    ODEINT_INTEGRATOR_DEFAULT_ATOL);

    size_t shape[2] = {self->priv->occupants->len, self->priv->occupants->len};
    sparse_matrix_init (self->priv->J,
                        SPARSE_MATRIX_STORAGE_CSR,
                        SPARSE_MATRIX_TYPE_DOUBLE,
                        shape,
                        nnz);

    // initialize the sparse slots beforehand because it is not thread-safe and
    // we want to keep the in order for fast access
    // TODO: find a way to remove the ordered clause
    #pragma omp parallel for collapse(2) ordered
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

                gsize colind[self->priv->J->shape[0]];
                gsize row_nnz = 0;

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

                        colind[row_nnz++] = _mirbooking_broker_get_occupant_index (self, other_occupant);
                    }
                }

                // substitute miRNAs
                GSList *occupant_list;
                for (occupant_list = target_site->occupants; occupant_list != NULL; occupant_list = occupant_list->next)
                {
                    MirbookingOccupant *other_occupant = occupant_list->data;
                    if (other_occupant->mirna != g_ptr_array_index (self->priv->mirnas, j))
                    {
                        colind[row_nnz++] = _mirbooking_broker_get_occupant_index (self, other_occupant);
                    }
                }

                // sort colind
                qsort (colind, row_nnz, sizeof (gsize), cmp_gsize);

                #pragma omp ordered
                sparse_matrix_reserve_range (self->priv->J,
                                             _mirbooking_broker_get_occupant_index (self, occupant),
                                             colind,
                                             row_nnz);
            }
        }
    }

    g_debug ("Prepared the first step in %lums", 1000 * (g_get_monotonic_time () - prepare_begin) / G_USEC_PER_SEC);

    return TRUE;
}

/*
 * Compute the system state.
 */
static void
_compute_F (double t, const double *y, double *F, void *user_data)
{
    MirbookingBroker *self = user_data;

    const gdouble *E  = y;
    const gdouble *S  = E  + self->priv->mirnas->len;
    const gdouble *ES = S  + self->priv->targets->len;
    const gdouble *P  = ES + self->priv->occupants->len;

    gdouble *dEdt  = F;
    gdouble *dSdt  = dEdt  + self->priv->mirnas->len;
    gdouble *dESdt = dSdt  + self->priv->targets->len;
    gdouble *dPdt  = dESdt + self->priv->occupants->len;

    gsize prime5_footprint = self->priv->prime5_footprint;
    gsize prime3_footprint = self->priv->prime3_footprint;

    memset (F, 0, sizeof (gdouble) * self->priv->y_len);

    gint i, j;
    #pragma omp parallel for collapse(2)
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

                gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

                g_assert (target_site->target == target);
                g_assert_cmpint (target_site->position, ==, seed_scores->positions[p]);
                g_assert (occupant->mirna == mirna);

                // The dissociation constant is derived from the duplex's Gibbs
                // free energy
                // We convert it to nanomolar and then using kappa, to the
                // concentration units which are typically FPKM from
                // high-throughput sequencing
                gdouble Kd = occupant->score;
                gdouble Km = occupant->enzymatic_score;

                gdouble kf   = self->priv->lambda / Kd;
                gdouble kr   = self->priv->lambda;
                gdouble kcat = kf * (Km - Kd);

                // FIXME: allow catalysis
                g_assert_cmpfloat (kcat, ==, 0.0);

                gdouble Stp = S[i] * _mirbooking_broker_get_target_site_vacancy (self,
                                                                                 target_site,
                                                                                 prime5_footprint,
                                                                                 prime3_footprint,
                                                                                 S[i],
                                                                                 ES);

                g_assert_cmpfloat (Stp, <=, S[i]);

                #pragma omp atomic
                dEdt[j] += -kf * E[j] * Stp + kr * ES[k] + kcat * ES[k];

                #pragma omp atomic
                dSdt[i] += -kcat * ES[k];

                dESdt[k] = kf * E[j] * Stp - kr * ES[k] - kcat * ES[k];

                #pragma omp atomic
                dPdt[i] += kcat * ES[k];
            }
        }
    }
}

/*
 * Compute the system Jacobian.
 */
static void
_compute_J (double t, const double *y, SparseMatrix *J, void *user_data)
{
    MirbookingBroker *self = user_data;

    const gdouble *E  = y;
    const gdouble *S  = y + self->priv->mirnas->len;
    const gdouble *ES = S + self->priv->targets->len;
    const gdouble *P  = ES + self->priv->occupants->len;

    gsize prime5_footprint = self->priv->prime5_footprint;
    gsize prime3_footprint = self->priv->prime3_footprint;

    gint i, j;
    #pragma omp parallel for collapse(2)
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

                gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

                g_assert (target_site->target == target);
                g_assert_cmpint (target_site->position, ==, seed_scores->positions[p]);
                g_assert (occupant->mirna == mirna);

                // The dissociation constant is derived from the duplex's Gibbs
                // free energy
                // We convert it to nanomolar and then using kappa, to the
                // concentration units which are typically FPKM from
                // high-throughput sequencing
                gdouble Kd = occupant->score;
                gdouble Km = occupant->enzymatic_score;

                g_assert_cmpfloat (Kd, >=, 0);
                g_assert_cmpfloat (Km, >=, 0);

                gdouble kf   = self->priv->lambda / Kd;
                gdouble kr   = self->priv->lambda;
                gdouble kcat = kf * (Km - Kd);

                g_assert_cmpfloat (kcat, ==, 0.0);

                gdouble Stp = S[i] * _mirbooking_broker_get_target_site_vacancy (self,
                                                                                 target_site,
                                                                                 prime5_footprint,
                                                                                 prime3_footprint,
                                                                                 S[i],
                                                                                 ES);

                // substitute target for the microRNA
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

                        gsize other_k = _mirbooking_broker_get_occupant_index (self, other_occupant);

                        gdouble dEdES  = -1; // always
                        gdouble dSdES  = (z == i && seed_scores->positions[p] == alternative_seed_scores->positions[w]) ? -Stp / (self->priv->S[i] - _mirbooking_broker_get_target_site_occupants_quantity (self, target_site, ES)) : 0;
                        gdouble dESdES = kf * (E[j] * dSdES + Stp * dEdES) - (kr + kcat) * (occupant == other_occupant ? 1 : 0);

                        sparse_matrix_set_double (self->priv->J,
                                                  k,
                                                  other_k,
                                                  -dESdES);
                    }
                }

                // substitute miRNA for the target
                GSList *other_occupant_list;
                for (other_occupant_list = target_site->occupants; other_occupant_list != NULL; other_occupant_list = other_occupant_list->next)
                {
                    MirbookingOccupant *other_occupant = other_occupant_list->data;

                    gsize other_k = _mirbooking_broker_get_occupant_index (self, other_occupant);

                    gdouble dEdES  = occupant->mirna == other_occupant->mirna ? -1 : 0;
                    gdouble dSdES  = -Stp / (self->priv->S[i] - _mirbooking_broker_get_target_site_occupants_quantity (self, target_site, ES));
                    gdouble dESdES = kf * (E[j] * dSdES + Stp * dEdES) - (kr + kcat) * (occupant == other_occupant ? 1 : 0);

                    if (mirna == other_occupant->mirna)
                    {
                        // it's already computed above
                        g_assert_cmpfloat (sparse_matrix_get_double (self->priv->J, k, other_k), ==, -dESdES);
                        continue;
                    }

                    sparse_matrix_set_double (self->priv->J,
                                              k,
                                              other_k,
                                              -dESdES);
                }
            }
        }
    }
}

/**
 * mirbooking_broker_evaluate:
 * @norm (out) (optional): The norm of the system.
 *
 * Evaluate the current state of the system.
 *
 * Returns: %TRUE if the evaluation is successful, otherwise %FALSE and @error
 * is set accordingly
 */
gboolean
mirbooking_broker_evaluate (MirbookingBroker          *self,
                            gdouble                   *norm,
                            GError                   **error)
{
    guint64 step_begin = g_get_monotonic_time ();

    if (g_once_init_enter (&self->priv->init))
    {
        g_return_val_if_fail (_mirbooking_broker_prepare_step (self),
                              FALSE);
        g_once_init_leave (&self->priv->init, 1);
    }

    _compute_F (self->priv->t,
                self->priv->y,
                self->priv->F,
                self);

    if (norm)
    {
        gdouble _norm = 0;
        gint i;
        #pragma omp parallel for reduction(+:_norm)
        for (i = 0; i < self->priv->y_len; i++)
        {
            _norm += pow (self->priv->F[i], 2);
        }
        *norm = sqrt (_norm);
    }

    return TRUE;
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
    guint64 step_begin = g_get_monotonic_time ();

    if (step_mode == MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE)
    {
        _compute_F (self->priv->t,
                    self->priv->y,
                    self->priv->F,
                    self);

        _compute_J (self->priv->t,
                    self->priv->y,
                    self->priv->J,
                    self);

        gboolean ret;
        ret = sparse_solver_solve (self->priv->solver,
                                   self->priv->J,
                                   self->priv->ES_delta,
                                   self->priv->dESdt);

        if (!ret)
        {
            g_set_error_literal (error,
                                 MIRBOOKING_ERROR,
                                 MIRBOOKING_ERROR_FAILED,
                                 "The solve step has failed.");
            return FALSE;
        }

        SparseSolverStatistics stats = sparse_solver_get_statistics (self->priv->solver);
        g_debug ("reorder-time: %fs factor-time: %fs solve-time: %fs flops: %f", stats.reorder_time, stats.factor_time, stats.solve_time, stats.flops);

        // apply the update
        gint i, j;
        #pragma omp parallel for collapse(2)
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

                    gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

                    #pragma omp atomic
                    self->priv->E[j]   -= step_size * self->priv->ES_delta[k];
                    self->priv->ES[k]  += step_size * self->priv->ES_delta[k];
                }
            }
        }
    }
    else if (step_mode == MIRBOOKING_BROKER_STEP_MODE_INTEGRATE)
    {
        gdouble t = self->priv->t;
        odeint_integrator_integrate (self->priv->integrator,
                                     _compute_F,
                                     self,
                                     self->priv->t + step_size);
        g_assert_cmpfloat (self->priv->t, ==, t + step_size);
    }
    else
    {
        g_assert_not_reached ();
    }

    g_debug ("Stepped in %lums", 1000 * (g_get_monotonic_time () - step_begin) / G_USEC_PER_SEC);

    return TRUE;
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
    return fmax (1 - pow (2, -0.0392 * k + 0.0054), 0);
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
        if (_mirbooking_broker_get_target_site_occupants_quantity (self, target_site, self->priv->ES))
        {
            gdouble proba = _mirbooking_broker_get_target_site_occupants_quantity (self, target_site, self->priv->ES) / target_quantity;
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
        g_assert_cmpfloat (pb_pmf (&pb, k), >=, 0);
        g_assert_cmpfloat (mean_silencing_by_number_of_sites (k), >=, 0);
        target_silencing += pb_pmf (&pb, k) * mean_silencing_by_number_of_sites (k);
    }

    pb_destroy (&pb);

    return target_silencing;
}
