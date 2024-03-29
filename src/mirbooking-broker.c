#include "mirbooking-broker.h"
#include "mirbooking-score-table-private.h"

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
#if HAVE_MKL_CBLAS
#include <mkl_cblas.h>
#elif HAVE_OPENBLAS
#include <openblas/cblas.h>
#else
#include <cblas.h>
#endif

#define RTOL 1e-6
#define ATOL 1e-8

#define COALESCE(x,d) (x == NULL ? (d) : (x))

typedef struct _MirbookingTargetPositions
{
    gsize     *positions;
    gsize      positions_len;
    GPtrArray *occupants;
} MirbookingTargetPositions;

static void
mirbooking_target_positions_clear (MirbookingTargetPositions *tss)
{
    g_free (tss->positions);
    g_ptr_array_unref (tss->occupants);
}

typedef struct
{
    gint                        rank;

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
    GArray     *target_positions; // (target, mirna) -> positions, scores and occupants in the target
    GArray     *occupants;

    /* transcription */
    gdouble    *targets_ktr;
    gdouble    *mirnas_ktr;

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
    MirbookingBrokerSparseSolver sparse_solver;
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
    PROP_RANK = 1,
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
        case PROP_RANK:
            self->priv->rank = g_value_get_int (value);
            break;
        case PROP_5PRIME_FOOTPRINT:
            self->priv->prime5_footprint = g_value_get_uint (value);
            break;
        case PROP_3PRIME_FOOTPRINT:
            self->priv->prime3_footprint = g_value_get_uint (value);
            break;
        case PROP_SCORE_TABLE:
            mirbooking_broker_set_score_table (self, g_value_get_object (value));
            break;
        case PROP_SPARSE_SOLVER:
            mirbooking_broker_set_sparse_solver (self, g_value_get_enum (value));
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static void
mirbooking_broker_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    MirbookingBroker *self = MIRBOOKING_BROKER (object);

    switch (property_id)
    {
        case PROP_RANK:
            g_value_set_int (value, self->priv->rank);
            break;
        case PROP_5PRIME_FOOTPRINT:
            g_value_set_uint (value, self->priv->prime5_footprint);
            break;
        case PROP_3PRIME_FOOTPRINT:
            g_value_set_uint (value, self->priv->prime3_footprint);
            break;
        case PROP_SCORE_TABLE:
            g_value_set_object (value, self->priv->score_table);
            break;
        default:
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
            break;
    }
}

static size_t
_mirbooking_broker_get_occupant_index (MirbookingBroker *self, const MirbookingOccupant *occupant)
{
    return occupant - &g_array_index (self->priv->occupants, MirbookingOccupant, 0);
}

gdouble
_mirbooking_broker_get_occupant_quantity (MirbookingBroker *self, const MirbookingOccupant *occupant, const gdouble *ES)
{
    gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

    g_assert_cmpint (k, >=, 0);
    g_assert_cmpint (k, <, self->priv->occupants->len);

    return ES[k];
}

gdouble
mirbooking_broker_get_occupant_quantity (MirbookingBroker *self, const MirbookingOccupant *occupant)
{
    g_return_val_if_fail (self->priv->init, 0.0);
    return _mirbooking_broker_get_occupant_quantity (self, occupant, self->priv->ES);
}

#if !GLIB_CHECK_VERSION(2,54,0)
static gboolean
g_ptr_array_find_with_equal_func (GPtrArray *array, gconstpointer elem, GEqualFunc equal_func, guint *index)
{
    guint i;
    for (i = 0; i < array->len; i++)
    {
        if (equal_func (g_ptr_array_index (array, i), elem))
        {
            if (index)
                *index = i;
            return TRUE;
        }
    }

    return FALSE;
}
#endif

/**
 * mirbooking_broker_set_occupant_quantity:
 * @self: A #MirbookingBroker
 * @occupant: A #MirbookingOccupant previously obtained from #mirbooking_broker_get_target_sites
 * @quantity: The new quantity of that occupant
 *
 * Set the concentration of an occupant to the given value.
 */
void
mirbooking_broker_set_occupant_quantity (MirbookingBroker *self, const MirbookingOccupant *occupant, gdouble quantity)
{
    g_return_if_fail (self->priv->init);
    g_return_if_fail (quantity >= 0);

    guint i;
    g_return_if_fail (g_ptr_array_find_with_equal_func (self->priv->mirnas, occupant->mirna, (GEqualFunc) mirbooking_sequence_equal, &i));

    gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

    self->priv->ES[k] = quantity;
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

    if (self->priv->target_sites)
    {
        g_array_unref (self->priv->target_sites);
        g_hash_table_unref (self->priv->target_sites_by_target);
        g_ptr_array_unref (self->priv->target_sites_by_target_index);
        g_array_unref (self->priv->target_positions);
        g_array_unref (self->priv->occupants);
    }

    if (self->priv->integrator)
    {
        g_free (self->priv->y);
        g_free (self->priv->F);
        g_free (self->priv->targets_ktr);
        g_free (self->priv->mirnas_ktr);
        odeint_integrator_free (self->priv->integrator);
    }

    if (self->priv->J)
    {
        sparse_matrix_clear (self->priv->J);
        g_free (self->priv->J);
    }

    if (self->priv->ES_delta)
        g_free (self->priv->ES_delta);

    sparse_solver_free (self->priv->solver);

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

    g_object_class_install_property (object_class, PROP_RANK,
                                     g_param_spec_int ("rank", "", "", 0, G_MAXINT, 0, G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
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

/**
 * mirbooking_broker_new_with_rank:
 *
 * Returns: (transfer full): A plain #Mirbooking instance with specified rank
 */
MirbookingBroker *
mirbooking_broker_new_with_rank (gint rank)
{
    return g_object_new (MIRBOOKING_BROKER_TYPE, "rank", rank, NULL);
}

/**
 * mirbooking_broker_get_rank:
 *
 * Obtain the rank of this broker in a distributed context.
 */
gint
mirbooking_broker_get_rank (MirbookingBroker *self)
{
    return self->priv->rank;
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
    if (self->priv->solver == NULL || self->priv->sparse_solver != sparse_solver)
    {
        if (self->priv->solver)
        {
            sparse_solver_free (self->priv->solver);
        }

        self->priv->sparse_solver = sparse_solver;
        self->priv->solver = sparse_solver_new ((SparseSolverMethod) sparse_solver);
        g_return_if_fail (self->priv->solver != NULL);

        g_object_notify (G_OBJECT (self), "sparse-solver");
    }
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
    g_return_if_fail (!self->priv->init);

    if (score_table != self->priv->score_table)
    {
        g_clear_object (&self->priv->score_table);
        self->priv->score_table = g_object_ref (score_table);
        g_object_notify (G_OBJECT (self), "score-table");
    }
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
 * mirbooking_broker_get_bound_mirna_quantity:
 *
 * Obtain the bound quantity of @mirna.
 */
gdouble
mirbooking_broker_get_bound_mirna_quantity (MirbookingBroker *self, MirbookingMirna *mirna)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    gdouble bound_quantity = 0;

    gsize k;
    #pragma omp parallel for reduction(+:bound_quantity)
    for (k = 0; k < self->priv->occupants->len; k++)
    {
        MirbookingOccupant *occupant = &g_array_index (self->priv->occupants, MirbookingOccupant, k);
        if (occupant->mirna == mirna)
        {
            bound_quantity += self->priv->ES[_mirbooking_broker_get_occupant_index (self, occupant)];
        }
    }

    return bound_quantity;
}

/**
 * mirbooking_broker_get_sequence_quantity:
 * @sequence: A #MirbookingSequence to retrieve quantity
 *
 * Retrieve the free quantity of @sequence.
 */
gdouble
mirbooking_broker_get_sequence_quantity (MirbookingBroker *self, MirbookingSequence *sequence)
{
    g_return_val_if_fail (g_hash_table_contains (self->priv->quantification, sequence), 0.0);

    if (self->priv->init)
    {
        if (MIRBOOKING_IS_MIRNA (sequence))
        {
            guint i;
            if (g_ptr_array_find_with_equal_func (self->priv->mirnas, sequence, (GEqualFunc) mirbooking_sequence_equal, &i))
            {
                return self->priv->E[i];
            }
        }
        else
        {
            guint i;
            if (g_ptr_array_find_with_equal_func (self->priv->targets, sequence, (GEqualFunc) mirbooking_sequence_equal, &i))
            {
                return self->priv->S[i];
            }
        }
    }
    else
    {
        return gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, sequence));
    }

    g_return_val_if_reached (0.0);
}

/**
 * mirbooking_broker_get_product_quantity:
 */
gdouble
mirbooking_broker_get_product_quantity (MirbookingBroker *self, MirbookingTarget *target)
{
    g_return_val_if_fail (self->priv->init, 0.0);
    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target, (GEqualFunc) mirbooking_sequence_equal, &i), 0.0);
    return self->priv->P[i];
}

/**
 * mirbooking_broker_set_product_quantity:
 */
void
mirbooking_broker_set_product_quantity (MirbookingBroker *self, MirbookingTarget *target, gdouble quantity)
{
    g_return_if_fail (self->priv->init);
    g_return_if_fail (quantity >= 0);
    guint i;
    g_return_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target, (GEqualFunc) mirbooking_sequence_equal, &i));
    self->priv->P[i] = quantity;
}

/**
 * mirbooking_set_sequence_quantity:
 * @sequence: A #MirbookingSequence being quantified for the
 * upcoming execution
 *
 * Set the free concentration of @sequence to @quantity. If @sequence is not
 * part of the system, it is added.
 *
 * The concentration of a #MirbookingMirna can be set to a negative value, as
 * long as the total amount present (see #mirbooking_broker_get_bound_mirna_quantity)
 * is positive.
 *
 * Note that no new sequence can be added this way once
 * #mirbooking_broker_evaluate has been called.
 */
void
mirbooking_broker_set_sequence_quantity (MirbookingBroker *self, MirbookingSequence *sequence, gdouble quantity)
{
    g_return_if_fail (MIRBOOKING_IS_MIRNA (sequence) || MIRBOOKING_IS_TARGET (sequence));
    g_return_if_fail (self->priv->init == 0 || g_hash_table_contains (self->priv->quantification, sequence));
    g_return_if_fail (MIRBOOKING_IS_MIRNA (sequence) || quantity >= 0);

    /* update the system */
    // TODO: have reverse-index for these use cases
    if (self->priv->init)
    {
        if (MIRBOOKING_IS_MIRNA (sequence))
        {
            g_return_if_fail (quantity + mirbooking_broker_get_bound_mirna_quantity (self, MIRBOOKING_MIRNA (sequence)) >= 0);
            guint i;
            if (g_ptr_array_find_with_equal_func (self->priv->mirnas, sequence, (GEqualFunc) mirbooking_sequence_equal, &i))
            {
                self->priv->E[i] = quantity;
            }
        }
        else
        {
            guint i;
            if (g_ptr_array_find_with_equal_func (self->priv->targets, sequence, (GEqualFunc) mirbooking_sequence_equal, &i))
            {
                self->priv->S[i] = quantity;
            }
        }
    }

    if (!g_hash_table_contains (self->priv->quantification, sequence))
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

    g_hash_table_insert (self->priv->quantification,
                         sequence,
                         gpointer_from_gfloat (quantity));
}

/**
 * mirbooking_broker_get_time:
 *
 * Get the time in seconds of the system.
 *
 * This is much more accurate to retrieve the time this way than to keep an
 * external counter because of some numerical errors that can accumulate when
 * stepping.
 */
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
    // find the lower target site
    *from_target_site = target_site - MIN (prime5_footprint, target_site->position);

    // find the upper target site
    *to_target_site = target_site + MIN (prime3_footprint,
                                         mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target_site->target)) - target_site->position - 1);

    g_assert ((*from_target_site)->target == target_site->target);
    g_assert ((*to_target_site)->target == target_site->target);
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
    gdouble vacancy = 1.0;

    const MirbookingTargetSite *from_target_site, *to_target_site;
    _mirbooking_broker_get_footprint_window (self,
                                             target_site,
                                             prime5_footprint,
                                             prime3_footprint,
                                             &from_target_site,
                                             &to_target_site);

    /*
     * Curiously, this corresponds to the following factorization of the joint
     * distribution of having all positions simultaneously unbound.
     *
     * Pr(p_1,p_2,...,p_n) = Pr(p_1) * Pr(p_2|p_1) * ... * Pr(p_n|p_1, p_2,... p_{n-1})
     *
     * Pr(p_1) is straightforward to calculate.
     *
     * For every other position, we compute the vacancy by assuming that every
     * preceding positions are unbounded as well.
     */

    const MirbookingTargetSite *last_tts = NULL;

    const MirbookingTargetSite *ts;
    for (ts = from_target_site; ts <= to_target_site; ts++)
    {
        /*
         * Here, we invert the footprint because we want to know where other
         * occupants can overlap this position.
         */
        const MirbookingTargetSite *fts, *tts;
        _mirbooking_broker_get_footprint_window (self,
                                                 ts,
                                                 prime3_footprint,
                                                 prime5_footprint,
                                                 &fts,
                                                 &tts);

        gdouble quantity = 0;
        const MirbookingTargetSite *nearby_target_site;
        for (nearby_target_site = fts; nearby_target_site <= tts; nearby_target_site++)
        {
            if (nearby_target_site > last_tts)
            {
                quantity += _mirbooking_broker_get_target_site_occupants_quantity (self, nearby_target_site, ES);
            }
        }

        vacancy *= (1 - (quantity / St));

        last_tts = tts;
    }

    return vacancy;
}

typedef struct _MirbookingScoredTargetSite
{
    MirbookingTargetSite *target_site;
    MirbookingScore       score;
} MirbookingScoredTargetSite;

static gint
cmp_gsize (gconstpointer a, gconstpointer b)
{
    return *(const gsize*)a - *(const gsize*)b;
}

static gboolean
_mirbooking_broker_prepare_step (MirbookingBroker *self, GError **error)
{
    guint64 prepare_begin = g_get_monotonic_time ();

    gsize target_sites_len = 0;

    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (self->priv->score_table != NULL, FALSE);
    g_return_val_if_fail (self->priv->solver != NULL, FALSE);

    guint i;
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
    self->priv->target_positions = g_array_sized_new (FALSE,
                                                         FALSE,
                                                         sizeof (MirbookingTargetPositions),
                                                         self->priv->targets->len * self->priv->mirnas->len);
    g_array_set_clear_func (self->priv->target_positions,
                            (GDestroyNotify) mirbooking_target_positions_clear);
    g_array_set_size (self->priv->target_positions,
                      self->priv->targets->len * self->priv->mirnas->len);

    // compute scores
    gsize occupants_len = 0;
    guint j;
    gboolean anyerror = FALSE;
    #pragma omp parallel for collapse(2) reduction(+:occupants_len)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTarget *target = g_ptr_array_index (self->priv->targets, i);
            MirbookingMirna *mirna   = g_ptr_array_index (self->priv->mirnas, j);
            MirbookingTargetPositions *seed_positions = &g_array_index (self->priv->target_positions, MirbookingTargetPositions, i * self->priv->mirnas->len + j);
            anyerror |= !mirbooking_score_table_compute_positions (MIRBOOKING_SCORE_TABLE (self->priv->score_table),
                                                                   mirna,
                                                                   target,
                                                                   &seed_positions->positions,
                                                                   &seed_positions->positions_len,
                                                                   error);

            occupants_len += seed_positions->positions_len;
        }
    }

    if (anyerror)
    {
        return FALSE;
    }

    g_debug ("Number of duplexes: %lu", occupants_len);

    // pre-allocate occupants in contiguous memory
    self->priv->occupants = g_array_sized_new (FALSE, FALSE, sizeof (MirbookingOccupant), occupants_len);
    g_array_set_size (self->priv->occupants, occupants_len);

    // intitialize occupants
    MirbookingOccupant *occupants = (MirbookingOccupant*) self->priv->occupants->data;
    gint k = 0;
    #pragma omp parallel for collapse(2) ordered
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingMirna *mirna = g_ptr_array_index (self->priv->mirnas, j);

            MirbookingTargetSite *target_sites = g_ptr_array_index (self->priv->target_sites_by_target_index,
                                                                    i);

            MirbookingTargetPositions *seed_positions = &g_array_index (self->priv->target_positions,
                                                                        MirbookingTargetPositions,
                                                                        i * self->priv->mirnas->len + j);

            gint _k;
            #pragma omp ordered
            {
                _k = k;
                k += seed_positions->positions_len;
            }

            seed_positions->occupants = g_ptr_array_sized_new (seed_positions->positions_len);

            guint p;
            for (p = 0; p < seed_positions->positions_len; p++)
            {
                MirbookingTargetSite *target_site = &target_sites[seed_positions->positions[p]];

                g_assert_cmpint (target_site->position, ==, seed_positions->positions[p]);

                MirbookingScore score = {0};
                anyerror |= !mirbooking_score_table_compute_score (self->priv->score_table,
                                                                   mirna,
                                                                   target_site->target,
                                                                   seed_positions->positions[p],
                                                                   &score,
                                                                   error);

                mirbooking_occupant_init (&occupants[_k + p],
                                          target_site->target,
                                          seed_positions->positions[p],
                                          mirna,
                                          score);

                /*
                 * Multiple microRNA might share this target site and prepended
                 * at once.
                 */
                #pragma omp critical
                target_site->occupants = g_slist_prepend (target_site->occupants, &occupants[_k + p]);
                g_ptr_array_add (seed_positions->occupants, &occupants[_k + p]);
            }
        }
    }

    if (anyerror)
    {
        return FALSE;
    }

    g_assert_cmpint (self->priv->occupants->len, ==, occupants_len);

    self->priv->targets_ktr  = g_new0 (gdouble, self->priv->targets->len);
    self->priv->mirnas_ktr  = g_new0 (gdouble, self->priv->mirnas->len);

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

    self->priv->ES_delta = g_new0 (gdouble, self->priv->occupants->len);

    // integrator
    self->priv->integrator = odeint_integrator_new (ODEINT_METHOD_DORMAND_PRINCE,
                                                    &self->priv->t,
                                                    self->priv->y,
                                                    self->priv->y_len,
                                                    ODEINT_INTEGRATOR_DEFAULT_RTOL,
                                                    ODEINT_INTEGRATOR_DEFAULT_ATOL);

    g_debug ("Prepared the first step in %lums", 1000 * (g_get_monotonic_time () - prepare_begin) / G_USEC_PER_SEC);

    return TRUE;
}

static gsize
absdiff (gsize a, gsize b)
{
    return MAX (a, b) - MIN (a, b);
}

static gdouble
_compute_kother (MirbookingBroker *self,
                 gsize             i,
                 gsize             position,
                 const gdouble    *ES,
                 gdouble           St)
{
    gdouble kother = 0;

    guint j;
    for (j = 0; j < self->priv->mirnas->len; j++)
    {
        // all the position of the other microrna
        MirbookingTargetPositions *tss = &g_array_index (self->priv->target_positions, MirbookingTargetPositions, self->priv->mirnas->len * i + j);

        guint q;
        for (q = 0; q < tss->positions_len; q++)
        {
            if (absdiff (tss->positions[q], position) > (self->priv->prime5_footprint + self->priv->prime3_footprint))
            {
                MirbookingOccupant *occupant_q = g_ptr_array_index (tss->occupants, q);
                kother += occupant_q->score.kcat * (_mirbooking_broker_get_occupant_quantity (self, occupant_q, ES) / St);
            }
        }
    }

    return kother;
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

    // basic transcription and degradation

    // dEdt = ktr - KDEGE * [E]
    cblas_dcopy (self->priv->mirnas->len,
                 self->priv->mirnas_ktr,
                 1,
                 dEdt,
                 1);
    cblas_daxpy (self->priv->mirnas->len,
                 -KDEGE,
                 E,
                 1,
                 dEdt,
                 1);

    // dSdt = ktr - KDEGS * [S]
    cblas_dcopy (self->priv->targets->len,
                 self->priv->targets_ktr,
                 1,
                 dSdt,
                 1);
    cblas_daxpy (self->priv->targets->len,
                 -KDEGS,
                 S,
                 1,
                 dSdt,
                 1);

    // dPdt = -KDEGP * P
    cblas_dcopy (self->priv->targets->len,
                 P,
                 1,
                 dPdt,
                 1);
    cblas_dscal (self->priv->targets->len,
                 -KDEGP,
                 dPdt,
                 1);

    guint i, j;
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
            MirbookingTargetPositions *seed_positions = &g_array_index (self->priv->target_positions,
                                                                        MirbookingTargetPositions,
                                                                        self->priv->mirnas->len * i + j);

            guint p;
            for (p = 0; p < seed_positions->positions_len; p++)
            {
                MirbookingOccupant *occupant = g_ptr_array_index (seed_positions->occupants, p);
                MirbookingTargetSite *target_site = &target_sites[seed_positions->positions[p]];

                gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

                g_assert (target_site->target == target);
                g_assert_cmpint (target_site->position, ==, seed_positions->positions[p]);
                g_assert (occupant->mirna == mirna);

                gdouble kf   = occupant->score.kf;
                gdouble kr   = occupant->score.kr;
                gdouble kcat = occupant->score.kcat;

                gdouble Stp = S[i] * _mirbooking_broker_get_target_site_vacancy (self,
                                                                                 target_site,
                                                                                 prime5_footprint,
                                                                                 prime3_footprint,
                                                                                 S[i],
                                                                                 ES);

                gdouble kother = _compute_kother (self, i, seed_positions->positions[p], ES, S[i]);

                #pragma omp atomic
                dEdt[j] += -kf * E[j] * Stp + kr * ES[k] + kcat * ES[k] + kother * ES[k];

                #pragma omp atomic
                dSdt[i] += - kcat * ES[k];

                dESdt[k] = kf * E[j] * Stp - kr * ES[k] - kcat * ES[k] - kother * ES[k];

                #pragma omp atomic
                dPdt[i] += kcat * ES[k];
            }
        }
    }
}

static void
_prepare_J (MirbookingBroker *self)
{
    self->priv->J = g_new0 (SparseMatrix, 1);

    // count nnz entries in the Jacobian
    gsize nnz = 0;
    guint i, j;
    #pragma omp parallel for collapse(2) reduction(+:nnz)
    for (i = 0; i < self->priv->targets->len; i++)
    {
        for (j = 0; j < self->priv->mirnas->len; j++)
        {
            MirbookingTargetSite *target_sites = g_ptr_array_index (self->priv->target_sites_by_target_index,
                                                                    i);
            MirbookingTargetPositions *seed_scores = &g_array_index (self->priv->target_positions,
                                                                       MirbookingTargetPositions,
                                                                       i * self->priv->mirnas->len + j);

            guint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                // substitute targets
                guint z;
                for (z = 0; z < self->priv->targets->len; z++)
                {
                    MirbookingTargetPositions *alternative_seed_scores = &g_array_index (self->priv->target_positions,
                                                                                           MirbookingTargetPositions,
                                                                                           z * self->priv->mirnas->len + j);
                    nnz += alternative_seed_scores->positions_len;
                }

                // substitute miRNAs (excluding this one as we consider it as a substitute target)
                nnz += g_slist_length (target_sites[seed_scores->positions[p]].occupants) - 1;
            }
        }
    }

    g_debug ("nnz: %lu, sparsity: %.2f%%", nnz,  100.0 * (1.0 - (gdouble) nnz / pow (self->priv->occupants->len, 2)));

    size_t shape[2] = {self->priv->occupants->len, self->priv->occupants->len};
    sparse_matrix_init (self->priv->J,
                        self->priv->sparse_solver == MIRBOOKING_BROKER_SPARSE_SOLVER_LAPACK ? SPARSE_MATRIX_STORAGE_DENSE : SPARSE_MATRIX_STORAGE_CSR,
                        SPARSE_MATRIX_TYPE_DOUBLE,
                        shape,
                        nnz);

    self->priv->J->hints |= SPARSE_MATRIX_HINT_SYMMETRIC_STRUCTURE;
    self->priv->J->hints |= SPARSE_MATRIX_HINT_POSITIVE_DEFINITE;

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
            MirbookingTargetPositions *seed_scores = &g_array_index (self->priv->target_positions,
                                                                       MirbookingTargetPositions,
                                                                       i * self->priv->mirnas->len + j);

            guint p;
            for (p = 0; p < seed_scores->positions_len; p++)
            {
                // footprint interactions
                MirbookingTargetSite *target_site = &target_sites[seed_scores->positions[p]];
                MirbookingOccupant *occupant = g_ptr_array_index (seed_scores->occupants, p);

                gsize colind[self->priv->J->shape[0]];
                gsize row_nnz = 0;

                // substitute target
                guint z;
                for (z = 0; z < self->priv->targets->len; z++)
                {
                    MirbookingTargetPositions *alternative_seed_scores = &g_array_index (self->priv->target_positions,
                                                                                           MirbookingTargetPositions,
                                                                                           z * self->priv->mirnas->len + j);

                    guint w;
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

    gsize prime5_footprint = self->priv->prime5_footprint;
    gsize prime3_footprint = self->priv->prime3_footprint;

    guint i, j;
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
            MirbookingTargetPositions *seed_positions = &g_array_index (self->priv->target_positions,
                                                                        MirbookingTargetPositions,
                                                                        self->priv->mirnas->len * i + j);

            guint p;
            for (p = 0; p < seed_positions->positions_len; p++)
            {
                MirbookingOccupant *occupant = g_ptr_array_index (seed_positions->occupants, p);
                MirbookingTargetSite *target_site = &target_sites[seed_positions->positions[p]];

                gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

                g_assert (target_site->target == target);
                g_assert_cmpint (target_site->position, ==, seed_positions->positions[p]);
                g_assert (occupant->mirna == mirna);

                gdouble kf   = occupant->score.kf;
                gdouble kr   = occupant->score.kr;
                gdouble kcat = occupant->score.kcat;

                gdouble Stp = S[i] * _mirbooking_broker_get_target_site_vacancy (self,
                                                                                 target_site,
                                                                                 prime5_footprint,
                                                                                 prime3_footprint,
                                                                                 S[i],
                                                                                 ES);

                // substitute target for the microRNA
                guint z;
                for (z = 0; z < self->priv->targets->len; z++)
                {
                    MirbookingTargetPositions *alternative_seed_positions = &g_array_index (self->priv->target_positions,
                                                                                            MirbookingTargetPositions,
                                                                                            z * self->priv->mirnas->len + j);
                    guint w;
                    for (w = 0; w < alternative_seed_positions->occupants->len; w++)
                    {
                        MirbookingOccupant *other_occupant = g_ptr_array_index (alternative_seed_positions->occupants, w);

                        g_assert (other_occupant->mirna == g_ptr_array_index (self->priv->mirnas, j));

                        gsize other_k = _mirbooking_broker_get_occupant_index (self, other_occupant);


                        gdouble kother = 0;
                        if (occupant == other_occupant)
                        {
                            kother = _compute_kother (self, i, seed_positions->positions[p], ES, S[i]);
                        }
                        else if (i == z && absdiff (seed_positions->positions[p], alternative_seed_positions->positions[w]) > (self->priv->prime5_footprint + self->priv->prime3_footprint))
                        {
                            /*
                             * Here, we account for the kother if a microRNA is
                             * shared for the pair of complexes because it's
                             * essentially free and speeds up the convergence.
                             *
                             * Ideally we would do if for all pair of
                             * complexes, but it has a combinatorial cost.
                             */
                            kother = occupant->score.kcat * (_mirbooking_broker_get_occupant_quantity (self, occupant, ES) / self->priv->S[i]);
                        }

                        gdouble dEdES  = -1; // always
                        gdouble dSdES  = (z == i && seed_positions->positions[p] == alternative_seed_positions->positions[w]) ? -1 : 0;
                        gdouble dESdES = kf * (E[j] * dSdES + Stp * dEdES) - (kr + kcat) * (occupant == other_occupant ? 1 : 0) - kother;

                        sparse_matrix_set_double (J,
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

                    gdouble kother = occupant == other_occupant ? _compute_kother (self, i, seed_positions->positions[p], ES, S[i]) : 0;

                    gdouble dEdES  = occupant->mirna == other_occupant->mirna ? -1 : 0;
                    gdouble dSdES  = -1;
                    gdouble dESdES = kf * (E[j] * dSdES + Stp * dEdES) - (kr + kcat) * (occupant == other_occupant ? 1 : 0) - kother;

                    sparse_matrix_set_double (J,
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
 * @self: A #MirbookingBroker
 * @error_ratio: (out) (optional): Error-to-tolerance ratio
 *
 * Evaluate the current state of the system.
 *
 * Returns: %TRUE if the evaluation is successful, otherwise %FALSE and @error
 * is set accordingly
 */
gboolean
mirbooking_broker_evaluate (MirbookingBroker  *self,
                            gdouble           *error_ratio,
                            GError           **error)
{
    if (g_once_init_enter (&self->priv->init))
    {
        g_return_val_if_fail (_mirbooking_broker_prepare_step (self, error),
                              FALSE);
        g_once_init_leave (&self->priv->init, 1);
    }

    _compute_F (self->priv->t,
                self->priv->y,
                self->priv->F,
                self);

    if (error_ratio)
    {
        gdouble _error_ratio = 0;
        gsize i;
        for (i = 0; i < self->priv->y_len; i++)
        {
            _error_ratio = fmax (_error_ratio, fabs (self->priv->F[i]) / (RTOL  * fabs (self->priv->y[i]) + ATOL));
        }

        *error_ratio = _error_ratio;
    }

    return TRUE;
}

/**
 * mirbooking_broker_step:
 * @step_mode: A #MirbookingBrokerStepMode
 * @step_size: A step size for integration or a fraction of the step for the
 * Newton-Raphson update
 *
 * Perform a step based on the current state of the system.
 *
 * Returns: %TRUE on success, otherwise @error is set
 */
gboolean
mirbooking_broker_step (MirbookingBroker         *self,
                        MirbookingBrokerStepMode  step_mode,
                        gdouble                   step_size,
                        GError                   **error)
{
    if (g_once_init_enter (&self->priv->init))
    {
        g_return_val_if_fail (_mirbooking_broker_prepare_step (self, error),
                              FALSE);
        g_once_init_leave (&self->priv->init, 1);
    }

    if (step_mode == MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE)
    {
        if (self->priv->rank == 0)
        {
            _compute_F (self->priv->t,
                        self->priv->y,
                        self->priv->F,
                        self);

            if (self->priv->J == NULL)
            {
                _prepare_J (self);
            }

            _compute_J (self->priv->t,
                        self->priv->y,
                        self->priv->J,
                        self);
        }

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
        guint i, j;
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
                MirbookingTargetPositions *seed_positions = &g_array_index (self->priv->target_positions,
                                                                            MirbookingTargetPositions,
                                                                            self->priv->mirnas->len * i + j);

                guint p;
                for (p = 0; p < seed_positions->positions_len; p++)
                {
                    MirbookingOccupant *occupant = g_ptr_array_index (seed_positions->occupants, p);

                    g_assert (occupant->mirna == mirna);

                    gsize k = _mirbooking_broker_get_occupant_index (self, occupant);

                    #pragma omp atomic
                    self->priv->E[j]   -= step_size * self->priv->ES_delta[k];
                    self->priv->ES[k]  += step_size * self->priv->ES_delta[k];
                }
            }
        }

        // Under the steady-state assumption, all substrate degradation is
        // compensated by transcription of new targets.
        // The system must be however reevaluated as we have applied an update.
        _compute_F (self->priv->t,
                    self->priv->y,
                    self->priv->F,
                    self);

        {
            // ktr = ktr - dEdt
            cblas_daxpy (self->priv->mirnas->len,
                         -1,
                         self->priv->dEdt,
                         1,
                         self->priv->mirnas_ktr,
                         1);

            // ktr = ktr - dSdt
            cblas_daxpy (self->priv->targets->len,
                         -1,
                         self->priv->dSdt,
                         1,
                         self->priv->targets_ktr,
                         1);

            // P = -(dSdt - ktr) / KDEGP
            // P = dSdt; P = -ktr + P; P = -1/KDEGP * P
            cblas_dcopy (self->priv->targets->len,
                         self->priv->dSdt,
                         1,
                         self->priv->P,
                         1);
            cblas_daxpy (self->priv->targets->len,
                         -1,
                         self->priv->targets_ktr,
                         1,
                         self->priv->P,
                         1);
            cblas_dscal (self->priv->targets->len,
                         -1.0 / KDEGP,
                         self->priv->P,
                         1);
        }
    }
    else if (step_mode == MIRBOOKING_BROKER_STEP_MODE_INTEGRATE)
    {
        odeint_integrator_integrate (self->priv->integrator,
                                     _compute_F,
                                     self,
                                     self->priv->t + step_size);
    }
    else
    {
        g_assert_not_reached ();
    }

    return TRUE;
}
/**
 * mirbooking_broker_get_mirna_transcription_rate:
 *
 * Get the rate of transcription of the given #MirbookingTarget.
 *
 * This is resolved when stepping with @MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE
 * using the steady-state assumption.
 */
gdouble
mirbooking_broker_get_mirna_transcription_rate (MirbookingBroker *self,
                                                MirbookingMirna  *mirna)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->mirnas, mirna, (GEqualFunc) mirbooking_sequence_equal, &i), 0);

    return self->priv->mirnas_ktr[i];
}

/**
 * mirbooking_broker_set_mirna_transcription_rate:
 *
 * Set the rate of transcription of @mirna to @transcription_rate.
 */
void
mirbooking_broker_set_mirna_transcription_rate (MirbookingBroker *self, MirbookingMirna *mirna, gdouble transcription_rate)
{
    g_return_if_fail (self->priv->init);
    g_return_if_fail (transcription_rate >= 0);

    guint i;
    g_return_if_fail (g_ptr_array_find_with_equal_func (self->priv->mirnas, mirna, (GEqualFunc) mirbooking_sequence_equal, &i));

    self->priv->mirnas_ktr[i] = transcription_rate;
}

/**
 * mirbooking_broker_get_mirna_degradation_rate:
 *
 * Obtain @mirna degradation rate.
 */
gdouble
mirbooking_broker_get_mirna_degradation_rate (MirbookingBroker *self,
                                              MirbookingMirna  *mirna)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->mirnas, mirna, (GEqualFunc) mirbooking_sequence_equal, &i), 0.0);

    return KDEGE * self->priv->E[i];
}

/**
 * mirbooking_broker_get_target_transcription_rate:
 *
 * Get the rate of transcription of the given #MirbookingTarget.
 *
 * This is resolved when stepping with @MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE
 * using the steady-state assumption.
 */
gdouble
mirbooking_broker_get_target_transcription_rate (MirbookingBroker *self,
                                                 MirbookingTarget *target)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target, (GEqualFunc) mirbooking_sequence_equal, &i), 0);

    return self->priv->targets_ktr[i];
}

/**
 * mirbooking_broker_get_target_degradation_rate:
 *
 * Get the rate of degradation of the given #MirbookingTarget.
 */
gdouble
mirbooking_broker_get_target_degradation_rate (MirbookingBroker *self,
                                               MirbookingTarget *target)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target, (GEqualFunc) mirbooking_sequence_equal, &i), 0);

    return KDEGS * self->priv->S[i];
}

/**
 * mirbooking_broker_set_target_transcription_rate:
 *
 * Set the rate of transcription of @target to @transcription_rate.
 */
void
mirbooking_broker_set_target_transcription_rate (MirbookingBroker *self, MirbookingTarget *target, gdouble transcription_rate)
{
    g_return_if_fail (self->priv->init);
    g_return_if_fail (transcription_rate >= 0);

    guint i;
    g_return_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target, (GEqualFunc) mirbooking_sequence_equal, &i));

    self->priv->targets_ktr[i] = transcription_rate;
}

/**
 * mirbooking_broker_get_product_degradation_rate:
 *
 * Get the rate of product degradation.
 *
 * This is resolved when stepping with @MIRBOOKING_BROKER_STEP_MODE_SOLVE_STEADY_STATE
 * using the steady-state assumption.
 */
gdouble
mirbooking_broker_get_product_degradation_rate (MirbookingBroker *self,
                                                MirbookingTarget *target)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target, (GEqualFunc) mirbooking_sequence_equal, &i), 0);

    return KDEGP * self->priv->P[i];
}

/**
 * mirbooking_broker_get_target_sites:
 *
 * Obtain the computed #MirbookingTargetSite array by this #MirbookingBroker.
 *
 * Returns: (element-type MirbookingTargetSite) (transfer none): A view of the
 * computed #MirbookingTargetSite
 */
const GArray *
mirbooking_broker_get_target_sites (MirbookingBroker *self)
{
    g_return_val_if_fail (self != NULL, NULL);
    g_return_val_if_fail (self->priv->init, NULL);

    return self->priv->target_sites;
}

/**
 * mirbooking_broker_get_target_site_quantity:
 *
 * Obtain the expected concentration of a #MirbookingTargetSite which account
 * for overlapping miRISC in the neighborhood.
 */
gdouble
mirbooking_broker_get_target_site_quantity (MirbookingBroker *self, const MirbookingTargetSite *target_site)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    gdouble St = mirbooking_broker_get_sequence_quantity (self, MIRBOOKING_SEQUENCE (target_site->target));

    return _mirbooking_broker_get_target_site_vacancy (self,
                                                       target_site,
                                                       self->priv->prime5_footprint,
                                                       self->priv->prime3_footprint,
                                                       St,
                                                       self->priv->ES) * St;
}

/**
 * mirbooking_broker_get_target_site_occupants_quantity:
 *
 * Obtain the total concentration of occupants occupying this target site
 * either directly or indirectly via the footprint.
 */
gdouble
mirbooking_broker_get_target_site_occupants_quantity (MirbookingBroker            *self,
                                                      const MirbookingTargetSite *target_site)
{
    gdouble total_occupants_quantity = 0;

    const MirbookingTargetSite *fts, *tts;
    _mirbooking_broker_get_footprint_window (self,
                                             target_site,
                                             self->priv->prime3_footprint,
                                             self->priv->prime5_footprint,
                                             &fts, &tts);

    const MirbookingTargetSite *ts;
    for (ts = fts; ts <= tts; ts++)
    {
        total_occupants_quantity += _mirbooking_broker_get_target_site_occupants_quantity (self, ts, self->priv->ES);
    }

    return total_occupants_quantity;
}

/**
 * mirbooking_broker_get_target_site_catalytic_rate:
 *
 * Obtain the target site catalytic rate.
 */
gdouble
mirbooking_broker_get_target_site_kother (MirbookingBroker *self,
                                          const MirbookingTargetSite *target_site)
{
    g_return_val_if_fail (self->priv->init, 0.0);

    guint i;
    g_return_val_if_fail (g_ptr_array_find_with_equal_func (self->priv->targets, target_site->target, (GEqualFunc) mirbooking_sequence_equal, &i), 0.0);

    return _compute_kother (self, i, target_site->position, self->priv->ES, self->priv->S[i]);
}

/**
 * mirbooking_broker_get_mirnas:
 *
 * Returns: (element-type MirbookingMirna) (transfer none):
 */
const GPtrArray *
mirbooking_broker_get_mirnas (MirbookingBroker *self)
{
    return self->priv->mirnas;
}

/**
 * mirbooking_broker_get_targets:
 *
 * Returns: (element-type MirbookingTarget) (transfer none):
 */
const GPtrArray *
mirbooking_broker_get_targets (MirbookingBroker *self)
{
    return self->priv->targets;
}

/**
 * mirbooking_broker_get_occupants:
 *
 * This is much faster to manipulate if the intent is to traverse all the
 * complexes regardless of their actual location.
 *
 * Returns: (element-type MirbookingOccupant) (transfer none): A view over the
 * occupants
 */
const GArray *
mirbooking_broker_get_occupants (MirbookingBroker *self)
{
    g_return_val_if_fail (self->priv->init, NULL);
    return self->priv->occupants;
}

/**
 * mirbooking_broker_get_target_occupants_pmf:
 * @target: The #MirbookingTarget for which we are retrieving the silencing
 * @pmf_len: (out): Length of the PMF which correspond to the number
 * of occupied sites plus one
 *
 * Compute the probability mass function of the number of occupied target sites
 * on a given target by modeling them with a Poisson-Binomial distribution.
 *
 * Returns: (array length=pmf_len): The probability mass function of the number
 * of bound miRISC complexes or %NULL if it cannot be computed
 */
gdouble *
mirbooking_broker_get_target_occupants_pmf (MirbookingBroker *self, MirbookingTarget *target, gsize *pmf_len)
{
    g_return_val_if_fail (self->priv->init, NULL);
    g_return_val_if_fail (g_hash_table_contains (self->priv->quantification, target), NULL);

    gfloat target_quantity = mirbooking_broker_get_sequence_quantity (self, MIRBOOKING_SEQUENCE (target));

    g_autoptr (GArray) probability_by_position = g_array_new (FALSE, FALSE, sizeof (gdouble));

    MirbookingTargetSite *target_site = g_hash_table_lookup (self->priv->target_sites_by_target, target);
    while (target_site < &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len) &&
           target_site->target == target)
    {
        gdouble occupants_quantity = _mirbooking_broker_get_target_site_occupants_quantity (self, target_site, self->priv->ES);
        if (occupants_quantity > 0)
        {
            gdouble p = occupants_quantity / target_quantity;
            g_array_append_val (probability_by_position, p);
        }

        ++target_site;
    }

    PoissonBinomial pb;
    pb_init (&pb,
             (gdouble*) probability_by_position->data,
             probability_by_position->len);

    gdouble* pmf = g_new (gdouble, 1 + probability_by_position->len);
    memcpy (pmf, pb.pmf, (1 + probability_by_position->len) * sizeof (gdouble));

    if (pmf_len)
        *pmf_len = 1 + probability_by_position->len;

    pb_destroy (&pb);

    return pmf;
}

/**
 * mirbooking_broker_get_target_expressed_fraction:
 *
 * Obtain the fraction of substrate @target that is expressed.
 */
gdouble
mirbooking_broker_get_target_expressed_fraction (MirbookingBroker *broker, MirbookingTarget *target)
{
    gdouble ef = 0;
    gdouble lambda = .27;

    gsize pmf_len;
    g_autofree gdouble *pmf = mirbooking_broker_get_target_occupants_pmf (broker, target, &pmf_len);

    gsize i;
    for (i = 0; i < pmf_len; i++)
    {
        ef += exp (-lambda * i) * pmf[i];
    }

    return ef;
}

static gboolean
write_output_to_tsv (MirbookingBroker  *mirbooking,
                     GOutputStream     *out,
                     GError           **error)
{
    g_autoptr (GDataOutputStream) output_f = g_data_output_stream_new (out);

    gchar *header = "gene_accession\t"
                    "gene_name\t"
                    "target_accession\t"
                    "target_name\t"
                    "target_quantity\t"
                    "position\t"
                    "mirna_accession\t"
                    "mirna_name\t"
                    "mirna_quantity\t"
                    "score\t"
                    "quantity\n";

    if (!g_data_output_stream_put_string (output_f, header, NULL, error))
    {
        return FALSE;
    }

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    gfloat target_quantity = 0;

    const MirbookingTargetSite *target_site;
    MirbookingTarget *cur_target = NULL;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        // recompute each time the target changes
        if (cur_target != target_site->target)
        {
            cur_target = target_site->target;
            target_quantity = mirbooking_broker_get_sequence_quantity (mirbooking,
                                                                       MIRBOOKING_SEQUENCE (target_site->target));
        }

        // report individual occupants
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            MirbookingOccupant *occupant = occupants->data;
            g_autofree gchar* line = g_strdup_printf ("%s\t%s\t%s\t%s\t%e\t%lu\t%s\t%s\t%e\t%e\t%e\n",
                                                      COALESCE (mirbooking_sequence_get_gene_accession (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                                                      COALESCE (mirbooking_sequence_get_gene_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                                                      COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                                                      target_quantity,
                                                      target_site->position + 1, // 1-based
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                                                      COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)), "N/A"),
                                                      mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (occupant->mirna)) + mirbooking_broker_get_bound_mirna_quantity (mirbooking, occupant->mirna),
                                                      MIRBOOKING_SCORE_KM (occupant->score) + (mirbooking_broker_get_target_site_kother (mirbooking, target_site) / occupant->score.kf),
                                                      mirbooking_broker_get_occupant_quantity (mirbooking, occupant));

            if (!g_data_output_stream_put_string (output_f, line, NULL, error))
            {
                return FALSE;
            }
        }
    }

    return TRUE;
}

static gboolean
write_output_to_tsv_detailed (MirbookingBroker  *mirbooking,
                              GOutputStream     *out,
                              GError           **error)
{
    g_autoptr (GDataOutputStream) output_f = g_data_output_stream_new (out);

    gchar *header = "gene_accession\t"
                    "gene_name\t"
                    "target_accession\t"
                    "target_name\t"
                    "target_quantity\t"
                    "position\t"
                    "mirna_accession\t"
                    "mirna_name\t"
                    "mirna_quantity\t"
                    "kf\t"
                    "kr\t"
                    "kcleave\t"
                    "krelease\t"
                    "kcat\t"
                    "kother\t"
                    "kd\t"
                    "km\t"
                    "quantity\n";

    if (!g_data_output_stream_put_string (output_f,
                                          header,
                                          NULL,
                                          error))
    {
        return FALSE;
    }

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    gfloat target_quantity = 0;

    const MirbookingTargetSite *target_site;
    MirbookingTarget *cur_target = NULL;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        // recompute each time the target changes
        if (cur_target != target_site->target)
        {
            cur_target = target_site->target;
            target_quantity = mirbooking_broker_get_sequence_quantity (mirbooking,
                                                                       MIRBOOKING_SEQUENCE (target_site->target));
        }

        // report individual occupants
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            MirbookingOccupant *occupant = occupants->data;
            g_autofree gchar *line = g_strdup_printf ("%s\t%s\t%s\t%s\t%e\t%lu\t%s\t%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                                                      COALESCE (mirbooking_sequence_get_gene_accession (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                                                      COALESCE (mirbooking_sequence_get_gene_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                                                      COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (target_site->target)), "N/A"),
                                                      target_quantity,
                                                      target_site->position + 1, // 1-based
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)),
                                                      COALESCE (mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)), "N/A"),
                                                      mirbooking_broker_get_sequence_quantity (mirbooking, MIRBOOKING_SEQUENCE (occupant->mirna)) + mirbooking_broker_get_bound_mirna_quantity (mirbooking, occupant->mirna),
                                                      occupant->score.kf,
                                                      occupant->score.kr,
                                                      occupant->score.kcleave,
                                                      occupant->score.krelease,
                                                      occupant->score.kcat,
                                                      mirbooking_broker_get_target_site_kother (mirbooking, target_site),
                                                      MIRBOOKING_SCORE_KD (occupant->score),
                                                      MIRBOOKING_SCORE_KM (occupant->score) + (mirbooking_broker_get_target_site_kother (mirbooking, target_site) / occupant->score.kf),
                                                      mirbooking_broker_get_occupant_quantity (mirbooking, occupant));
            if (!g_data_output_stream_put_string (output_f, line, NULL, error))
            {
                return FALSE;
            }
        }
    }

    return TRUE;
}

static gboolean
write_output_to_gff3 (MirbookingBroker *mirbooking, GOutputStream *out, GError **error)
{
    g_autoptr (GDataOutputStream) output_f = g_data_output_stream_new (out);

    if (!g_data_output_stream_put_string (output_f, "##gff-version 3\n", NULL, error))
    {
        return FALSE;
    }

    const GArray *target_sites = mirbooking_broker_get_target_sites (mirbooking);

    gint i = 1;

    const MirbookingTargetSite *target_site;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        // report individual occupants
        GSList *occupants;
        for (occupants = target_site->occupants; occupants != NULL; occupants = occupants->next)
        {
            // the sequence ontology for 'miRNA_target_site' is 'SO:0000934'
            MirbookingOccupant *occupant = occupants->data;
            g_autofree gchar *line = g_strdup_printf ("%s\tmiRBooking\tmiRNA_target_site\t%lu\t%lu\t%e\t.\t.\tID=%d;Name=%s;Alias=%s\n",
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)),
                                                      (gsize) MAX (1, (gssize) target_site->position + 1 - (gssize) mirbooking->priv->prime5_footprint),
                                                      MIN (mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target_site->target)), target_site->position + 1 + mirbooking->priv->prime3_footprint),
                                                      mirbooking_broker_get_occupant_quantity (mirbooking, occupant),
                                                      i++,
                                                      mirbooking_sequence_get_name (MIRBOOKING_SEQUENCE (occupant->mirna)),
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (occupant->mirna)));
            if (!g_data_output_stream_put_string (output_f,
                                                  line,
                                                  NULL,
                                                  error))
            {
                return FALSE;
            }
        }
    }

    return TRUE;
}

static gboolean
write_output_to_wiggle (MirbookingBroker *broker, GOutputStream *out, GError **error)
{
    g_autoptr (GDataOutputStream) output_f = g_data_output_stream_new (out);

    const GArray *target_sites = mirbooking_broker_get_target_sites (broker);

    if (!g_data_output_stream_put_string (output_f, "track type=wiggle_0\n", NULL, error))
    {
        return FALSE;
    }

    MirbookingTarget *target = NULL;
    const MirbookingTargetSite *target_site;
    for (target_site = &g_array_index (target_sites, MirbookingTargetSite, 0);
         target_site < &g_array_index (target_sites, MirbookingTargetSite, target_sites->len);
         target_site++)
    {
        if (target != target_site->target)
        {
            g_autofree gchar *line = g_strdup_printf ("variableStep chrom=%s\n",
                                                      mirbooking_sequence_get_accession (MIRBOOKING_SEQUENCE (target_site->target)));
            if (!g_data_output_stream_put_string (output_f,
                                                  line,
                                                  NULL,
                                                  error))
            {
                return FALSE;
            }
            target = target_site->target;
        }

        gdouble St  = mirbooking_broker_get_sequence_quantity (broker, MIRBOOKING_SEQUENCE (target_site->target));
        gdouble Stp = mirbooking_broker_get_target_site_occupants_quantity (broker, target_site);

        // only report positions with activity
        if (Stp > 0)
        {
            g_autofree gchar *line = g_strdup_printf ("%lu %f\n",
                                                      target_site->position + 1,
                                                      Stp / St);
            if (!g_data_output_stream_put_string (output_f,
                                                  line,
                                                  NULL,
                                                  error))
            {
                return FALSE;
            }
        }
    }

    return TRUE;
}

typedef struct _MirbookingBrokerOutputFormatMeta
{
    MirbookingBrokerOutputFormat  output_format;
    const gchar            *nick;
    gboolean (*write) (MirbookingBroker *broker, GOutputStream *output, GError **error);
} MirbookingBrokerOutputFormatMeta;

static const MirbookingBrokerOutputFormatMeta OUTPUT_FORMAT_META[] =
{
    {MIRBOOKING_BROKER_OUTPUT_FORMAT_TSV,          "tsv",          write_output_to_tsv},
    {MIRBOOKING_BROKER_OUTPUT_FORMAT_TSV_DETAILED, "tsv-detailed", write_output_to_tsv_detailed},
    {MIRBOOKING_BROKER_OUTPUT_FORMAT_GFF3,         "gff3",         write_output_to_gff3},
    {MIRBOOKING_BROKER_OUTPUT_FORMAT_WIG,          "wig",          write_output_to_wiggle}
};

gboolean
mirbooking_broker_write_output_to_stream (MirbookingBroker              *self,
                                          GOutputStream                 *out,
                                          MirbookingBrokerOutputFormat   output_format,
                                          GError                       **error)
{
    g_return_val_if_fail (self->priv->init, FALSE);

    gsize i;
    for (i = 0; i < sizeof (OUTPUT_FORMAT_META); i++)
    {
        if (OUTPUT_FORMAT_META[i].output_format == output_format)
        {
            return OUTPUT_FORMAT_META[i].write (self,
                                                out,
                                                error);
        }
    }

    g_return_val_if_reached (FALSE);
}

gboolean
mirbooking_broker_write_output_to_file (MirbookingBroker              *self,
                                        GFile                         *output_file,
                                        MirbookingBrokerOutputFormat   output_format,
                                        GError                       **error)
{
    g_autoptr (GOutputStream) out = G_OUTPUT_STREAM (g_file_replace (output_file,
                                                                     NULL,
                                                                     FALSE,
                                                                     G_FILE_CREATE_NONE,
                                                                     NULL,
                                                                     error));

    g_return_val_if_fail (out != NULL, FALSE);

    return mirbooking_broker_write_output_to_stream (self,
                                                     out,
                                                     output_format,
                                                     error);
}
