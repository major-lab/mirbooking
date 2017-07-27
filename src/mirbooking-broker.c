#include "mirbooking-broker.h"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>

typedef struct
{
    gfloat                threshold;
    gfloat                log_base;
    gsize                 prime5_footprint;
    gsize                 prime3_footprint;
    MirbookingScoreTable *score_table;

    GSList               *targets;
    GSList               *mirnas;

    GHashTable           *quantification; // #MirbookingSequence -> #gfloat (initial quantity)

    /* all the target sites, stored contiguously */
    GArray               *target_sites;
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
    self->priv->quantification = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                   (GEqualFunc) mirbooking_sequence_equal);
}

enum
{
    PROP_THRESHOLD = 1,
    PROP_LOG_BASE,
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
        case PROP_THRESHOLD:
            self->priv->threshold = g_value_get_float (value);
            break;
        case PROP_LOG_BASE:
            self->priv->log_base = g_value_get_float (value);
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
        case PROP_THRESHOLD:
            g_value_set_float (value, self->priv->threshold);
            break;
        case PROP_LOG_BASE:
            g_value_set_float (value, self->priv->log_base);
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
mirbooking_occupant_new (MirbookingMirna *mirna, guint quantity)
{
    MirbookingOccupant *ret = g_slice_new (MirbookingOccupant);
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

    g_slist_free_full (self->priv->targets, g_object_unref);
    g_slist_free_full (self->priv->mirnas, g_object_unref);

    if (self->priv->target_sites)
    {
        g_array_free (self->priv->target_sites, TRUE);
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

    g_object_class_install_property (object_class, PROP_THRESHOLD,
                                     g_param_spec_float ("threshold", "", "", 0, 1, MIRBOOKING_BROKER_DEFAULT_THRESHOLD, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_LOG_BASE,
                                     g_param_spec_float ("log-base", "", "", 1, G_MAXFLOAT, MIRBOOKING_BROKER_DEFAULT_LOG_BASE, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
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
mirbooking_broker_set_threshold (MirbookingBroker *self,
                                 gfloat            threshold)
{
    self->priv->threshold = threshold;
}

void
mirbooking_broker_set_log_base (MirbookingBroker *self,
                                gfloat            log_base)
{
    self->priv->log_base = log_base;
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

union guintptr
{
    guint u;
    gpointer p;
};

static gpointer
gpointer_from_guint (guint u)
{
    union guintptr ptr = { .u = u };
    return ptr.p;
}

static guint
guint_from_gpointer (gpointer p)
{
    union guintptr ptr = { .p = p };
    return ptr.u;
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
            self->priv->mirnas = g_slist_prepend (self->priv->mirnas, g_object_ref (sequence));
        }
        else
        {
            self->priv->targets = g_slist_prepend (self->priv->targets, g_object_ref (sequence));
        }
    }
}

static int
sequence_cmp_desc (MirbookingSequence *a, MirbookingSequence *b, MirbookingBroker *mirbooking)
{
    gfloat a_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, a));
    gfloat b_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, b));
    return (a_quantity < b_quantity) - (a_quantity > b_quantity);
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
    self->priv->targets = g_slist_sort_with_data (self->priv->targets, (GCompareDataFunc) sequence_cmp_desc, self);
    self->priv->mirnas  = g_slist_sort_with_data (self->priv->mirnas, (GCompareDataFunc) sequence_cmp_desc, self);

    GSList *target_list;
    for (target_list = self->priv->targets; target_list != NULL; target_list = target_list->next)
    {
        target_sites_len += mirbooking_sequence_get_sequence_length (target_list->data) - 7;
    }

    // prepare an contiguous array
    self->priv->target_sites = g_array_sized_new (FALSE,
                                                  FALSE,
                                                  sizeof (MirbookingTargetSite),
                                                  target_sites_len);

    // automatically clear the target sites
    g_array_set_clear_func (self->priv->target_sites,
                            (GDestroyNotify) mirbooking_target_site_clear);

    // intialize sites
    for (target_list = self->priv->targets; target_list != NULL; target_list = target_list->next)
    {
        MirbookingTarget *target = target_list->data;

        gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target));

        gsize position;
        for (position = 0; position < seq_len - 7; position++)
        {
            MirbookingTargetSite target_site;
            target_site.target    = g_object_ref (target);
            target_site.position  = position;
            target_site.occupants = NULL;
            target_site.occupancy = 0;
            g_array_append_val (self->priv->target_sites, target_site);
        }
    }

    // memoize assigned quantities
    // otherwise we would have to lookup the whole #MirbookingTargetSite array
    g_autoptr (GHashTable) assigned_mirna_quantities = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                                         (GEqualFunc) mirbooking_sequence_equal);

    MirbookingTargetSite *target_site = &g_array_index (self->priv->target_sites,
                                                        MirbookingTargetSite,
                                                        0);

    while (target_site < &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len))
    {
        MirbookingTarget *target = target_site->target;

        gfloat available_target_quantity = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, target));

        guint vacancy = available_target_quantity;

        // find the lower target site
        MirbookingTargetSite *from_target_site = MAX (target_site - self->priv->prime5_footprint, &g_array_index (self->priv->target_sites, MirbookingTargetSite, 0));

        // if we overlap the preceding target, move right
        while (from_target_site->target != target) from_target_site++;

        // find the upper target site
        MirbookingTargetSite *to_target_site = MIN (target_site + self->priv->prime3_footprint, &g_array_index (self->priv->target_sites, MirbookingTargetSite, self->priv->target_sites->len));

        // if we overlap the next target, move left
        while (to_target_site->target != target) to_target_site--;

        // minimize vacancy around the footprint
        MirbookingTargetSite *nearby_target_site;
        for (nearby_target_site = from_target_site; nearby_target_site < to_target_site; nearby_target_site++)
        {
            vacancy = MIN (vacancy, available_target_quantity - nearby_target_site->occupancy);
        }

        GSList *prev_mirna_list = NULL;
        GSList *mirna_list;
        for (mirna_list = self->priv->mirnas; mirna_list != NULL; prev_mirna_list = mirna_list, mirna_list = mirna_list->next)
        {
            MirbookingMirna *mirna = mirna_list->data;

            if (vacancy == 0)
            {
                break; // the site has been filled
            }

            guint available_mirna_quantity = floorf (gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, mirna)));
            guint assigned_mirna_quantity  = guint_from_gpointer (g_hash_table_lookup (assigned_mirna_quantities, mirna));
            guint unassigned_mirna_quantity = available_mirna_quantity - assigned_mirna_quantity;

            if (unassigned_mirna_quantity == 0)
            {
                // mirna is depleted

                if (prev_mirna_list == NULL)
                {
                    self->priv->mirnas = mirna_list->next;
                }
                else
                {
                    prev_mirna_list->next = mirna_list->next;
                }

                mirna_list->next = NULL;
                g_slist_free_full (mirna_list, g_object_unref);

                continue;
            }

            g_autoptr (GError) error = NULL;
            gfloat seed_score = mirbooking_score_table_compute_score (self->priv->score_table,
                                                                      MIRBOOKING_SEQUENCE (mirna),
                                                                      1,
                                                                      MIRBOOKING_SEQUENCE (target_site->target),
                                                                      target_site->position,
                                                                      7,
                                                                      &error);

            if (seed_score >= self->priv->threshold)
            {
                guint occupants = floorf (available_target_quantity * logf (unassigned_mirna_quantity) / logf (self->priv->log_base) * seed_score) + 1;

                occupants = MIN (occupants, unassigned_mirna_quantity);
                occupants = MIN (occupants, vacancy);

                g_assert_cmpint (vacancy, >, 0);
                g_assert_cmpint (occupants, >, 0);
                g_assert_cmpint (vacancy, >=, occupants);

                // occupy the site
                MirbookingOccupant *occupant = mirbooking_occupant_new (mirna, occupants);
                target_site->occupants = g_slist_prepend (target_site->occupants, occupant);
                target_site->occupancy += occupants;

                // update the vacancy of the current site
                vacancy -= occupants;

                // update the assigned quantity for the mirna
                g_hash_table_insert (assigned_mirna_quantities,
                                     mirna,
                                     gpointer_from_guint (assigned_mirna_quantity + occupants));
            }
        }

        ++target_site;
    }

    // clear internal quantification & mirnas
    // references of used sequences are kept in the target sites
    g_hash_table_remove_all (self->priv->quantification);

    g_slist_free_full (self->priv->targets, g_object_unref);
    self->priv->targets = NULL;

    g_slist_free_full (self->priv->mirnas, g_object_unref);
    self->priv->mirnas = NULL;

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