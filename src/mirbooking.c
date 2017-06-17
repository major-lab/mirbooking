#include "mirbooking.h"

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
    MirbookingTargetSite *target_sites;
    gsize                 target_sites_len;
} MirbookingPrivate;

struct _Mirbooking
{
    GObject            parent_instance;
    MirbookingPrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (Mirbooking, mirbooking, G_TYPE_OBJECT)

static void
mirbooking_init (Mirbooking *self)
{
    self->priv = g_new0 (MirbookingPrivate, 1);
    self->priv->quantification = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                   (GEqualFunc) mirbooking_sequence_equal);
}

enum
{
    PROP_THRESHOLD = 1,
    PROP_LOG_BASE,
    PROP_SCORES_TABLE,
    PROP_5PRIME_FOOTPRINT,
    PROP_3PRIME_FOOTPRINT
};

static void
mirbooking_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    Mirbooking *self = MIRBOOKING_MIRBOOKING (object);

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
        default:
            g_assert_not_reached ();
    }
}
static void
mirbooking_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    Mirbooking *self = MIRBOOKING_MIRBOOKING (object);

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
        default:
            g_assert_not_reached ();
    }
}

static MirbookingOccupant *
mirbooking_occupant_new (MirbookingMirna *mirna, guint quantity)
{
    MirbookingOccupant *ret = g_slice_new (MirbookingOccupant);
    ret->mirna              = mirna;
    ret->quantity           = quantity;
    return ret;
}

static void
mirbooking_occupant_free (MirbookingOccupant *self, gpointer user_data)
{
    g_slice_free (MirbookingOccupant, self);
}

static void
mirbooking_finalize (GObject *object)
{
    Mirbooking *self = MIRBOOKING_MIRBOOKING (object);

    if (self->priv->score_table)
    {
        g_object_unref (self->priv->score_table);
    }

    g_hash_table_unref (self->priv->quantification);

    gint i;
    MirbookingTargetSite *target_site = self->priv->target_sites;
    for (i = 0; i < self->priv->target_sites_len; i++)
    {
        g_slist_free_full (target_site->occupants, (GDestroyNotify) mirbooking_occupant_free);
        ++target_site;
    }
    g_free (self->priv->target_sites);

    g_slist_free_full (self->priv->targets, g_object_unref);
    g_slist_free_full (self->priv->mirnas, g_object_unref);

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_parent_class)->finalize (object);
}

static void
mirbooking_class_init (MirbookingClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->set_property = mirbooking_set_property;
    object_class->get_property = mirbooking_get_property;
    object_class->finalize     = mirbooking_finalize;

    g_object_class_install_property (object_class, PROP_THRESHOLD,
                                     g_param_spec_float ("threshold", "", "", 0, 1, MIRBOOKING_DEFAULT_THRESHOLD, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_LOG_BASE,
                                     g_param_spec_float ("log-base", "", "", 1, G_MAXFLOAT, MIRBOOKING_DEFAULT_LOG_BASE, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_5PRIME_FOOTPRINT,
                                     g_param_spec_uint ("prime5-footprint", "", "", 0, G_MAXUINT, MIRBOOKING_DEFAULT_5PRIME_FOOTPRINT, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_3PRIME_FOOTPRINT,
                                     g_param_spec_uint ("prime3-footprint", "", "", 0, G_MAXUINT, MIRBOOKING_DEFAULT_3PRIME_FOOTPRINT, G_PARAM_CONSTRUCT | G_PARAM_READWRITE));
}

Mirbooking *
mirbooking_new ()
{
    return g_object_new (MIRBOOKING_TYPE, NULL);
}

void
mirbooking_set_threshold (Mirbooking *self, gfloat threshold)
{
    self->priv->threshold = threshold;
}

void
mirbooking_set_log_base (Mirbooking *self, gfloat log_base)
{
    self->priv->log_base = log_base;
}

void
mirbooking_set_5prime_footprint (Mirbooking *self,
                                 gsize       footprint)
{
    self->priv->prime5_footprint = footprint;
}

void
mirbooking_set_3prime_footprint (Mirbooking *self,
                                 gsize       footprint)
{
    self->priv->prime3_footprint = footprint;
}

void
mirbooking_set_score_table (Mirbooking *self, MirbookingScoreTable *score_table)
{
    self->priv->score_table = g_object_ref (score_table);
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

void
mirbooking_set_sequence_quantity (Mirbooking *self, MirbookingSequence *sequence, gfloat quantity)
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
sequence_cmp_desc (MirbookingSequence *a, MirbookingSequence *b, Mirbooking *mirbooking)
{
    gfloat a_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, a));
    gfloat b_quantity = gfloat_from_gpointer (g_hash_table_lookup (mirbooking->priv->quantification, b));
    return (a_quantity < b_quantity) - (a_quantity > b_quantity);
}

gboolean
mirbooking_run (Mirbooking *self, GError **error)
{
    gsize target_sites_len = 0;

    g_return_val_if_fail (self != NULL, FALSE);
    g_return_val_if_fail (self->priv->score_table != NULL, FALSE);

    // sort internal targets and mirnas in descending quantity
    self->priv->targets = g_slist_sort_with_data (self->priv->targets, (GCompareDataFunc) sequence_cmp_desc, self);
    self->priv->mirnas  = g_slist_sort_with_data (self->priv->mirnas, (GCompareDataFunc) sequence_cmp_desc, self);

    GSList *target_list;
    for (target_list = self->priv->targets; target_list != NULL; target_list = target_list->next)
    {
        target_sites_len += mirbooking_sequence_get_sequence_length (target_list->data) - 7;
    }

    // allocate *all* the target sites in a single & contiguous chunk
    self->priv->target_sites     = g_new (MirbookingTargetSite, target_sites_len);
    self->priv->target_sites_len = target_sites_len;

    MirbookingTargetSite *target_site = self->priv->target_sites;

    // intialize sites
    for (target_list = self->priv->targets; target_list != NULL; target_list = target_list->next)
    {
        MirbookingTarget *target = target_list->data;

        gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (target));

        gsize position;
        for (position = 0; position < seq_len - 7; position++)
        {
            target_site->target           = target;
            target_site->position      = position;
            target_site->occupants = NULL;
            ++target_site;
        }
    }

    // memoize assigned quantities
    // otherwise we would have to lookup the whole #MirbookingTargetSite array
    g_autoptr (GHashTable) assigned_mirna_quantities = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                                         (GEqualFunc) mirbooking_sequence_equal);

    target_site = self->priv->target_sites;
    while (target_site < self->priv->target_sites + self->priv->target_sites_len)
    {
        MirbookingTarget *target = target_site->target;

        gfloat available_target_quantity = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, target));

        guint vacancy = available_target_quantity;

        // find the lower target site
        MirbookingTargetSite *from_target_site = target_site;
        while (from_target_site > self->priv->target_sites           &&
               (from_target_site - 1)->target == target_site->target &&
               (target_site - from_target_site) < self->priv->prime5_footprint)
        {
            --from_target_site;
        }

        // find the upper target site
        MirbookingTargetSite *to_target_site = target_site;
        while (to_target_site < self->priv->target_sites + self->priv->target_sites_len &&
               (to_target_site + 1)->target == target_site->target                      &&
               (to_target_site - target_site) < self->priv->prime3_footprint)
        {
            ++to_target_site;
        }

        // minimize vacancy around the footprint
        MirbookingTargetSite *nearby_target_site;
        for (nearby_target_site = from_target_site; nearby_target_site < to_target_site; nearby_target_site++)
        {
            vacancy = MIN (vacancy, available_target_quantity - nearby_target_site->occupancy);
        }

        GSList *mirna_list;
        for (mirna_list = self->priv->mirnas; mirna_list != NULL; mirna_list = mirna_list->next)
        {
            MirbookingMirna *mirna = mirna_list->data;

            if (vacancy == 0)
            {
                break; // the site has been filled
            }

            gfloat available_mirna_quantity = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, mirna));
            guint assigned_mirna_quantity  = guint_from_gpointer (g_hash_table_lookup (assigned_mirna_quantities, mirna));

            if (available_mirna_quantity - assigned_mirna_quantity < 1.0f)
            {
                continue; // mirna is depleted
            }

            gfloat seed_score = mirbooking_score_table_compute_score (self->priv->score_table,
                                                                      MIRBOOKING_SEQUENCE (target_site->target),
                                                                      target_site->position,
                                                                      MIRBOOKING_SEQUENCE (mirna),
                                                                      1,
                                                                      7);

            if (seed_score >= self->priv->threshold)
            {
                gfloat unassigned_mirna_quantity = available_mirna_quantity - assigned_mirna_quantity;

                guint occupants = floorf (available_target_quantity * logf (unassigned_mirna_quantity) / logf (self->priv->log_base) * seed_score) + 1;

                occupants = MIN (occupants, floorf (unassigned_mirna_quantity));
                occupants = MIN (occupants, vacancy);

                g_assert_cmpint (vacancy, >, 0);
                g_assert_cmpint (occupants, >, 0);
                g_assert_cmpint (vacancy, >=, occupants);

                // occupy the site
                MirbookingOccupant *occupant = mirbooking_occupant_new (mirna, occupants);
                occupant->mirna = mirna;
                occupant->quantity = occupants;
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

    return TRUE;
}

MirbookingTargetSite *
mirbooking_get_target_sites (Mirbooking *self, gsize *len)
{
    g_return_val_if_fail (self != NULL, NULL);
    g_return_val_if_fail (self->priv->target_sites != NULL, NULL);

    *len = self->priv->target_sites_len;
    return self->priv->target_sites;
}
