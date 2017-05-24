#include "mirbooking.h"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>

typedef struct
{
    gfloat               threshold;
    gfloat               log_base;
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
    PROP_SCORES_TABLE
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
    }
}

static MirbookingMirnaQuantity *
mirbooking_mirna_quantity_new (MirbookingMirna *mirna, gfloat quantity)
{
    MirbookingMirnaQuantity *ret = g_slice_new (MirbookingMirnaQuantity);
    ret->mirna    = mirna;
    ret->quantity = quantity;
    return ret;
}

static void
mirbooking_mirna_quantity_free (MirbookingMirnaQuantity *self, gpointer user_data)
{
    g_slice_free (MirbookingMirnaQuantity, self);
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
        g_slist_free_full (target_site->mirna_quantities, (GDestroyNotify) mirbooking_mirna_quantity_free);
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

void
mirbooking_set_target_quantity (Mirbooking *self, MirbookingTarget *target, gfloat quantity)
{
    if (g_hash_table_insert (self->priv->quantification,
                             target,
                             gpointer_from_gfloat (quantity)));
    {
        self->priv->targets = g_slist_prepend (self->priv->targets, g_object_ref (target));
    }
}

void
mirbooking_set_mirna_quantity (Mirbooking *self, MirbookingMirna *mirna, gfloat quantity)
{
    if (g_hash_table_insert (self->priv->quantification,
                             mirna,
                             gpointer_from_gfloat (quantity)))
    {
        self->priv->mirnas = g_slist_prepend (self->priv->mirnas, g_object_ref (mirna));
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

        gsize site_offset;
        for (site_offset = 0; site_offset < seq_len - 7; site_offset++)
        {
            target_site->target           = target;
            target_site->site_offset      = site_offset;
            target_site->mirna_quantities = NULL;
            ++target_site;
        }
    }

    // memoize assigned quantities
    g_autoptr (GHashTable) assigned_quantities = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                                   (GEqualFunc) mirbooking_sequence_equal);

    target_site = self->priv->target_sites;
    while (target_site < self->priv->target_sites + self->priv->target_sites_len)
    {
        MirbookingTarget *target = target_site->target;

        gfloat available_target_quantity = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, target));
        gfloat assigned_target_quantity = gfloat_from_gpointer (g_hash_table_lookup (assigned_quantities, target));

        // TODO: compute actual vacancy (i.e. scan nearby sites) for the target
        guint vacancy = available_target_quantity - assigned_target_quantity;

        if (vacancy == 0)
        {
            ++target_site; // TODO: jump to the next target
            continue;
        }

        GSList *mirna_list;
        for (mirna_list = self->priv->mirnas; mirna_list != NULL; mirna_list = mirna_list->next)
        {
            MirbookingMirna *mirna = mirna_list->data;

            if (vacancy == 0)
            {
                break; // the site has been filled
            }

            gfloat seed_score = mirbooking_score_table_compute_score (self->priv->score_table,
                                                                      MIRBOOKING_SEQUENCE (target_site->target),
                                                                      target_site->site_offset,
                                                                      MIRBOOKING_SEQUENCE (mirna),
                                                                      1,
                                                                      7);

            if (seed_score >= self->priv->threshold)
            {
                gfloat available_mirna_quantity = gfloat_from_gpointer (g_hash_table_lookup (self->priv->quantification, mirna));
                gfloat assigned_mirna_quantity  = gfloat_from_gpointer (g_hash_table_lookup (assigned_quantities, mirna));

                guint unassigned_mirna_quantity = available_mirna_quantity - assigned_mirna_quantity;

                // the mirna has been used up
                if (unassigned_mirna_quantity == 0)
                {
                    continue;
                }

                guint occupants = (guint) (available_target_quantity * logf (unassigned_mirna_quantity) / logf (self->priv->log_base) * seed_score) + 1;

                occupants = MIN (occupants, unassigned_mirna_quantity);
                occupants = MIN (occupants, vacancy);

                // occupy the site
                MirbookingMirnaQuantity *mirna_quantity = mirbooking_mirna_quantity_new (mirna, occupants);
                mirna_quantity->mirna = mirna;
                mirna_quantity->quantity = occupants;
                target_site->mirna_quantities = g_slist_prepend (target_site->mirna_quantities, mirna_quantity);
                vacancy -= occupants;

                // update the assigned quantity for the mirna and the targe
                g_hash_table_insert (assigned_quantities, mirna, gpointer_from_gfloat (assigned_mirna_quantity + occupants));
                assigned_target_quantity += occupants;
            }
        }

        // update assigned quantity for the target
        g_hash_table_insert (assigned_quantities, target, gpointer_from_gfloat (assigned_target_quantity));

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
