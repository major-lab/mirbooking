#include "mirbooking.h"

#include <stdlib.h>

typedef struct
{
    gdouble               threshold;
    gdouble               log_base;
    MirbookingScoreTable *score_table;

    GSList               *targets;
    GSList               *mirnas;

    GHashTable           *quantification; // #MirbookingSequence -> #gfloat

    GHashTable           *targets_sites;  //  #MirbookingTarget -> #GSList of #MirbookingTargetSite
} MirbookingPrivate;

struct _Mirbooking
{
    GObject            parent_instance;
    MirbookingPrivate *priv;
};

struct _MirbookingTargetSite
{
    gsize   site_offset;    // offset in #MirbookingTarget
    GSList *mirnas;         // all attached #MirbookingMirna
    GHashTable *quantification; // #MirbookingMirna -> gfloat (in the pointer) or %NULL if unused
};

G_DEFINE_TYPE_WITH_PRIVATE (Mirbooking, mirbooking, G_TYPE_OBJECT)

static void
mirbooking_init (Mirbooking *self)
{
    self->priv = g_new0 (MirbookingPrivate, 1);
    self->priv->quantification = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
                                                   (GEqualFunc) mirbooking_sequence_equal);
    self->priv->targets_sites = g_hash_table_new ((GHashFunc) mirbooking_sequence_hash,
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
            self->priv->threshold = g_value_get_double (value);
            break;
        case PROP_LOG_BASE:
            self->priv->log_base = g_value_get_double (value);
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
            g_value_set_double (value, self->priv->threshold);
            break;
        case PROP_LOG_BASE:
            g_value_set_double (value, self->priv->log_base);
            break;
    }
}

static void
mirbooking_target_site_free (MirbookingTargetSite *self)
{
    g_hash_table_unref (self->quantification);
    g_slist_free_full (self->mirnas, g_object_unref);
    g_free (self);
}

static void
free_slist_of_site_target (gpointer key, GSList *list, gpointer user_data)
{
    g_slist_free_full (list, (GDestroyNotify) mirbooking_target_site_free);
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
    g_slist_free_full (self->priv->targets, g_object_unref);
    g_slist_free_full (self->priv->mirnas, g_object_unref);
    g_hash_table_foreach (self->priv->targets_sites, (GHFunc) free_slist_of_site_target, NULL);
    g_hash_table_unref (self->priv->targets_sites);

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
                                     g_param_spec_double ("threshold", "", "", 0, 1, MIRBOOKING_DEFAULT_THRESHOLD, G_PARAM_READWRITE));
    g_object_class_install_property (object_class, PROP_LOG_BASE,
                                     g_param_spec_double ("log-base", "", "", 1, G_MAXDOUBLE, MIRBOOKING_DEFAULT_LOG_BASE, G_PARAM_READWRITE));
}

Mirbooking *
mirbooking_new ()
{
    return g_object_new (MIRBOOKING_TYPE, NULL);
}

void
mirbooking_set_threshold (Mirbooking *self, gdouble threshold)
{
    self->priv->threshold = threshold;
}

void
mirbooking_set_log_base (Mirbooking *self, gdouble log_base)
{
    self->priv->log_base = log_base;
}

void
mirbooking_set_score_table (Mirbooking *self, MirbookingScoreTable *score_table)
{
    self->priv->score_table = g_object_ref (score_table);
}

void
mirbooking_set_target_quantity (Mirbooking *self, MirbookingTarget *target, gfloat quantity)
{
    if (g_hash_table_insert (self->priv->quantification,
                             target,
                             (gpointer) (gintptr) quantity))
    {
        self->priv->targets = g_slist_prepend (self->priv->targets, g_object_ref (target));
    }
}

void
mirbooking_set_mirna_quantity (Mirbooking *self, MirbookingMirna *mirna, gfloat quantity)
{
    if (g_hash_table_insert (self->priv->quantification,
                             mirna,
                             (gpointer) (gintptr) quantity))
    {
        self->priv->mirnas = g_slist_prepend (self->priv->mirnas, g_object_ref (mirna));
    }
}

static gint
mirna_hash (MirbookingSequence *a)
{
    return g_str_hash (mirbooking_sequence_get_accession (a));
}

static gboolean
mirna_equal (MirbookingSequence *a, MirbookingSequence *b)
{
    return g_str_equal (mirbooking_sequence_get_accession (a), mirbooking_sequence_get_accession (b));
}

static void
mirbooking_target_site_set_mirna_quantity (MirbookingTargetSite *self,
                                           MirbookingMirna      *mirna,
                                           gfloat                quantity)
{
    if (g_hash_table_insert (self->quantification,
                             mirna,
                             (gpointer) (gintptr) quantity))
    {
        self->mirnas = g_slist_prepend (self->mirnas, g_object_ref (mirna));
    }
}

gboolean
mirbooking_run (Mirbooking *self, GError **error)
{
    g_return_val_if_fail (MIRBOOKING_IS_MIRBOOKING (self), FALSE);

    GSList *target;
    for (target = self->priv->targets; target != NULL; target = target->next)
    {
        gsize seq_len = mirbooking_sequence_get_sequence_length (target->data);

        MirbookingTargetSite *target_sites = g_new0 (MirbookingTargetSite,
                                                     seq_len - 7);

        // build a *fake* linked-list even though they are all stored contigously
        GSList *target_sites_list = NULL;

        gsize site_offset;
        for (site_offset = 0; site_offset < seq_len - 7; site_offset++)
        {
            target_sites[site_offset].site_offset = site_offset;
            g_slist_prepend (target_sites_list, &target_sites[site_offset]);
        }

        g_hash_table_insert (self->priv->targets_sites,
                             target->data,
                             target_sites_list);

        /*
        gsize site_offset;
        for (site_offset = 0; site_offset < ; site_offset++)
        {
            // some random offset on the target for now
            gsize site_offset = random () % (mirbooking_sequence_get_sequence_length (target->data) - 7);

            MirbookingTargetSite *target_site = target_sites + offset;

            target_site->offset = offset;
            target_site->quantification = g_hash_table_new ((GHashFunc) mirna_hash,
                                                            (GEqualFunc) mirna_equal);

            GSList *mirna;
            for (mirna = self->priv->mirnas; mirna != NULL; mirna = mirna->next)
            {
                gfloat score;

                score = mirbooking_score_table_compute_score (self->priv->score_table,
                                                              target->data,
                                                              site_offset,
                                                              mirna->data,
                                                              1,
                                                              7);

                if (score >= self->priv->threshold)
                {
                    gfloat quantity = (gfloat) (gintptr) g_hash_table_lookup (self->priv->quantification, mirna->data);
                    mirbooking_target_site_set_mirna_quantity (target_site, mirna->data, quantity);
                }
            }

            GSList *target_sites = g_hash_table_lookup (self->priv->targets_sites,
                                                        target->data);

            g_hash_table_insert (self->priv->targets_sites,
                                 target->data,
                                 g_slist_prepend (target_sites, target_site));
        }
        */
    }

    return TRUE;
}

GSList *
mirbooking_get_targets (Mirbooking *self)
{
    return self->priv->targets;
}

GSList *
mirbooking_get_target_sites (Mirbooking *self, MirbookingTarget *target)
{
    return g_hash_table_lookup (self->priv->targets_sites, target);
}

gsize
mirbooking_target_site_get_site_offset (MirbookingTargetSite *self)
{
    return self->site_offset;
}

GSList *
mirbooking_target_site_get_mirnas (MirbookingTargetSite *self)
{
    return self->mirnas;
}

gfloat
mirbooking_target_site_get_mirna_quantity (MirbookingTargetSite *self, MirbookingMirna *mirna)
{
    return (gfloat) (gintptr) g_hash_table_lookup (self->quantification, mirna);
}
