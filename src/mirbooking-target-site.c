#include "mirbooking-target-site.h"

typedef struct
{
    MirbookingTarget *target;
    gsize             offset;
    GSList           *mirnas;
    GHashTable       *quantification; // #MirbookingMirna -> gfloat
} MirbookingTargetSitePrivate;

struct _MirbookingTargetSite
{
    GObject                      parent_instance;
    MirbookingTargetSitePrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingTargetSite, mirbooking_target_site, G_TYPE_OBJECT)

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
mirbooking_target_site_init (MirbookingTargetSite *self)
{
    self->priv = g_new0 (MirbookingTargetSitePrivate, 1);
    self->priv->quantification = g_hash_table_new ((GHashFunc) mirna_hash,
                                                   (GEqualFunc) mirna_equal);
}

enum
{
    PROP_TARGET = 1
};

static void
mirbooking_target_site_set_property (GObject *object, guint property_id, const GValue *value, GParamSpec *pspec)
{
    MirbookingTargetSite *self = MIRBOOKING_TARGET_SITE (object);

    switch (property_id)
    {
        case PROP_TARGET:
            self->priv->target = g_value_get_object (value);
            break;
    }
}

static void
mirbooking_target_site_get_property (GObject *object, guint property_id, GValue *value, GParamSpec *pspec)
{
    MirbookingTargetSite *self = MIRBOOKING_TARGET_SITE (object);

    switch (property_id)
    {
        case PROP_TARGET:
            g_value_set_object (value, self->priv->target);
            break;
    }
}

static void
mirbooking_target_site_finalize (GObject *object)
{
    MirbookingTargetSite *self = MIRBOOKING_TARGET_SITE (object);

    g_hash_table_unref (self->priv->quantification);
    g_slist_free_full (self->priv->mirnas, g_object_unref);
    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_target_site_parent_class)->finalize (object);
}

static void
mirbooking_target_site_class_init (MirbookingTargetSiteClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->set_property = mirbooking_target_site_set_property;
    object_class->get_property = mirbooking_target_site_get_property;
    object_class->finalize     = mirbooking_target_site_finalize;

    g_object_class_install_property (object_class,
                                     PROP_TARGET,
                                     g_param_spec_object ("target", "", "", MIRBOOKING_TYPE_TARGET,
                                                          G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE));
}

MirbookingTargetSite *
mirbooking_target_site_new (MirbookingTarget *target, gsize offset)
{
    return g_object_new (MIRBOOKING_TYPE_TARGET_SITE, "target", target, NULL);
}

gsize
mirbooking_target_site_get_offset (MirbookingTargetSite *self)
{
    return 0;
}

void
mirbooking_target_site_set_mirna_quantity (MirbookingTargetSite *self,
                                           MirbookingMirna      *mirna,
                                           gfloat                quantity)
{
    if (g_hash_table_insert (self->priv->quantification,
                             mirna,
                             (gpointer) (gintptr) quantity))
    {
        self->priv->mirnas = g_slist_prepend (self->priv->mirnas, g_object_ref (mirna));
    }
}

gfloat
mirbooking_target_site_get_mirna_quantity (MirbookingTargetSite *self, MirbookingMirna *mirna)
{
    return (gfloat) (gintptr) g_hash_table_lookup (self->priv->quantification, mirna);
}

GSList *
mirbooking_target_site_get_mirnas (MirbookingTargetSite *self)
{
    return self->priv->mirnas;
}
