#include "mirbooking-target.h"

typedef struct
{
    gsize cds_offset;
    gsize cds_len;
} MirbookingTargetPrivate;

struct _MirbookingTarget
{
    MirbookingSequence       parent_instance;
    MirbookingTargetPrivate *priv;
};

G_DEFINE_TYPE (MirbookingTarget, mirbooking_target, MIRBOOKING_TYPE_SEQUENCE)

static void
mirbooking_target_init (MirbookingTarget *self)
{
    self->priv = g_new0 (MirbookingTargetPrivate, 1);
}

static void
mirbooking_target_finalize (GObject *object)
{
    MirbookingTarget *self = MIRBOOKING_TARGET (object);

    g_free (self->priv);

    G_OBJECT_CLASS (mirbooking_target_parent_class)->finalize (object);
}

static void
mirbooking_target_class_init (MirbookingTargetClass *klass)
{
    GObjectClass *object_class = G_OBJECT_CLASS (klass);

    object_class->finalize = mirbooking_target_finalize;
}

MirbookingTarget *
mirbooking_target_new (const gchar *accession)
{
    return g_object_new (MIRBOOKING_TYPE_TARGET, "accession", accession, NULL);
}

void
mirbooking_target_set_cds (MirbookingTarget *self, gsize offset, gsize len)
{
    gsize seq_len = mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (self));

    g_return_if_fail (offset + len <= seq_len);

    self->priv->cds_offset = offset;
    self->priv->cds_len    = len;
}

const gsize
mirbooking_target_get_cds (MirbookingTarget *self, gsize *len)
{
    gsize ret;

    ret = self->priv->cds_offset;
    *len = self->priv->cds_len;

    return ret;
}
