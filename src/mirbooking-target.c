#include "mirbooking-target.h"

typedef struct _MirbookingTargetPrivate
{
    gfloat *accessibility_scores; /* accessibility of each position */
} MirbookingTargetPrivate;

struct _MirbookingTarget
{
    MirbookingSequence       parent_instance;
    MirbookingTargetPrivate *priv;
};

G_DEFINE_TYPE_WITH_PRIVATE (MirbookingTarget, mirbooking_target, MIRBOOKING_TYPE_SEQUENCE)

static void
mirbooking_target_init (MirbookingTarget *self)
{
    self->priv = g_new0 (MirbookingTargetPrivate, 1);
}

static void
mirbooking_target_finalize (GObject *object)
{
    MirbookingTarget *self = MIRBOOKING_TARGET (object);

    if (self->priv->accessibility_scores)
    {
        g_free (self->priv->accessibility_scores);
    }
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

MirbookingTarget *
mirbooking_target_new_with_name (const gchar *accession, const gchar *name)
{
    return g_object_new (MIRBOOKING_TYPE_TARGET, "accession", accession, "name", name, NULL);
}

gfloat
mirbooking_target_get_accessibility_score (MirbookingTarget *self,
                                           gsize             position)
{
    return self->priv->accessibility_scores ? self->priv->accessibility_scores[position] : 0;
}

void
mirbooking_target_set_accessibility_score (MirbookingTarget *self,
                                           gsize             position,
                                           gfloat            score)
{
    if (self->priv->accessibility_scores == NULL)
    {
        self->priv->accessibility_scores = g_new0 (gfloat, mirbooking_sequence_get_sequence_length (MIRBOOKING_SEQUENCE (self)));
    }
    self->priv->accessibility_scores[position] = score;
}
