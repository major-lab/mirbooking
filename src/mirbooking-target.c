#include "mirbooking-target.h"

struct _MirbookingTarget
{
    MirbookingSequence parent_instance;
};

G_DEFINE_TYPE (MirbookingTarget, mirbooking_target, MIRBOOKING_TYPE_SEQUENCE)

static void
mirbooking_target_init (MirbookingTarget *self)
{
}

static void
mirbooking_target_class_init (MirbookingTargetClass *klass)
{

}

MirbookingTarget *
mirbooking_target_new (const gchar *accession)
{
    return g_object_new (MIRBOOKING_TYPE_TARGET, "accession", accession, NULL);
}
