#include "mirbooking-mirna.h"

struct _MirbookingMirna
{
    MirbookingSequence parent_instance;
};

G_DEFINE_TYPE (MirbookingMirna, mirbooking_mirna, MIRBOOKING_TYPE_SEQUENCE)

static void
mirbooking_mirna_init (MirbookingMirna *self)
{

}

static void
mirbooking_mirna_class_init (MirbookingMirnaClass *klass)
{

}

MirbookingMirna *
mirbooking_mirna_new (const gchar *accession)
{
    return g_object_new (MIRBOOKING_TYPE_MIRNA, "accession", accession, NULL);
}

MirbookingMirna *
mirbooking_mirna_new_with_name (const gchar *accession, const gchar *name)
{
    return g_object_new (MIRBOOKING_TYPE_MIRNA, "accession", accession, "name", name, NULL);
}
