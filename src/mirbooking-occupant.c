#include "mirbooking-occupant.h"

void
mirbooking_occupant_init (MirbookingOccupant* self, MirbookingTarget *target, gsize position, MirbookingMirna *mirna, MirbookingScore score)
{
    self->target = g_object_ref (target);
    self->position = position;
    self->mirna = g_object_ref (mirna);
    self->score = score;
}

void
mirbooking_occupant_clear (MirbookingOccupant *self)
{
    g_object_unref (self->target);
    g_object_unref (self->mirna);
}

guint
mirbooking_occupant_hash (const MirbookingOccupant *a)
{
    // this is the djb2 hashing described at http://www.cse.yorku.ca/~oz/hash.html
    // and also used by g_str_hash
    guint hash = 5381;
    hash = hash * 33 + mirbooking_sequence_hash (MIRBOOKING_SEQUENCE (a->target));
    hash = hash * 33 + g_int_hash (&(a->position));
    hash = hash * 33 + mirbooking_sequence_hash (MIRBOOKING_SEQUENCE (a->mirna));
    return hash;
}

    gboolean
mirbooking_occupant_equal (const MirbookingOccupant *a,
                           const MirbookingOccupant *b)
{
    return a->position == b->position &&
        mirbooking_sequence_equal (MIRBOOKING_SEQUENCE (a->mirna), MIRBOOKING_SEQUENCE (b->mirna)) &&
        mirbooking_sequence_equal (MIRBOOKING_SEQUENCE (a->target), MIRBOOKING_SEQUENCE (b->target));
}
