#ifndef __MIRBOOKING_OCCUPANT_H__
#define __MIRBOOKING_OCCUPANT_H__

#include "mirbooking-mirna.h"
#include "mirbooking-target.h"
#include "mirbooking-score-table.h"

#include <glib.h>

G_BEGIN_DECLS

typedef struct _MirbookingOccupant
{
    MirbookingTarget *target;
    gsize             position;
    MirbookingMirna  *mirna;
    MirbookingScore   score;
} MirbookingOccupant;

void     mirbooking_occupant_init  (MirbookingOccupant *self,
                                    MirbookingTarget   *target,
                                    gsize               position,
                                    MirbookingMirna    *mirna,
                                    MirbookingScore     score);
void     mirbooking_occupant_clear (MirbookingOccupant *self);
gboolean mirbooking_occupant_equal (const MirbookingOccupant *a,
                                    const MirbookingOccupant *b);
guint    mirbooking_occupant_hash  (const MirbookingOccupant *a);

G_END_DECLS

#endif /* __MIRBOOKING_OCCUPANT_H__ */
