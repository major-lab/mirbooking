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

G_END_DECLS

#endif /* __MIRBOOKING_OCCUPANT_H__ */
