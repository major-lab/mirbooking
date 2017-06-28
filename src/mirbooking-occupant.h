#ifndef __MIRBOOKING_OCCUPANT_H__
#define __MIRBOOKING_OCCUPANT_H__

#include "mirbooking-mirna.h"

#include <glib.h>

G_BEGIN_DECLS

typedef struct _MirbookingOccupant
{
    MirbookingMirna *mirna;
    guint            quantity;
} MirbookingOccupant;

G_END_DECLS

#endif /* __MIRBOOKING_OCCUPANT_H__ */
