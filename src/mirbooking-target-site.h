#ifndef __MIRBOOKING_TARGET_SITE_H__
#define __MIRBOOKING_TARGET_SITE_H__

#include "mirbooking-target.h"

#include <glib.h>

G_BEGIN_DECLS

typedef struct _MirbookingTargetSite MirbookingTargetSite;

/**
 * MirbookingTargetSite:
 * @occupants: (element-type MirbookingOccupant)
 */
struct _MirbookingTargetSite
{
    MirbookingTarget *target;
    gsize             position;
    GSList           *occupants; // #GSList of #MirbookingOccupant
    guint             occupancy; // sum of all occupants' quantities
};

#endif /* __MIRBOOKING_TARGET_SITE_H__ */
