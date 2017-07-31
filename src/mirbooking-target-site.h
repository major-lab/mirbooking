#ifndef __MIRBOOKING_TARGET_SITE_H__
#define __MIRBOOKING_TARGET_SITE_H__

#include "mirbooking-target.h"

#include <glib.h>

G_BEGIN_DECLS

typedef struct _MirbookingTargetSite MirbookingTargetSite;

/**
 * MirbookingTargetSite:
 * @occupants: (element-type MirbookingOccupant) A #GSList of occupants
 *
 * The @occupants entry could have been a #GArray, but in practice, most target
 * sites are unoccupied so it's much more compact to use the default %NULL
 * initializer for the linked list.
 */
struct _MirbookingTargetSite
{
    MirbookingTarget *target;
    gsize             position;
    GSList           *occupants; // #GSList of #MirbookingOccupant
};

#endif /* __MIRBOOKING_TARGET_SITE_H__ */
