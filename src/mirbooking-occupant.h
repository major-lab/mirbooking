#ifndef __MIRBOOKING_OCCUPANT_H__
#define __MIRBOOKING_OCCUPANT_H__

#include "mirbooking-mirna.h"

#include <glib.h>

G_BEGIN_DECLS

typedef struct _MirbookingOccupant
{
    MirbookingMirna *mirna;
    gdouble          quantity;
    gdouble          cleaved_quantity;
    /* component of the Jacbobian for each step */
    gdouble ES_jac[4];
    gdouble E_jac[4];
    gdouble S_jac[4];
    gdouble P_jac[4];
} MirbookingOccupant;

G_END_DECLS

#endif /* __MIRBOOKING_OCCUPANT_H__ */
