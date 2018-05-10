#ifndef __MIRBOOKING_BROKER_INTEGRATOR_H__
#define __MIRBOOKING_BROKER_INTEGRATOR_H__

#include <glib-object.h>

G_BEGIN_DECLS

GType mirbooking_broker_integrator_get_type (void) G_GNUC_CONST;
#define MIRBOOKING_BROKER_INTEGRATOR_ENUM mirbooking_broker_integrator_get_type ()

typedef enum _MirbookingBrokerIntegrator
{
    MIRBOOKING_BROKER_INTEGRATOR_EULER = 0,
    MIRBOOKING_BROKER_INTEGRATOR_HEUNS,
    MIRBOOKING_BROKER_INTEGRATOR_RUNGE_KUTTA
} MirbookingBrokerIntegrator;

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_INTEGRATOR_H__ */
