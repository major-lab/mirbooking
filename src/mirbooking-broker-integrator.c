#include "mirbooking-broker-integrator.h"

GType
mirbooking_broker_integrator_get_type (void)
{
    static gsize type_id_once;

    if (g_once_init_enter (&type_id_once))
    {
        static const GEnumValue values[] =
        {
            {0, "EULER",       "euler"},
            {1, "HEUNS",       "heuns"},
            {2, "RUNGE-KUTTA", "runge-kutta"},
            {0, NULL, NULL}
        };

        GType type = g_enum_register_static ("MirbookingBrokerIntegrator",
                                             values);
        g_once_init_leave (&type_id_once, type);
    }

    return type_id_once;
}
