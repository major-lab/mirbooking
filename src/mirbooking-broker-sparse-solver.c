#include "mirbooking-broker-sparse-solver.h"

GType
mirbooking_broker_sparse_solver_get_type (void)
{
    static gsize type_id_once;

    if (g_once_init_enter (&type_id_once))
    {
        static const GEnumValue values[] =
        {
            {0, "SUPERLU",  "superlu"},
            {2, "UMFPACK",  "umfpack"},
            {3, "MKL_DSS",  "mkl-dss"},
            {4, "MKL_CLUSTER", "mkl-cluster"},
            {5, "CUSOLVER", "cusolver"},
            {0, NULL,       NULL}
        };

        GType type = g_enum_register_static ("MirbookingBrokerSparseSolver",
                                             values);
        g_once_init_leave (&type_id_once, type);
    }

    return type_id_once;
}
