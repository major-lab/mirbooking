#include "mirbooking-broker-sparse-solver.h"
#include <sparse.h>

G_DEFINE_AUTOPTR_CLEANUP_FUNC (SparseSolver, sparse_solver_free);

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
            {6, "PARDISO",  "pardiso"},
            {0, NULL,       NULL}
        };

        GType type = g_enum_register_static ("MirbookingBrokerSparseSolver",
                                             values);
        g_once_init_leave (&type_id_once, type);
    }

    return type_id_once;
}

static MirbookingBrokerSparseSolver MIRBOOKING_BROKER_SPARSE_SOLVER_WITH_PRIORITY[] =
{
    MIRBOOKING_BROKER_SPARSE_SOLVER_MKL_DSS,
    MIRBOOKING_BROKER_SPARSE_SOLVER_PARDISO,
    MIRBOOKING_BROKER_SPARSE_SOLVER_UMFPACK,
    MIRBOOKING_BROKER_SPARSE_SOLVER_SUPERLU
};

MirbookingBrokerSparseSolver
mirbooking_broker_sparse_solver_get_default (void)
{
    gint i;
    for (i = 0; i < sizeof (MIRBOOKING_BROKER_SPARSE_SOLVER_WITH_PRIORITY) / sizeof (MIRBOOKING_BROKER_SPARSE_SOLVER_WITH_PRIORITY[0]); i++)
    {
        if (mirbooking_broker_sparse_solver_is_available (MIRBOOKING_BROKER_SPARSE_SOLVER_WITH_PRIORITY[i]))
        {
            return MIRBOOKING_BROKER_SPARSE_SOLVER_WITH_PRIORITY[i];
        }
    }

    /*
     * This is the fallback solver even if it's not available.
     */
    return MIRBOOKING_BROKER_SPARSE_SOLVER_SUPERLU;
}

gboolean
mirbooking_broker_sparse_solver_is_available (MirbookingBrokerSparseSolver sparse_solver)
{
    g_autoptr (SparseSolver) ss = sparse_solver_new ((SparseSolverMethod) sparse_solver);
    return ss != NULL;
}
