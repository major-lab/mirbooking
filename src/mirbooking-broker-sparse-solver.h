#ifndef __MIRBOOKING_BROKER_SPARSE_SOLVER_H__
#define __MIRBOOKING_BROKER_SPARSE_SOLVER_H__

#include <glib-object.h>

G_BEGIN_DECLS

GType mirbooking_broker_sparse_solver_get_type (void) G_GNUC_CONST;
#define MIRBOOKING_BROKER_SPARSE_SOLVER_ENUM mirbooking_broker_sparse_solver_get_type ()
typedef enum _MirbookingBrokerSparseSolver
{
    MIRBOOKING_BROKER_SPARSE_SOLVER_SUPERLU  = 0,
    MIRBOOKING_BROKER_SPARSE_SOLVER_UMFPACK  = 2,
    MIRBOOKING_BROKER_SPARSE_SOLVER_MKL_DSS  = 3,
    MIRBOOKING_BROKER_SPARSE_SOLVER_MKL_CLUSTER = 4,
    MIRBOOKING_BROKER_SPARSE_SOLVER_CUSOLVER = 5,
    MIRBOOKING_BROKER_SPARSE_SOLVER_PARDISO = 6
} MirbookingBrokerSparseSolver;

MirbookingBrokerSparseSolver mirbooking_broker_sparse_solver_get_default (void);

gboolean mirbooking_broker_sparse_solver_is_available (MirbookingBrokerSparseSolver sparse_solver);

G_END_DECLS

#endif /* __MIRBOOKING_BROKER_SPARSE_SOLVER_H__ */
