#include "sparse-matrix.h"
#include "sparse-solver.h"
#include "sparse-solver-private.h"

#include <stdlib.h>
#include <string.h>

/**
 * sparse_solver_new:
 *
 * Returns: A #SparseSolver, or %NULL if the solver is not available.
 */
SparseSolver *
sparse_solver_new (SparseSolverMethod method)
{
    SparseSolver *ret = malloc (sizeof (SparseSolver));

    memset (&ret->statistics, 0, sizeof (SparseSolverStatistics));
    ret->verbose = 0;

    #define PREPARE_SOLVER(solver_uc,solver)  \
    case SPARSE_SOLVER_METHOD_##solver_uc:     \
        ret->init  = sparse_##solver##_init;  \
        ret->clear = sparse_##solver##_clear; \
        ret->solve = sparse_##solver##_solve; \
        break;

    switch (method)
    {
            PREPARE_SOLVER(LAPACK,lapack)
#if HAVE_SUPERLU
            PREPARE_SOLVER(SUPERLU,superlu)
#endif
#if HAVE_SUPERLU_MT
            PREPARE_SOLVER(SUPERLU_MT,superlu_mt)
#endif
#if HAVE_UMFPACK
            PREPARE_SOLVER(UMFPACK,umfpack)
#endif
#if HAVE_MKL_DSS
            PREPARE_SOLVER(MKL_DSS,mkl_dss)
#endif
#if HAVE_MKL_CLUSTER
            PREPARE_SOLVER(MKL_CLUSTER,mkl_cluster)
#endif
#if HAVE_CUSOLVER
            PREPARE_SOLVER(CUSOLVER,cusolver)
#endif
#if HAVE_PARDISO
            PREPARE_SOLVER(PARDISO,pardiso)
#endif
        default:
            return NULL;
    }

    ret->init (ret);

    return ret;
}

/**
 * sparse_solver_set_verbose:
 * @verbose: A verbosity level where 0 is no message (the default) and >=1
 * print resource usage statistics.
 *
 * Set the verbosity level of the solver.
 */
void
sparse_solver_set_verbose (SparseSolver *solver, int verbose)
{
    solver->verbose = verbose;
}

/**
 * sparse_solver_solve:
 * @solver: A #SparseSolver
 * @A: A square sparse matrix
 * @x: A vector of shape A->shape[1] whose element type match A->type
 * @b: A vector of size A->shape[1] whose element type match A->type
 *
 * Solve linear system of the form Ax=b where @A is sparse.
 */
int
sparse_solver_solve (SparseSolver *solver, SparseMatrix *A, void *x, void *b)
{
    return solver->solve (solver, A, x, b);
}

SparseSolverStatistics
sparse_solver_get_statistics (SparseSolver *solver)
{
    return solver->statistics;
}


/**
 * sparse_solver_free:
 */
void
sparse_solver_free (SparseSolver *solver)
{
    solver->clear (solver);
    free (solver);
}
