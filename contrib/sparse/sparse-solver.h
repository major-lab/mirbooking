#ifndef __SPARSE_SOLVER_H__
#define __SPARSE_SOLVER_H__

typedef struct _SparseSolver SparseSolver;

typedef enum _SparseSolverMethod
{
    SPARSE_SOLVER_METHOD_LAPACK,
    SPARSE_SOLVER_METHOD_SUPERLU,
    SPARSE_SOLVER_METHOD_SUPERLU_MT,
    SPARSE_SOLVER_METHOD_UMFPACK,
    SPARSE_SOLVER_METHOD_MKL_DSS,
    SPARSE_SOLVER_METHOD_MKL_CLUSTER,
    SPARSE_SOLVER_METHOD_CUSOLVER,
    SPARSE_SOLVER_METHOD_PARDISO
} SparseSolverMethod;

typedef struct _SparseSolverStatistics
{
    double reorder_time;
    double factor_time;
    double solve_time;
    double flops;
} SparseSolverStatistics;

SparseSolver *         sparse_solver_new            (SparseSolverMethod solver_method);
int                    sparse_solver_solve          (SparseSolver *solver, SparseMatrix *A, void *x, void *b);
SparseSolverStatistics sparse_solver_get_statistics (SparseSolver *solver);
void                   sparse_solver_set_verbose    (SparseSolver *solver, int verbose);
void                   sparse_solver_free           (SparseSolver *solver);

#endif /* __SPARSE_SOLVER_H__ */
