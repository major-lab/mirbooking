#ifndef __SPARSE_SOLVER_PRIVATE_H__
#define __SPARSE_SOLVER_PRIVATE_H__

#define DECLARE_SOLVER(solver)                  \
extern void                                     \
sparse_##solver##_init  (SparseSolver *solver); \
extern void                                     \
sparse_##solver##_clear (SparseSolver *solver); \
extern int                                      \
sparse_##solver##_solve (SparseSolver *solver,  \
                         SparseMatrix *A,       \
                         void         *x,       \
                         const void   *b);

DECLARE_SOLVER(lapack)
DECLARE_SOLVER(superlu)
DECLARE_SOLVER(superlu_mt)
DECLARE_SOLVER(umfpack)
DECLARE_SOLVER(mkl_dss)
DECLARE_SOLVER(mkl_cluster)
DECLARE_SOLVER(cusolver)
DECLARE_SOLVER(pardiso)

#define memcpy_loop(to, from, n) \
{                                \
    size_t i;                    \
    for (i = 0; i < (n); i++)    \
    {                            \
        (to)[i] = (from)[i];     \
                                 \
    }                            \
}

/*
 * Solver implementations have access to the sparse solver struct internals.
 */

struct _SparseSolver
{
    int verbose;
    /* virtual slots */
    void (*init)  (SparseSolver *solver);
    int  (*solve) (SparseSolver *solver, SparseMatrix *A, void *x, const void *b);
    void (*clear) (SparseSolver *solver);
    /* solver-specific private storage */
    void *storage;
    SparseSolverStatistics statistics;
};

static inline void
sparse_matrix_set_solver_storage (SparseMatrix *M, void *solver_storage, void (*solver_storage_destroy) (void*), void *owner)
{
    if (M->solver_storage_destroy)
    {
        M->solver_storage_destroy (M->solver_storage);
    }

    M->solver_storage_owner = owner;
    M->solver_storage = solver_storage;
    M->solver_storage_destroy = solver_storage_destroy;
}

#endif /* __SPARSE_SOLVER_PRIVATE_H__ */
