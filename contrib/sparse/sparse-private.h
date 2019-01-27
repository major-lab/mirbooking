#ifndef __SPARSE_PRIVATE_H__
#define __SPARSE_PRIVATE_H__

/*
 * Solver implementations have access to the sparse solver struct internals.
 */

struct _SparseSolver
{
    int verbose;
    /* virtual slots */
    void (*init)  (SparseSolver *solver);
    int  (*solve) (SparseSolver *solver, SparseMatrix *A, void *x, void *b);
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

#endif /* __SPARSE_PRIVATE_H__ */
