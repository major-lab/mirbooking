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
};

#endif /* __SPARSE_PRIVATE_H__ */
