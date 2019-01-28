#include "sparse.h"
#include "sparse-private.h"
#include <suitesparse/umfpack.h>
#include <math.h>
#include <assert.h>

#define memcpy_loop(to, from, n) \
{                                \
    int i;                       \
    for (i = 0; i < (n); i++)    \
    {                            \
        (to)[i] = (from)[i];     \
                                 \
    }                            \
}

void
sparse_umfpack_init (SparseSolver *solver)
{

}

void
sparse_umfpack_clear (SparseSolver *solver)
{

}

typedef struct _SolverStorage {
    void *symbolic;
} SolverStorage;

static void
free_solver_storage (void* ptr)
{
    SolverStorage *storage = ptr;
    if (storage->symbolic)
        umfpack_dl_free_symbolic (&storage->symbolic);
    free (ptr);
}

int
sparse_umfpack_solve (SparseSolver *solver, SparseMatrix *A, void *x, void *b)
{
    double control[UMFPACK_CONTROL];
    double info[UMFPACK_INFO];

    umfpack_dl_defaults (control);

    control[UMFPACK_PRL] = solver->verbose ? 2 : UMFPACK_DEFAULT_PRL;

    assert (A->type == SPARSE_MATRIX_TYPE_DOUBLE);
    assert (A->shape[0] < LONG_MAX);
    assert (A->s.csr.nnz < LONG_MAX);

    if (A->shape[0] == 0 || A->shape[1] == 0)
    {
        return 1;
    }

    long *rowptr, *colind;

    if (sizeof (long) == sizeof (size_t))
    {
        rowptr = (long*) A->s.csr.rowptr;
        colind = (long*) A->s.csr.colind;
    }
    else
    {
        rowptr = malloc ((A->shape[0] + 1) * sizeof (size_t));
        colind = malloc ((A->shape[0] + 1) * sizeof (size_t));
        memcpy_loop (rowptr, A->s.csr.rowptr, A->shape[0] + 1);
        memcpy_loop (colind, A->s.csr.colind, A->s.csr.nnz);
    }

    SolverStorage *solver_storage = A->solver_storage;
    if (A->solver_storage_owner != solver)
    {
        solver_storage = calloc (1, sizeof (SolverStorage));

        umfpack_dl_symbolic (A->shape[0],
                             A->shape[1],
                             rowptr,
                             colind,
                             A->data,
                             &solver_storage->symbolic,
                             control,
                             info);

        if (info[UMFPACK_STATUS] != UMFPACK_OK)
        {
            free_solver_storage (solver_storage);
            goto cleanup;
        }

        sparse_matrix_set_solver_storage (A,
                                          solver_storage,
                                          free_solver_storage,
                                          solver);
    }

    void *numeric;
    umfpack_dl_numeric (rowptr,
                        colind,
                        A->data,
                        solver_storage->symbolic,
                        &numeric,
                        control,
                        info);

    if (info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
    {
        goto cleanup_numeric;
    }
    else if (info[UMFPACK_STATUS] != UMFPACK_OK)
    {
        goto cleanup;
    }

    umfpack_dl_solve (UMFPACK_At,
                      rowptr,
                      colind,
                      A->data,
                      x,
                      b,
                      numeric,
                      control,
                      info);

    if (info[UMFPACK_STATUS] != UMFPACK_OK)
    {
        goto cleanup_numeric;
    }

    solver->statistics.reorder_time = info[UMFPACK_SYMBOLIC_WALLTIME];
    solver->statistics.factor_time  = info[UMFPACK_NUMERIC_WALLTIME];
    solver->statistics.solve_time   = info[UMFPACK_SOLVE_WALLTIME];
    solver->statistics.flops        = info[UMFPACK_FLOPS];

cleanup_numeric:
    umfpack_dl_free_numeric (&numeric);

cleanup:
    umfpack_dl_report_info (control, info);
    umfpack_dl_report_status (control, info[UMFPACK_STATUS]);

    if (sizeof (long) != sizeof (size_t))
    {
        free (rowptr);
        free (colind);
    }

    return info[UMFPACK_STATUS] == UMFPACK_OK;
}
