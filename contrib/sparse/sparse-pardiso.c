#include "sparse-matrix.h"
#include "sparse-solver.h"
#include "sparse-solver-private.h"

#include <stdlib.h>
#include <assert.h>
#include <omp.h>

extern void pardisoinit        (void *, int *, int *, int *, double *, int *);
extern void pardiso            (void *, int *, int *, int *, int *, int *,
                                double *, int *, int *, int *, int *, int *,
                                int *, const double *, double *, int *, double *);
extern void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern void pardiso_chkvec     (int *, int *, const double *, int *);
extern void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                                double *, int *);

typedef struct _SolverStorage
{
    void *pt[64];
    int  *rowptr;
    int  *colind;
} SolverStorage;

static void
free_solver_storage (void *ptr)
{
    SolverStorage *solver_storage = ptr;

    int iparm[64] = {0};
    double dparm[64] = {0};
    int phase = -1;
    int cleanup_error;
    int msglvl = 0;
    int maxfct = 1;
    pardiso (solver_storage->pt, &maxfct, NULL, NULL, &phase, NULL, NULL, NULL,
             NULL, NULL,NULL, iparm, &msglvl, NULL, NULL, &cleanup_error, dparm);
    assert (cleanup_error == 0);

    free (solver_storage->rowptr);
    free (solver_storage->colind);
    free (ptr);
}

void
sparse_pardiso_init (SparseSolver *solver)
{

}

void
sparse_pardiso_clear (SparseSolver *solver)
{

}

int
sparse_pardiso_solve (SparseSolver *solver, SparseMatrix *A, void *x, const void *b)
{
    int maxfct = 1;
    int mtype;
    int phase;
    int n = A->shape[0];
    int perm;
    int nrhs = 1;
    int iparm[64] = {0};
    int msglvl = solver->verbose ? 1 : 0;
    int error = 0;
    int mnum = 1;
    double dparm[64] = {0};

    if (A->hints & SPARSE_MATRIX_HINT_SYMMETRIC && A->hints & SPARSE_MATRIX_HINT_POSITIVE_DEFINITE)
    {
        mtype = 2;
    }
    else if (A->hints & SPARSE_MATRIX_HINT_SYMMETRIC)
    {
        mtype = -2;
    }
    else if (A->hints & SPARSE_MATRIX_HINT_SYMMETRIC_STRUCTURE)
    {
        mtype = 1;
    }
    else
    {
        mtype = 11;
    }

    iparm[0] = 0; /* use defaults */
    iparm[2] = omp_get_num_threads (); /* nproc */

    SolverStorage *solver_storage = A->solver_storage;
    if (A->solver_storage_owner != solver)
    {
        solver_storage = malloc (sizeof (SolverStorage));

        int solver_ = 0;
        pardisoinit (solver_storage->pt, &mtype, &solver_, iparm, dparm, &error);
        if (error != 0)
        {
            free (solver_storage);
            goto cleanup;
        }

        solver_storage->rowptr = malloc ((A->shape[0] + 1) * sizeof (int));
        solver_storage->colind = malloc (A->s.csr.nnz * sizeof (int));
        memcpy_loop (solver_storage->rowptr, A->s.csr.rowptr, A->shape[0] + 1);
        memcpy_loop (solver_storage->colind, A->s.csr.colind, A->s.csr.nnz);

        int i;
        for (i = 0; i < n + 1; i++)
        {
            solver_storage->rowptr[i] += 1;
        }

        for (i = 0; i < A->s.csr.nnz; i++)
        {
            solver_storage->colind[i] += 1;
        }

        pardiso_chkmatrix (&mtype, &n, A->data, solver_storage->rowptr, solver_storage->colind, &error);
        assert (error == 0);

        pardiso_chkvec (&n, &nrhs, b, &error);
        assert (error == 0);

        // reorder
        phase = 11;
        pardiso (solver_storage->pt, &maxfct, &mnum, &mtype, &phase,
                 &n, A->data, solver_storage->rowptr, solver_storage->colind, &perm, &nrhs,
                 iparm, &msglvl, b, x, &error, dparm);
        if (error != 0)
        {
            free_solver_storage (solver_storage);
            goto cleanup;
        }

        sparse_matrix_set_solver_storage (A, solver_storage, free_solver_storage, solver);
    }

    // factor and solve
    phase = 23;
    pardiso (solver_storage->pt, &maxfct, &mnum, &mtype, &phase,
             &n, A->data, solver_storage->rowptr, solver_storage->colind, &perm, &nrhs,
             iparm, &msglvl, b, x, &error, dparm);

cleanup:
    return error == 0;
}
