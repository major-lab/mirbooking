#include "sparse.h"
#include "sparse-private.h"

#include <assert.h>
#include <mkl_dss.h>
#include <stdio.h>
#include <stdlib.h>

void
sparse_mkl_dss_init (SparseSolver *solver)
{

}

void
sparse_mkl_dss_clear (SparseSolver *solver)
{

}

typedef struct _SolverStorage
{
    _MKL_DSS_HANDLE_t handle;
    MKL_INT *rowptr;
    MKL_INT *colind;
} SolverStorage;

void
free_solver_storage (void *ptr)
{
    SolverStorage *solver_storage = ptr;
    MKL_INT delete_opt = 0;
    dss_delete (solver_storage->handle, delete_opt);
    if (sizeof (MKL_INT) != sizeof (size_t))
    {
        free (solver_storage->rowptr);
        free (solver_storage->colind);
    }

    free (solver_storage);
}

int
sparse_mkl_dss_solve (SparseSolver *solver, SparseMatrix *A, void *x, void *b)
{
    int ret = 0;

    if (A->shape[0] == 0 || A->shape[1] == 0)
    {
        return 1;
    }

    SolverStorage *solver_storage = A->solver_storage;
    if (A->solver_storage_owner != solver)
    {
        solver_storage = malloc (sizeof (SolverStorage));

        if (sizeof (MKL_INT) == sizeof (size_t))
        {
            solver_storage->rowptr  = (MKL_INT*) A->s.csr.rowptr;
            solver_storage->colind  = (MKL_INT*) A->s.csr.colind;
        }
        else
        {
            solver_storage->rowptr  = malloc ((A->shape[0] + 1) * sizeof (MKL_INT));
            solver_storage->colind  = malloc (A->s.csr.nnz * sizeof (MKL_INT));
            memcpy_loop (solver_storage->rowptr, A->s.csr.rowptr, A->shape[0] + 1);
            memcpy_loop (solver_storage->colind, A->s.csr.colind, A->s.csr.nnz);
        }

        MKL_INT create_opt = MKL_DSS_ZERO_BASED_INDEXING;

        if (A->type == SPARSE_MATRIX_TYPE_FLOAT)
        {
            create_opt |= MKL_DSS_SINGLE_PRECISION;
        }
        else
        {
            assert (A->type == SPARSE_MATRIX_TYPE_DOUBLE);
        }

        ret = dss_create (solver_storage->handle, create_opt);

        if (ret != MKL_DSS_SUCCESS)
        {
            free_solver_storage (solver_storage);
            return 0;
        }

        MKL_INT define_opt = 0;
        MKL_INT m = A->shape[0];
        MKL_INT n = A->shape[1];
        MKL_INT nnz = A->s.csr.nnz;

        assert (m == A->shape[0]);
        assert (n == A->shape[1]);
        assert (nnz == A->s.csr.nnz);

        if (A->hints & SPARSE_MATRIX_HINT_SYMMETRIC)
        {
            define_opt |= MKL_DSS_SYMMETRIC;
        }
        else
        {
            define_opt |= MKL_DSS_NON_SYMMETRIC;
        }

        ret = dss_define_structure (solver_storage->handle,
                                    define_opt,
                                    solver_storage->rowptr,
                                    m,
                                    n,
                                    solver_storage->colind,
                                    nnz);
        if (ret != MKL_DSS_SUCCESS)
        {
            free_solver_storage (solver_storage);
            goto cleanup;
        }

        MKL_INT reorder_opt = 0;
        ret = dss_reorder (solver_storage->handle, reorder_opt, NULL);
        if (ret != MKL_DSS_SUCCESS)
        {
            free_solver_storage (solver_storage);
            goto cleanup;
        }

        sparse_matrix_set_solver_storage (A,
                                          solver_storage,
                                          free_solver_storage,
                                          solver);
    }

    MKL_INT factor_opt = 0;

    if (A->hints & SPARSE_MATRIX_HINT_POSITIVE_DEFINITE)
    {
        factor_opt |= MKL_DSS_POSITIVE_DEFINITE;
    }
    else
    {
        factor_opt |= MKL_DSS_INDEFINITE;
    }

    ret = dss_factor_real (solver_storage->handle, factor_opt, A->data);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    MKL_INT solve_opt = 0;
    MKL_INT nrhs = 1;
    ret = dss_solve_real (solver_storage->handle, solve_opt, b, nrhs, x);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    MKL_INT statistics_opt = 0;
    double ret_values[7];
    ret = dss_statistics (solver_storage->handle,
                          statistics_opt,
                          "ReorderTime,FactorTime,SolveTime,Flops",
                          ret_values);

    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    solver->statistics.reorder_time = ret_values[0];
    solver->statistics.factor_time  = ret_values[1];
    solver->statistics.solve_time   = ret_values[2];
    solver->statistics.flops        = ret_values[3];

cleanup:

    return ret == MKL_DSS_SUCCESS;
}
