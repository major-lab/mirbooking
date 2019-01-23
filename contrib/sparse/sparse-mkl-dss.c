#include "sparse.h"
#include "sparse-private.h"

#include <assert.h>
#include <mkl_dss.h>
#include <stdio.h>
#include <stdlib.h>

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
sparse_mkl_dss_init (SparseSolver *solver)
{

}

void
sparse_mkl_dss_clear (SparseSolver *solver)
{

}

int
sparse_mkl_dss_solve (SparseSolver *solver, SparseMatrix *A, void *x, void *b)
{
    _MKL_DSS_HANDLE_t handle;
    int ret = 0;

    if (A->shape[0] == 0 || A->shape[1] == 0)
    {
        return 1;
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

    ret = dss_create (handle, create_opt);

    if (ret != MKL_DSS_SUCCESS)
    {
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

    MKL_INT *rowptr, *colind, *rowperm;
    if (sizeof (MKL_INT) == sizeof (size_t))
    {
        rowptr  = (MKL_INT*) A->s.csr.rowptr;
        colind  = (MKL_INT*) A->s.csr.colind;
    }
    else
    {
        rowptr  = malloc ((A->shape[0] + 1) * sizeof (MKL_INT));
        colind  = malloc (A->s.csr.nnz * sizeof (MKL_INT));
        memcpy_loop (rowptr, A->s.csr.rowptr, A->shape[0] + 1);
        memcpy_loop (colind, A->s.csr.colind, A->s.csr.nnz);
    }

    ret = dss_define_structure (handle,
                                define_opt,
                                rowptr,
                                m,
                                n,
                                colind,
                                nnz);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    if (A->solver_storage == NULL)
    {
        A->solver_storage = calloc (A->shape[0], sizeof (MKL_INT));
    }

    rowperm = A->solver_storage;

    MKL_INT reorder_opt = (rowperm[0] != rowperm[A->shape[0] - 1]) ? MKL_DSS_MY_ORDER : MKL_DSS_GET_ORDER;
    ret = dss_reorder (handle, reorder_opt, rowperm);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
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

    ret = dss_factor_real (handle, factor_opt, A->data);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    MKL_INT solve_opt = 0;
    MKL_INT nrhs = 1;
    ret = dss_solve_real (handle, solve_opt, b, nrhs, x);

    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    MKL_INT statistics_opt = 0;
    double ret_values[7];
    ret = dss_statistics (handle,
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
    if (sizeof (MKL_INT) != sizeof (size_t))
    {
        free (rowptr);
        free (colind);
    }

    MKL_INT delete_opt = 0;
    dss_delete (handle, delete_opt);

    return ret == MKL_DSS_SUCCESS;
}
