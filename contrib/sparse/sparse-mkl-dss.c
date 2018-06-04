#include "sparse.h"

#include <mkl_dss.h>

int
sparse_mkl_dss_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b)
{
    _MKL_DSS_HANDLE_t handle;
    int ret = 0;

    MKL_INT create_opt = MKL_DSS_ZERO_BASED_INDEXING;
    ret = dss_create (handle, create_opt);

    if (ret != MKL_DSS_SUCCESS)
    {
        return 0;
    }

    MKL_INT define_opt = MKL_DSS_NON_SYMMETRIC;
    MKL_INT m = A->shape[0];
    MKL_INT n = A->shape[1];
    MKL_INT nnz = A->s.csr.nnz;

    ret = dss_define_structure (handle,
                                define_opt,
                                (MKL_INT*) A->s.csr.rowptr,
                                m,
                                n,
                                (MKL_INT*) A->s.csr.colind,
                                nnz);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    MKL_INT reorder_opt = (A->rowperm[0] != A->rowperm[A->shape[0] - 1]) ? MKL_DSS_MY_ORDER : MKL_DSS_GET_ORDER;
    ret = dss_reorder (handle, reorder_opt, (MKL_INT*) A->rowperm);
    if (ret != MKL_DSS_SUCCESS)
    {
        goto cleanup;
    }

    MKL_INT factor_opt = MKL_DSS_INDEFINITE;
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

    MKL_INT delete_opt = 0;
cleanup:
    dss_delete (handle, delete_opt);

    return ret == MKL_DSS_SUCCESS;
}
