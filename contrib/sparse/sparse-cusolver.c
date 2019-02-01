#include "sparse.h"
#include "sparse-private.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cusolverSp.h>
#include <cuda_device_runtime_api.h>
#include <cuda_runtime_api.h>
#include <limits.h>

void
sparse_cusolver_init (SparseSolver *solver)
{

}

void
sparse_cusolver_clear (SparseSolver *solver)
{

}

int
sparse_cusolver_solve (SparseSolver *solver,
                       SparseMatrix *A,
                       void         *x,
                       void         *b)
{
    cusolverSpHandle_t handle;
    cusparseMatDescr_t descrA;
    cusolverStatus_t status;
    cudaError_t err;

    assert (A->shape[0] < INT_MAX);
    assert (A->s.csr.nnz < INT_MAX);

    status = cusolverSpCreate (&handle);

    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        goto cleanup;
    }

    status = cusparseCreateMatDescr (&descrA);

    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        goto cleanup_handle;
    }

    int *rowptr = malloc ((A->shape[0]+1) * sizeof (int));
    int *colind = malloc (A->s.csr.nnz * sizeof (int));

    memcpy_loop (rowptr, A->s.csr.rowptr, A->shape[0] + 1);
    memcpy_loop (colind, A->s.csr.colind, A->s.csr.nnz);

    int singularity;
    if (A->type == SPARSE_MATRIX_TYPE_FLOAT)
    {
        status = cusolverSpScsrlsvluHost (handle,
                                          A->shape[0],
                                          A->s.csr.nnz,
                                          descrA,
                                          A->data,
                                          rowptr,
                                          colind,
                                          b,
                                          0,
                                          1,
                                          x,
                                          &singularity);
    }
    else
    {
        status = cusolverSpDcsrlsvluHost (handle,
                                          A->shape[0],
                                          A->s.csr.nnz,
                                          descrA,
                                          A->data,
                                          rowptr,
                                          colind,
                                          b,
                                          0,
                                          1,
                                          x,
                                          &singularity);
    }

    free (rowptr);
    free (colind);

cleanup_handle:
    cusolverSpDestroy (handle);

cleanup:

    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        err = cudaGetLastError ();
        fprintf (stderr, "%s\n", cudaGetErrorString (err));
    }

    return status == CUSOLVER_STATUS_SUCCESS;
}
