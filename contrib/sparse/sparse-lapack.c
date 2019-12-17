#include "sparse.h"
#include "sparse-private.h"

#include <stdlib.h>
#include <string.h>
#if HAVE_MKL_LAPACK
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#include <assert.h>

void
sparse_lapack_init (SparseSolver *solver)
{

}

void
sparse_lapack_clear (SparseSolver *solver)
{

}

int
sparse_lapack_solve (SparseSolver *solver,
                      SparseMatrix *A,
                      void         *x,
                      const void   *b)
{
    lapack_int ret;

    assert (A->storage == SPARSE_MATRIX_STORAGE_DENSE);
    assert (A->shape[0] == A->shape[1]);

    lapack_int *pivot = malloc (A->shape[0] * sizeof (lapack_int));

    if (A->type == SPARSE_MATRIX_TYPE_DOUBLE)
    {
        memcpy (x, b, A->shape[1] * sizeof (double));
        ret = LAPACKE_dgesv (LAPACK_ROW_MAJOR,
                             A->shape[0],
                             1,
                             A->data,
                             A->shape[1],
                             pivot,
                             (double*) x,
                             1);
    }
    else if (A->type == SPARSE_MATRIX_TYPE_FLOAT)
    {
        memcpy (x, b, A->shape[1] * sizeof (float));
        ret = LAPACKE_sgesv (LAPACK_ROW_MAJOR,
                             A->shape[0],
                             1,
                             A->data,
                             A->shape[1],
                             pivot,
                             (float*) x,
                             1);
    }
    else
    {
        assert (0);
    }

    free (pivot);

    return ret == 0;
}

