#include "sparse.h"
#include "sparse-superlu-private.h"
#include "sparse-superlu-mt-private.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct _SparseSolver
{
    int (*solve) (SparseSolver *solver, SparseMatrix *A, double *x, double *b);
};

void
sparse_matrix_init (SparseMatrix        *matrix,
                    SparseMatrixStorage  storage,
                    size_t               shape[2],
                    size_t               nnz)
{
    matrix->storage      = storage;
    memcpy (matrix->shape, shape, 2 * sizeof (size_t));

    /* CSR */
    assert (matrix->storage == SPARSE_MATRIX_STORAGE_CSR);
    matrix->s.csr.nnz    = nnz;
    matrix->s.csr.colind = calloc (nnz, sizeof (size_t));
    matrix->s.csr.rowptr = calloc (shape[0] + 1, sizeof (size_t));

    matrix->data         = calloc (nnz, sizeof (double));
}

void
sparse_matrix_clear (SparseMatrix *matrix)
{
    free (matrix->s.csr.colind);
    free (matrix->data);
    free (matrix->s.csr.rowptr);
}

/**
 * sparse_matrix_get_value:
 */
double
sparse_matrix_get_value (SparseMatrix *matrix, int i, int j)
{
    int k = -1;

    assert (i < matrix->shape[0]);
    assert (j < matrix->shape[1]);

    int z;
    for (z = matrix->s.csr.rowptr[i]; z < matrix->s.csr.rowptr[i+1]; z++)
    {
        if (matrix->s.csr.colind[z] == j)
        {
            k = z;
            break;
        }
    }

    if (k == -1)
    {
        return 0;
    }

    return matrix->data[k];
}

/**
 * sparse_matrix_set_value:
 *
 * Set a value at given (i, j) coordinate.
 */
inline void
sparse_matrix_set_value (SparseMatrix *matrix, int i, int j, double v)
{
    assert (matrix->s.csr.rowptr);
    assert (i < matrix->shape[0]);
    assert (j < matrix->shape[1]);

    int k = -1;
    int z;

    for (z = matrix->s.csr.rowptr[i]; z < matrix->s.csr.rowptr[i+1]; z++)
    {
        if (matrix->s.csr.colind[z] == j)
        {
            k = z;
            break;
        }
    }

    if (k == -1)
    {
        // ensure we're not overwriting the next row unless this is already the
        // last one
        assert ((i == matrix->shape[0] - 1) || matrix->s.csr.rowptr[i+2] == 0);

        // ensure we have some memory left
        assert (matrix->s.csr.rowptr[i+1] < matrix->s.csr.nnz);

        // unused row
        if (matrix->s.csr.rowptr[i + 1] == 0)
        {
            matrix->s.csr.rowptr[i + 1] = matrix->s.csr.rowptr[i];
        }

        // allocate a slot
        k = matrix->s.csr.rowptr[i+1]++;
        matrix->s.csr.colind[k] = j;
    }

    matrix->data[k] = v;
}

SparseSolver *
sparse_solver_new (SparseSolverMethod method)
{
    SparseSolver *ret = malloc (sizeof (SparseSolver));
    switch (method)
    {
        case SPARSE_SOLVER_METHOD_SUPERLU:
            ret->solve = superlu_solve;
            break;
        case SPARSE_SOLVER_METHOD_SUPERLU_MT:
            ret->solve = superlu_mt_solve;
            break;
        default:
            assert (0);
    }
    return ret;
}

/**
 * sparse_solver_solve:
 * @solver: A #SparseSolver
 * @A: A square sparse matrix
 * @x: A vector of shape A->shape[1]
 * @b: A vector of size A->shape[1]
 *
 * Solve linear system of the form Ax=b where @A is sparse.
 */
int
sparse_solver_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b)
{
    assert (A->shape[0] == A->shape[1]);
    return solver->solve (solver, A, x, b);
}

/**
 * sparse_solver_free:
 */
void
sparse_solver_free (SparseSolver *solver)
{
    free (solver);
}
