#include "sparse.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

extern int
sparse_superlu_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b);

extern int
sparse_superlu_mt_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b);

extern int
sparse_umfpack_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b);

extern int
sparse_mkl_dss_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b);

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

    matrix->colperm = calloc (shape[1], sizeof (size_t));
    matrix->rowperm = calloc (shape[0], sizeof (size_t));;
}

void
sparse_matrix_clear (SparseMatrix *matrix)
{
    free (matrix->s.csr.colind);
    free (matrix->data);
    free (matrix->s.csr.rowptr);
    free (matrix->colperm);
    free (matrix->rowperm);
}

/**
 * sparse_matrix_get_value:
 */
double
sparse_matrix_get_value (SparseMatrix *matrix, size_t i, size_t j)
{
    ssize_t k = -1;

    assert (i < matrix->shape[0]);
    assert (j < matrix->shape[1]);

    // FIXME: use bsearch
    size_t z;
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
void
sparse_matrix_set_value (SparseMatrix *matrix, size_t i, size_t j, double v)
{
    assert (matrix->s.csr.rowptr);
    assert (i < matrix->shape[0]);
    assert (j < matrix->shape[1]);

    ssize_t k = -1;
    size_t z;

    // FIXME: use bsearch
    for (z = matrix->s.csr.rowptr[i]; z < matrix->s.csr.rowptr[i+1]; z++)
    {
        if (matrix->s.csr.colind[z] >= j)
        {
            k = z;
            break;
        }
    }

    if (k == -1 || matrix->s.csr.colind[k] != j)
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
        ++matrix->s.csr.rowptr[i+1];

        if (k == -1)
        {
            // append
            k = matrix->s.csr.rowptr[i+1] - 1;
        }
        else
        {
            // insertion
            memmove (&matrix->s.csr.colind[k + 1], &matrix->s.csr.colind[k], (matrix->s.csr.rowptr[i+1] - 1 - k) * sizeof (size_t));
            memmove (&matrix->data[k + 1], &matrix->data[k], (matrix->s.csr.rowptr[i+1] - 1 - k) * sizeof (double));
        }

        matrix->s.csr.colind[k] = j;
    }

    matrix->data[k] = v;
}

/**
 * sparse_solver_new:
 *
 * Returns: A #SparseSolver, or %NULL if the solver is not available.
 */
SparseSolver *
sparse_solver_new (SparseSolverMethod method)
{
    SparseSolver *ret = malloc (sizeof (SparseSolver));
    switch (method)
    {
        case SPARSE_SOLVER_METHOD_SUPERLU:
#if HAVE_SUPERLU
            ret->solve = sparse_superlu_solve;
#else
            return NULL;
#endif
            break;
        case SPARSE_SOLVER_METHOD_SUPERLU_MT:
#if HAVE_SUPERLU_MT
            ret->solve = sparse_superlu_mt_solve;
#else
            return NULL;
#endif
            break;
        case SPARSE_SOLVER_METHOD_UMFPACK:
#if HAVE_UMFPACK
            ret->solve = sparse_umfpack_solve;
#else
            return NULL;
#endif
            break;
        case SPARSE_SOLVER_METHOD_MKL_DSS:
#if HAVE_MKL
            ret->solve = sparse_mkl_dss_solve;
#else
            return NULL;
#endif
            break;
        case SPARSE_SOLVER_METHOD_MKL_CLUSTER:
#if HAVE_MKL & HAVE_MPI
            ret->solve = sparse_mkl_cluster_solve;
#else
            return NULL;
#endif
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
