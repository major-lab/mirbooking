#include "sparse.h"
#include "sparse-private.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define DECLARE_SOLVER(solver)                  \
extern void                                     \
sparse_##solver##_init  (SparseSolver *solver); \
extern void                                     \
sparse_##solver##_clear (SparseSolver *solver); \
extern int                                      \
sparse_##solver##_solve (SparseSolver *solver,  \
                         SparseMatrix *A,       \
                         void         *x,       \
                         void         *b);

DECLARE_SOLVER(superlu)
DECLARE_SOLVER(superlu_mt)
DECLARE_SOLVER(umfpack)
DECLARE_SOLVER(mkl_dss)
DECLARE_SOLVER(mkl_cluster)
DECLARE_SOLVER(cusolver)

static size_t
_size_for_type (SparseMatrixType type)
{
    switch (type)
    {
        case SPARSE_MATRIX_TYPE_FLOAT:
            return sizeof (float);
        case SPARSE_MATRIX_TYPE_DOUBLE:
            return sizeof (double);
        default:
            assert (0);
    }
}

void
sparse_matrix_init (SparseMatrix        *matrix,
                    SparseMatrixStorage  storage,
                    SparseMatrixType     type,
                    size_t               shape[2],
                    size_t               nnz)
{
    matrix->storage      = storage;
    matrix->type         = type;
    memcpy (matrix->shape, shape, 2 * sizeof (size_t));

    /* CSR */
    assert (matrix->storage == SPARSE_MATRIX_STORAGE_CSR);
    matrix->s.csr.nnz    = nnz;
    matrix->s.csr.colind = calloc (nnz, sizeof (size_t));
    matrix->s.csr.rowptr = calloc (shape[0] + 1, sizeof (size_t));

    /* default is zero */
    memset (&matrix->default_data, 0, sizeof (matrix->default_data));

    /* storage */
    matrix->data = calloc (nnz, _size_for_type (matrix->type));

    matrix->colperm = calloc (shape[1], sizeof (size_t));
    matrix->rowperm = calloc (shape[0], sizeof (size_t));
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
 * sparse_matrix_reserve_range:
 *
 * Pre-allocate a range of non-zero entries along a row.
 */
void
sparse_matrix_reserve_range (SparseMatrix *matrix, size_t i, size_t *colind, size_t n)
{
    assert ((i == matrix->shape[0] - 1) || matrix->s.csr.rowptr[i+1] == 0);
    assert (matrix->s.csr.rowptr[i] + n <= matrix->s.csr.nnz);

    matrix->s.csr.rowptr[i + 1] = matrix->s.csr.rowptr[i] + n;
    memcpy (matrix->s.csr.colind + matrix->s.csr.rowptr[i], colind, n * sizeof (size_t));
}

static int
cmp_size_t (const void *a, const void *b)
{
    return *(const size_t*)a - *(const size_t*)b;
}

static ssize_t
_sparse_matrix_get_index (SparseMatrix *matrix, size_t i, size_t j)
{
    size_t *k = NULL;

    assert (i < matrix->shape[0]);
    assert (j < matrix->shape[1]);

    k = bsearch (&j, matrix->s.csr.colind + matrix->s.csr.rowptr[i],
                 matrix->s.csr.rowptr[i+1] - matrix->s.csr.rowptr[i],
                 sizeof (size_t),
                 cmp_size_t);

    return k ? (k - matrix->s.csr.colind) : -1;
}

/**
 * sparse_matrix_get_float:
 */
float
sparse_matrix_get_float (SparseMatrix *matrix, size_t i, size_t j)
{
    ssize_t k = _sparse_matrix_get_index (matrix, i, j);

    if (k == -1)
    {
        return matrix->default_data.f;
    }

    return ((float*)matrix->data)[k];
}

/**
 * sparse_matrix_get_double:
 */
double
sparse_matrix_get_double (SparseMatrix *matrix, size_t i, size_t j)
{
    ssize_t k = _sparse_matrix_get_index (matrix, i, j);

    if (k == -1)
    {
        return matrix->default_data.d;
    }

    return ((double*)matrix->data)[k];
}


ssize_t
_sparse_matrix_reserve_index (SparseMatrix *matrix, size_t i, size_t j)
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
            memmove (matrix->data + (k + 1) * _size_for_type (matrix->type),
                     matrix->data + k * _size_for_type (matrix->type),
                     (matrix->s.csr.rowptr[i+1] - 1 - k) * _size_for_type (matrix->type));
        }

        matrix->s.csr.colind[k] = j;
    }

    return k;
}

/**
 * sparse_matrix_set_value:
 *
 * Set a value at given (i, j) coordinate.
 */
void
sparse_matrix_set_float (SparseMatrix *matrix, size_t i, size_t j, float v)
{
    ssize_t k = _sparse_matrix_reserve_index (matrix, i, j);

    assert (matrix->type == SPARSE_MATRIX_TYPE_FLOAT);

    ((float*)matrix->data)[k] = v;
}

/**
 * sparse_matrix_set_value:
 *
 * Set a value at given (i, j) coordinate.
 */
void
sparse_matrix_set_double (SparseMatrix *matrix, size_t i, size_t j, double v)
{
    ssize_t k = _sparse_matrix_reserve_index (matrix, i, j);

    assert (matrix->type == SPARSE_MATRIX_TYPE_DOUBLE);

    ((double*)matrix->data)[k] = v;
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

    memset (&ret->statistics, 0, sizeof (SparseSolverStatistics));
    ret->verbose = 0;

    #define PREPARE_SOLVER(solver_uc,solver)  \
    case SPARSE_SOLVER_METHOD_##solver_uc:     \
        ret->init  = sparse_##solver##_init;  \
        ret->clear = sparse_##solver##_clear; \
        ret->solve = sparse_##solver##_solve; \
        break;

    switch (method)
    {
#if HAVE_SUPERLU
            PREPARE_SOLVER(SUPERLU,superlu)
#endif
#if HAVE_SUPERLU_MT
            PREPARE_SOLVER(SUPERLU_MT,superlu_mt)
#endif
#if HAVE_UMFPACK
            PREPARE_SOLVER(UMFPACK,umfpack);
#endif
#if HAVE_MKL_DSS
            PREPARE_SOLVER(MKL_DSS,mkl_dss);
#endif
#if HAVE_MKL_CLUSTER
            PREPARE_SOLVER(MKL_CLUSTER,mkl_cluster);
#endif
#if HAVE_CUSOLVER
            PREPARE_SOLVER(CUSOLVER,cusolver)
#endif
        default:
            return NULL;
    }

    ret->init (ret);

    return ret;
}

/**
 * sparse_solver_set_verbose:
 * @verbose: A verbosity level where 0 is no message (the default) and >=1
 * print resource usage statistics.
 *
 * Set the verbosity level of the solver.
 */
void
sparse_solver_set_verbose (SparseSolver *solver, int verbose)
{
    solver->verbose = verbose;
}

/**
 * sparse_solver_solve:
 * @solver: A #SparseSolver
 * @A: A square sparse matrix
 * @x: A vector of shape A->shape[1] whose element type match A->type
 * @b: A vector of size A->shape[1] whose element type match A->type
 *
 * Solve linear system of the form Ax=b where @A is sparse.
 */
int
sparse_solver_solve (SparseSolver *solver, SparseMatrix *A, void *x, void *b)
{
    assert (A->shape[0] == A->shape[1]);

    return solver->solve (solver, A, x, b);
}

SparseSolverStatistics
sparse_solver_get_statistics (SparseSolver *solver)
{
    return solver->statistics;
}


/**
 * sparse_solver_free:
 */
void
sparse_solver_free (SparseSolver *solver)
{
    solver->clear (solver);
    free (solver);
}
