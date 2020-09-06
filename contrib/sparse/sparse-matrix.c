#include "sparse-matrix.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

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
    if (matrix->storage == SPARSE_MATRIX_STORAGE_CSR)
    {
        matrix->s.csr.nnz    = nnz;
        matrix->s.csr.colind = calloc (nnz, sizeof (size_t));
        matrix->s.csr.rowptr = calloc (shape[0] + 1, sizeof (size_t));
    }
    else if (matrix->storage == SPARSE_MATRIX_STORAGE_DENSE)
    {
        // no structure to initialize
    }
    else
    {
        assert (0);
    }

    /* default is zero */
    memset (&matrix->default_data, 0, sizeof (matrix->default_data));

    /* storage */
    matrix->data = calloc (nnz, _size_for_type (matrix->type));

    /* solver-specific storage */
    matrix->solver_storage_owner = NULL;
    matrix->solver_storage = NULL;
    matrix->solver_storage_destroy = NULL;
}

void
sparse_matrix_clear (SparseMatrix *matrix)
{
    if (matrix->storage == SPARSE_MATRIX_STORAGE_CSR)
    {
        free (matrix->s.csr.colind);
        free (matrix->s.csr.rowptr);
    }
    free (matrix->data);
    if (matrix->solver_storage_destroy)
    {
        matrix->solver_storage_destroy (matrix->solver_storage);
    }
}

/**
 * sparse_matrix_reserve_range:
 *
 * Pre-allocate a range of non-zero entries along a row.
 */
void
sparse_matrix_reserve_range (SparseMatrix *matrix, size_t i, size_t *colind, size_t n)
{
    if (matrix->storage == SPARSE_MATRIX_STORAGE_DENSE)
        return;

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

    if (matrix->storage == SPARSE_MATRIX_STORAGE_DENSE)
        return matrix->shape[0] * i + j;

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
    assert (i < matrix->shape[0]);
    assert (j < matrix->shape[1]);

    if (matrix->storage == SPARSE_MATRIX_STORAGE_DENSE)
        return matrix->shape[0] * i + j;

    assert (matrix->s.csr.rowptr);

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
            memmove ((char*) matrix->data + (k + 1) * _size_for_type (matrix->type),
                     (char*) matrix->data + k * _size_for_type (matrix->type),
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
