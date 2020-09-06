#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

#include <stddef.h>

typedef enum _SparseMatrixStorage
{
    SPARSE_MATRIX_STORAGE_CSR,
    SPARSE_MATRIX_STORAGE_DENSE
} SparseMatrixStorage;

typedef enum _SparseMatrixType
{
    SPARSE_MATRIX_TYPE_FLOAT,
    SPARSE_MATRIX_TYPE_DOUBLE
} SparseMatrixType;

typedef enum _SparseMatrixHint
{
    SPARSE_MATRIX_HINT_SYMMETRIC_STRUCTURE = 1 << 0, /* full matrix is necessary */
    SPARSE_MATRIX_HINT_SYMMETRIC           = 1 << 1, /* only the upper triangle is necessary */
    SPARSE_MATRIX_HINT_POSITIVE_DEFINITE   = 1 << 2
} SparseMatrixHint;

typedef struct _SparseMatrix
{
    SparseMatrixStorage storage;
    SparseMatrixType    type;
    /* hints for solvers */
    SparseMatrixHint hints;
    size_t shape[2];
    union
    {
        struct
        {
            size_t  nnz;
            size_t *colind;
            size_t *rowptr;
        } csr;
    } s;
    union
    {
        float f;
        double d;
    } default_data;
    void *data;
    /* solver-specific data */
    void *solver_storage_owner;
    void *solver_storage;
    void (*solver_storage_destroy) (void*);
} SparseMatrix;

void   sparse_matrix_init      (SparseMatrix *matrix, SparseMatrixStorage storage, SparseMatrixType type, size_t shape[2], size_t nnz);
void   sparse_matrix_clear     (SparseMatrix *matrix);

void sparse_matrix_reserve_range (SparseMatrix *matrix, size_t i, size_t *colind, size_t n);

float sparse_matrix_get_float (SparseMatrix *matrix, size_t i, size_t j);
void  sparse_matrix_set_float (SparseMatrix *matrix, size_t i, size_t j, float v);

double sparse_matrix_get_double (SparseMatrix *matrix, size_t i, size_t j);
void   sparse_matrix_set_double (SparseMatrix *matrix, size_t i, size_t j, double v);

double sparse_matrix_get_double (SparseMatrix *matrix, size_t i, size_t j);
void   sparse_matrix_set_double (SparseMatrix *matrix, size_t i, size_t j, double v);

#endif /* __SPARSE_MATRIX_H__ */
