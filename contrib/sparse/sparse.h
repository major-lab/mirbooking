#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <stddef.h>

typedef enum _SparseMatrixStorage
{
    SPARSE_MATRIX_STORAGE_CSR
} SparseMatrixStorage;

typedef struct _SparseMatrix
{
    int     storage;
    size_t  shape[2];
    union
    {
        struct
        {
            size_t  nnz;
            int    *colind;
            int    *rowptr;
        } csr;
    } s;
    double *data;
} SparseMatrix;

void   sparse_matrix_init      (SparseMatrix *matrix, SparseMatrixStorage storage, size_t shape[2], size_t nnz);
void   sparse_matrix_clear     (SparseMatrix *matrix);
double sparse_matrix_get_value (SparseMatrix *matrix, int i, int j);
void   sparse_matrix_set_value (SparseMatrix *matrix, int i, int j, double v);

typedef enum _SparseSolverMethod
{
    SPARSE_SOLVER_METHOD_SUPERLU,
    SPARSE_SOLVER_METHOD_SUPERLU_MT
} SparseSolverMethod;

typedef struct _SparseSolver SparseSolver;

SparseSolver * sparse_solver_new   (SparseSolverMethod solver_method);
int            sparse_solver_solve (SparseSolver *solver, SparseMatrix *A, double *x, double *b);
void           sparse_solver_free  (SparseSolver *solver);

#endif /* __SPARSE_H__ */
