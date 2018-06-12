#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <stddef.h>

typedef enum _SparseMatrixStorage
{
    SPARSE_MATRIX_STORAGE_CSR
} SparseMatrixStorage;

typedef enum _SparseMatrixHint
{
    SPARSE_MATRIX_HINT_SYMMETRIC         = 1 << 0, /* only the upper triangle is necessary */
    SPARSE_MATRIX_HINT_POSITIVE_DEFINITE = 1 << 1
} SparseMatrixHint;

typedef struct _SparseMatrix
{
    SparseMatrixStorage storage;
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
    double *data;
    /* optimal row and col permutations */
    size_t *colperm;
    size_t *rowperm;
    /* hints for solvers */
    SparseMatrixHint hints;
} SparseMatrix;

void   sparse_matrix_init      (SparseMatrix *matrix, SparseMatrixStorage storage, size_t shape[2], size_t nnz);
void   sparse_matrix_clear     (SparseMatrix *matrix);
double sparse_matrix_get_value (SparseMatrix *matrix, size_t i, size_t j);
void   sparse_matrix_set_value (SparseMatrix *matrix, size_t i, size_t j, double v);

typedef enum _SparseSolverMethod
{
    SPARSE_SOLVER_METHOD_SUPERLU,
    SPARSE_SOLVER_METHOD_SUPERLU_MT,
    SPARSE_SOLVER_METHOD_UMFPACK,
    SPARSE_SOLVER_METHOD_MKL_DSS,
    SPARSE_SOLVER_METHOD_MKL_CLUSTER
} SparseSolverMethod;

typedef struct _SparseSolver SparseSolver;

SparseSolver * sparse_solver_new         (SparseSolverMethod solver_method);
void           sparse_solver_set_verbose (SparseSolver *solver, int verbose);
int            sparse_solver_solve       (SparseSolver *solver, SparseMatrix *A, double *x, double *b);
void           sparse_solver_free        (SparseSolver *solver);

#endif /* __SPARSE_H__ */
