#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <stddef.h>

typedef enum _SparseMatrixStorage
{
    SPARSE_MATRIX_STORAGE_CSR
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

typedef struct _SparseSolver SparseSolver;

typedef enum _SparseSolverMethod
{
    SPARSE_SOLVER_METHOD_SUPERLU,
    SPARSE_SOLVER_METHOD_SUPERLU_MT,
    SPARSE_SOLVER_METHOD_UMFPACK,
    SPARSE_SOLVER_METHOD_MKL_DSS,
    SPARSE_SOLVER_METHOD_MKL_CLUSTER,
    SPARSE_SOLVER_METHOD_CUSOLVER
} SparseSolverMethod;

typedef struct _SparseSolverStatistics
{
    double reorder_time;
    double factor_time;
    double solve_time;
    double flops;
} SparseSolverStatistics;

SparseSolver *         sparse_solver_new            (SparseSolverMethod solver_method);
int                    sparse_solver_solve          (SparseSolver *solver, SparseMatrix *A, void *x, void *b);
SparseSolverStatistics sparse_solver_get_statistics (SparseSolver *solver);
void                   sparse_solver_set_verbose    (SparseSolver *solver, int verbose);
void                   sparse_solver_free           (SparseSolver *solver);

#endif /* __SPARSE_H__ */
