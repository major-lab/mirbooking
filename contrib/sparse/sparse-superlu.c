#include "sparse.h"
#include "sparse-private.h"

#if HAVE_SUPERLU_LOWERCASE_INCDIR
#include <superlu/slu_ddefs.h>
#else
#include <SuperLU/slu_ddefs.h>
#endif
#include <assert.h>
#include <limits.h>

void
sparse_superlu_init (SparseSolver *solver)
{

}

void
sparse_superlu_clear (SparseSolver *solver)
{

}

int
sparse_superlu_solve (SparseSolver *solver,
                      SparseMatrix *A,
                      void         *x,
                      const void   *b)
{
    superlu_options_t options = {0};
    SuperMatrix AA, L, U, BB;
    SuperLUStat_t stat;
    int_t info;

    assert (A->storage == SPARSE_MATRIX_STORAGE_CSR);
    assert (A->type == SPARSE_MATRIX_TYPE_DOUBLE);
    assert (A->shape[0] < INT_MAX);
    assert (A->shape[1] < INT_MAX);
    assert (A->s.csr.nnz < INT_MAX);

    set_default_options (&options);

    StatInit (&stat);

    if (A->solver_storage_owner != solver)
    {
        sparse_matrix_set_solver_storage (A,
                                          calloc (A->shape[0] + A->shape[1], sizeof (int)),
                                          free,
                                          solver);
    }

    int *rowperm = A->solver_storage;
    int *colperm = rowperm + A->shape[0];

    if (A->shape[0] && rowperm[0] != rowperm[A->shape[0] - 1])
    {
        options.RowPerm = MY_PERMR;
    }

    if (A->shape[1] && colperm[0] != colperm[A->shape[1] - 1])
    {
        options.ColPerm = MY_PERMC;
    }

    int *colind = intMalloc (A->s.csr.nnz);
    int *rowptr = intMalloc ((A->shape[0] + 1));

    memcpy_loop (colind, A->s.csr.colind, A->s.csr.nnz);
    memcpy_loop (rowptr, A->s.csr.rowptr, A->shape[0] + 1);
    memcpy (x, b, A->shape[0] * sizeof (double));

    dCreate_CompRow_Matrix (&AA,
                            A->shape[0],
                            A->shape[1],
                            A->s.csr.nnz,
                            A->data,
                            colind,
                            rowptr,
                            SLU_NR,
                            SLU_D,
                            SLU_GE);

    dCreate_Dense_Matrix (&BB,
                          A->shape[0],
                          1,
                          x,
                          A->shape[0],
                          SLU_DN,
                          SLU_D,
                          SLU_GE);

    dgssv (&options,
           &AA,
           colperm,
           rowperm,
           &L,
           &U,
           &BB,
           &stat,
           &info);

    if (solver->verbose)
    {
        StatPrint (&stat);
    }

    SUPERLU_FREE (colind);
    SUPERLU_FREE (rowptr);

    Destroy_SuperMatrix_Store (&BB);
    Destroy_SuperMatrix_Store (&AA);
    Destroy_SuperNode_Matrix (&L);
    Destroy_CompRow_Matrix (&U);
    StatFree (&stat);

    return info == 0;
}
