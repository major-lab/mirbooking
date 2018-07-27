#include "sparse.h"

#include <SuperLU/slu_ddefs.h>
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
                      void         *b)
{
    superlu_options_t options = {0};
    SuperMatrix AA, L, U, BB;
    SuperLUStat_t stat;
    int_t info;

    assert (A->s.csr.nnz <= INT_MAX);

    set_default_options (&options);

    StatInit (&stat);

    if (A->colperm[0] != A->colperm[A->shape[0] - 1])
    {
        options.ColPerm = MY_PERMC;
    }

    if (A->rowperm[0] != A->rowperm[A->shape[0] - 1])
    {
        options.RowPerm = MY_PERMR;
    }

    assert (A->shape[0] < INT_MAX);
    assert (A->s.csr.nnz < INT_MAX);

    int *colind = malloc (A->s.csr.nnz * sizeof (int));
    int *rowptr = malloc ((A->shape[0] + 1) * sizeof (int));
    int *colperm = malloc (A->shape[1] * sizeof (int));
    int *rowperm = malloc (A->shape[0] * sizeof (int));

    int i;
    for (i = 0; i < A->s.csr.nnz; i++)
    {
        colind[i] = A->s.csr.colind[i];
    }

    for (i = 0; i < A->shape[0] + 1; i++)
    {
        rowptr[i] = A->s.csr.rowptr[i];
    }

    for (i = 0; i < A->shape[0]; i++)
    {
        rowperm[i] = A->rowperm[i];
    }

    for (i = 0; i < A->shape[1]; i++)
    {
        colperm[i] = A->colperm[i];
    }

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

    memcpy (x, b, A->shape[0] * sizeof (double));

    dCreate_Dense_Matrix (&BB,
                          A->shape[1],
                          1,
                          x,
                          A->shape[1],
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

    if (options.RowPerm != MY_PERMR)
    {
        for (i = 0; i < A->shape[0]; i++)
        {
            A->rowperm[i] = rowperm[i];
        }
    }

    if (options.ColPerm != MY_PERMC)
    {
        for (i = 0; i < A->shape[1]; i++)
        {
            A->colperm[i] = colperm[i];
        }
    }

    Destroy_CompRow_Matrix (&L);
    Destroy_CompRow_Matrix (&U);
    StatFree (&stat);
    free (colind);
    free (rowptr);
    free (colperm);
    free (rowperm);

    return info == 0;
}
