#include "sparse.h"

#include <SuperLUMT/slu_mt_ddefs.h>
#include <assert.h>
#include <limits.h>

/* this API is not exposed in SuperLUMT */
extern void
dCreate_CompRow_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, double *nzval,
                                  int_t *colind, int_t *rowptr,
                                  Stype_t stype, Dtype_t dtype, Mtype_t mtype);

void
sparse_superlu_mt_init (SparseSolver *solver)
{

}

void
sparse_superlu_mt_clear (SparseSolver *solver)
{

}

int
sparse_superlu_mt_solve (SparseSolver *solver,
                         SparseMatrix *A,
                         void         *x,
                         void         *b)
{
    SuperMatrix AA, L, U, BB;
    int_t info;

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

    pdgssv (8,
            &AA,
            colperm,
            rowperm,
            &L,
            &U,
            &BB,
            &info);

    Destroy_SuperMatrix_Store (&L);
    Destroy_SuperMatrix_Store (&U);

    return info == 0;
}


