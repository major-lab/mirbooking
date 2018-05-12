#include "sparse-superlu-private.h"

#include <SuperLU/slu_ddefs.h>

int
superlu_solve (SparseSolver *solver,
               SparseMatrix *A,
               double       *x,
               double       *b)
{
    superlu_options_t options = {0};
    SuperMatrix AA, L, U, BB;
    SuperLUStat_t stat;
    int info;

    set_default_options (&options);

    StatInit (&stat);

    int *perm_r = malloc (A->shape[0] * sizeof (int));
    int *perm_c = malloc (A->shape[1] * sizeof (int));

    dCreate_CompRow_Matrix (&AA,
                            A->shape[0],
                            A->shape[1],
                            A->s.csr.nnz,
                            A->data,
                            A->s.csr.colind,
                            A->s.csr.rowptr,
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
           perm_c,
           perm_r,
           &L,
           &U,
           &BB,
           &stat,
           &info);

    Destroy_CompRow_Matrix (&L);
    Destroy_CompRow_Matrix (&U);
    StatFree (&stat);
    free (perm_r);
    free (perm_c);

    return info == 0;
}
