#include "sparse-superlu-mt-private.h"

#include <SuperLUMT/slu_mt_ddefs.h>

int
superlu_mt_solve (SparseSolver *solver,
                  SparseMatrix *A,
                  double       *x,
                  double       *b)
{
    SuperMatrix AA, L, U, BB;
    int info;

    int *perm_r = malloc (A->shape[0] * sizeof (int));
    int *perm_c = malloc (A->shape[1] * sizeof (int));

    // FIXME: this is not working
    AA.Stype = SLU_NR;
    AA.Dtype = SLU_D;
    AA.Mtype = SLU_GE;
    AA.nrow = A->shape[0];
    AA.ncol = A->shape[1];
    AA.Store = malloc (sizeof (NRformat));
    NRformat *AA_store = AA.Store;
    AA_store->nnz = A->s.csr.nnz;
    AA_store->nzval = A->data;
    AA_store->colind = A->s.csr.colind;
    AA_store->rowptr = A->s.csr.rowptr;

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
            perm_c,
            perm_r,
            &L,
            &U,
            &BB,
            &info);

    free (AA.Store);
    Destroy_SuperMatrix_Store (&L);
    Destroy_SuperMatrix_Store (&U);
    free (perm_r);
    free (perm_c);

    return info == 0;
}


