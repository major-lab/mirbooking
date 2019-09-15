#include "sparse.h"
#include "sparse-private.h"

#include <SuperLUMT/slu_mt_ddefs.h>
#include <assert.h>
#include <limits.h>
#include <omp.h>

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
    int_t *rowperm;
    int_t info;

    assert (A->storage == SPARSE_MATRIX_STORAGE_CSR);
    assert (A->type == SPARSE_MATRIX_TYPE_DOUBLE);
    assert (A->shape[0] < INT_MAX);
    assert (A->shape[1] < INT_MAX);
    assert (A->s.csr.nnz < INT_MAX);

    int_t *colind = intMalloc (A->s.csr.nnz * sizeof (int_t));
    int_t *rowptr = intMalloc ((A->shape[0] + 1) * sizeof (int_t));

    memcpy_loop (colind, A->s.csr.colind, A->s.csr.nnz);
    memcpy_loop (rowptr, A->s.csr.rowptr, A->shape[0] + 1)
    memcpy (x, b, A->shape[0] * sizeof (double));

    dCreate_CompCol_Matrix (&AA,
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
                          A->shape[1],
                          1,
                          x,
                          A->shape[1],
                          SLU_DN,
                          SLU_D,
                          SLU_GE);

    if (A->solver_storage_owner != solver)
    {
        rowperm = intMalloc (A->shape[0]);

        int_t permc_spec = 1;
        get_perm_c (permc_spec, &AA, rowperm);

        sparse_matrix_set_solver_storage (A,
                                          rowperm,
                                          free,
                                          solver);
    }

    rowperm = A->solver_storage;

    pdgssv (omp_get_num_threads (),
            &AA,
            rowperm,
            rowperm,
            &L,
            &U,
            &BB,
            &info);

    SUPERLU_FREE (colind);
    SUPERLU_FREE (rowptr);
    Destroy_SuperMatrix_Store (&BB);
    Destroy_SuperMatrix_Store (&AA);
    Destroy_SuperNode_Matrix (&L);
    Destroy_CompCol_Matrix (&U);

    return info == 0;
}


