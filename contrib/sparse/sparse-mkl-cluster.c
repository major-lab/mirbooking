#include "sparse.h"
#include "sparse-private.h"

#include <mkl_dss.h>
#include <mkl_cluster_sparse_solver.h>
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>

#ifdef MKL_ILP64
#define cluster_sparse_solver cluster_sparse_solver_64
#endif

void
sparse_mkl_cluster_init (SparseSolver *solver)
{

}

void
sparse_mkl_cluster_clear (SparseSolver *solver)
{

}

int
sparse_mkl_cluster_solve (SparseSolver *solver,
                          SparseMatrix *A,
                          void         *x,
                          void         *b)
{
    void *handle[64] = {0};
    int ret = 0;
    double reorder_begin, reorder_end, factor_begin, factor_end, solve_begin, solve_end;

    int flag;
    MPI_Initialized (&flag);
    assert (flag); /* make sure that MPI is already initialized */

    MPI_Fint comm = MPI_Comm_c2f (MPI_COMM_WORLD);
    int rank;
    assert (MPI_Comm_rank (MPI_COMM_WORLD, &rank) == MPI_SUCCESS);

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 11;
    MKL_INT phase;
    MKL_INT n = A != NULL ? A->shape[0] : 0;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT iparm[64] = {0};
    MKL_INT error = 0;

    // set defaults
    iparm[0]  =  1;
    iparm[1]  =  2;
    iparm[9]  =  13;
    iparm[10] =  1;
    iparm[12] =  1;
    iparm[17] = -1;
    iparm[20] =  1;

    // extra flags
    iparm[26] =  1; // additional checks for CSR
    iparm[34] =  1; // zero-based indexing

    assert (rank != 0 || (A != NULL && b != NULL));
    assert (x != NULL);

    MKL_INT *rowptr = NULL, *colind = NULL, *rowperm = NULL;
    double *data = NULL;
    if (A != NULL)
    {
        assert (A->shape[0] == A->shape[1]);

        // empty matrices cannot work in this setup
        assert (A->shape[0] > 0);
        assert (A->shape[1] > 0);

        if (sizeof (MKL_INT) == sizeof (size_t))
        {
            rowptr  = (MKL_INT*) A->s.csr.rowptr;
            colind  = (MKL_INT*) A->s.csr.colind;
        }
        else
        {
            rowptr = malloc ((A->shape[0] + 1)  * sizeof (MKL_INT));
            colind = malloc (A->s.csr.nnz * sizeof (MKL_INT));
            memcpy_loop (rowptr, A->s.csr.rowptr, A->shape[0] + 1);
            memcpy_loop (colind, A->s.csr.colind, A->s.csr.nnz);
        }

        data = A->data;
    }

    // TODO
    reorder_begin = 0;
    reorder_end = 0;
    factor_begin = 0;
    factor_end = 0;

    reorder_begin = MPI_Wtime ();

    phase = 11;
    cluster_sparse_solver (&handle,
                           &maxfct,
                           &mnum,
                           &mtype,
                           &phase,
                           &n,
                           data,
                           rowptr,
                           colind,
                           rowperm,
                           &nrhs,
                           iparm,
                           &msglvl,
                           b,
                           x,
                           &comm,
                           &error);

    reorder_end = MPI_Wtime ();

    ret = (error == 0);

    if (error != 0)
    {
        goto cleanup;
    }

    factor_begin = MPI_Wtime ();

    phase = 22;
    cluster_sparse_solver (&handle,
                           &maxfct,
                           &mnum,
                           &mtype,
                           &phase,
                           &n,
                           data,
                           rowptr,
                           colind,
                           rowperm,
                           &nrhs,
                           iparm,
                           &msglvl,
                           b,
                           x,
                           &comm,
                           &error);

    factor_end = MPI_Wtime ();

    ret = (error == 0);

    if (error != 0)
    {
        goto cleanup;
    }

    solve_begin = MPI_Wtime ();

    phase = 33;
    cluster_sparse_solver (&handle,
                           &maxfct,
                           &mnum,
                           &mtype,
                           &phase,
                           &n,
                           data,
                           rowptr,
                           colind,
                           rowperm,
                           &nrhs,
                           iparm,
                           &msglvl,
                           b,
                           x,
                           &comm,
                           &error);

    solve_end = MPI_Wtime ();

    ret = (error == 0);

    if (error != 0)
    {
        goto cleanup;
    }

    // broadcast the solution vector
    #ifdef MKL_ILP64
    MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    #else
    MPI_Bcast (&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    #endif
    MPI_Bcast (x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    solver->statistics.reorder_time = reorder_end - reorder_begin;
    solver->statistics.factor_time  = factor_end  - factor_begin;
    solver->statistics.solve_time   = solve_end   - solve_begin;

cleanup:
    /* cleanup */
    phase = -1;
    cluster_sparse_solver (handle,
                           &maxfct,
                           &mnum,
                           &mtype,
                           &phase,
                           &n,
                           data,
                           rowptr,
                           colind,
                           rowperm,
                           &nrhs,
                           iparm,
                           &msglvl,
                           b,
                           x,
                           &comm,
                           &error);

    assert (error == 0);

    if (A != NULL && sizeof (MKL_INT) != sizeof (size_t))
    {
        free (rowptr);
        free (colind);
    }

    return ret;
}
