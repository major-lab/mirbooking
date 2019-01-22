#include <stdlib.h>
#include <sparse.h>
#include <assert.h>
#include <mpi.h>

int main (int argc, char **argv)
{
    int provided;
    assert (MPI_Init_thread (&argc, &argv, MPI_THREAD_FUNNELED, &provided) == MPI_SUCCESS);
    assert (provided == MPI_THREAD_FUNNELED);

    int rank;
    assert (MPI_Comm_rank (MPI_COMM_WORLD, &rank) == MPI_SUCCESS);

    SparseMatrix J;
    size_t shape[2] = {2,2};
    sparse_matrix_init (&J, SPARSE_MATRIX_STORAGE_CSR, SPARSE_MATRIX_TYPE_DOUBLE, shape, 2);

    sparse_matrix_set_double (&J, 0, 1, 1);
    sparse_matrix_set_double (&J, 1, 0, 2);

    assert (sparse_matrix_get_double (&J, 0, 0) == 0);
    assert (sparse_matrix_get_double (&J, 1, 1) == 0);

    assert (sparse_matrix_get_double (&J, 0, 1) == 1);
    assert (sparse_matrix_get_double (&J, 1, 0) == 2);

    double x[2] = {0, 0};
    double b[2] = {1, 1};

    SparseSolver *solver = sparse_solver_new (SPARSE_SOLVER_METHOD_MKL_CLUSTER);
    assert (solver);
    if (rank == 0)
    {
        assert (sparse_solver_solve (solver,
                                     &J,
                                     x,
                                     b));
    }
    else
    {
        assert (sparse_solver_solve (solver, NULL, x, NULL));
    }

    assert (x[0] == 0.5);
    assert (x[1] == 1);

    if (rank == 0)
    {
        assert (sparse_solver_solve (solver,
                                     &J,
                                     x,
                                     b));
    }
    else
    {
        assert (sparse_solver_solve (solver, NULL, x, NULL));
    }

    // printf ("%f %f %d %d\n", x[0], x[1], J.rowperm[0], J.rowperm[1]);
    assert (x[0] == 0.5);
    assert (x[1] == 1);

    sparse_solver_free (solver);
    sparse_matrix_clear (&J);

    MPI_Finalize ();

    return EXIT_SUCCESS;
}
