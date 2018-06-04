#include <sparse.h>
#include <assert.h>
#include <stdlib.h>

int
main (void)
{
    SparseMatrix J;
    size_t shape[2] = {2,2};
    sparse_matrix_init (&J, SPARSE_MATRIX_STORAGE_CSR, shape, 2);

    sparse_matrix_set_value (&J, 0, 1, 1);
    sparse_matrix_set_value (&J, 1, 0, 2);

    assert (sparse_matrix_get_value (&J, 0, 0) == 0);
    assert (sparse_matrix_get_value (&J, 1, 1) == 0);

    assert (sparse_matrix_get_value (&J, 0, 1) == 1);
    assert (sparse_matrix_get_value (&J, 1, 0) == 2);

    double x[2] = {0, 0};
    double b[2] = {1, 1};

    assert (J.rowperm[0] == J.rowperm[J.shape[1] - 1]);

    SparseSolver *solver = sparse_solver_new (SPARSE_SOLVER_METHOD_SUPERLU);
    assert (solver);
    assert (sparse_solver_solve (solver,
                                 &J,
                                 x,
                                 b));

    assert (x[0] == 0.5);
    assert (x[1] == 1);
    assert (J.rowperm[0] != J.rowperm[J.shape[1] - 1]);

    // resolving with permcol/permrow
    assert (sparse_solver_solve (solver,
                                 &J,
                                 x,
                                 b));

    assert (x[0] == 0.5);
    assert (x[1] == 1);

    sparse_solver_free (solver);
    sparse_matrix_clear (&J);

    return EXIT_SUCCESS;
}
