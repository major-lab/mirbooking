#include "sparse.h"

int
superlu_solve (SparseSolver *solver,
               SparseMatrix *A,
               double       *x,
               double       *b);
