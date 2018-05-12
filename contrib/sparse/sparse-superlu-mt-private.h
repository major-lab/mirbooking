#include "sparse.h"

int
superlu_mt_solve (SparseSolver *solver,
                  SparseMatrix *A,
                  double       *x,
                  double       *b);
