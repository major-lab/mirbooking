#include <stddef.h>

typedef struct _PoissonBinomial {
    double *p;
    size_t n;
    double *pmf;
} PoissonBinomial;

void   pb_init    (PoissonBinomial *pb, double *p, size_t n);
void   pb_destroy (PoissonBinomial *pb);
double pb_pmf     (PoissonBinomial *pb, unsigned int k);
