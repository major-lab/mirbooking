#include "pb.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#if HAVE_FFTW3
#include <complex.h>
#include <fftw3.h>
#endif

static void
pb_compute_pmf (PoissonBinomial *pb)
{
#if HAVE_FFTW3
    fftw_complex *in;
    fftw_plan plan;
    size_t N = pb->n + 1;

    in = fftw_malloc (sizeof (fftw_complex) * N);

    plan = fftw_plan_dft_c2r_1d (N,
                                 in,
                                 pb->pmf,
                                 FFTW_ESTIMATE);

    // initialize
    size_t n;
    for (n = 0; n < N; n++)
    {
        in[n] = 1.0;

        size_t k;
        for (k = 0; k < pb->n; k++)
        {
            in[n] *= pb->p[k] * cexp (-I * 2.0 * M_PI * n / N) + (1 - pb->p[k]);
        }
    }

    fftw_execute (plan);

    for (n = 0; n < N; n++)
    {
        pb->pmf[n] /= N;
        pb->pmf[n] = fmax (pb->pmf[n], 0);
    }

    fftw_destroy_plan (plan);
    fftw_free (in);
#else
    pb->pmf[0] = 1;
    int i;
    for (i = 0; i < pb->n; i++)
    {
        pb->pmf[0] *= (1 - pb->p[i]);
    }

    int k;
    for (k = 1; k < pb->n + 1; k++)
    {
        double pk = 0;
        int i;
        for (i = 0; i < k; i++)
        {
            double ti = 0;
            int j;
            for (j = 0; j < pb->n; j++)
            {
                ti += pow (pb->p[j] / (1 - pb->p[j]), i + 1);
            }
            pk += (i % 2 == 0 ? 1 : -1) * pb->pmf[k - i - 1] * ti;
        }
        pb->pmf[k] = pk / k;
    }
#endif
}

/**
 *
 */
void
pb_init (PoissonBinomial *pb, double *p, size_t n)
{
    pb->p = memcpy (malloc (sizeof (double) * n), p, sizeof (double) * n);
    pb->n = n;
    pb->pmf = malloc (sizeof (double) * (n + 1));
    pb_compute_pmf (pb);
}

/**
 *
 */
double
pb_pmf (PoissonBinomial *pb, unsigned int k)
{
    return pb->pmf[k];
}

void
pb_destroy (PoissonBinomial *pb)
{
    free (pb->p);
    free (pb->pmf);
}
