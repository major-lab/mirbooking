#include "odeint.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* TODO:
 * - implicit methods
 * - off-GPU
 */

/**
 * The maximum number of steps.
 */
#define MAX_STEP 32

/**
 * @name: Name of the method
 * @steps: Number of steps
 * @last_step_is_update: Indicate if the last step correspond to the update, in
 * which case no final update is computed and the last transient 'y' is used.
 * @c: Whole step weight
 * @a: Previous steps weights
 * @b: Final step weights
 * @order: Order of the final update used to adjust the step size
 * @e: Error estimate
 * @error_estimate_order: Order of the error estimate (unused)
 */
typedef struct _OdeIntMeta
{
    const char *name;
    int         steps;
    int         last_step_is_update;
    int         error_estimate;
    double      c[MAX_STEP];
    double      a[(MAX_STEP*(MAX_STEP-1))/2];
    double      b[MAX_STEP];
    int         order;
    double      e[MAX_STEP];
    int         error_estimate_order;
} OdeIntMeta;

/*
 * Butcher table
 *
 * For method that provide no error estimate, we reuse the last step, which
 * will yield no error.
 */
const OdeIntMeta INTEGRATOR_META[] =
{
    {"euler", 1, 0, 0,
        {0.0},
        {},
        {1.0}, 1},
    {"heuns", 2, 0, 0,
        {0.0,     1.0},
        {1.0},
        {1.0/2.0, 1.0/2.0}, 2},
    // FIXME: this method yield very unaccurate results
    {"heuns-euler", 2, 0, 1,
        {0.0,     1.0},
        {1.0},
        {1.0/2.0, 1.0/2.0}, 2,
        {1.0,     0.0},     1},
    {"bogacki-shampine", 4, 1, 1,
        {0.0,      1.0/2.0, 3.0/4.0, 1.0},
        {1.0/2.0,
         0.0,      3.0/4.0,
         2.0/9.0,  1.0/3.0, 4.0/9.0},
        {2.0/9.0,  1.0/3.0, 4.0/9.0, 0.0},     3,
        {7.0/24.0, 1.0/4.0, 1.0/3.0, 1.0/8.0}, 2},
    {"runge-kutta", 4, 0, 0,
        {0.0,     1.0/2.0, 1.0/2.0, 1.0},
        {1.0/2.0,
         0.0,     1.0/2.0,
         0.0,     0.0,     1.0},
        {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}, 4},
    {"runge-kutta-fehlberg", 6, 0, 1,
        {0.0,           1.0/4.0,        3.0/8.0,        12.0/13.0,       1.0,       1.0/2.0},
        {1.0/4.0,
         3.0/32.0,      9.0/32.0,
         1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,
         439.0/216.0,   -8.0,           3680.0/513.0,   -845/4104.0,
         -8.0/27.0,     2.0,            -3544.0/2565.0, 1859.0/4104.0,   -11.0/40.0},
        {16.0/135.0,    0.0,            6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0}, 5,
        {25.0/216.0,    0.0,            1408.0/2565.0,  2197.0/4104.0,   -1.0/5.0,  0.0},      4},
    {"cash-karp", 6, 0, 1,
        {0.0,            1.0/5.0,     3.0/10.0,        3.0/5.0,          1.0,           7.0/8.0},
        {1.0/5.0,
         3.0/40.0,       9.0/40.0,
         3.0/10.0,       -9.0/10.0,   6.0/5.0,
         -11.0/54.0,     5.0/2.0,     -70.0/27.0,      35.0/27.0,
         1631.0/55296.0, 175.0/512.0, 575.0/13824.0,   44275.0/110592.0, 253.0/4096.0},
        {37.0/378.0,     0.0,         250.0/621.0,     125.0/594.0,      0.0,           512.0/1771.0}, 5,
        {2825.0/27648.0, 0.0,         18575.0/48384.0, 13525.0/55296.0,  277.0/14336.0, 1.0/4.0},      4},
    {"dormand-prince", 7, 1, 1,
        {0.0,            1.0/5.0,         3.0/10.0,       4.0/5.0,     8.0/9.0,           1.0,          1.0},
        {1.0/5.0,
         3.0/40.0,       9.0/40.0,
         44.0/45.0,      -56.0/15.0,      32.0/9.0,
         19372.0/6561,   -25360.0/2187.0, 64448.0/6561.0, 212.0/729.0,
         9017.0/3168,    -355.0/33.0,     46732.0/5247.0, 49.0/176.0,  -5103.0/18656.0,
         35.0/384.0,     0.0,             500.0/1113.0,   125.0/192.0, -2187.0/6784.0,    11.0/84.0},
        {35.0/384.0,     0.0,             500.0/1113.0,   125.0/192.0, -2187.0/6784.0,    11.0/84.0,    0.0},      5,
        {5179.0/57600.0, 0.0,             7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0}, 4}
};

struct _OdeIntIntegrator
{
    const OdeIntMeta *integrator_meta;
    double *t;
    double *y;
    size_t  n;
    double  rtol;
    double  atol;
    double *transient_y;
    double *transient_ye;
    double *F; // dy/dt
};

/**
 * odeint_integrator_new:
 * @t0: Initial timestep
 * @y0: Initial value for the system which will be updated throughout the
 * integration
 * @n: Number of equations and variables (each variable has a differential equation w.r.t. time 't')
 * @rtol: Relative tolerance (i.e. number of significant digits)
 * @atol: Absolute tolerance
 */
OdeIntIntegrator *
odeint_integrator_new (OdeIntMethod  method,
                       double       *t0,
                       double       *y0,
                       size_t        n,
                       double        rtol,
                       double        atol)
{
    OdeIntIntegrator *ret = malloc (sizeof (OdeIntIntegrator));

    ret->integrator_meta = &INTEGRATOR_META[method];

    ret->t = t0;
    ret->y = y0;
    ret->n = n;
    ret->F = calloc (ret->integrator_meta->steps * n, sizeof (double));

    /* transient states for multi-step methods */
    ret->transient_y = malloc ( n * sizeof (double) );
    ret->transient_ye = malloc ( n * sizeof (double) );

    ret->rtol = rtol;
    ret->atol = atol;

    assert (ret->F != NULL);

    return ret;
}

/**
 * odeint_integrator_integrate:
 * @func: The function to integrate
 * @user_data:
 * @w: The domain width for which the integration is performed
 *
 * Integrate dy/dt = @func(t, y) over [t, t + @w].
 */
void
odeint_integrator_integrate (OdeIntIntegrator *self,
                             OdeIntFunc        func,
                             void             *user_data,
                             double            w)
{
    double t0 = *self->t;
    // reference: Watts, “Starting Step Size for an ODE Solver.”
    // TODO: We can use an estimate of the second-order derivatives to scale
    // the step size
    double h  = sqrt (2.0) * pow (self->atol, 1.0 / (self->integrator_meta->order + 1));

    /* transient state for the multi-step method */
    double *y  = self->transient_y;
    double *ye = self->transient_ye;

    while (self->rtol * fabs (w - (*self->t - t0)) >= self->atol)
    {
        // ensure we always land exactly on the upper integration bound
        // if we went too far, the step size will be negative and still point
        // toward the integration bound
        h = fmin (h, w - (*self->t - t0));

        int step;
        for (step = 0; step < self->integrator_meta->steps; step++)
        {
            /* restore state at the beginning of the step */
            memcpy (y, self->y, self->n * sizeof (double));

            int i, prev_step;
            for (prev_step = 0; prev_step < step; prev_step++)
            {
                #pragma omp parallel for
                for (i = 0; i < self->n; i++)
                {
                    y[i] += h * self->integrator_meta->a[(step * step + step)/2 + prev_step] * self->F[prev_step * self->n + i];
                }
            }

            /* update from previous steps */
            func (*self->t + h * self->integrator_meta->c[step], y, self->F + (step * self->n), user_data);
        }

        // final update
        // For some methods, the last step (stored in y) is the same as the
        // update
        // FIXME: reusing the last step is broken?
        // if (!self->integrator_meta->last_step_is_update)
        {
            memcpy (y, self->y, self->n * sizeof (double));
            int i, step;
            for (step = 0; step < self->integrator_meta->steps; step++)
            {
                #pragma omp parallel for
                for (i = 0; i < self->n; i++)
                {
                    y[i] += h * self->integrator_meta->b[step] * self->F[step * self->n + i];
                }
            }
        }

        // error estimate
        if (self->integrator_meta->error_estimate)
        {
            memcpy (ye, self->y, self->n * sizeof (double));
            int i, step;
            for (step = 0; step < self->integrator_meta->steps; step++)
            {
                #pragma omp parallel for
                for (i = 0; i < self->n; i++)
                {
                    ye[i] += h * self->integrator_meta->e[step] * self->F[step * self->n + i];
                }
            }

            double error = 0;
            #pragma omp parallel for reduction(max:error)
            for (i = 0; i < self->n; i++)
            {
                error = fmax (error, fabs (y[i] - ye[i]));
            }

            if (self->rtol * error <= self->atol)
            {
                *self->t += h;
                memcpy (self->y, y, self->n * sizeof (double));
            }

            // adjust step size
            h *= fmin (fmax (0.84 * pow (self->atol / (self->rtol * error), 1.0 / self->integrator_meta->order), 0.1), 4.0);
        }
        else
        {
            *self->t += h;
            memcpy (self->y, y, self->n * sizeof (double));
        }
    }
}

void odeint_integrator_free (OdeIntIntegrator *self)
{
    free (self->transient_y);
    free (self->transient_ye);
    free (self->F);
    free (self);
}
