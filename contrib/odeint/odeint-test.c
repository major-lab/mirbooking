#include <odeint.h>
#include <assert.h>
#include <math.h>

/* dE/dt = -kf[E][S] + kr[ES]
 * dS/dt = -kf[E][S] + kr[ES]
 * dES/dt = kf[E][S] - kr[ES]
 */
static void
f (double t, const double *y, double *F, void *user_data)
{
    F[0] = - y[0] * y[1] + y[2];
    F[1] = - y[0] * y[1] + y[2];
    F[2] =   y[0] * y[1] - y[2];
}

int
main (void)
{
    double output[ODEINT_METHOD_DORMAND_PRINCE + 1] = {0};

    OdeIntMethod method;
    for (method = 0; method <= ODEINT_METHOD_DORMAND_PRINCE; method++)
    {
        double t = 0;
        double y[4] = {50, 5, 0};

        OdeIntIntegrator *integrator = odeint_integrator_new (method,
                                                              &t,
                                                              y,
                                                              3,
                                                              ODEINT_INTEGRATOR_DEFAULT_RTOL,
                                                              ODEINT_INTEGRATOR_DEFAULT_ATOL);

        int i;
        for (i = 1; i < 10; i++)
        {
            double expected_t = i * 1e-3;
            odeint_integrator_integrate (integrator,
                                         f,
                                         NULL,
                                         expected_t);
            assert (fabs (t - expected_t) < ODEINT_INTEGRATOR_DEFAULT_RTOL * expected_t + ODEINT_INTEGRATOR_DEFAULT_ATOL);
        }

        // at equilibrium, we have [E][S]/[ES] == kr/kf
        output[method] = y[0];

        odeint_integrator_free (integrator);
    }

    // all methods should give the same output
    int i, j;
    for (i = 0; i <= ODEINT_METHOD_DORMAND_PRINCE; i++)
    {
        for (j = 0; j <= ODEINT_METHOD_DORMAND_PRINCE; j++)
        {
            assert (fabs (output[i] - output[j]) < 10 * (ODEINT_INTEGRATOR_DEFAULT_RTOL * fabs (output[j]) + ODEINT_INTEGRATOR_DEFAULT_ATOL));
        }
    }


    return 0;
}
