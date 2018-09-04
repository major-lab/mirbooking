# OdeInt

Ordinary Differential Equation integrator

This package provide a simple and heavily parallelized ODE integrator that
implements most explicit methods via Butcher tables.

It is designed to work off-GPU given that the function callback resides there.

It solves systems in the form of $dy/dt = f(t, y)$.

```c
#include <stdio.h>
#include <odeint.h>

static void
f (double t, double *y, double *F, void *user_data)
{
    F[0] = - y[0] * y[1] + y[2]; // d[E]/dt
    F[1] = - y[0] * y[1] + y[2]; // d[S]/dt
    F[2] =   y[0] * y[1] - y[2]; // d[ES]/dt
}

int main (void)
{
  double t = 0;
  double y[4] = {10, 5, 0}; // initial [E], [S] and [ES] concentrations

  OdeIntIntegrator *integrator = odeint_integrator_new (ODEINT_METHOD_DORMAND_PRINCE,
                                                        &t, /* t is updated with exact time */
                                                        y,
                                                        3,
                                                        ODEINT_INTEGRATOR_DEFAULT_RTOL,
                                                        ODEINT_INTEGRATOR_DEFAULT_ATOL);

  // every 5 timestep
  int i;
  for (i = 0; i < 100; i++)
  {
    odeint_integrator_integrate (integrator, f, NULL, 1e-2);
    printf ("t: %f [E]:%f [S]: %f [ES]: %f", t, y[0], y[1], y[2]);
  }
}
```

