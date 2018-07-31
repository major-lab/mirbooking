# OdeInt

Ordinary Differential Equation integrator

This package provide a simple and heavily parallelized ODE integrator that
implements most explicit methods via Butcher tables.

It is designed to work off-GPU given that the function callback resides there.

It solves systems in the form of $dy/dt = f(t, y)$.

```c
double t = 0;
double y[4] = {0, 0, 0, 0};

OdeIntIntegrator *integrator = odeint_integrator_new (ODEINT_METHOD_DORMAND_PRINCE,
                                                      t,
                                                      y,
                                                      4,
                                                      1e-8);

// every 5 timestep
int i;
for (i = 0; i < 100; i++)
{
  odeint_integrator_integrate (integrator, f, NULL, 5);
  // y contains the state of 't + 5 * i'
}
```

