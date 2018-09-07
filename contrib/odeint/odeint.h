#include <stddef.h>

#define ODEINT_INTEGRATOR_DEFAULT_RTOL 1e-6
#define ODEINT_INTEGRATOR_DEFAULT_ATOL 1e-12

typedef enum _OdeIntMethod
{
    ODEINT_METHOD_EULER,
    ODEINT_METHOD_HEUNS,
    ODEINT_METHOD_HEUNS_EULER,
    ODEINT_METHOD_BOGACKI_SHAMPINE,
    ODEINT_METHOD_RUNGE_KUTTA,
    ODEINT_METHOD_RUNGE_KUTTA_FELHBERG,
    ODEINT_METHOD_CASH_KARP,
    ODEINT_METHOD_DORMAND_PRINCE
} OdeIntMethod;

typedef struct _OdeIntIntegrator OdeIntIntegrator;

/**
 * @t: The timestep
 * @y: The state of the ODE
 * @F: dy/dt
 */
typedef void (*OdeIntFunc) (double t, const double *y, double *F, void *user_data);

OdeIntIntegrator * odeint_integrator_new       (OdeIntMethod method, double *t0, double *y0, size_t n, double rtol, double atol);
void               odeint_integrator_integrate (OdeIntIntegrator *self, OdeIntFunc func, void* user_data, double w);
void               odeint_integrator_free      (OdeIntIntegrator *self);
