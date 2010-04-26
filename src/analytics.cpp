#include "integrate.h"

double integrate(function f, void * params, double lower_bound, double upper_bound)
{
	gsl_integration_workspace * w 
	 = gsl_integration_workspace_alloc (mem);

	double result, error;
	double bounds[2] = {lower_bound, upper_bound};

	gsl_function F;
	F.function = f;
	F.params = params;

	gsl_integration_qags (&F, bounds[0], bounds[1], abserr, relerr,
						mem, w, &result, &error); 
	gsl_integration_workspace_free (w);
	return result;
}


