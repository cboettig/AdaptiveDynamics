#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#define abserr 1e-2
#define relerr 1e-2
#define mem 10000
typedef double (* function) (double x, void * params); 

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


/** Mutational Kernel */
double M(double y, double x, double sigma_mu)
{ 
	exp( -gsl_pow_2( (x-y)/(2*sigma_mu)))/(sqrt(2*M_PI)*sigma_mu);
}

/** Invasion fitness of y in monomorphic population x*/
double S(double y, double x, par
