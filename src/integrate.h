#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#define DIM 1
#define NPTS 100

#define abserr 1e-5
#define relerr 1e-5
#define mem 10000
typedef double (* function) (double x, void * params); 
double integrate(function f, void * params, double lower_bound, double upper_bound);

//int func (double t, const double y[], double f[], void *params);
//void gsl_ode(gsl_odeiv_system sys, double Xo, double maxtime);
