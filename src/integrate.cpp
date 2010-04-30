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

void gsl_ode(gsl_odeiv_system sys, double Xo, double maxtime){

	/* Range of time to simulate*/	
	double t = 0.0, t1 = maxtime;
	/* Initial step size, will be modified as needed by adaptive alogorithm */
	double h = 1e-6;
	/* initial conditions */
	double y[DIM] = { Xo };


	/* Define method as Embedded Runge-Kutta Prince-Dormand (8,9) method */
	const gsl_odeiv_step_type * T
	 = gsl_odeiv_step_rk4;
//	 = gsl_odeiv_step_rk8pd;
	/* allocate stepper for our method of correct dimension*/
	gsl_odeiv_step * s
	 = gsl_odeiv_step_alloc (T, DIM);
	/* control will maintain absolute and relative error */
	gsl_odeiv_control * c
	 = gsl_odeiv_control_y_new (1e-6, 0.0);
	/* allocate the evolver */
	gsl_odeiv_evolve * e
	 = gsl_odeiv_evolve_alloc (DIM);

	/*dummy to make easy to switch to regular printing */
	double ti = t1; 
	int i;

	/* Uncomment the outer for loop to print *
	 * 100 sample pts at regular intervals   */
	for (i = 1; i <= NPTS; i++){   ti = i * t1 / NPTS;
		while (t < ti)
		{
			int status = gsl_odeiv_evolve_apply (e, c, s,
												&sys,
												&t, t1,
												&h, y);
			if (status != GSL_SUCCESS)
			   break;
		}
		printf("%.6e %.6e\n",t,y[0]);  
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
}




