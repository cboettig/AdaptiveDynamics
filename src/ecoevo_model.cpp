#include "branch_vec.hh"

/*
double C(double x, double y, par_list * pars){return exp( -2*gsl_pow_4(x-y)*gsl_pow_2(pars->mc) );}
double K(double x, par_list * pars){return pars->ko*exp( -2*gsl_pow_4(x)*gsl_pow_2(pars->mk) );}
*/

/** Competition kernel */
double C(double x, double y, par_list * pars)
{
	return exp( -gsl_pow_2(x-y)*pars->mc );
}

/** Carrying capacity*/
double K(double x, par_list * pars)
{
	return pars->ko*exp( -gsl_pow_2(x)*pars->mk) ;
}

/** Mutational Kernel */
double M(double y, double x, par_list * pars)
{ 
	double sigma_mu = pars->sigma_mu;
	return exp( -0.5*gsl_pow_2( (x-y)/(sigma_mu)))/(sqrt(2*M_PI)*sigma_mu);
}

/** Invasion fitness of y in monomorphic population x*/
double S(double y, double x, par_list * pars)
{
	return 1 - C(y,x, pars)*K(x, pars)/K(y,pars);
}

/** an integrand: product of M(y,x) and S(y,x) */
double M_x_S(double x, void * params)
{
	par_list * pars = (par_list *) params;
	double y = pars->y;
	return M(y,x,pars)*S(y,x,pars);
}
/** Probability of mutant arising and establishing in coexistence region*/
double P(double x, par_list * pars)
{
	double MIN = -10, MAX = 10;
	double ratio = pars->mc/pars->mk;
	double phi = x*(ratio + 1) / (ratio - 1);
	double psi = x*(ratio - 1) / (ratio + 1);
	return pars->mu*K(x,pars)*(integrate(M_x_S, pars, MIN, phi) +  integrate(M_x_S, pars, psi, MAX) );
}

/* // ode solver for mean path not implemented yet!
int func (double t, const double y[], double f[], void *params)
{
	par_list * pars = (par_list *) params;
	f[0] = -pars->mu*gsl_pow_2(pars->sigma_mu)*K(y[0], pars)*y[0]*pars->mk;
}
double mean_path(double t, par_list * pars)
{
	int *jac = NULL;
	gsl_odeiv_system sys = {func, jac, DIM, pars};
	gsl_ode(sys, pars->xo, 1e4);

}
*/

/** Approximation of the canonical equation path for constant K(x) */
double mean_path_approx(double t, par_list * pars)
{
	return pars->xo*exp( - t*pars->mu*gsl_pow_2(pars->sigma_mu)*pars->ko*2*pars->mk);
}

double Pintegrand(double t, void * params)
{
	par_list * pars = (par_list *) params;
	return P(mean_path_approx(t, pars), pars);
}

/** First approximation of waiting time distribution for phase 1 */
double waiting_time_1(double T, par_list * pars)
{
	return P(mean_path_approx(T, pars), pars) * exp(-integrate(Pintegrand, pars, 0, T));
}

double mean_integrand(double t, void * params)
{
	par_list * pars = (par_list *) params;
	return t*waiting_time_1(t, pars);
}

double mean_waiting_time_1(par_list * pars)
{
	return integrate(mean_integrand, pars, 0, 1e5);
}

