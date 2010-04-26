#include "branch_vec.hh"
/*
double C(double x, double y, par_list * pars){return exp( -2*gsl_pow_4(x-y)*gsl_pow_2(pars->mc) );}
double K(double x, par_list * pars){return pars->ko*exp( -2*gsl_pow_4(x)*gsl_pow_2(pars->mk) );}
double Kinv(double x, par_list * pars){return pars->ik*exp( 2*gsl_pow_4(x)*gsl_pow_2(pars->mk) );}
*/

double C(double x, double y, par_list * pars)
{
	return exp( -gsl_pow_2(x-y)*pars->mc );
}
double K(double x, par_list * pars)
{
	return pars->ko*exp( -gsl_pow_2(x)*pars->mk) ;
}

//double Kinv(double x, par_list * pars){return pars->ik*exp( gsl_pow_2(x)*pars->mk) ;}




