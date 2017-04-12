#ifndef DCONSTANTS_h
#define DCONSTANTS_h

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace dlib
{

	const double dc_h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
	const double dc_hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
	const double dc_ech = GSL_CONST_MKSA_ELECTRON_CHARGE;
	const double dc_c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	const double dc_eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	const double dc_mu0 = GSL_CONST_MKSA_VACUUM_PERMEABILITY;
	const double dc_kB = GSL_CONST_MKSA_BOLTZMANN;
	const double dc_me = GSL_CONST_MKSA_MASS_ELECTRON;

	const gsl_complex GSL_COMPLEX_I = gsl_complex_rect(0.0,1.0) ;
	const gsl_complex GSL_COMPLEX_NEGI = gsl_complex_rect(0.0,-1.0) ;
}

#endif
