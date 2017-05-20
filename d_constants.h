/*! \file
 * Contains shorthand definitions of various physical constants and
 * complex numbers
 */

#ifndef DCONSTANTS_h
#define DCONSTANTS_h

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace dlib
{

	/*! \brief Planck's constant (MKSA) */
	const double dc_h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;
	/*! \brief Reduced Planck's constant (MKSA) */
	const double dc_hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR; 
	/*! \brief Electron charge (MKSA) */
	const double dc_ech = GSL_CONST_MKSA_ELECTRON_CHARGE;
	/*! \brief Speed of light (MKSA) */
	const double dc_c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	/*! \brief Vacuum permittivity (MKSA) */
	const double dc_eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
	/*! \brief Vacuum permeability (MKSA) */
	const double dc_mu0 = GSL_CONST_MKSA_VACUUM_PERMEABILITY;
	/*! \brief Boltzmann constant (MKSA) */
	const double dc_kB = GSL_CONST_MKSA_BOLTZMANN;
	/*! \brief Electron mass (MKSA) */
	const double dc_me = GSL_CONST_MKSA_MASS_ELECTRON;

	/*! \brief gsl_complex i */
	const gsl_complex gslc_I = gsl_complex_rect(0.0,1.0);
	/*! \brief gsl_complex -i */
	const gsl_complex gslc_negI = gsl_complex_rect(0.0,-1.0);
	/*! \brief gsl_complex 1 */
	const gsl_complex gslc_one = gsl_complex_rect(1.0,0.0);
	/*! \brief gsl_complex -1 */
	const gsl_complex gslc_negone = gsl_complex_rect(-1.0,0.0);
}

#endif
