#ifndef D_FUNCTIONS_h

#include <gsl/gsl_complex.h>

#include <iostream>
#include <iomanip>
#include <cmath>

void print_gsl_complex( 
		gsl_complex z, int precision=3, double zerocut=1.0e-8 );

#endif
