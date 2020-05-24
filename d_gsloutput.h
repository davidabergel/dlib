/*! \file
 * Contains output routines for gsl_complex numbers
 */
#ifndef D_GSLOUTPUT_h
#define D_GSLOUTPUT_h

#include <gsl/gsl_complex.h>

#include <iostream>
#include <iomanip>
#include <cmath>

namespace dlib
{
    /*! \brief Prints a gsl_complex z to cout with specified precision.
    * If z<zerocut then zero is displayed */
    void print_gslc( 
            gsl_complex z, int precision=3, double zerocut=1.0e-8 );
}

#endif
