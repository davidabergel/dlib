#include "./d_gsloutput.h"

using namespace dlib;

void dlib::print_gslc( gsl_complex z, int precision, double zerocut )
{
    // Save cout state to restore later
    std::ios state(NULL);
    state.copyfmt(std::cout);

    std::cout << std::setprecision(precision);

    // Check special cases:
    // Close to 0
    if( fabs( GSL_REAL(z) ) < zerocut && fabs( GSL_IMAG(z) ) < zerocut )
        std::cout << "0";
    // Close to i
    else if( fabs( GSL_REAL(z) ) < zerocut && fabs(GSL_IMAG(z)-1.0) < zerocut )
        std::cout << "i";
    // Close to -i
    else if( fabs( GSL_REAL(z) ) < zerocut && fabs(GSL_IMAG(z)+1.0) < zerocut )
        std::cout << "-i";
    // Imaginary
    else if( fabs( GSL_REAL(z) ) < zerocut )
        std::cout << std::showpos << GSL_IMAG(z) << "i";
    // Real
    else if( fabs( GSL_IMAG(z) ) < zerocut )
        std::cout << GSL_REAL(z);
    // General output
    else
        std::cout << GSL_REAL(z) << std::showpos << GSL_IMAG(z)
            << std::noshowpos << "i";

    // Restore cout state
    std::cout.copyfmt(state);
}

