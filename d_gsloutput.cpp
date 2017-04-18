#include "./d_gsloutput.h"

using namespace dlib;

void dlib::print_gsl_complex( gsl_complex z, int precision, double zerocut )
{
	std::ios state(NULL);
	state.copyfmt(std::cout);

	std::cout << std::setprecision(precision);

	if( fabs( GSL_REAL(z) ) < zerocut && fabs( GSL_IMAG(z) ) < zerocut )
		std::cout << "0";
	else if( fabs( GSL_REAL(z) ) < zerocut )
		std::cout << std::showpos << GSL_IMAG(z) << "i";
	else if( fabs( GSL_IMAG(z) ) < zerocut )
		std::cout << GSL_REAL(z);
	else
		std::cout << GSL_REAL(z) << std::showpos << GSL_IMAG(z)
			<< std::noshowpos << "i";
	
	std::cout.copyfmt(state);
}

