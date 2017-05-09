#ifndef D_FUNCTIONS_h

// Need this for command line output
#include <iostream>

// Need these for gslc functions
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

namespace dlib
{
	double stepsize( double xmin, double xmax, int xpts );

	gsl_complex gslc_sum( gsl_complex a, gsl_complex b );
	gsl_complex gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c );
	gsl_complex gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c, 
			gsl_complex d );

	gsl_complex gslc_prod( gsl_complex a, gsl_complex b );
	gsl_complex gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c );
	gsl_complex gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c, 
			gsl_complex d );

	struct dIntParams { double xmin; double xmax; int xpts; };

	class dIntegrator
	{
		public :
			dIntegrator() {}

			double Trapezoid(
					double (*f)(double,void*), void *params,
					dlib::dIntParams *ip );
			double Simpson( 
					double (*f)(double,void*), void *params, 
					dlib::dIntParams *ip );

			~dIntegrator() {}

		private :
	};
}

#endif
