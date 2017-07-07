/*! \file d_functions.h contains definitions of stepsize, gslc_sum,
 * gslc_prod, etc
 */
#ifndef D_FUNCTIONS_h

// Need this for command line output
#include <iostream>
#include <fstream>

// Need this for sqrt etc in TwoDVector
#include <cmath>

// Need these for gslc functions
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/*! \brief Operators for arithmetic with gsl_complex. Moved outside
 * the dlib namespace for visibility */
gsl_complex operator+( const gsl_complex &a, const gsl_complex &b );
gsl_complex operator+( const gsl_complex &a, const double &b );
gsl_complex operator+( const double &a, const gsl_complex &b );

gsl_complex operator-( const gsl_complex &a, const gsl_complex &b );
gsl_complex operator-( const gsl_complex &a, const double &b );
gsl_complex operator-( const double &a, const gsl_complex &b );

gsl_complex operator*( const gsl_complex &a, const gsl_complex &b );
gsl_complex operator*( const gsl_complex &a, const double &b );
gsl_complex operator*( const double &a, const gsl_complex &b );

gsl_complex operator/( const gsl_complex &a, const gsl_complex &b );
gsl_complex operator/( const gsl_complex &a, const double &b );
gsl_complex operator/( const double &a, const gsl_complex &b );

namespace dlib
{
	/*! \brief stepsize returns step size for a range [xmin,xmax] with
	 * xpts points */
	double stepsize( double xmin, double xmax, int xpts );

	/*! \brief The various legacy gslc_sum functions return the sum of
	 * gsl_complex numbers */
	gsl_complex gslc_sum( gsl_complex a, gsl_complex b );
	gsl_complex gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c );
	gsl_complex gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c, 
			gsl_complex d );

	/*! \brief The various gslc_prod functions return the product of
	 * gsl_complex numbers */
	gsl_complex gslc_prod( gsl_complex a, gsl_complex b );
	gsl_complex gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c );
	gsl_complex gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c, 
			gsl_complex d );

	/*! \brief A structure containing integration range and number of
	 * points to be used with the dIntegrator class */
	// TODO Convert this to a class?
	struct dIntParams { double xmin; double xmax; int xpts; };

	/*! \brief A structure containing integration ranges and number of
	 * points for 2D integration to be used with the d2DIntegrator class
	 */
	// TODO Covert this to a class?
	struct d2DIntParams { double xmin; double xmax; int xpts;
		double ymin; double ymax; int ypts; };

	/*!  \brief A class representing two dimensional real-valued
	 * vectors.  */
	class TwoDVector
	{
		public :
			TwoDVector( double x, double y );

			/// Return the x component of the vector
			double getx() { return _x; }
			/// Return the y component of the vector
			double gety() { return _y; }

			/// Return the length of the vector
			double getmod() { return sqrt(_x*_x + _y*_y); }
			/// Return the squared modulus of the vector
			double getmodsq() { return _x*_x + _y*_y; }

			/// Return the dot product x.y of this vector with another
			double dotproduct( TwoDVector *vec );
			/// Return |x-y| between this and another vector
			double modvecdiff( TwoDVector *vec );

		private :
			double _x, _y;
	};

	/*! \brief A class for doing one-dimensional real-valued
	 * integrations.  */
	class dIntegrator
	{
		public :
			dIntegrator( double (*f)(double,void*) ) :
				_f(f), _iginit(false)
			{}

			int make_integrand( void *params, dlib::dIntParams *ip );

			/*!
			 * Trapezoid uses the trapezoid rule.
			 * f is a pointer to the function to integrate, params are
			 * the parameters (must be typecast within f), and ip are
			 * the integration parameters
			 */
			double Trapezoid( dlib::dIntParams *ip );
			/*!
			 * Simpson uses the vanilla Simpson rule.
			 * f is a pointer to the function to integrate, params are
			 * the parameters (must be typecast within f), and ip are
			 * the integration parameters
			 */
			double Simpson( dlib::dIntParams *ip );

			~dIntegrator() 
			{
				if( _iginit == true ) delete[] _ig;
			}


		private :
			double (*_f)(double,void*);
			double *_ig;
			bool _iginit;
	};

	class d2DIntegrator
	{
		public :
			d2DIntegrator( double (*f)(double,double,void*)  ) :
				_f(f), _iginit(false)
			{}

			/*!
			 * Simpson uses the vanilla 2D Simpson rule.
			 * f is a pointer to the function to integrate, params are
			 * the parameters (must be typecast within f), and ip are
			 * the integration parameters. xpts and ypts must be odd
			 * integers greater than or equal to 5.
			 */
			int make_integrand( void *params, dlib::d2DIntParams *ip );

			/*! \brief Once the integrand is constructed, it can be
			 * dumped to a file with this function */
			int print_integrand( 
					std::ofstream *fout, dlib::d2DIntParams *ip );

			double Simpson( dlib::d2DIntParams *ip );

			~d2DIntegrator() 
			{
				if( _iginit == true ) delete[] _ig;
			}

		private :
			double (*_f)(double,double,void*);
			double *_ig;
			bool _iginit;
	};

}

#endif
