#include "./d_functions.h"

gsl_complex operator+( const gsl_complex &a, const gsl_complex &b )
{
	return gsl_complex_add( a, b );
}

gsl_complex operator+( const gsl_complex &a, const double &b )
{
	return gsl_complex_add_real( a, b );
}

gsl_complex operator+( const double &a, const gsl_complex &b )
{
	return gsl_complex_add_real( b, a );
}

gsl_complex operator-( const gsl_complex &a, const gsl_complex &b )
{
	return gsl_complex_sub( a, b );
}

gsl_complex operator-( const gsl_complex &a, const double &b )
{
	return gsl_complex_sub_real( a, b );
}

gsl_complex operator-( const double &a, const gsl_complex &b )
{
	return gsl_complex_add_real( gsl_complex_negative( b ), a );
}

gsl_complex operator*( const gsl_complex &a, const gsl_complex &b )
{
	return gsl_complex_mul( a, b );
}

gsl_complex operator*( const gsl_complex &a, const double &b )
{
	return gsl_complex_mul_real( a, b );
}

gsl_complex operator*( const double &a, const gsl_complex &b )
{
	return gsl_complex_mul_real( b, a );
}

gsl_complex operator/( const gsl_complex &a, const gsl_complex &b )
{
	return gsl_complex_div( a, b );
}

gsl_complex operator/( const gsl_complex &a, const double &b )
{
	return gsl_complex_div_real( a, b );
}

gsl_complex operator/( const double &a, const gsl_complex &b )
{
	return gsl_complex_mul_real( gsl_complex_inverse(b), a );
}

namespace dlib
{
	double stepsize( double xmin, double xmax, int xpts )
	{
		return (xmax-xmin)/(double)(xpts-1);
	}


	gsl_complex gslc_sum( gsl_complex a, gsl_complex b )
	{
		return gsl_complex_add( a, b );
	}

	gsl_complex gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c )
	{
		return gsl_complex_add( a, gsl_complex_add( b, c ) );
	}

	gsl_complex gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c, 
				gsl_complex d )
	{
		return gsl_complex_add( gsl_complex_add(a,b), gsl_complex_add(c,d) );
	}

	gsl_complex gslc_prod( gsl_complex a, gsl_complex b )
	{
		return gsl_complex_mul( a, b );
	}

	gsl_complex gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c )
	{
		return gsl_complex_mul( a, gsl_complex_mul( b, c ) );
	}

	gsl_complex gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c, 
				gsl_complex d )
	{
		return gsl_complex_mul( gsl_complex_mul(a,b), gsl_complex_mul(c,d) );
	}

	TwoDVector::TwoDVector( double x, double y ) : _x(x), _y(y) {}

	double dlib::TwoDVector::dotproduct( TwoDVector *vec )
	{
		return _x*vec->getx() + _y*vec->gety();
	}

	double TwoDVector::modvecdiff( TwoDVector *vec )
	{
		double xdiff = _x - vec->getx();
		double ydiff = _y - vec->gety();
		return sqrt( xdiff*xdiff + ydiff*ydiff );
	}

	double dIntegrator::Trapezoid(
			double (*f)(double,void*), void *params, dlib::dIntParams *ip )
	{
		double *ig = new double[ip->xpts];
		double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );

#pragma omp parallel for schedule(dynamic)
		for( int jj=0; jj< ip->xpts-1; jj++ )
		{
			double thisx = ip->xmin + (double)jj*xstep;
			ig[jj] = f(thisx,params);
		}

		double intsum = 0.5*( ig[0] + ig[ip->xpts-1] );
		for( int jj=1; jj<ip->xpts-1; jj++ )
			intsum += ig[jj];

		delete[] ig;
		return xstep*intsum;
	}

	double dIntegrator::Simpson( 
			double (*f)(double,void*), void *params, dlib::dIntParams *ip )
	{
		if( ip->xpts % 2 == 0 )
		{
			std::cout << "ERROR: xpts in dIntegrator::Simpson is "
				<< ip->xpts << ": must be odd" << std::endl;
			return 0;
		}

		double *ig = new double[ip->xpts];
		double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );

#pragma omp parallel for schedule(dynamic)
		for( int jj=0; jj<ip->xpts; jj++ )
		{
			double thisx = ip->xmin + (double)jj*xstep;
			ig[jj] = f(thisx,params);
		}

		double intsum = ig[0] + ig[ip->xpts-1] + 4*ig[ip->xpts-2];
		for( int jj=1; jj<ip->xpts-2; jj+=2 )
		{
			//std::cout << jj << " ";
			intsum += 4.0*ig[jj];
			intsum += 2.0*ig[jj+1];
		}
		//std::cout << std::endl;
		delete[] ig;
		return intsum*xstep/3.0;
	}
}
