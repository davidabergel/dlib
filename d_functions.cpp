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

	int d2DIntegrator::make_integrand( double (*f)(double,double,void*),
					void *params, dlib::d2DIntParams *ip )
	{
		double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );
		double ystep = dlib::stepsize( ip->ymin, ip->ymax, ip->ypts );

		// Allocate memory
		_ig = new double[ip->xpts*ip->ypts];

		for( int iy=0; iy< ip->ypts; iy++ )
		{
			double thisy = ip->ymin + (double)iy*ystep;
			int igoffset = iy*ip->xpts;

			//TODO Parallelize this
			for( int ix=0; ix< ip->xpts; ix++ )
			{
				double thisx = ip->xmin + (double)ix*xstep;
				_ig[ igoffset + ix ] = f( thisx, thisy, params );
			}
		}

		_iginit = true;
		
		return 0;
	}
		
	double d2DIntegrator::Simpson( dlib::d2DIntParams *ip )
	{
		if( _iginit == false )
		{
			std::cout << "WARNING: Attempting to integrate when integrand not set" 
				<< std::endl << "   returning 0" << std::endl;
			return 0.0;
		}

		int xpts = ip->xpts;
		int ypts = ip->ypts;

		double xstep = dlib::stepsize( ip->xmin, ip->xmax, xpts );
		double ystep = dlib::stepsize( ip->ymin, ip->ymax, ypts );
		// Sum integral
		double runsum=0.0;

		// Corners
		runsum += _ig[0] + _ig[xpts-1] + _ig[xpts*(ypts-1)] + _ig[xpts*ypts-1];

		// First row
		runsum += 4.0*_ig[1];
		for( int ix=2; ix<xpts-2; ix+=2 )
		{
			runsum += 2.0*_ig[ix] + 4.0*_ig[ix+1];
		}

		// Second row
		int secoffset = xpts;
		runsum += 4.0*_ig[secoffset] + 16.0*_ig[secoffset+1];
		for( int ix=2; ix<xpts-2; ix+=2 )
		{
			runsum += 8.0*_ig[secoffset+ix] + 16.0*_ig[secoffset+ix+1];
		}
		runsum += 4.0*_ig[secoffset+xpts-1];


		// Middle rows
		for( int iy=2; iy<ypts-2; iy+=2 )
		{
			int thisoffset = iy*xpts;
			// First column
			runsum += 2.0*_ig[thisoffset] + 4.0*_ig[thisoffset+xpts];
			// Second column
			runsum += 8.0*_ig[thisoffset+1] + 16.0*_ig[thisoffset+xpts+1];
			// Last column
			runsum += 2.0*_ig[thisoffset+xpts-1]
				+ 4.0*_ig[thisoffset+2*xpts-1];
			for( int ix=2; ix<xpts-2; ix+=2 )
			{
				runsum += 4.0*_ig[thisoffset+ix] + 8.0*_ig[thisoffset+ix+1];
				runsum += 8.0*_ig[thisoffset+xpts+ix]
					+ 16.0*_ig[thisoffset+xpts+ix+1];
			}
		}

		// Last row
		int lastoffset = xpts*(ypts-1);
		runsum += 4.0*_ig[lastoffset+1];
		for( int ix=2; ix<xpts-2; ix+=2 )
		{
			runsum += 2.0*_ig[lastoffset+ix] + 4.0*_ig[lastoffset+ix+1];
		}
		
		return runsum*xstep*ystep/9.0;
	}

}
