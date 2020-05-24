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

    int dIntegrator::make_integrand( void *params, dlib::dIntParams *ip )
    {
        double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );

        // Check if memory already allocated. If it is, delete
        if( _iginit == true ) delete[] _ig;

        // Allocate memory
        _ig = new double[ip->xpts];

    #pragma omp parallel for
        for( int ix=0; ix< ip->xpts; ix++ )
        {
            double thisx = ip->xmin + (double)ix*xstep;
            _ig[ ix ] = _f( thisx, params );
        }

        _iginit = true;
        
        return 0;
    }

    double dIntegrator::Trapezoid( dlib::dIntParams *ip )
    {
        if( _iginit == false )
        {
            std::cout << "WARNING: _ig not populated in dIntegrator::Trapezoid. Returning 0." << std::endl;
            return 0.0;
        }

        double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );

        double intsum = 0.5*( _ig[0] + _ig[ip->xpts-1] );
        for( int jj=1; jj<ip->xpts-1; jj++ )
            intsum += _ig[jj];

        return xstep*intsum;
    }

    double dIntegrator::Simpson( dlib::dIntParams *ip )
    {
        if( _iginit == false )
        {
            std::cout << "WARNING: _ig not populated in dIntegrator::Simpson. Returning 0." << std::endl;
            return 0.0;
        }

        double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );

        double intsum = _ig[0] + _ig[ip->xpts-1] + 4*_ig[ip->xpts-2];
        for( int jj=1; jj<ip->xpts-2; jj+=2 )
        {
            //std::cout << jj << " ";
            intsum += 4.0*_ig[jj];
            intsum += 2.0*_ig[jj+1];
        }
        //std::cout << std::endl;
        return intsum*xstep/3.0;
    }

    int d2DIntegrator::make_integrand( 
                    void *params, dlib::d2DIntParams *ip )
    {
        double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );
        double ystep = dlib::stepsize( ip->ymin, ip->ymax, ip->ypts );

        // Allocate memory
        if( _iginit == true ) delete[] _ig;

        _ig = new double[ip->xpts*ip->ypts];

        for( int iy=0; iy< ip->ypts; iy++ )
        {
            double thisy = ip->ymin + (double)iy*ystep;
            int igoffset = iy*ip->xpts;

    #pragma omp parallel for
            for( int ix=0; ix< ip->xpts; ix++ )
            {
                double thisx = ip->xmin + (double)ix*xstep;
                _ig[ igoffset + ix ] = _f( thisx, thisy, params );
            }
        }

        _iginit = true;
        
        return 0;
    }
        
    int d2DIntegrator::print_integrand( 
            std::ofstream *fout, dlib::d2DIntParams *ip )
    {
        double xstep = dlib::stepsize( ip->xmin, ip->xmax, ip->xpts );
        double ystep = dlib::stepsize( ip->ymin, ip->ymax, ip->ypts );

        for( int iy = 0; iy < ip->ypts; iy++ )
        {
            double thisy = ip->ymin + (double)iy * ystep;
            int iyoffset = ip->xpts*iy;
            for( int ix = 0; ix < ip->xpts; ix++ )
            {
                double thisx = ip->xmin + (double)ix * xstep;
                *fout << thisx << " " << thisy << " " 
                    << _ig[ iyoffset + ix ] << std::endl;
            }
            *fout << std::endl;
        }
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
