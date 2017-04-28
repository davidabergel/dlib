#include "./d_functions.h"

double dlib::stepsize( double xmin, double xmax, int xpts )
{
	return (xmax-xmin)/(double)(xpts-1);
}

gsl_complex dlib::gslc_sum( gsl_complex a, gsl_complex b )
{
	return gsl_complex_add( a, b );
}

gsl_complex dlib::gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c )
{
	return gsl_complex_add( a, gsl_complex_add( b, c ) );
}

gsl_complex dlib::gslc_sum( gsl_complex a, gsl_complex b, gsl_complex c, 
			gsl_complex d )
{
	return gsl_complex_add( gsl_complex_add(a,b), gsl_complex_add(c,d) );
}

gsl_complex dlib::gslc_prod( gsl_complex a, gsl_complex b )
{
	return gsl_complex_mul( a, b );
}

gsl_complex dlib::gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c )
{
	return gsl_complex_mul( a, gsl_complex_mul( b, c ) );
}

gsl_complex dlib::gslc_prod( gsl_complex a, gsl_complex b, gsl_complex c, 
			gsl_complex d )
{
	return gsl_complex_mul( gsl_complex_mul(a,b), gsl_complex_mul(c,d) );
}
