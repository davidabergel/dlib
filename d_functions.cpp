#include "./d_functions.h"

double dlib::stepsize( double xmin, double xmax, int xpts )
{
	return (xmax-xmin)/(double)(xpts-1);
}
