#include "./d_functions.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace dlib;

double testintegrand( double x, void * );

int main()
{
	cout << "Running testing..." << endl;

	dIntParams dip = { 0.0, M_PI, 11 };

	dIntegrator inter = dIntegrator();
	double intTrap = inter.Trapezoid( &testintegrand, NULL, &dip );
	double intSimp = inter.Simpson( &testintegrand, NULL, &dip );

	const double mapleres= 0.0;
	cout << "Got intTrap = " << intTrap << "  (error " << intTrap-mapleres 
		<< ")   intSimp = " << intSimp << "  (error " << intSimp-mapleres
		<< ")" << endl;

	return 0;
}

double testintegrand( double x, void * )
{
	return cos(x) + cos(3*x) + cos(5*x) + cos(6*x);
}

