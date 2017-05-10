#include "./d_functions.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace dlib;

double testintegrand( double x, void * );

int main()
{
	cout << "Running testing..." << endl;

	dIntParams dip = { 0.2, M_PI+0.2, 11 };

	dIntegrator inter = dIntegrator();
	cout << "Trapezoid:" << endl;
	double intTrap = inter.Trapezoid( &testintegrand, NULL, &dip );
	cout << endl;
	cout << "Simpson:";
	double intSimp = inter.Simpson( &testintegrand, NULL, &dip );
	cout << endl;

	//const double mapleres= M_PI*0.5;
	const double mapleres= 1.173457665;
	cout << "Got intTrap = " << intTrap << "  (error " << intTrap-mapleres 
		<< ")   intSimp = " << intSimp << "  (error " << intSimp-mapleres
		<< ")" << endl;

	return 0;
}

double testintegrand( double x, void * )
{
	cout << "  testintegrand(" << x << ")" << endl;
	return cos(x) - sin(2*x) + cos(3*x)*cos(3*x) - sin(6*x);
}

