#include "./d_functions.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace dlib;

double testintegrand( double x, void * );

int main()
{
	cout << "Running testing..." << endl;

	dIntParams dip = { 0.0, 1.0, 1001 };

	dIntegrator inter = dIntegrator();
	double intres = inter.Simpson( &testintegrand, NULL, &dip );

	const double mapleres= 0.2166666667;
	cout << "Got intres = " << intres 
		<< "   intres-maplere = " << intres-mapleres
		<< endl;

	return 0;
}

double testintegrand( double x, void * )
{
	return 0.1*pow(x,5) + 0.8*pow(x,3);
}

