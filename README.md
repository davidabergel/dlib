dlib

A small library of C++ functions that I regularly use.

Requires standard C++ libraries. (cmath, iostream)
Requires gsl: https://www.gnu.org/software/gsl/
Requires compilation with -fopenmp flag

d_constants: 
Definitions of physical constants and complex numbers.

d_functions:
stepsize: Returns stepsize given min, max, and npoints.
overloaded aritmetic operators for gsl_complex
gslc_sum and gslc_prod: sums and products of gsl_complex numbers
dIntegrator: A class for integrating functions that currently contains
trapezium and Simpson methods.
d2DIntegrator: A class for integrating functions of two variables.
Currently contains Simpson method.
TwoDVector: A class for two-dimensional real vectors.


d_gsloutput: 
Various pretty printing functions, especially of gsl_complex types.
