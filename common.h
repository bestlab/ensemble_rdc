
#include <string>
#include <iostream>
#include <cmath>

#ifndef _COMMON_H
#define _COMMON_H

using namespace std;

void print_err( const char * function_name, const char * message );

void print_err( const string & function_name, const string & message);

namespace linalg
{
	// left matrix*vector mult.
	void multiply_Ax ( const double mat[3][3], double & ax, double & ay, 
			double & az, double & bx, double & by, double & bz );

	// jacobi matrix inversion: stolen from vmd code - heh,heh!
	int jacobi( double a[3][3], double d[3], double v[3][3]);
	
	void normalise( double x[3] );
	void normalise( double & x, double & y, double & z );
	double modx ( double x[3] );
	double dotp ( double x[3], double y[3] );
	void crossp (double a[3], double b[3], double c[3]);
	double angle (double x[3], double y[3] );
	double det( double M[3][3] );
	void zero_mat( double M[3][3] );
	void ident_mat( double M[3][3] );
}

#endif

