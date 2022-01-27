/*
   file - swapbytes.h
   author - Robert Best
   purpose - routines to do swapping of int's long's float's and doubles
 */

#ifndef _SWAPBYTES_H
#define _SWAPBYTES_H

#include <cstdio>
//#include "config.h"
//
using namespace std;

void reverse_int(int * N);
void reverse_long(long * N);
void reverse_float(float * N);
void reverse_double(double * N);
void reverse_int_array( int array[], int size );
void reverse_long_array( long array[], int size );
void reverse_float_array( float array[], int size );
void reverse_double_array( double array[], int size );
int sw_float_fread( FILE *, float * val, bool swap_bytes );
int sw_int_fread( FILE *, int * val, bool swap_bytes );

#endif			// _SWAPBYTES_H


