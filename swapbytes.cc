/*
   file - swapbytes.cc
   author - Robert Best
   purpose - routines to do swapping of int's long's float's and doubles
 */
#include "swapbytes.h"

// shamelessly stolen from VMD
void reverse_int(int * N)
{
	char byteArray[4] = {'0','0','0','0'};
	char * bytePointer = (char*)N;

	byteArray[0]  =  *bytePointer;
	byteArray[1]  =  *(bytePointer+1);
	byteArray[2]  =  *(bytePointer+2);
	byteArray[3]  =  *(bytePointer+3);

	*bytePointer     = byteArray[3];
	*(bytePointer+1) = byteArray[2];
	*(bytePointer+2) = byteArray[1];
	*(bytePointer+3) = byteArray[0];
}

void reverse_long(long * N)
{
	char byteArray[8] = {'0','0','0','0','0','0','0','0'};
	char * bytePointer = (char*)N;

	byteArray[0]  =  *bytePointer;
	byteArray[1]  =  *(bytePointer+1);
	byteArray[2]  =  *(bytePointer+2);
	byteArray[3]  =  *(bytePointer+3);
	byteArray[4]  =  *(bytePointer+4);
	byteArray[5]  =  *(bytePointer+5);
	byteArray[6]  =  *(bytePointer+6);
	byteArray[7]  =  *(bytePointer+7);

	*bytePointer     = byteArray[7];
	*(bytePointer+1) = byteArray[6];
	*(bytePointer+2) = byteArray[5];
	*(bytePointer+3) = byteArray[4];
	*(bytePointer+4) = byteArray[3];
	*(bytePointer+5) = byteArray[2];
	*(bytePointer+6) = byteArray[1];
	*(bytePointer+7) = byteArray[0];
}

void reverse_float(float * N)
{
	char byteArray[4] = {'0','0','0','0'};
	char * bytePointer = (char*)N;

	byteArray[0]  =  *bytePointer;
	byteArray[1]  =  *(bytePointer+1);
	byteArray[2]  =  *(bytePointer+2);
	byteArray[3]  =  *(bytePointer+3);

	*bytePointer     = byteArray[3];
	*(bytePointer+1) = byteArray[2];
	*(bytePointer+2) = byteArray[1];
	*(bytePointer+3) = byteArray[0];
}

// also shamelessly stolen from VMD
void reverse_double(double * N)
{
	char byteArray[8] = {'0','0','0','0','0','0','0','0'};
	char * bytePointer = (char*)N;

	byteArray[0]  =  *bytePointer;
	byteArray[1]  =  *(bytePointer+1);
	byteArray[2]  =  *(bytePointer+2);
	byteArray[3]  =  *(bytePointer+3);
	byteArray[4]  =  *(bytePointer+4);
	byteArray[5]  =  *(bytePointer+5);
	byteArray[6]  =  *(bytePointer+6);
	byteArray[7]  =  *(bytePointer+7);

	*bytePointer     = byteArray[7];
	*(bytePointer+1) = byteArray[6];
	*(bytePointer+2) = byteArray[5];
	*(bytePointer+3) = byteArray[4];
	*(bytePointer+4) = byteArray[3];
	*(bytePointer+5) = byteArray[2];
	*(bytePointer+6) = byteArray[1];
	*(bytePointer+7) = byteArray[0];
}

void reverse_int_array( int array[], int size )
{
	for (int i=0; i<size; i++) {
		reverse_int( & array[i] );
	}
}

void reverse_long_array( long array[], int size )
{
	for (int i=0; i<size; i++) {
		reverse_long( & array[i] );
	}
}

void reverse_float_array( float array[], int size )
{
	for (int i=0; i<size; i++) {
		reverse_float( & array[i] );
	}
}

void reverse_double_array( double array[], int size )
{
	for (int i=0; i<size; i++) {
		reverse_double( & array[i] );
	}
}


int sw_float_fread( FILE * fname, float * val, bool swap_bytes )
{
	fread ( val, sizeof(float), 1, fname );
	if (swap_bytes)
		reverse_float(val);
	return 0;
}
	
int sw_int_fread( FILE * fname, int * val, bool swap_bytes )
{
	fread ( val, sizeof(int), 1, fname );
	if (swap_bytes)
		reverse_int(val);
	return 0;
}
