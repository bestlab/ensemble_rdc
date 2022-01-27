/*
Author: Robert Best
Date: Tue Oct 19 14:21:24 SAST 1999
Purpose: Implementation of trajectory file i/o using interface defined
	by abstract class.
*/

/*
   bugs - output does not write the header you specify.

   changelog -

   Thu Feb 10 23:09:14 SAST 2000
   RB - added goto_frame function

   Thu Jan 27 09:19:12 SAST 2000
   RB - added automatic byte-flipping for DCD file reading. 
   data still written in native format

   Mon Jul  2 08:54:41 BST 2001
   RB - added TrajSet to allow a group of trajectory files
   to be grouped and handled as one.
 */

#include <cstdio>
#include <iostream>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include "traj.h"
#include "common.h"

#undef DEBUG		// we hope !
//#define DEBUG

#define XTC_MAGIC 1995
#define MAXABS INT_MAX-2
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x):(y))
#endif
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
static int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))


static int sizeofint(const int size) {
    unsigned int num = 1;
    int num_of_bits = 0;
    
    while (size >= num && num_of_bits < 32) {
	num_of_bits++;
	num <<= 1;
    }
    return num_of_bits;
}

static int sizeofints( const int num_of_ints, unsigned int sizes[]) {
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++) {	
	tmp = 0;
	for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
	    tmp = bytes[bytecnt] * sizes[i] + tmp;
	    bytes[bytecnt] = tmp & 0xff;
	    tmp >>= 8;
	}
	while (tmp != 0) {
	    bytes[bytecnt++] = tmp & 0xff;
	    tmp >>= 8;
	}
	num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
	num_of_bits++;
	num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;

}
    
static void sendbits(int buf[], int num_of_bits, int num) {
    
    unsigned int cnt, lastbyte;
    int lastbits;
    unsigned char * cbuf;
    
    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte =(unsigned int) buf[2];
    while (num_of_bits >= 8) {
	lastbyte = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
	cbuf[cnt++] = lastbyte >> lastbits;
	num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
	lastbyte = (lastbyte << num_of_bits) | num;
	lastbits += num_of_bits;
	if (lastbits >= 8) {
	    lastbits -= 8;
	    cbuf[cnt++] = lastbyte >> lastbits;
	}
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits>0) {
	cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}

static void sendints(int buf[], const int num_of_ints, const int num_of_bits,
	unsigned int sizes[], unsigned int nums[]) {

    int i;
    unsigned int bytes[32], num_of_bytes, bytecnt, tmp;

    tmp = nums[0];
    num_of_bytes = 0;
    do {
	bytes[num_of_bytes++] = tmp & 0xff;
	tmp >>= 8;
    } while (tmp != 0);

    for (i = 1; i < num_of_ints; i++) {
	if (nums[i] >= sizes[i]) {
	    fprintf(stderr,"major breakdown in sendints num %u doesn't "
		    "match size %u\n", nums[i], sizes[i]);
	    exit(1);
	}
	/* use one step multiply */    
	tmp = nums[i];
	for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
	    tmp = bytes[bytecnt] * sizes[i] + tmp;
	    bytes[bytecnt] = tmp & 0xff;
	    tmp >>= 8;
	}
	while (tmp != 0) {
	    bytes[bytecnt++] = tmp & 0xff;
	    tmp >>= 8;
	}
	num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8) {
	for (i = 0; i < num_of_bytes; i++) {
	    sendbits(buf, 8, bytes[i]);
	}
	sendbits(buf, num_of_bits - num_of_bytes * 8, 0);
    } else {
	for (i = 0; i < num_of_bytes-1; i++) {
	    sendbits(buf, 8, bytes[i]);
	}
	sendbits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}

static int receivebits(int buf[], int num_of_bits) {

    int cnt, num; 
    unsigned int lastbits, lastbyte;
    unsigned char * cbuf;
    int mask = (1 << num_of_bits) -1;

    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];
    
    num = 0;
    while (num_of_bits >= 8) {
	lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
	num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
	num_of_bits -=8;
    }
    if (num_of_bits > 0) {
	if (lastbits < num_of_bits) {
	    lastbits += 8;
	    lastbyte = (lastbyte << 8) | cbuf[cnt++];
	}
	lastbits -= num_of_bits;
	num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num; 
}

static void receiveints(int buf[], const int num_of_ints, int num_of_bits,
	unsigned int sizes[], int nums[]) {
    int bytes[32];
    int i, j, num_of_bytes, p, num;
    
    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
	bytes[num_of_bytes++] = receivebits(buf, 8);
	num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
	bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints-1; i > 0; i--) {
	num = 0;
	for (j = num_of_bytes-1; j >=0; j--) {
	    num = (num << 8) | bytes[j];
	    p = num / sizes[i];
	    bytes[j] = p;
	    num = num - p * sizes[i];
	}
	nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

int read_int(FILE *f, bool sfbit) {
	int32_t i32;
	int64_t i64;
	int x;
	size_t junk;
	if (sfbit) {
		junk=fread(&i64, sizeof(int64_t), 1, f);
		x = i64;
	} else {
		junk=fread(&i32, sizeof(int32_t), 1, f);
		x = i32;
	}
	return x;
}


float dotpn(float *a, float *b, int n)
{
	float d = 0.0;
	for (int i=0; i<n; i++) {
		d+=a[i]*b[i];
	}
	return d;
}

DCDITrajFile::DCDITrajFile()
{
	sfbit = false;
	N=0;
	file_open=0;
	header="";
	curr_frame=0;
	swap_bytes = 0;
}

DCDITrajFile::DCDITrajFile( const char * filename )
{
	sfbit = false;
	file_open = 0;
	swap_bytes = 0;
	open( filename );	
}

DCDITrajFile::~DCDITrajFile()
{
	close();
}

/*
   DCDITrajFile::open						
   ____________________________________________________________________ 
   function to open and initialize dcd file for input		
 */
int DCDITrajFile::open( const char * filename )
{
	int icount, all_atom_rec, free_atom_rec, CrystLen;
	int64_t *freeatoms64;
	const int charmm_string_length = 80;
	char tmpstr[ charmm_string_length + 1 ];
	char typestr[5];
	int tmp, ntitl;
	long int tmp_l;
	//bool sfbita;
//	int * tmp_ptr = & tmp;
	int32_t icntrl[20];		// mimic CHARMM data structure
	int64_t icntrl64[20];		// mimic CHARMM data structure
	sfbit = false;
	sfbita = false;
	if (file_open) {
		print_err( "DCDITrajFile::open", 
				"A trajectory is already open");
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	
	swap_bytes = 0;
	dcdfile = fopen ( filename, "r" );
	if ( dcdfile == NULL ) {
		print_err( "DCDITrajFile::open", 
				"could not open trajectory file");
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	
	//icount = fread( &tmp_l, sizeof(long int), 1, dcdfile ); 
	//fprintf(stdout,"tmp_l = %12li\n", tmp_l);
	// read mysterious first block
	icount = fread( &tmp, sizeof(int32_t), 1, dcdfile ); 
	//fprintf(stdout,"tmp = %12i\n", tmp);
	if (tmp == 4 + 20 * sizeof(int64_t)) {
		sfbit = true;
		sfbita = true;
	}
#ifdef DEBUG
	if (tmp == 4 + 20 * sizeof(int32_t)) {
		fprintf(stdout,"DCD file claims to have 32 bit integers\n");
	} else if (tmp == 4 + 20 * sizeof(int64_t)) {
		//sfbit = true;
		fprintf(stdout,"DCD file claims to have 64 bit integers\n");
	} else {
		fprintf(stdout,"DCD file has unknown integers!!!\n");
	}
#endif
	
	if ( (tmp != ( 4 + 20 * sizeof(int))) && !sfbit ) { 
			reverse_int( & tmp );
			if ( tmp == ( 4 + 20 * sizeof(int) ) ) {
#ifdef DEBUG
				cout << "Enabling on-the-fly byte swapping ..." 
					<< endl;
#endif
				swap_bytes = 1;
			} else {
				print_err( "DCDITrajFile::open", 
"A Mysterious Header error has occurred!\n\
This dcd file could be in the wrong binary format (i.e. byte-order).\n\
In this case you could use a utility like flipdcd to convert it.\n\
However, since this library is meant to handle this automatically,\n\
it is MOST LIKELY that the file you have specified is not actually\n\
a dcd file, or has been corrupted (e.g. by ascii ftp)."); 
#ifndef SWIG
				throw BinFormatErr();
#else
				exit(1);
#endif
			}
	}

	icount = fread( typestr, sizeof(char), 4, dcdfile ); // get file type 
	typestr[5] = '\0';
	if ( strncmp( typestr, "CORD", 4) == 0 ) {
		//fprintf(stdout,"COORDINATES\n");
		type = 'c';
	} else if ( strncmp( typestr, "VELD", 4) == 0 ) {
		//fprintf(stdout,"VELOCITIES\n");
		type = 'v';
	} else {
		// maybe lying about int size?
		fseek(dcdfile, 0L, SEEK_SET);
	        icount = fread( &tmp, sizeof(int64_t), 1, dcdfile ); 
	        icount = fread( typestr, sizeof(char), 4, dcdfile ); // get file type 
		if ( strncmp( typestr, "CORD", 4) == 0 ) {
			fprintf(stdout,"COORDINATES\n");
			type = 'c';
			sfbita = true;
			fprintf(stdout,"DCD file actually has 64 bit integers\n");
		} else if ( strncmp( typestr, "VELD", 4) == 0 ) {
			fprintf(stdout,"VELOCITIES\n");
			type = 'v';
			sfbita = true;
			fprintf(stdout,"DCD file actually has 64 bit integers\n");
		} else {
			fprintf(stderr,"Unknown file type: %s\n", typestr);
			print_err( "DCDITrajFile::open", 
					"Could not determine traj. file type:"
					" coordinate or velocity");
#ifndef SWIG
			throw FileFormatErr();
#else
			exit(1);
#endif
		}
	}

	if (sfbit) {
		icount = fread( icntrl64, sizeof(int64_t), 20, dcdfile ); // read 
	} else {
		icount = fread( icntrl, sizeof(int32_t), 20, dcdfile ); // read 
	}
	if ( swap_bytes ) 
		reverse_int_array( icntrl, 20 );
	if (sfbit) {
		for (int i=0; i<20; i++) {
			icntrl[i] = icntrl64[i];
		}
	}
	//CHARMM icntrl block

	// milk icntrl for all it's worth!______________________________
	nframes = icntrl[0];
	itime =	icntrl[1];
	if ( type == 'c' )
		nsave = icntrl[2];
	else if ( type == 'v' )
		nsave = icntrl[4];
	nsteps = icntrl[3];
	dof = icntrl[7];
	nfixed = icntrl[8];
	delta = static_cast<float>(icntrl[9]);	// A hack: beware!
	if (icntrl[10])
		qcrystal = 1;
	else 
		qcrystal = 0;
	version = icntrl[19];

#ifdef DEBUG
	cout << "nframes " << nframes << '\n' 
		<< "itime " << itime << '\n'
		<< "nsave " << nsave << '\n'
		<< "nsteps " << nsteps << '\n'
		<< "dof " << dof << '\n'
		<< "nfixed " << nfixed << '\n'
		<< "delta " << delta << '\n'
		<< "qcrystal " << qcrystal << '\n'
		<< "version " << version << endl;
#endif
	//___________________________________________________________________

	read_int(dcdfile,sfbita);
	read_int(dcdfile,sfbita);
	ntitl = read_int(dcdfile,false);
	//cout << " ntitl = " << ntitl << endl;
	if ( swap_bytes ) {
		reverse_int( & ntitl );
	}
	
	header = "";
	for ( int i = 0; i < ntitl; i++)
	{
		icount =  fread( tmpstr, charmm_string_length, 1, dcdfile );
		tmpstr[charmm_string_length] = '\0';
		header += tmpstr;
	}
	read_int(dcdfile,sfbita);
	read_int(dcdfile,sfbita);
	//read_int(dcdfile,false);
	N = read_int(dcdfile,sfbit);
	//cout << " natom = " << N << endl;
	if ( swap_bytes )
		reverse_int( & N );
	read_int(dcdfile,sfbita);

	//-------------- MEMORY ALLOCATION -------------------------
	/*
#ifdef DEBUG
	cout << "Trying to allocate coordinate space ..." << endl;
#endif
*/
	X = new float[N];
	Y = new float[N];
	Z = new float[N];
	/*
	cout << "Successfully allocated coordinate arrays!" << endl;
	*/
	if (nfixed != 0)
	{
		/*
#ifdef DEBUG
	cout << "Trying to allocate free atom space ..." << endl;
#endif
*/
		freeatoms = new int[N - nfixed]; // no free atoms = moron
		if (sfbit) {
			freeatoms64 = new int64_t[N-nfixed];
		}
		freex = new float[N - nfixed];
		freey = new float[N - nfixed];
		freez = new float[N - nfixed];
		/*
#ifdef DEBUG
	cout << "Successfully allocated free atom arrays!" << endl;
#endif
*/
	}
	//-----------------------------------------------------------
	
	// read in free-atom list if necessary
	if (nfixed != 0)
	{
		read_int(dcdfile,false);
		if (sfbit) {
			icount = fread(freeatoms64,sizeof(int64_t),N-nfixed,dcdfile);
			for (int i=0; i<N-nfixed; i++) {
				freeatoms[i] = freeatoms64[i];
			}
		} else {
			icount = fread(freeatoms,sizeof(int32_t),N-nfixed,dcdfile);
		}
		read_int(dcdfile,false);
		if ( swap_bytes )
			reverse_int_array( freeatoms, N - nfixed );
	}
	
	//--------------- set up dimensions -----------------------------------
	//if (sfbit) {
	//	all_atom_rec = 6 * sizeof(int64_t) + 3 * N * sizeof(float);
	//	free_atom_rec = 6 * sizeof(int64_t) + 3 * (N - nfixed) * sizeof(float);
	//} else {
		all_atom_rec = 6 * sizeof(int32_t) + 3 * N * sizeof(float);
		free_atom_rec = 6 * sizeof(int32_t) + 3 * (N - nfixed) * sizeof(float);
	//}
	if (qcrystal) 
	{
	//	if (sfbit) {
	//		CrystLen = 2 * sizeof(int64_t) + 6 * sizeof(double);
	//	} else {
			CrystLen = 2 * sizeof(int32_t) + 6 * sizeof(double);
	//	}
		all_atom_rec += CrystLen;
		free_atom_rec += CrystLen;
	}
	//---------------------------------------------------------------------

	start_record = ftell( dcdfile ); 	// current pos = beginning of records
	
	/* aside: options for seek start posn.
	   SEEK_SET = 0;
	   SEEK_CUR = 1;
	   SEEK_END = 2;
	 */
	// find the end of file = "0 offset from end"
	long EoF = fseek( dcdfile, 0L, SEEK_END );
	if (EoF) 
	{
		print_err( "DCDITrajFile::open", 
				"Can't read trajectory file: ");
		print_err( "DCDITrajFile::open", 
				filename);
#ifndef SWIG
		throw EofErr();
#else
		exit(1);
#endif
	}
#ifdef DEBUG
	cout << "\nstart_record " << start_record << '\n'
		<< "all_atom_rec " << all_atom_rec << '\n'
		<< "free_atom_rec " << free_atom_rec << '\n'
		<< "ftell(dcdfile) " << ftell(dcdfile) << endl;
#endif
	actual_nframes = (ftell( dcdfile ) - start_record - all_atom_rec)/free_atom_rec + 1;
	/*
#ifdef DEBUG
	cout << coorsets << nframes << endl;
#endif */
	/*
	if ( actual_nframes != nframes )
	{
		print_err( "DCDITrajFile::open",
			"Num of frames in file != num frames in header!");
		cerr << "continuing anyway...!" << endl;
		cerr << "Actual number of frames approx = " << coorsets << endl;
#ifndef SWIG
		//throw FileFormatErr();
#else

//		exit(1);
#endif
	}
	*/
	
	fseek( dcdfile, start_record, SEEK_SET );  // seek records again
	
	curr_frame=0;
	file_open = 1;

	return 0;
}

/*
   DCDITrajFile::read_frame()					
   ____________________________________________________________________ 
   pop next frame out of trajectory				
 */
void DCDITrajFile::read_frame( double * X_ret, double * Y_ret, double * Z_ret, int num_at) 
{
	/*
	int tmp;
	
	if (! flopen() )
	{
		print_err( "DCDITrajFile::read_frame", 
				"traj. file not open for reading");
#ifndef SWIG
		throw ReadErr();
#else
		exit(1);
#endif
	}
	
	if ( curr_frame >= nframes || curr_frame < 0)
	{
		print_err( "DCDITrajFile::read_frame", 
				"requested frame out of range");
#ifndef SWIG
		throw EofErr();
#else
		exit(1);
#endif
	}

	if ( num_at != N )
	{
		print_err( "DCDITrajFile::read_frame", 
			"number of atoms in trajectory != number of atoms asked for");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}

	if (qcrystal)
	{
		read_int(dcdfile,false);
		fread( crystal_data, sizeof(double), 6, dcdfile);
		if ( swap_bytes )
			reverse_double_array( crystal_data, 6 );
		read_int(dcdfile,false);
	}

	// first frame always has all coords; if no fixed atoms, then all
	// coords present.
	if ( curr_frame == 0 || nfixed == 0 ) 
	{
		read_int(dcdfile,false);
		fread( X, sizeof(float), N, dcdfile );
		read_int(dcdfile,false);
		read_int(dcdfile,false);
		fread( Y, sizeof(float), N, dcdfile );
		read_int(dcdfile,false);
		read_int(dcdfile,false);
		fread( Z, sizeof(float), N, dcdfile );
		read_int(dcdfile,false);
		if ( swap_bytes ) {
			reverse_float_array( X, N );
			reverse_float_array( Y, N );
			reverse_float_array( Z, N );
		}
	}
	else if ( nfixed > 0 ) // only update fixed atoms
	{
		read_int(dcdfile,false);
		fread( freex, sizeof(float), N-nfixed, dcdfile );
		read_int(dcdfile,false);
		read_int(dcdfile,false);
		fread( freey, sizeof(float), N-nfixed, dcdfile );
		read_int(dcdfile,false);
		read_int(dcdfile,false);
		fread( freez, sizeof(float), N-nfixed, dcdfile );
		read_int(dcdfile,false);
		if ( swap_bytes ) {
			reverse_float_array( freex, N-nfixed );
			reverse_float_array( freey, N-nfixed );
			reverse_float_array( freez, N-nfixed );
		}

		for (int i=0; i<(N-nfixed); i++)
		{
			X[freeatoms[i]-1] = freex[i]; // NB FORTRAN array correction
			Y[freeatoms[i]-1] = freey[i];
			Z[freeatoms[i]-1] = freez[i]; 
		}
	}
	*/

	// read into internal X,Y,Z array
	read_frame(X,Y,Z,num_at);
	// convert: float -> double for return
	for (int i=0; i<N; i++)
	{
		X_ret[i] = X[i]; 
		Y_ret[i] = Y[i]; 
		Z_ret[i] = Z[i];
	}

	//curr_frame++;

}

/*
   DCDITrajFile::read_frame()					
   ____________________________________________________________________ 
   pop next frame out of trajectory				
 */
void DCDITrajFile::read_frame( float *XX, float *YY, float *ZZ, int num_at) 
{
	int tmp;
	size_t junk;
	
	if (! flopen() )
	{
		print_err( "DCDITrajFile::read_frame", 
				"traj. file not open for reading");
#ifndef SWIG
		throw ReadErr();
#else
		exit(1);
#endif
	}
	
	if ( curr_frame >= nframes || curr_frame < 0)
	{
		print_err( "DCDITrajFile::read_frame", 
				"requested frame out of range");
#ifndef SWIG
		throw EofErr();
#else
		exit(1);
#endif
	}

	if ( num_at != N )
	{
		print_err( "DCDITrajFile::read_frame", 
			"number of atoms in trajectory != number of atoms asked for");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}

	if (qcrystal)
	{
		read_int(dcdfile,sfbita);
		junk=fread( crystal_data, sizeof(double), 6, dcdfile);
		if ( swap_bytes )
			reverse_double_array( crystal_data, 6 );
		read_int(dcdfile,sfbita);
	}

	// first frame always has all coords; if no fixed atoms, then all
	// coords present.
	if ( curr_frame == 0 || nfixed == 0 ) 
	{
		read_int(dcdfile,sfbita);
		junk=fread( XX, sizeof(float), N, dcdfile );
		read_int(dcdfile,sfbita);
		read_int(dcdfile,sfbita);
		junk=fread( YY, sizeof(float), N, dcdfile );
		read_int(dcdfile,sfbita);
		read_int(dcdfile,sfbita);
		junk=fread( ZZ, sizeof(float), N, dcdfile );
		read_int(dcdfile,sfbita);
		if ( swap_bytes ) {
			reverse_float_array( XX, N );
			reverse_float_array( YY, N );
			reverse_float_array( ZZ, N );
		}
	}
	else if ( nfixed > 0 ) // only update fixed atoms
	{
		read_int(dcdfile,sfbita);
		junk=fread( freex, sizeof(float), N-nfixed, dcdfile );
		read_int(dcdfile,sfbita);
		read_int(dcdfile,sfbita);
		junk=fread( freey, sizeof(float), N-nfixed, dcdfile );
		read_int(dcdfile,sfbita);
		read_int(dcdfile,sfbita);
		junk=fread( freez, sizeof(float), N-nfixed, dcdfile );
		read_int(dcdfile,sfbita);
		if ( swap_bytes ) {
			reverse_float_array( freex, N-nfixed );
			reverse_float_array( freey, N-nfixed );
			reverse_float_array( freez, N-nfixed );
		}

		for (int i=0; i<(N-nfixed); i++)
		{
			XX[freeatoms[i]-1] = freex[i]; // NB FORTRAN array correction
			YY[freeatoms[i]-1] = freey[i];
			ZZ[freeatoms[i]-1] = freez[i]; 
		}
	}
	//fprintf(stdout,"X[0] = %12.6f\n",XX[0]);
	curr_frame++;
}




void DCDITrajFile::close()
{
	if (file_open == 1)
	{
		fclose( dcdfile );
		/*
#ifdef DEBUG
		cout << "Trying to remove coordinate memory alloc ..." << endl;
#endif
*/
		delete [] X;
		delete [] Y;
		delete [] Z;
		/*
#ifdef DEBUG
		cout << "Successfully removed coordinate memory alloc!" << endl;
#endif */
		if ( nfixed > 0 )
		{
			/*
#ifdef DEBUG
			cout << "Trying to remove free atom memory alloc ..." << endl;
#endif
*/
			delete [] freex;
			delete [] freey;
			delete [] freez;
			delete [] freeatoms;
			/*
#ifdef DEBUG
			cout << "Successfully removed free atom memory alloc!" << endl;
#endif
*/
		}
		file_open = 0;
	}
}

/*
   DCDITrajFile::frames_left()
   ____________________________________________________________________
   test whether there are frames left in the file (handy for loops)
 */
bool DCDITrajFile::frames_left() const
{
	if ( curr_frame >= nframes || curr_frame < 0)
		return 0;
	else 
		return 1;
}

/*
   DCDITrajFile::skip()						
   ____________________________________________________________________ 
   skip i frames forwards (i positive)				
 */
void DCDITrajFile::skip( int i )
{
	int all_atom_rec = 6 * sizeof(int) + 3 * N * sizeof(float);
	int free_atom_rec = 6 * sizeof(int) + 3 * (N - nfixed) * sizeof(float);
	if (qcrystal) 
	{
		int CrystLen = 2 * sizeof(int) + 6 * sizeof(double);
		all_atom_rec += CrystLen;
		free_atom_rec += CrystLen;
	}
	long int offset;
	if (i < 0)
	{
		print_err( "DCDITrajFile::skip", 
				"Error: can't skip backwards(yet!)");
#ifndef SWIG
		throw IllegalRequestErr();
#else
		exit(1);
#endif
	}
	else if (i==0)
		return; 	// our work is done!

	if ((curr_frame + i) >= nframes)
	{
		print_err( "DCDITrajFile::skip", 
				"Error: trying to skip beyond end of file");
#ifndef SWIG
		throw EofErr();
#else
		exit(1);
#endif
	}

	if (curr_frame == 0)
	offset = all_atom_rec + (i-1)*free_atom_rec;
	else
	offset = free_atom_rec * i;
	
	fseek ( dcdfile, offset, SEEK_CUR );
	curr_frame += i;
}

/************************************************************************
 *	DCDITrajFile::rewind()						*
 * ____________________________________________________________________ *
 * 	go to beginning of file						*
 ************************************************************************/
void DCDITrajFile::rewind()
{
	fseek( dcdfile, start_record, SEEK_SET );  // seek records again
	curr_frame = 0;
}

void DCDITrajFile::goto_frame( int i )
{
	rewind();
	skip(i);
}

/************************************************************************
 *	DCDITrajFile::get_crystal_data()				*
 * ____________________________________________________________________ *
 * 	get the "XTLABC" 'shape matrix' written by charmm		*
 ************************************************************************/
void DCDITrajFile::get_crystal_data( double rdata[6] ) const
{
	for (int i=0; i<6; i++)
		rdata[i] = crystal_data[i];
}

/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------
  									
  			ITRAJSET STUFF				
  								
  -----------------------------------------------------------------------
  ------------------------------------------------------------------------*/

/*
ITrajSet::ITrajSet( int n )
{
	N=0;
	ntrajfile=0;
	header="";
	trajfiles.resize(n);
	file_open=0;
	currentfile=-1; currentframe=-1;
}

ITrajSet::ITrajSet( BaseITrajFile * tfiles, int ntfile )
{
	lbounds = new int[ntfile];
	flengths = new int[ntfile];
	setup( BaseITrajFile * tfiles, int ntfile );
}

void ITrajSet::setup( BaseITrajFile * tfiles, int ntfile )
{
	file_open = 0;
	trajfiles = tfiles;
	ntrajfiles = ntfile;
	lbounds = new int[ntfile];
	flengths = new int[ntfile];
	// check files all OK
	int itime = tfiles[0].initial_step();
	delta = tfiles[0].frame_freq();
	nfixed = tfiles[0].fixed_atoms();

	for (int i=0; i<ntfile; i++) {
		if ( ! tfiles[i].flopen() ) {
			throw ReadErr();
		} else {
			lbounds[i] = 
		}
	}
	currentfile = 0;
	currentframe = 0;

}


DCDITrajFile::~DCDITrajFile()
{
	close();
}

//------------------------------------------------------------------------
//	DCDITrajFile::open
//
// 	function to open and initialize dcd file for input
//	apologies to leif laaksonen!
//------------------------------------------------------------------------
int DCDITrajFile::open( const char * filename )
{
	int icount;
	const int charmm_string_length = 80;
	char tmpstr[ charmm_string_length + 1 ];
	char typestr[5];
	int tmp, ntitl;
	int icntrl[20];		// mimic CHARMM data structure
	if (file_open)
	{
		cerr << "\nMessage from DCDOTrajFile::open() : "
			<< "\nA trajectory is already open" << endl;
		throw OpenErr();
	}
	
	dcdfile = fopen ( filename, "r" );
	if ( dcdfile == NULL )
	{
		cerr << "\nMessage from DCDOTrajFile::open() : "
			<< "\ncould not open trajectory file" << filename
			<< " for input\n" << endl;
		throw OpenErr();
	}
	
	// read mysterious first block
	icount = fread( &tmp, sizeof(int), 1, dcdfile ); 
	
	if ( tmp != ( 4 + 20 * sizeof(int) ) )  // usu. = 84
	{
		cerr << "\nA Mysterious Header error has occurred!"
			<< "\nThis dcd file is probably in the wrong"
			<< "\nbinary format (i.e. byte-order)."
			<< "\nUse a utility like flipdcd to convert it." 
			<< endl;
		exit(1);
	}
	
	icount = fread( typestr, sizeof(char), 4, dcdfile );// get file type 
	typestr[5] = '\0';
	if ( strncmp( typestr, "CORD", 4) == 0 )
		type = 'c';
	else if ( strncmp( typestr, "VELD", 4) == 0 )
		type = 'v';
	else
	{
		cerr << "\nCould not determine traj. file type:"
		       	"coordinate or velocity\n"
			<< endl;
		exit(1);
	}
			
	icount = fread( icntrl, sizeof(int), 20, dcdfile ); // read 
							//CHARMM icntrl block
	
	// milk icntrl for all it's worth!______________________________
	nframes = icntrl[0];
	itime =	icntrl[1];
	if ( type == 'c' )
		nsave = icntrl[2];
	else if ( type == 'v' )
		nsave = icntrl[4];
	nsteps = icntrl[3];
	dof = icntrl[7];
	nfixed = icntrl[8];
	delta = static_cast<float>(icntrl[9]);	// A hack: beware!
	if (icntrl[10])
		qcrystal = 1;
	else 
		qcrystal = 0;
	version = icntrl[19];
	
#ifdef DEBUG
	cout << "nframes " << nframes << '\n' 
		<< "itime " << itime << '\n'
		<< "nsave " << nsave << '\n'
		<< "nsteps " << nsteps << '\n'
		<< "dof " << dof << '\n'
		<< "nfixed " << nfixed << '\n'
		<< "delta " << delta << '\n'
		<< "qcrystal " << qcrystal << '\n'
		<< "version " << version << endl;
#endif
	//___________________________________________________________________

	icount = fread( &tmp, sizeof(int), 1, dcdfile );
	icount = fread( &tmp, sizeof(int), 1, dcdfile );
	icount = fread( &ntitl, sizeof(int), 1, dcdfile );
	
	header = "";
	for ( int i = 0; i < ntitl; i++)
	{
		icount =  fread( tmpstr, charmm_string_length, 1, dcdfile );
		tmpstr[charmm_string_length] = '\0';
		header += tmpstr;
	}
	
	icount = fread( &tmp, sizeof(int), 1, dcdfile );
	icount = fread( &tmp, sizeof(int), 1, dcdfile );

	icount = fread( &N, sizeof(int), 1, dcdfile );
	icount = fread( &tmp, sizeof(int), 1, dcdfile );

	//-------------- MEMORY ALLOCATION -------------------------
#ifdef DEBUG
	cout << "About to allocate coordinate space ..." << endl;
#endif
	X = new float[N];
	Y = new float[N];
	Z = new float[N];
	cout << "Successfully allocated coordinate arrays!" << endl;
	if (nfixed != 0)
	{
#ifdef DEBUG
	cout << "About to allocate free atom space ..." << endl;
#endif
		freeatoms = new int[N - nfixed]; // no free atoms = moron
		freex = new float[N - nfixed];
		freey = new float[N - nfixed];
		freez = new float[N - nfixed];
#ifdef DEBUG
	cout << "Successfully allocated free atom arrays!" << endl;
#endif
	}
	//-----------------------------------------------------------
	
	// read in free-atom list if necessary
	if (nfixed != 0)
	{
		icount = fread( &tmp, sizeof(int), 1, dcdfile );
		icount = fread( freeatoms, sizeof(int), N - nfixed, dcdfile );
		icount = fread( &tmp, sizeof(int), 1, dcdfile );
	}
	
	//--------------- set up dimensions -----------------------------------
	int all_atom_rec = 6 * sizeof(int) + 3 * N * sizeof(float);
	int free_atom_rec = 6 * sizeof(int) + 3 * (N - nfixed) * sizeof(float);
	if (qcrystal) 
	{
		int CrystLen = 2 * sizeof(int) + 6 * sizeof(double);
		all_atom_rec += CrystLen;
		free_atom_rec += CrystLen;
	}
	//---------------------------------------------------------------------

	start_record = ftell( dcdfile ); 	// current pos = beginning of records
	
	// aside: options for seek start posn.
	//   SEEK_SET = 0;
	//   SEEK_CUR = 1;
	//   SEEK_END = 2;
	// find the end of file = "0 offset from end"
	long EoF = fseek( dcdfile, 0L, SEEK_END );
	if (EoF) 
	{
		cerr << "\nCan't read trajectory file: " << filename << endl;
		exit(1);
	}
#ifdef DEBUG
	cout << "\nstart_record " << start_record << '\n'
		<< "all_atom_rec " << all_atom_rec << '\n'
		<< "free_atom_rec " << free_atom_rec << '\n'
		<< "ftell(dcdfile) " << ftell(dcdfile) << endl;
#endif
	int coorsets = (ftell( dcdfile ) - start_record - all_atom_rec)/free_atom_rec + 1;
#ifdef DEBUG
	cout << coorsets << nframes << endl;
#endif
	if ( coorsets != nframes )
	{
		cerr << "\nNum of frames in file != num frames in header!\n" << endl;
		exit(1);
	}
	
	fseek( dcdfile, start_record, SEEK_SET );  // seek records again
	
	curr_frame=0;
	file_open = 1;
}

//------------------------------------------------------------------------
//	DCDITrajFile::read_frame
//
// 	pop next frame out of trajectory
//------------------------------------------------------------------------
void DCDITrajFile::read_frame( double * X_ret, double * Y_ret, double * Z_ret, int num_at) 
{
	int tmp;
	
#ifdef DEBUG
	cout << "About to read frame " << curr_frame << endl;
#endif
	if (! flopen() )
	{
		cerr << "\ntraj. file not open for reading\n" << endl;
		exit(1);
	}
	
	if ( curr_frame >= nframes || curr_frame < 0)
	{
		cerr << "\nrequested frame out of range\n" << endl;
		exit(1);
	}

	if ( num_at != N )
	{
		cerr << "\nnumber of atoms in trajectory != number of atoms asked for\n"
			<< endl;
		exit(1);
	}

	if (qcrystal)
	{
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( crystal_data, sizeof(double), 6, dcdfile);
		fread( &tmp, sizeof(int), 1, dcdfile );
#ifdef DEBUG
		cout << "Crystal Data for this frame is: " << endl;
		for (int i = 0; i<6; i++)
		{
			cout << crystal_data[i] << endl;
		}
#endif
	}

	// first frame always has all coords; if no fixed atoms, then all
	// coords present.
	if ( curr_frame == 0 || nfixed == 0 ) 
	{
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( X, sizeof(float), N, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( Y, sizeof(float), N, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( Z, sizeof(float), N, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
	}
	else if ( nfixed > 0 ) // only update fixed atoms
	{
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( freex, sizeof(float), N-nfixed, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( freey, sizeof(float), N-nfixed, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );
		fread( freez, sizeof(float), N-nfixed, dcdfile );
		fread( &tmp, sizeof(int), 1, dcdfile );

		for (int i=0; i<(N-nfixed); i++)
		{
			X[freeatoms[i]-1] = freex[i]; // NB FORTRAN array correction
			Y[freeatoms[i]-1] = freey[i];
			Z[freeatoms[i]-1] = freez[i]; 
		}
	}
		

	// convert: float -> double for return
	for (int i=0; i<N; i++)
	{
		X_ret[i] = X[i]; 
		Y_ret[i] = Y[i]; 
		Z_ret[i] = Z[i];
	}

#ifdef DEBUG
	cout << "Successfully read frame " << curr_frame << endl;
#endif
	curr_frame++;

}

void DCDITrajFile::close()
{
	if (file_open == 1)
	{
		fclose( dcdfile );
#ifdef DEBUG
		cout << "About to remove coordinate memory alloc ..." << endl;
#endif
		delete [] X;
		delete [] Y;
		delete [] Z;
#ifdef DEBUG
		cout << "Successfully removed coordinate memory alloc!" << endl;
#endif
		if ( nfixed > 0 )
		{
#ifdef DEBUG
			cout << "About to remove free atom memory alloc ..." << endl;
#endif
			delete [] freex;
			delete [] freey;
			delete [] freez;
			delete [] freeatoms;
#ifdef DEBUG
			cout << "Successfully removed free atom memory alloc!" << endl;
#endif
		}
		file_open = 0;
	}
}

void DCDITrajFile::skip( int i )
{
	int all_atom_rec = 6 * sizeof(int) + 3 * N * sizeof(float);
	int free_atom_rec = 6 * sizeof(int) + 3 * (N - nfixed) * sizeof(float);
	if (qcrystal) 
	{
		int CrystLen = 2 * sizeof(int) + 6 * sizeof(double);
		all_atom_rec += CrystLen;
		free_atom_rec += CrystLen;
	}
	long int offset;
	if (i < 0)
	{
		cerr << "\nError: can't skip backwards(yet!)\n" << endl;
		exit(1);
	}
	else if (i==0)
		return; 	// our work is done!

	if ((curr_frame + i) >= nframes)
	{
		cerr << "\nError: trying to skip beyond end of file\n" << endl;
		exit(1);
	}

	if (curr_frame == 0)
	offset = all_atom_rec + (i-1)*free_atom_rec;
	else
	offset = free_atom_rec * i;
	
	fseek ( dcdfile, offset, SEEK_CUR );
	curr_frame += i;
}

void DCDITrajFile::rewind()
{
	fseek( dcdfile, start_record, SEEK_SET );  // seek records again
	curr_frame = 0;
}

// hack to see what this "XTLABC" 'shape matrix' is all about
void DCDITrajFile::get_crystal_data( double rdata[6] ) const
{
	for (int i=0; i<6; i++)
		rdata[i] = crystal_data[i];
}

*/

/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------
  									
 			DCDOTRAJFILE STUFF
  								
  -----------------------------------------------------------------------
  ------------------------------------------------------------------------*/

void pad(char *s, int len)
{
	int curlen;
	int i;

	curlen=strlen(s);

	if (curlen>len)
	{
		s[len]=0;
		return;
	}

	for (i=curlen; i<len; i++)
	{
		s[i]=' ';
	}

	s[i]=0;
}

/************************************************************************
 *	DCDOTrajFile::DCDOTrajFile()					*
 * ____________________________________________________________________ *
 * 	Default constructor						*
 * 	This will set all header values to sensible defaults, so under	*
 *	most circumstances, the user will only have to specify		*
 *	the number of frames, the number of atoms and the name of the	*
 *	file to open. The main exception to this is if there is crystal	* 
 *	information							*
 ************************************************************************/
DCDOTrajFile::DCDOTrajFile()
{
	type = 'c';
	nframes = 0;
	itime = 0;
	nsave = 50;
	nsteps = 0;
	dof = 0;
	nfixed = 0;
	delta = 0.2045;		// 1 fs in AKMA units (see usage.doc)
	qcrystal = 0;
	version = 27;
	N = 0;
	dcdfile = NULL;
	start_record = 0;
	file_open = 0;
	header = "";
	curr_frame = 0;
	X = Y = Z = freex = freey = freez = NULL;
	freeatoms = NULL;
}

/************************************************************************
 *	DCDOTrajFile::DCDOTrajFile()					*
 * ____________________________________________________________________ *
 * 	open dcd file trajectory for output				*
 ************************************************************************/
DCDOTrajFile::DCDOTrajFile( const char * filename )
{
	type = 'c';
	nframes = 0;
	itime = 0;
	nsave = 50;
	nsteps = 0;
	dof = 0;
	nfixed = 0;
	delta = 0.2045;		// 1 fs in AKMA units (see usage.doc)
	qcrystal = 0;
	version = 27;
	N = 0;
	dcdfile = NULL;
	start_record = 0;
	file_open = 0;
	header = "";
	curr_frame = 0;
	X = Y = Z = freex = freey = freez = NULL;
	freeatoms = NULL;
	open( filename );	
}

/************************************************************************
 *	DCDOTrajFile::DCDOTrajFile()					*
 * ____________________________________________________________________ *
 * 	construct output trajectory based on input			*
 ************************************************************************/
DCDOTrajFile::DCDOTrajFile( const DCDITrajFile & inptraj )
{
	setup( inptraj );
}

/************************************************************************
 *	DCDOTrajFile::~DCDOTrajFile()					*
 * ____________________________________________________________________ *
 * 	Die								*
 ************************************************************************/
DCDOTrajFile::~DCDOTrajFile()
{
	close();
}

/************************************************************************
 *	DCDOTrajFile::setup()						*
 * ____________________________________________________________________ *
 * 	setup parameters for output trajectory based on input		*
 ************************************************************************/
void DCDOTrajFile::setup( const DCDITrajFile & inptraj)
{
	type = inptraj.traj_type();
	nframes = inptraj.total_frames();
	itime = inptraj.initial_step();
	nsave = inptraj.frame_freq();
	nsteps = inptraj.total_steps();
	dof = inptraj.deg_free();
	nfixed = inptraj.fixed_atoms();
	delta = inptraj.step_size();
	qcrystal = inptraj.crystal();
	version = 27;
	N = inptraj.num_atoms();
	inptraj.read_title(header);
	X = new float[N];
	Y = new float[N];
	Z = new float[N];
	if (nfixed != 0)
	{
		freex = new float[N-nfixed];
		freey = new float[N-nfixed];
		freez = new float[N-nfixed];
		freeatoms = new int[N-nfixed];
	}
	
}

/************************************************************************
 *	DCDOTrajFile::open()						*
 * ____________________________________________________________________ *
 * 	open a dcd file for output					*
 ************************************************************************/
int DCDOTrajFile::open( const char * filename )
{
	if (file_open)
	{
		print_err( "DCDOTrajFile::open", 
				"A output trajectory is already open"); 
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	
	dcdfile = fopen ( filename, "w" );
	if ( dcdfile == NULL )
	{
		print_err( "DCDOTrajFile::open", 
				"Could not open trajectory file for output:");
		print_err( "DCDOTrajFile::open", filename);
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	file_open = 1;
	header_written=0;
	curr_frame = 0;
	return 0;
}

void DCDOTrajFile::close()
{
	if (file_open == 1)
	{
		if (curr_frame != nframes) {
			rewrite_frames_and_steps( curr_frame );
		}
		fclose( dcdfile );
#ifdef DEBUG
		cout << "Trying to remove coordinate memory alloc ..." << endl;
#endif
		delete [] X;
		delete [] Y;
		delete [] Z;
#ifdef DEBUG
		cout << "Successfully removed coordinate memory alloc!" << endl;
#endif
		if ( nfixed > 0 )
		{
#ifdef DEBUG
			cout << "Trying to remove free atom memory alloc ..." << endl;
#endif
			delete [] freex;
			delete [] freey;
			delete [] freez;
			delete [] freeatoms;
#ifdef DEBUG
			cout << "Successfully removed free atom memory alloc!" << endl;
#endif
		}
		file_open = 0;
	}
}

/************************************************************************
 *	DCDOTrajFile::write_header()					*
 * ____________________________________________________________________ *
 * 	write dcd header; requires that necessary info have previously	*
 *	been initialised.						*
 ************************************************************************/
void DCDOTrajFile::write_header()
{
	int icount;
	const int charmm_string_length = 80;
	char tmpstr[ charmm_string_length + 1 ] = "* A HEADER WRITTEN BY ROB'S PROGRAM";
	char typestr[5];
	int tmp, ntitl;
	int icntrl[20];		// mimic CHARMM data structure
	if (! file_open)
	{
		print_err( "DCDOTrajFile::write_header()", 
		"Error: trying to write to unopened file");
#ifndef SWIG
		throw WriteErr();
#else
		exit(1);
#endif
	}
	
	// write mysterious first block
	tmp = 4 + 20 * sizeof(int);
	fwrite( &tmp, sizeof(int), 1, dcdfile ); 
	
	if (type == 'c')
		strcpy (typestr, "CORD");
	else if (type == 'v')
		strcpy (typestr, "VELD");
	else 
	{
		print_err( "DCDOTrajFile::write_header()", 
				"Invalid trajectory file type set:\
				must be coordinate or velocity");
#ifndef SWIG
		throw InternalDataErr();
#else
		exit(1);
#endif
	} 
	fwrite( typestr, sizeof(char), 4, dcdfile );

	// initialize icntrl
	for ( int i = 0; i< 20; i++)	
		icntrl[i] = 0;
	icntrl[0] = nframes;
	icntrl[1] = itime;
	if ( type == 'c' )
		icntrl[2] = nsave;
	else
		icntrl[4] = nsave;
	icntrl[3] = nsteps;
	icntrl[7] = dof;
	icntrl[8] = nfixed;
	icntrl[9] = static_cast<int>(delta);
	icntrl[10] = qcrystal;
	icntrl[19] = version;
	// write icntrl
	fwrite( icntrl, sizeof(int), 20, dcdfile ); 
	
	//___________________________________________________________________

	tmp = 84;
	fwrite( &tmp, sizeof(int), 1, dcdfile );
	fwrite( &tmp, sizeof(int), 1, dcdfile );
	tmp = 1;
	fwrite( &tmp, sizeof(int), 1, dcdfile );
	pad ( tmpstr, 80 );
	fwrite( tmpstr, charmm_string_length, 1, dcdfile );
	
	tmp = 84;
	fwrite( &tmp, sizeof(int), 1, dcdfile );
	tmp = 4;
	fwrite( &tmp, sizeof(int), 1, dcdfile );
	tmp = N;
	fwrite( &tmp, sizeof(int), 1, dcdfile );
	tmp = 4;
	fwrite( &tmp, sizeof(int), 1, dcdfile );

	// write out free-atom list if necessary
	// not complete yet
	if (nfixed != 0)
	{
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( freeatoms, sizeof(int), N - nfixed, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
	}
	
	//--------------- set up dimensions -----------------------------------
	int all_atom_rec = 6 * sizeof(int) + 3 * N * sizeof(float);
	int free_atom_rec = 6 * sizeof(int) + 3 * (N - nfixed) * sizeof(float);
	if (qcrystal) 
	{
		int CrystLen = 2 * sizeof(int) + 6 * sizeof(double);
		all_atom_rec += CrystLen;
		free_atom_rec += CrystLen;
	}
	//---------------------------------------------------------------------

	start_record = ftell( dcdfile ); 	// current pos = beginning of records
	header_written = 1;
	
}

/************************************************************************
 *	DCDOTrajFile::write_frame					*
 * ____________________________________________________________________ *
 * 	pop next frame out of trajectory				*
 ************************************************************************/
void DCDOTrajFile::write_frame( double * XX, double * YY, double * ZZ, int num_at) 
{
	// convert: double -> float for writing
	for (int i=0; i<N; i++)
	{
		X[i] = XX[i]; 
		Y[i] = YY[i]; 
		Z[i] = ZZ[i];
	}
	write_frame(X,Y,Z,num_at);
}


void DCDOTrajFile::write_frame( float * XX, float * YY, float * ZZ, int num_at) 
{
	int tmp;
	
	if (! file_open )
	{
		print_err( "DCDOTrajFile::write_frame()",
				"Traj. file not open for writing");
#ifndef SWIG
		throw WriteErr();
#else
		exit(1);
#endif
	}
	if (! header_written )
	{
		print_err( "DCDOTrajFile::write_frame()",
				"header not yet written for trajectory");
#ifndef SWIG
		throw NoHeaderWrittenErr();
#else
		exit(1);
#endif
	}
	if ( curr_frame >= nframes || curr_frame < 0)
	{
		print_err( "DCDOTrajFile::write_frame()",
		"You've already written quite enough frames!");
#ifndef SWIG
		throw EofErr();
#else
		exit(1);
#endif
	}
	if ( num_at != N )
	{
		print_err( "DCDOTrajFile::write_frame()",
		"number of atoms in trajectory != number of atoms asked for");
#ifndef SWIG
		throw InternalDataErr();
#else
		exit(1);
#endif
	}

	if (qcrystal)
	{
		tmp = 84;
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( crystal_data, sizeof(double), 6, dcdfile);
		fwrite( &tmp, sizeof(int), 1, dcdfile );
	}

	// convert: double -> float for writing
	/*
	for (int i=0; i<N; i++)
	{
		X[i] = X_val[i]; 
		Y[i] = Y_val[i]; 
		Z[i] = Z_val[i];
	}
	*/
	// first frame always has all coords; if no fixed atoms, then all
	// coords present.
	if ( curr_frame == 0 || nfixed == 0 ) 
	{
		//tmp = 80 + sizeof(float) * N;
		tmp = sizeof(float) * N;
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( XX, sizeof(float), N, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( YY, sizeof(float), N, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( ZZ, sizeof(float), N, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
	}
	else if ( nfixed > 0 ) // only update fixed atoms
	{
		for (int i=0; i<(N-nfixed); i++)
		{
			freex[i] = XX[freeatoms[i]-1]; // NB FORTRAN array correction
			freey[i] = YY[freeatoms[i]-1];
			freez[i] = ZZ[freeatoms[i]-1]; 
		}
		//tmp = 80 + sizeof(float) * (N-nfixed);
		tmp = sizeof(float) * (N-nfixed);
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( freex, sizeof(float), N-nfixed, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( freey, sizeof(float), N-nfixed, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );
		fwrite( freez, sizeof(float), N-nfixed, dcdfile );
		fwrite( &tmp, sizeof(int), 1, dcdfile );

	}
		
	curr_frame++;

}

/*
 * set crystal data on DCDOTrajFile to specified array
 */
void DCDOTrajFile::set_crystal_data(double *xtaldata)
{
	for (int i=0; i<6; i++) 
		crystal_data[i]=xtaldata[i];
}

/************************************************************************
 *		REWRITE NUMBER OF FRAMES IN HEADER			*
 * -------------------------------------------------------------------- *
 *	Allows total number of frames to be rewritten in header,	*
 *	even *after* frames have been written in the trajectory.	*
 ************************************************************************/
void DCDOTrajFile::rewrite_frames_and_steps( int f ) 
{
	int tmp =f;
	if (curr_frame > f)
	{
		print_err( "DCDOTrajFile::write_frame()",
				"Already more frames than that in the file");
#ifndef SWIG
		throw EofErr();
#else
		exit(1);
#endif
	}

	long curr_offset = ftell( dcdfile );
	long offset = sizeof(int) + 4*sizeof(char);
	fseek ( dcdfile, offset, SEEK_SET );	// get to start of icntrl
	fwrite ( &tmp, sizeof(int), 1, dcdfile );	// write f
	fseek ( dcdfile, 2*sizeof(int), SEEK_CUR);
	tmp *= nsave;
	fwrite ( &tmp, sizeof(int), 1, dcdfile );	// write f
	fseek ( dcdfile, curr_offset, SEEK_SET );
	nframes = f;
}



XTCITrajFile::XTCITrajFile()
{
	file_open=0;
	header="";
	curr_frame=0;
	qcrystal = false;
}

XTCITrajFile::XTCITrajFile(const char * file)
{
	file_open = 0;
	qcrystal = false;
	open(file);
}

XTCITrajFile::~XTCITrajFile()
{
}

int XTCITrajFile::open(const char * file)
{
        int magic;
	qcrystal = false;
        XDR *xd;
	if (file_open) {
		print_err( "XTCITrajFile::open", 
				"A trajectory is already open");
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	
	xtcfile = fopen ( file, "r" );
	if ( xtcfile == NULL ) {
		print_err( "XTCITrajFile::open", 
				"could not open trajectory file");
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	
	xdr = new XDR;
	int xx;
	xdrstdio_create(xdr,xtcfile,XDR_DECODE);
	xx = xdr_int(xdr,&magic);
	if ( magic != XTC_MAGIC ) {
		print_err( "XTCITrajFile::open", 
				"wrong magic - is it really an XTC file?");
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}
	fprintf(stdout,"magic = %i\n",magic);
	xx = xdr_int(xdr,&N);
	fprintf(stdout,"N = %i\n",N);
	xx = xdr_int(xdr,&istep);
	fprintf(stdout,"step = %i\n",istep);
	xx = xdr_float(xdr,&itime);
	fprintf(stdout,"itime = %12.6e\n",itime);
	//fprintf(stdout,"xdr_int = %i\n",xx);
	//fprintf(stdout,"magic = %i\n",magic);
	
	XYZ = new float[3*N];
	X =  new float[N];
	Y =  new float[N];
	Z =  new float[N];
	first = true;
	moredata = true;
  
	return 0;
}

void XTCITrajFile::close()
{
	if (file_open) 
		fclose(xtcfile);
	delete xdr;
	delete [] XYZ;
	delete [] X;
	delete [] Y;
	delete [] Z;
}

void XTCITrajFile::rewind()
{
}

void XTCITrajFile::skip(int s)
{
}

int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision, char mode) {
    

    static int *ip = NULL;
    static int oldsize;
    static int *buf;

    int minint[3], maxint[3], mindiff, *lip, diff;
    int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int minidx, maxidx;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int flag, k;
    int smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    float *lfp, lf;
    int tmp, *thiscoord,  prevcoord[3];
    unsigned int tmpcoord[30];

    int bufsize, xdrid, lsize;
    unsigned int bitsize;
    float inv_precision;
    int errval = 1;

    bitsizeint[0] = bitsizeint[1] = bitsizeint[2] = 0;
    prevcoord[0]  = prevcoord[1]  = prevcoord[2]  = 0;
    
    /* find out if xdrs is opened for reading or for writing */
    xdrid = 0;
    /*
    while (xdridptr[xdrid] != xdrs) {
	xdrid++;
	if (xdrid >= MAXID) {
	    fprintf(stderr, "xdr error. no open xdr stream\n");
	    exit (1);
	}
    }
    */
    if ((mode == 'w') || (mode == 'a')) {

	/* xdrs is open for writing */

	if (xdr_int(xdrs, size) == 0)
	    return 0;
	size3 = *size * 3;
	/* when the number of coordinates is small, don't try to compress; just
	 * write them as floats using xdr_vector
	 */
	if (*size <= 9 ) {
            return (xdr_vector(xdrs, (char *) fp, (unsigned int)size3, 
                    (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
	}
	
	xdr_float(xdrs, precision);
	if (ip == NULL) {
	    ip = (int *)malloc((size_t)(size3 * sizeof(*ip)));
	    if (ip == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    bufsize = int(size3 * 1.2);
	    buf = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
	    if (buf == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    oldsize = *size;
	} else if (*size > oldsize) {
	    ip = (int *)realloc(ip, (size_t)(size3 * sizeof(*ip)));
	    if (ip == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    bufsize = int(size3 * 1.2);
	    buf = (int *)realloc(buf, (size_t)(bufsize * sizeof(*buf)));
	    if (buf == NULL) {
		fprintf(stderr,"realloc failed\n");
		exit(1);
	    }
	    oldsize = *size;
	}
	/* buf[0-2] are special and do not contain actual data */
	buf[0] = buf[1] = buf[2] = 0;
	minint[0] = minint[1] = minint[2] = INT_MAX;
	maxint[0] = maxint[1] = maxint[2] = INT_MIN;
	prevrun = -1;
	lfp = fp;
	lip = ip;
	mindiff = INT_MAX;
	oldlint1 = oldlint2 = oldlint3 = 0;
	while(lfp < fp + size3 ) {
	    /* find nearest integer */
	    if (*lfp >= 0.0)
		lf = *lfp * *precision + 0.5;
	    else
		lf = *lfp * *precision - 0.5;
	    if (fabs(lf) > MAXABS) {
		/* scaling would cause overflow */
		errval = 0;
	    }
	    lint1 = int(lf);
	    if (lint1 < minint[0]) minint[0] = lint1;
	    if (lint1 > maxint[0]) maxint[0] = lint1;
	    *lip++ = lint1;
	    lfp++;
	    if (*lfp >= 0.0)
		lf = *lfp * *precision + 0.5;
	    else
		lf = *lfp * *precision - 0.5;
	    if (fabs(lf) > MAXABS) {
		/* scaling would cause overflow */
		errval = 0;
	    }
	    lint2 = int(lf);
	    if (lint2 < minint[1]) minint[1] = lint2;
	    if (lint2 > maxint[1]) maxint[1] = lint2;
	    *lip++ = lint2;
	    lfp++;
	    if (*lfp >= 0.0)
		lf = *lfp * *precision + 0.5;
	    else
		lf = *lfp * *precision - 0.5;
	    if (fabs(lf) > MAXABS) {
		/* scaling would cause overflow */
		errval = 0;
	    }
	    lint3 = int(lf);
	    if (lint3 < minint[2]) minint[2] = lint3;
	    if (lint3 > maxint[2]) maxint[2] = lint3;
	    *lip++ = lint3;
	    lfp++;
	    diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
	    if (diff < mindiff && lfp > fp + 3)
		mindiff = diff;
	    oldlint1 = lint1;
	    oldlint2 = lint2;
	    oldlint3 = lint3;
	}
	xdr_int(xdrs, &(minint[0]));
	xdr_int(xdrs, &(minint[1]));
	xdr_int(xdrs, &(minint[2]));
	
	xdr_int(xdrs, &(maxint[0]));
	xdr_int(xdrs, &(maxint[1]));
	xdr_int(xdrs, &(maxint[2]));
	
	if ((float)maxint[0] - (float)minint[0] >= MAXABS ||
		(float)maxint[1] - (float)minint[1] >= MAXABS ||
		(float)maxint[2] - (float)minint[2] >= MAXABS) {
	    /* turning value in unsigned by subtracting minint
	     * would cause overflow
	     */
	    errval = 0;
	}
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
	    bitsizeint[0] = sizeofint(sizeint[0]);
	    bitsizeint[1] = sizeofint(sizeint[1]);
	    bitsizeint[2] = sizeofint(sizeint[2]);
	    bitsize = 0; /* flag the use of large sizes */
	} else {
	    bitsize = sizeofints(3, sizeint);
	}
	lip = ip;
	luip = (unsigned int *) ip;
	smallidx = FIRSTIDX;
	while (smallidx < LASTIDX && magicints[smallidx] < mindiff) {
	    smallidx++;
	}
	xdr_int(xdrs, &smallidx);
	maxidx = MIN(LASTIDX, smallidx + 8) ;
	minidx = maxidx - 8; /* often this equal smallidx */
	smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
	larger = magicints[maxidx] / 2;
	i = 0;
	while (i < *size) {
	    is_small = 0;
	    thiscoord = (int *)(luip) + i * 3;
	    if (smallidx < maxidx && i >= 1 &&
		    abs(thiscoord[0] - prevcoord[0]) < larger &&
		    abs(thiscoord[1] - prevcoord[1]) < larger &&
		    abs(thiscoord[2] - prevcoord[2]) < larger) {
		is_smaller = 1;
	    } else if (smallidx > minidx) {
		is_smaller = -1;
	    } else {
		is_smaller = 0;
	    }
	    if (i + 1 < *size) {
		if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
			abs(thiscoord[1] - thiscoord[4]) < smallnum &&
			abs(thiscoord[2] - thiscoord[5]) < smallnum) {
		    /* interchange first with second atom for better
		     * compression of water molecules
		     */
		    tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
			thiscoord[3] = tmp;
		    tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
			thiscoord[4] = tmp;
		    tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
			thiscoord[5] = tmp;
		    is_small = 1;
		}
    
	    }
	    tmpcoord[0] = thiscoord[0] - minint[0];
	    tmpcoord[1] = thiscoord[1] - minint[1];
	    tmpcoord[2] = thiscoord[2] - minint[2];
	    if (bitsize == 0) {
		sendbits(buf, bitsizeint[0], tmpcoord[0]);
		sendbits(buf, bitsizeint[1], tmpcoord[1]);
		sendbits(buf, bitsizeint[2], tmpcoord[2]);
	    } else {
		sendints(buf, 3, bitsize, sizeint, tmpcoord);
	    }
	    prevcoord[0] = thiscoord[0];
	    prevcoord[1] = thiscoord[1];
	    prevcoord[2] = thiscoord[2];
	    thiscoord = thiscoord + 3;
	    i++;
	    
	    run = 0;
	    if (is_small == 0 && is_smaller == -1)
		is_smaller = 0;
	    while (is_small && run < 8*3) {
		if (is_smaller == -1 && (
			SQR(thiscoord[0] - prevcoord[0]) +
			SQR(thiscoord[1] - prevcoord[1]) +
			SQR(thiscoord[2] - prevcoord[2]) >= smaller * smaller)) {
		    is_smaller = 0;
		}

		tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
		tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
		tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;
		
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];

		i++;
		thiscoord = thiscoord + 3;
		is_small = 0;
		if (i < *size &&
			abs(thiscoord[0] - prevcoord[0]) < smallnum &&
			abs(thiscoord[1] - prevcoord[1]) < smallnum &&
			abs(thiscoord[2] - prevcoord[2]) < smallnum) {
		    is_small = 1;
		}
	    }
	    if (run != prevrun || is_smaller != 0) {
		prevrun = run;
		sendbits(buf, 1, 1); /* flag the change in run-length */
		sendbits(buf, 5, run+is_smaller+1);
	    } else {
		sendbits(buf, 1, 0); /* flag the fact that runlength did not change */
	    }
	    for (k=0; k < run; k+=3) {
		sendints(buf, 3, smallidx, sizesmall, &tmpcoord[k]);	
	    }
	    if (is_smaller != 0) {
		smallidx += is_smaller;
		if (is_smaller < 0) {
		    smallnum = smaller;
		    smaller = magicints[smallidx-1] / 2;
		} else {
		    smaller = smallnum;
		    smallnum = magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
	    }
	}
	if (buf[1] != 0) buf[0]++;
	xdr_int(xdrs, &(buf[0])); /* buf[0] holds the length in bytes */
        return errval * (xdr_opaque(xdrs, (char *)&(buf[3]), (unsigned int)buf[0]));
    } else {
	
	/* xdrs is open for reading */
	
	if (xdr_int(xdrs, &lsize) == 0) 
	    return 0;
	if (*size != 0 && lsize != *size) {
	    fprintf(stderr, "wrong number of coordinates in xdr3dfcoord; "
		    "%d arg vs %d in file", *size, lsize);
	}
	*size = lsize;
	size3 = *size * 3;
	if (*size <= 9) {
            return (xdr_vector(xdrs, (char *) fp, (unsigned int)size3, 
                    (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
	}
	xdr_float(xdrs, precision);
	if (ip == NULL) {
	    ip = (int *)malloc((size_t)(size3 * sizeof(*ip)));
	    if (ip == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    bufsize = int(size3 * 1.2);
	    buf = (int *)malloc((size_t)(bufsize * sizeof(*buf)));
	    if (buf == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    oldsize = *size;
	} else if (*size > oldsize) {
	    ip = (int *)realloc(ip, (size_t)(size3 * sizeof(*ip)));
	    if (ip == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    bufsize = int(size3 * 1.2);
	    buf = (int *)realloc(buf, (size_t)(bufsize * sizeof(*buf)));
	    if (buf == NULL) {
		fprintf(stderr,"malloc failed\n");
		exit(1);
	    }
	    oldsize = *size;
	}
	buf[0] = buf[1] = buf[2] = 0;
	
	xdr_int(xdrs, &(minint[0]));
	xdr_int(xdrs, &(minint[1]));
	xdr_int(xdrs, &(minint[2]));

	xdr_int(xdrs, &(maxint[0]));
	xdr_int(xdrs, &(maxint[1]));
	xdr_int(xdrs, &(maxint[2]));
		
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
	    bitsizeint[0] = sizeofint(sizeint[0]);
	    bitsizeint[1] = sizeofint(sizeint[1]);
	    bitsizeint[2] = sizeofint(sizeint[2]);
	    bitsize = 0; /* flag the use of large sizes */
	} else {
	    bitsize = sizeofints(3, sizeint);
	}
	
	if (xdr_int(xdrs, &smallidx) == 0)	
	    return 0;
	maxidx = MIN(LASTIDX, smallidx + 8) ;
	minidx = maxidx - 8; /* often this equal smallidx */
	smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	larger = magicints[maxidx];

    	/* buf[0] holds the length in bytes */

	if (xdr_int(xdrs, &(buf[0])) == 0)
	    return 0;
        if (xdr_opaque(xdrs, (char *)&(buf[3]), (unsigned int)buf[0]) == 0)
	    return 0;
	buf[0] = buf[1] = buf[2] = 0;
	
	lfp = fp;
	inv_precision = 1.0 / * precision;
	run = 0;
	i = 0;
	lip = ip;
	while ( i < lsize ) {
	    thiscoord = (int *)(lip) + i * 3;

	    if (bitsize == 0) {
		thiscoord[0] = receivebits(buf, bitsizeint[0]);
		thiscoord[1] = receivebits(buf, bitsizeint[1]);
		thiscoord[2] = receivebits(buf, bitsizeint[2]);
	    } else {
		receiveints(buf, 3, bitsize, sizeint, thiscoord);
	    }
	    
	    i++;
	    thiscoord[0] += minint[0];
	    thiscoord[1] += minint[1];
	    thiscoord[2] += minint[2];
	    
	    prevcoord[0] = thiscoord[0];
	    prevcoord[1] = thiscoord[1];
	    prevcoord[2] = thiscoord[2];
	    
	   
	    flag = receivebits(buf, 1);
	    is_smaller = 0;
	    if (flag == 1) {
		run = receivebits(buf, 5);
		is_smaller = run % 3;
		run -= is_smaller;
		is_smaller--;
	    }
	    if (run > 0) {
		thiscoord += 3;
		for (k = 0; k < run; k+=3) {
		    receiveints(buf, 3, smallidx, sizesmall, thiscoord);
		    i++;
		    thiscoord[0] += prevcoord[0] - smallnum;
		    thiscoord[1] += prevcoord[1] - smallnum;
		    thiscoord[2] += prevcoord[2] - smallnum;
		    if (k == 0) {
			/* interchange first with second atom for better
			 * compression of water molecules
			 */
			tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
				prevcoord[0] = tmp;
			tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
				prevcoord[1] = tmp;
			tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
				prevcoord[2] = tmp;
			*lfp++ = prevcoord[0] * inv_precision;
			*lfp++ = prevcoord[1] * inv_precision;
			*lfp++ = prevcoord[2] * inv_precision;
		    } else {
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
		    }
		    *lfp++ = thiscoord[0] * inv_precision;
		    *lfp++ = thiscoord[1] * inv_precision;
		    *lfp++ = thiscoord[2] * inv_precision;
		}
	    } else {
		*lfp++ = thiscoord[0] * inv_precision;
		*lfp++ = thiscoord[1] * inv_precision;
		*lfp++ = thiscoord[2] * inv_precision;		
	    }
	    smallidx += is_smaller;
	    if (is_smaller < 0) {
		smallnum = smaller;
		if (smallidx > FIRSTIDX) {
		    smaller = magicints[smallidx - 1] /2;
		} else {
		    smaller = 0;
		}
	    } else if (is_smaller > 0) {
		smaller = smallnum;
		smallnum = magicints[smallidx] / 2;
	    }
	    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	}
    }
    return 1;
}

void XTCITrajFile::read_frame(float *XX, float *YY, float *ZZ, int n)
{
	float a[3], b[3], c[3];
	float A,B,C,alpha,beta,gamma,precision,tmpf;
	int tmpi,xx,m;
	//fprintf(stderr," in read_frame now...1\n"); fflush(stderr);
	for (int j=0; j<3; j++) {
		xdr_float(xdr,&a[j]);
		//fprintf(stdout,"a[%i] = %12.6e\n",j,a[j]);
	}
	for (int j=0; j<3; j++) {
		xdr_float(xdr,&b[j]);
	//	fprintf(stdout,"b[%i] = %12.6e\n",j,b[j]);
	}
	for (int j=0; j<3; j++) {
		xdr_float(xdr,&c[j]);
	//	fprintf(stdout,"c[%i] = %12.6e\n",j,c[j]);
	}
	//fprintf(stderr," in read_frame now...1.5\n"); fflush(stderr);
	A = sqrt(dotpn(a,a,3));
	B = sqrt(dotpn(a,a,3));
	C = sqrt(dotpn(a,a,3));
	alpha = acos(dotpn(b,c,3)/(B*C));
	beta  = acos(dotpn(a,c,3)/(A*C));
	gamma = acos(dotpn(a,b,3)/(A*B));
	//fprintf(stderr," in read_frame now...2\n"); fflush(stderr);
	// A,B,C in nm; alpha,beta,gamma in rad;
	// convert to Angstroms and degrees
	crystal_data[0] = A*10.0;
	crystal_data[1] = B*10.0;
	crystal_data[2] = C*10.0;
	crystal_data[3] = alpha*180./M_PI;
	crystal_data[4] = beta*180./M_PI;
	crystal_data[5] = gamma*180./M_PI;
	//fprintf(stderr," in read_frame now...3\n"); fflush(stderr);
	int ii = xdr3dfcoord(xdr,XYZ,&N,&precision,'r');
	if (!ii) {
		moredata = false;
	}
	//fprintf(stderr," in read_frame now...4\n"); fflush(stderr);
	for (int i=0; i<N; i++) {
		XX[i] = XYZ[3*i]*10.;
		YY[i] = XYZ[3*i+1]*10.;
		ZZ[i] = XYZ[3*i+2]*10.;
	}
		xx = xdr_int(xdr,&m);
		xx = xdr_int(xdr,&m);
		xx = xdr_int(xdr,&m);
		xx = xdr_float(xdr,&tmpf);
	if (!xx) {
		moredata = false;
	}
	//fprintf(stderr," in read_frame now...5\n"); fflush(stderr);
	return;
}



void XTCITrajFile::read_frame(double *XX, double *YY, double *ZZ, int n)
{
	read_frame(X,Y,Z,n);
	for (int i=0;i<n;i++) {
		XX[i] = X[i];
		YY[i] = Y[i];
		ZZ[i] = Z[i];
	}
}

/*
bool XTCITrajFile::frames_left() const
{
	return false;
}
*/

void XTCITrajFile::get_crystal_data(double *rdata) const
{ 
	for (int i=0; i<6; i++)
		rdata[i] = crystal_data[i];
}


