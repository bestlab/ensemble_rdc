/*
	Trajectory file handles, with functions for reading
	& writing frames, skipping frames, etc.
*/

#ifndef _TRAJFILE_H
#define _TRAJFILE_H

#include <cstdio>
#include <cstdlib>
//#include <unistd.h>
#include <string>
#include <cstring>
#include <vector>
//#include "genfunlib.h"
#include "swapbytes.h"
#include <rpc/rpc.h>
#include <rpc/xdr.h>
//#include "config.h"

using namespace std;

/* 
Notes:
   In the classes below, there is an important distinction between dynamics
   steps, which are the timesteps at which the equations of motion are
   integrated, and trajectory frames, which are the coordinate sets written
   out at a frequency >= dynamic integration frequency.
 */ 

/*
   Change Log:
   !! put your changes here !!
   28 nov 1999: Robert Best: Added output facilities
   
   Thu Feb 24 23:16:05 SAST 2000
   RB - added on-the-fly byte-swapping 
   */

/***************************************************************
 *         BaseFile:                                           *
 *         Basic functions such as open/close, which should    *
 *         be provided by all trajectory classes               *
 *         - gives consistent interface                        *
 ***************************************************************/
class BaseFile
{
	public:
		/* virtual destructor */
		//virtual ~BaseFile() {}
		/* standard file operations */
		virtual int open( const char * ) = 0;
		virtual void close() = 0;
		virtual int flopen() const = 0;
		/* Classes to be thrown by exceptions */
#ifndef SWIG
		class FileErr{};
		class OpenErr : public FileErr{};
		class ReadErr : public FileErr{};
		class WriteErr : public FileErr{};
		class EofErr : public ReadErr{};
#endif
};

/***************************************************************
 *         BaseITrajFile:                                      *
 *         Base class for input trajectory files; all          *
 *         functions should use this class rather than         *
 *         the derivate DCDITrajFile, etc.                     *
 ***************************************************************/
class BaseITrajFile : public BaseFile
{
	public:
		virtual void rewind() =0;
		virtual void skip(int) = 0;		// skip n frames
		virtual int initial_step() const = 0; 	// initial dynamics step
		virtual int current_frame() const= 0;	// current dynamics step
		virtual int total_frames() const = 0;	// total # of frames
		virtual int frame_freq() const = 0;   	// frequency in dynsteps
							// for writing frames
		virtual double step_size() const = 0;	// step sz in AKMA units
		virtual int num_atoms() const = 0; 	// total number of atoms.
		virtual int fixed_atoms() const = 0;	// number of fixed atoms
		virtual void read_frame(double *, double *, double *, int) = 0; 
		virtual void read_frame(float *, float *, float *, int) = 0; 
		virtual bool crystal() const = 0;
		virtual void get_crystal_data( double[6] ) const =0;
		virtual bool frames_left() const = 0;
#ifndef SWIG
		class InpTrajErr{};
		class WrongNumAtomsErr : public InpTrajErr{};
		class BinFormatErr : public InpTrajErr{};
		class FileFormatErr : public InpTrajErr{};
		class IllegalRequestErr : public InpTrajErr{};
#endif

};
	

/*
   ITrajSet
   Class representing a set of input trajectory        
   files. Acts like any other Trajectory (hopefully).
   Need to assume (at present) that timestep between 
   frames is identical in all trajectories. Checks are
   performed to ensure that the files "fit together"
 */
/*
class ITrajSet : public BaseITrajFile
{
	private:
		vector< BaseITrajFile > trajfiles;
		int ntrajfile;
		vector< int > lbounds;	// position of first frame in file
					// relative to first in set.
		vector< int > flengths;	// number of frames in each file.
		int currentfile;
		int currentframe;	// relative to first
		int nframes;		// number of frames in file set 
		int itime;		// cumulative dynamics step number 
					//   of first frame
		int nsave;		// frequency (in dyn. steps ) 
					// for saving coordinates or velocities
		int nsteps;		// total dynamics steps in file
		int dof;		// number of degrees of freedom
		int nfixed;		// number of fixed atoms
		double delta;		// timestep size
		bool qcrystal; 		// crystal ?
		int N;			// total number of atoms 
		// file objects, status
		int file_open;		// checks whether ALL files are open
		string header;		// header of first file
		int curr_frame;

	public:
		// constructors and destructors
		ITrajSet( int n = 60 );
		ITrajSet( const vector< BaseITrajFile > & tfiles,int ntrajfile);
		~ITrajSet();
		// stuff to satisfy BaseFile
		int open( const char * ) {;}	// does nothing
		void setup( const vector< BaseITrajFile > & tfiles, 
				int ntrajfile);
		void append( const vector< BaseITrajFile > & tfiles );
		void close();
		int flopen() const;
		// stuff to satisfy BaseITrajFile
		void rewind();
		void skip(int);			// skip n frames
		int initial_step() const; 	// initial dynamics step
		int current_frame() const;	// current dynamics step
		int total_frames() const;	// total number of frames
		int frame_freq() const;		// frequency 
						// (in dynamics steps)
						// for writing frames
		double step_size() const;	// step size in AKMA 
						// time units
		int num_atoms() const;	 	// total number of atoms.
		int fixed_atoms() const;	// number of fixed atoms
		void read_frame( double *, double *, double *, int ); 
		bool crystal() const;
		void get_crystal_data( double[6] ) const; 
		bool frames_left();
};
*/

/***************************************************************
 *         BaseOTrajFile:                                      *
 *         Base class for output trajectory files; all         *
 *         functions should use this class rather than         *
 *         the derivate DCDOTrajFile, etc.                     *
 ***************************************************************/
class BaseOTrajFile : public BaseFile
{
	public:
		virtual void initial_step(int) =0;
		virtual void num_atoms(int) = 0;
		virtual void write_header() =0;
		virtual void write_frame( double *, double *, double *, int ) =0;
#ifndef SWIG
		class OutTrajErr{};
		class WrongNumAtomsErr : public OutTrajErr{};
		class InternalDataErr : public OutTrajErr{};
		class IllegalRequestErr : public OutTrajErr{};
		class NoHeaderWrittenErr : public OutTrajErr{};
#endif
};


/***************************************************************
 *         DCDITrajFile:                                       *
 *         class representing traditional charmm binary        *
 *         dcd trajectory format (input only)                  *
 ***************************************************************/
class DCDITrajFile : public BaseITrajFile
{
	private:
		// dcd file stats
		char type;		// "CORD" -> c, "VELD" -> v
		bool swap_bytes;	// whether to automagically swap bytes
		int nframes;		// number of frames in file -- according to header
		int actual_nframes;	// actual number of frames written
		int itime;		// cumulative dynamics step number of first frame
		int nsave;		/* frequency (in dyn. steps ) 
					   for saving coordinates or velocities */
		int nsteps;		// total dynamics step in file
		int dof;		// number of degrees of freedom
		int nfixed;		// number of fixed atoms
		double delta;		// timestep size = static_cast<float>(icntrl[9])
		bool qcrystal; 		// crystal ?
		int version;		// CHARMM version
		int N;			// total number of atoms 
		double crystal_data[6];	// what's this?
		// file objects, status
		FILE * dcdfile;
		long int start_record;
		int file_open;
		string header;
		int curr_frame;
		float * X;
		float * Y;
		float * Z;
		float * freex;
		float * freey;
		float * freez;
		int * freeatoms; 	// free atom list
		bool sfbit;		// 64 bit integers?
		bool sfbita;		// 64 bit aligned

	public:
		// constructors, destructors, etc.
		DCDITrajFile();
		DCDITrajFile( const char * );
		~DCDITrajFile();
		// stuff to satisfy BaseFile
		int open( const char * );
		void close();
		int flopen() const { return file_open; }
		// stuff to satisfy BaseITrajFile
		void rewind();
		void skip(int);			// skip n frames
		int initial_step() const { return itime; }
		int current_frame() const { return curr_frame; }
		int total_frames() const { return nframes; }
		int actual_frames() const { return actual_nframes; }
		int frame_freq() const { return nsave; }
		double step_size() const { return delta; }
		int num_atoms() const { return N; }
		int fixed_atoms() const { return nfixed; }
		void read_frame( double *, double *, double *, int ); 
		void read_frame( float *, float *, float *, int ); 
		bool frames_left() const;
		// other stuff
		bool swapping_bytes() const {return swap_bytes; }
		void get_crystal_data( double [6] ) const;
		int charmm_version() const { return version; }
		void goto_frame( int );
		char traj_type() const { return type; }
		int total_steps() const { return nsteps; }
		int deg_free() const { return dof; }
		bool crystal() const { return qcrystal; }
		void read_title( string & returnstring ) const
		{
			returnstring = header;
		}
		bool sixtyfourbit() const { return sfbit; } 
};


/***************************************************************
 *         DCDOTrajFile:                                       *
 *         class representing traditional charmm binary        *
 *         dcd trajectory format (output only)                 *
 ***************************************************************/
class DCDOTrajFile : public BaseOTrajFile
{
	private:
		// dcd file stats
		char type;		// "CORD" -> c, "VELD" -> v
		int nframes;		// number of frames in file
		int itime;		// cumulative dynamics step number of first frame
		int nsave;		/* frequency (in dyn. steps ) 
					   for saving coordinates or velocities */
		int nsteps;		// total dynamics step in file
		int dof;		// number of degrees of freedom
		int nfixed;		// number of fixed atoms
		double delta;		// timestep size = static_cast<float>(icntrl[9])
		bool qcrystal; 		// crystal ?
		int version;		// CHARMM version
		int N;			// total number of atoms 
		double crystal_data[6];	// what's this?
		// file objects, status
		FILE * dcdfile;
		long int start_record;
		int file_open;
		int header_written;
		string header;
		int curr_frame;
		float * X;
		float * Y;
		float * Z;
		float * freex;
		float * freey;
		float * freez;
		int * freeatoms; 	// free atom list

	public:
		// constructors, destructors, etc.
		DCDOTrajFile();
		DCDOTrajFile( const char * );
		DCDOTrajFile( const DCDITrajFile & );
		~DCDOTrajFile();
		// stuff to satisfy BaseFile
		int open( const char * );
		void close();
		int flopen() const { return file_open; }
		// stuff to satisfy BaseOTrajFile
		void initial_step(int s) { itime = s; } 
		void num_atoms(int num) { N = num; }
		void write_header();
		void write_frame( double *, double *, double *, int );
		void write_frame( float *, float *, float *, int );
		// other stuff
		void set_crystal( bool iscrystal ) { qcrystal=iscrystal; }
		bool get_crystal() const { return qcrystal; }
		/*
		void set_icntrl( int, int, int, int, int, int, int, int,
				double );
		void set_freeat ( int, double *, double *, double * );
		void set_type ( char );
		*/
		void setup ( const DCDITrajFile & );
		//void set_dof(int d) { dof = d; }
		void set_frames( int f ) { nframes = f; }
		void set_steps( int s ) { nsteps = s; }
		void set_crystal_data( double * xtaldata );
		void rewrite_frames_and_steps( int f );
};

class XTCITrajFile : public BaseITrajFile
{
	private:
		// xtc file stats
		int nframes;		// number of frames in file -- according to header
		//int actual_nframes;	// actual number of frames written
		int dof;		// number of degrees of freedom
		int nfixed;		// number of fixed atoms
		int N;			// total number of atoms 
		int nsave;
		double delta;
		int istep;
		float itime;
		bool qcrystal;
		double crystal_data[6];	// what's this?
		// file objects, status
		FILE * xtcfile;
		XDR * xdr;
		long int start_record;
		int file_open;
		string header;
		int curr_frame;
		bool first;
		float * X;
		float * Y;
		float * Z;
		float *XYZ;
		float box[3][3];
		bool moredata;

	public:
		// constructors, destructors, etc.
		XTCITrajFile();
		XTCITrajFile( const char * );
		~XTCITrajFile();
		// stuff to satisfy BaseFile
		int open( const char * );
		void close();
		int flopen() const { return file_open; }
		// stuff to satisfy BaseITrajFile
		void rewind();
		void skip(int);			// skip n frames
		int initial_step() const { return istep; }
		int current_frame() const { return curr_frame; }
		int total_frames() const { return nframes; }
		int frame_freq() const { return 0; }
		double step_size() const { return 0.0; }
		int num_atoms() const { return N; }
		int fixed_atoms() const { return 0; }
		//void read_frame( double *, double *, double *, int ); 
		void read_frame( float *, float *, float *, int ); 
		void read_frame( double *, double *, double *, int ); 
		bool frames_left() const { return moredata; };
		void get_crystal_data( double [6] ) const;
		bool crystal() const { return qcrystal; }
		// other stuff
		/*
		int charmm_version() const { return version; }
		void goto_frame( int );
		char traj_type() const { return type; }
		int total_steps() const { return nsteps; }
		int deg_free() const { return dof; }
		void read_title( string & returnstring ) const
		{
			returnstring = header;
		}
		bool sixtyfourbit() const { return sfbit; } 
		*/
};

#endif 		// _TRAJFILE_H

