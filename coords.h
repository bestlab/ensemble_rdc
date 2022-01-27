/*
   File - coords.h
   Author - Robert Best
   Date - April 1999
   Description - Declarations and Class Definitions for coords.cc which 
   implements the concept of a Charmm coordinate set

 */

/* CHANGE LOG:
   21 nov 1999: Included boxlength, Bound and Xtal_data in coords;
   changed functions which previously took these as arguments to use the
   local data.
   Input from trajectories revamped to used new hierarchy of traj file
   classes in TrajFile.h, which now allow input from "crystal"-type trajectories
   and allow much more info to be gleaned from the dcd file. dcd >> coords
   changed to automagically change box size at each step based on crystal info
   in dcd file. NB crystal data contains much more info, which could be used for
   arbitrary parallelopiped boxes!
 */

#ifndef _COORDS_H
#define _COORDS_H 

#include <math.h>
#include <string>
#include <cstring>
#include "traj.h"

#include <iostream>
#include <iomanip>

#define DEBUG


typedef ios_base::fmtflags FMTFLAGS;
using namespace std;

//******************** Globals ****************************//

const long MAX_ATOMS = 60120;
const int CHARMM_STR = 5;
enum BoundCond { CUBIC, OCTAHEDRAL, NONE };


// TODO: replace AtomType by simple string from libstdc++
class AtomType
{
	public:
		AtomType() { strcpy(Type,""); }
		AtomType(char * at_str);
		void Set(const char * at_str);
		void Set(AtomType at_typ);
		void Get(char * at_str) const;
		char* Get() const { char *tmp=new char[5]; strcpy(tmp,Type); return tmp; }

		//-------------Operators--------------------
		const AtomType & operator=(const AtomType & Rhs);
		int operator==(const AtomType & Rhs) const;
		friend ostream& operator<< (ostream&, const AtomType &);

	private:
		static const int MAX=5; // 4 characters + null terminator
		char Type[MAX];
};

//end added
class coords
{

	public:

		//--------Constructors and Destructors------------//
		coords();        // allocates only pointer space
		coords(int N);  // allocates correct free store 
		coords(const char*, const char * type= "crd");    
		// from CRD file
		coords(const coords&);  // Copy Constructor
		coords(const coords&, int * selection, int nselect,
				bool resrenum=1);  
		// Selection/Copy Constructor
		~coords();

		//************* General Functions *****************//

		int Coorwrite(const char *outfile = "output.crd", 
				const char *comment = "DEFAULT COMMENT");
		int Coorread(const char *infile, const char * type = "crd");    
		int Get_N() const { return N; }
		int Get_Atomno(int i) const { return Atomno[i]; }
		int Get_Resno(int i) const { return Resno[i]; }
		void Get_Resname(int i, char* res) const { strcpy(res, 
				Resname[i]); }
		void Get_Atomtype(int i, char* typ) const { Type[i].Get(typ); }
		void Get_Atomtype(int i, AtomType & typ) const {typ = Type[i]; }
		const string & Get_Atomtype(int i) const { return type[i]; }
		double Get_X(int i) const { return X[i]; }
		double Get_Y(int i) const { return Y[i]; }
		double Get_Z(int i) const { return Z[i]; }
		void Get_XYZ(int i, double * crd) const { 
			crd[0] = X[i]; crd[1] = Y[i]; crd[2] = Z[i];
		}
		void Get_Segid(int i, char* seg) const {strcpy(seg, Segid[i]); }
		void Get_Resid(int i, char* res) const {strcpy(res, Resid[i]); }
		double Get_Weight (int i) const { return Weight[i]; }
		int Get_Bound () const { return Bound; }
		bool Get_octcharmm() const {return octcharmm; }
		double Get_boxlength () const { return boxlength; }
		double Get_Volume() const;
		void copy( coords & rhs );
//		void add( const coords & rhs );
//		void add( double rhs );
//		void multiply( double cnst );
		void Get_Selection( int * selection, int nselect, 
				float * X, float * Y, float * Z) const;
		void Get_Selection( int * selection, int nselect, double * X, 
				double * Y, double * Z) const;
		void Get_Selection( int * selection, int nselect, 
				vector<double> & XYZ) const;
		void Set_Selection( int * selection, int nselect, 
				vector<double> & XYZ);
		void Set_XYZ( double *IX, double *IY, double *IZ, int n ) {
			for (int i=0; i<n; i++) {
				X[i] = IX[i];
				Y[i] = IY[i];
				Z[i] = IZ[i];
			}
		}

		// TODO: This is a hack: remove ASAP!
		void Get_Subset( const coords & source );
		void copy_selection(int * selection, int nselect, 
				coords & B);
		void Set_N(int n)  { N=n; }
		void Set_Atomno(int i, int n) { Atomno[i] = n; }
		void Set_Resno(int i, int  n) { Resno[i]=n; }
		void Set_Resname(int i, const char* res) 
		{ strcpy(Resname[i],res); }
		void Set_Atomtype(int i, char* typ) 
		{ Type[i].Set(typ); type[i] = string(typ); }
		void Set_Atomtype(int i, AtomType typ) { Type[i].Set(typ); }
		int Set_Atomtype(int i, string s); 
		void Set_X(int i, double d) { X[i]=d; }
		void Set_Y(int i, double d) { Y[i]=d; }
		void Set_Z(int i, double d) { Z[i]=d; }
		void Set_Segid(int i, const char* seg) {strcpy(Segid[i], seg);}
		void Set_Resid(int i, const char* res) {strcpy(Resid[i], res);}
		void Set_Weight (int i, double w) { Weight[i] = w; }
		void Set_Bound (BoundCond bnd) { Bound = bnd; }
		void Set_octcharmm(bool val) { octcharmm = val; 
						if (val) orot_init(); }
		void Set_boxlength (double l) { boxlength = l; }
		void add( int i, double x, double y, double z );

		// General utilities

		// selection utilities; n = size of array
		//			  sel = num of atoms selected
		/*
		   void select_atoms( const int * resnum,
		   const string * atomtypes, int * select_array, int n, 
		   int & sel ); 
		   void select_segid( const char * segid, int * select_array, 
		   int n, int & sel ); 
		   void select_resnum( int resnum, int * select_array, int n, 
		   int & sel ); 
		 */
		// periodic boundary conditions - still need to do diff for 
		// octahedral
		double xdiff(int i, int j) ;
		double ydiff(int i, int j); 
		double zdiff(int i, int j); 
		// the next is more effecient if all three diffs need to be 
		// calc'd
		void diffvec(int i, int j, double & xdiff, 
				double & ydiff,	double & zdiff );
		//m.kuttel changes robert's distance command to take inot 
		// account
		// periodic boundary conditions
		double dist(int i, int j);
		double dist_from_point( double cm[3],int i);
		void vector_from_point( double cm[3],int i, double & dx, 
				double & dy, double & dz ) const;
		double dihedral_angle(int at1, int at2, int a3, int a4);
		double angle(int at1, int at2, int a3);
		//later fix so atom mass is assumed
		void center_of_mass(int atom_list[], int no_atoms,
				double mass[], double CM[3]);
		void center_of_mass(double CM[3]);
		void center_of_mass(double CM[3], int * selection, int n);
		void rotate( double rotmat[3][3], double &x, double &y, 
				double &z );
		void rotate( double rotmat[3][3] );
		void rotate(double rotmat[3][3], int * selection, int nselect );
		double rgyr();
		double rgyr( int * selection, int n );
		void translate( double x, double y, double z );
		void translate( double x, double y, double z, 
				int * selection, int nselect );
		//void dipole( int * selection, int n, double dip[3] );
		//void quadrupole ( int * selection, int n, double quad[3][3]);
		void assign_masses ();
		void assign_masses ( int * massarray, int n );
		void assign_masses ( int * massarray, int * selection,
				int nselect );
		void set_uniform_weights();
		int count_residues( const string & segid ) const;
		// fit coordinates using rotation by kabsch method
		double find_kabsch_fit( coords & comp, int * selection, 
				int n, double rot[3][3] );
		double find_kabsch_fit( coords & comp, int * selection, 
				int * selection_comp,
				int n, double rot[3][3] );
		void cat ( const coords & set2);
		static void set_defaults( const coords & new_def );
		friend double rms ( const coords & set1, 
				const coords & set2 );
		friend double rms ( const coords & set1, 
				const coords & set2,
				int * selection, int selected );
		friend double sel_rms2 ( const coords & set1, 
				const coords & set2,
				int * selection, int selected, 
				double * rmssel );

		friend class PDBOutput; 
		//**************** Operators **********************//
		friend ostream& operator<< (ostream&, const coords&);
		friend istream& operator>> (istream&, coords&);
		friend void operator>> (BaseITrajFile&, coords&);
		friend void operator<< (coords&, BaseOTrajFile&);
//		friend ostream& operator+ (ostream&, const coords&);
		coords &  operator= (const coords & rhs);
		coords &  operator+= (const coords & rhs);
		coords &  operator*= (const double & rhs);

		//****************** Friends **********************//

		//*************** Exceptions **********************//
#ifndef SWIG
		class FileErr{};
		class OpenErr : public FileErr{};
		class ReadErr : public FileErr{};
		class WriteErr : public FileErr{};
		class MathErr{};
		class DivZeroErr : public MathErr{};
		class coordsErr{};
		class WrongNumAtomsErr : public coordsErr{};
		class BoundaryErr : public coordsErr{};
		class StringLenErr : public coordsErr {};
#endif

	private:

		static coords default_crd;
		BoundCond Bound;
		bool octcharmm;
		int N;                            // number of atoms in set
		int *Atomno;
		int *Resno;
		char (*Resname)[CHARMM_STR];
		AtomType * Type;
		string * type;		// will replace AtomType in future
		double *X;                    // xcoor array
		double *Y;
		double *Z;
		char (*Segid)[CHARMM_STR];                
		char (*Resid)[CHARMM_STR];
		double *Weight;
		double Xtal_data[6];	// holds upper triangle of charmm 
					// "shape matrix" from dcd file
		double boxlength;
		double orot[3][3];
		void orot_init();
};

/*
 * for writing coordinate data in pdb format; allows multiple chains
 * and multiple models
 */
class PDBOutput {
	public:
		PDBOutput(const string & pdbname, 
				char * title = "PDB Output");
		~PDBOutput();
		void append_chain(const coords & coords, char chain=' ');
		void append_model(const coords & coords);
		void close();
#ifndef SWIG
		class FileErr{};
#endif
	private:
		string outfilename;
		//ofstream outfile;
		FILE * outfile;
		int natom, model, mdl_natom;
};

#endif     // _COORDS_H

