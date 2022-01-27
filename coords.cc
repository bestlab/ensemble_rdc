/*
File: coords.cc
Author: Robert Best
Date: April 1999
Description: Declarations and Class Definitions for coords.C which 
implements the concept of a Charmm coordinate set
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <unistd.h>              // provides exit() system call  
#include <iomanip>             // provides setw()
#include <ctype.h>
#include <cstdio>
#include "coords.h"
#include "common.h"

#undef DEBUG		// we hope !

/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------
  									
 			ATOMTYPE CLASS DEFINITIONS
  								
  -----------------------------------------------------------------------
  ------------------------------------------------------------------------*/

void coords::orot_init() {
	orot[0][0]=-1.0/3.0; orot[0][1]=2.0/3.0; orot[0][2]=2.0/3.0; 
	orot[1][0]=2.0/3.0; orot[1][1]=-1.0/3.0; orot[1][2]=2.0/3.0; 
	orot[2][0]=2.0/3.0; orot[2][1]=2.0/3.0; orot[2][2]=-1.0/3.0; 
}

//functions for atomtype class
AtomType::AtomType(char * at_str)
{
	strcpy(Type, at_str);
}

void AtomType::Set(const char * at_str)
{
	strcpy(Type, at_str);
}

void AtomType::Set(AtomType at_typ)
{
  for(int i=0; i<MAX; i++)
    Type[i]=at_typ.Type[i];
}

void AtomType::Get(char * at_str) const
{
	strcpy(at_str,Type);
}

const AtomType & AtomType::operator=(const AtomType & Rhs)
{
	strcpy(Type, Rhs.Type);
  return *this;
}

int AtomType::operator==(const AtomType & Rhs) const
{
	int i =0;
	while(i<MAX) {
		if (Type[i]=='*' || Rhs.Type[i]=='*') // allow terminal wildcards
			return 1;
		else if (Type[i] != Rhs.Type[i])
			return 0;
		else if (Type[i] == '\0')
			return 1;
		i++;
	}
	return 1;	// necessary??

}

inline double anint( double x )
{
	return ceil( x - 0.50 );
}

/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------
  									
 			COORCHARMM CLASS DEFINITIONS
  								
  -----------------------------------------------------------------------
  ------------------------------------------------------------------------*/

//************** Constructors and Destructors ***************// 


// Do nothing constructor.
coords::coords()
{
	*this = default_crd;
	/*
	N=0;
	Atomno=Resno=0;                  // NULL pointers all
	Resname=Segid=Resid=0;
	X=Y=Z=0;
	Weight=0;
	Type=0;
	type = 0;
	boxlength=0;
	Bound = NONE;
	*/
}


// Default constructor (can take a number of atoms)
// not recommended for use, but included for completeness (max atoms
// held by a "large" build of Charmm.
coords::coords(int NUM)
{
	Atomno = new int[NUM];
	Resno = new int[NUM];
	Resname = new char[NUM][CHARMM_STR];
	Type = new AtomType[NUM];
	type = new string[NUM];
	X = new double[NUM];
	Y = new double[NUM];
	Z = new double[NUM];
	Segid = new char[NUM][CHARMM_STR];
	Resid = new char[NUM][CHARMM_STR];
	Weight = new double[NUM];
	Bound = NONE;
	octcharmm = false;
	N = NUM;
	boxlength = 0;
}

coords coords::default_crd = coords(100);

/*
 * count number of ATOM records in pdb file
 */
int pdb_count( istream & inp ) 
{
	const int slen = 200;
	char inpline[slen];
	int n=0;
	streampos recall = inp.tellg();
	while ( ! inp.eof() ) {
		inp.getline(inpline, slen, '\n');
		if ( ! strncmp(inpline, "ATOM", 4) && inpline[16] != 'B')
				n++; 
	}
	inp.seekg(recall);	// ??
	return n;
}

// Initialization using crdfile or pdb file ...
coords::coords(const char *infile, const char * filetype)
{
	ifstream crdfile(infile);
	if (!crdfile.good()) {
		print_err( "coords::coords()", 
				"Could not open input file");
#ifndef SWIG
		throw OpenErr();
#else
		exit(1);
#endif
	}

	const int slen = 200;
	char tmpstr[slen];
	string s;

	if (strcmp(filetype,"crd")==0) {
		//fixing, as problem if 1st character is part of 
		// the number! (large crd files)
		while (crdfile.peek() == '*') {
			crdfile.getline(tmpstr,slen,'\n');
		}

		crdfile >> N;

		Bound = NONE;
		octcharmm = false;
		boxlength = 0;
		Atomno = new int[N];
		Resno = new int[N];
		Resname = new char[N][CHARMM_STR];
		Type = new AtomType[N];
		type = new string[N];
		X = new double[N];
		Y = new double[N];
		Z = new double[N];
		Segid = new char[N][CHARMM_STR];
		Resid = new char[N][CHARMM_STR];
		Weight = new double[N];
		char atomtype_str[CHARMM_STR];
//		crdfile.getline(tmpstr,slen,'\n');
		crdfile.setf(ios::fixed | ios::right);
		crdfile.precision(5);
		for (int atom = 0; atom < N; atom++) {
			// needs to be able to handle abutting columns
			// somehow ...
			char tmpat[10], tmpres[10], tmps[10];
			crdfile >> setw(6) >> tmpat
				>> setw(6) >> tmpres;
			Atomno[atom] = atoi(tmpat);
			Resno[atom] = atoi(tmpres);
			crdfile >> Resname[atom] 
				>> type[atom]
				>> X[atom] 
				>> Y[atom] 
				>> Z[atom] 
				>> setw(5) >> Segid[atom]
				>> setw(5) >> Resid[atom] 
				>> setw(10) >> Weight[atom];
			Type[atom].Set(type[atom].c_str()); // FIXME
		}
	} else if (strcmp(filetype,"pdb")==0) {
		// LIMITED pdb support - pdb is a complex format
		// for now only reads major occupancy of partial 
		// occupancy
		N = pdb_count( crdfile );
		Bound = NONE;
		octcharmm = false;
		boxlength = 0;
		Atomno = new int[N];
		Resno = new int[N];
		Resname = new char[N][CHARMM_STR];
		Type = new AtomType[N];
		type = new string[N];
		X = new double[N];
		Y = new double[N];
		Z = new double[N];
		Segid = new char[N][CHARMM_STR];
		Resid = new char[N][CHARMM_STR];
		Weight = new double[N];
		int i = 0;
		int cresno = 0; char cresid[CHARMM_STR];
		cresid[0] = '\0';
		crdfile.close();
		crdfile.open(infile);
		crdfile.clear();
		crdfile.setf(ios::fixed | ios::right);
		crdfile.precision(5);
		int b;
		float occup;
		const char * format =
			"%2*%5i%*%4s%*%3s%2*%4i%4*%8.3f%8.3f%8.3f%6.2f%6.2f%6*%4s";
		char tmptype[5];  char alt;
		//cout << "got here at least" << endl;
		while ( ! crdfile.eof() ) {
			crdfile.getline(tmpstr, slen, '\n');
			float ftmp;
			if (strncmp(tmpstr, "ATOM", 4) == 0) {
				sscanf(tmpstr+6,"%5i", &Atomno[i]);
				strncpy(tmptype, tmpstr+12, 4 );
				tmptype[4] = '\0';
				alt = tmpstr[16];
				sscanf(tmpstr+17,"%4s", Resname[i]);
				sscanf(tmpstr+22,"%4s", Resid[i]);
				sscanf(tmpstr+30,"%8f", &ftmp); X[i] = ftmp;
				sscanf(tmpstr+38,"%8f", &ftmp); Y[i] = ftmp;
				sscanf(tmpstr+46,"%8f", &ftmp); Z[i] = ftmp;
				sscanf(tmpstr+54,"%6f", &occup);
				sscanf(tmpstr+60,"%6f", &ftmp); Weight[i]=ftmp;
				sscanf(tmpstr+72,"%4s", Segid[i]);
				/*sscanf( tmpstr, format, Atomno[i], tmptype,
						Resname[i], Resid[i],
						X[i], Y[i], Z[i],
						occup, Weight[i], Segid[i]);
						*/
				/*
				crdfile >> Atomno[i] >> setw(4) >> type[i] 
					>> Resname[i]
					>> Resid[i] >> X[i] >> Y[i] >> Z[i]
					>> s >> Weight[i] 
					>> setw(5) >> Segid[i];
				if ( ! strcmp(Segid[i],"1") ) 
					crdfile >> setw(5) >> Segid[i];
					*/
				if ( strcmp(cresid, Resid[i]) ) {
					strcpy(cresid, Resid[i]);
					cresno++;
				}
				Resno[i] = cresno;
				type[i] = string(tmptype);
				Type[i].Set(tmptype);
				//crdfile.getline(tmpstr,slen,'\n');
				if ( occup > 0.5 )
					i++;
				else if (occup == 0.5 && alt == 'A')
					i++;
			}
		}
	} else {
		print_err("coords::coords()", 
				"Unrecognised coordinate file type");
	}
	crdfile.close();
}

/* 
 *  efficient coordinate copy if the two coordinate sets have
 *  the same size
 */
void coords::copy( coords & Rhs )
{
	if (N != Rhs.N) {
		print_err( "coords::coords()", 
			"Coordinate sets have different number of atoms");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	Bound = Rhs.Bound;
	octcharmm = Rhs.octcharmm;
	boxlength = Rhs.boxlength;

	for (int atom =0; atom<N; atom++) {
		Atomno[atom] = Rhs.Atomno[atom];
		Resno[atom] = Rhs.Resno[atom];
		strcpy(Resname[atom], Rhs.Resname[atom]);
		Type[atom]=Rhs.Type[atom];		// FIXME
		type[atom] = Rhs.type[atom];
		X[atom] = Rhs.X[atom];
		Y[atom] = Rhs.Y[atom];
		Z[atom] = Rhs.Z[atom];
		strcpy(Segid[atom], Rhs.Segid[atom]);
		strcpy(Resid[atom], Rhs.Resid[atom]);
		Weight[atom] = Rhs.Weight[atom];
	}
}

/* 
 * copy a selection from set B to this crd set. nselect must = *this->N
 */
void coords::copy_selection( int * sel, int nsel, coords & Rhs )
{
	if (N != nsel) {
		print_err( "coords::coords()", 
			"Coordinate size mismatch");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	Bound = Rhs.Bound;
	octcharmm = Rhs.octcharmm;
	boxlength = Rhs.boxlength;

	for (int i =0; i<nsel; i++) {
		Atomno[i] = Rhs.Atomno[sel[i]];
		Resno[i] = Rhs.Resno[sel[i]];
		strcpy(Resname[i], Rhs.Resname[sel[i]]);
		Type[i]=Rhs.Type[sel[i]];		// FIXME
		type[i] = Rhs.type[sel[i]];
		X[i] = Rhs.X[sel[i]];
		Y[i] = Rhs.Y[sel[i]];
		Z[i] = Rhs.Z[sel[i]];
		strcpy(Segid[i], Rhs.Segid[sel[i]]);
		strcpy(Resid[i], Rhs.Resid[sel[i]]);
		Weight[i] = Rhs.Weight[sel[i]];
	}
}

// Copy Constructor
coords::coords(const coords& Rhs)
{
	*this = Rhs;		// see defn below
}

/*
coords::coords(const coords& Rhs)
{
	N=Rhs.N;
	Bound = Rhs.Bound;
	boxlength = Rhs.boxlength;

	Atomno = new int[N];
	Resno = new int[N];
	Resname = new char[N][CHARMM_STR];
	Type = new AtomType[N];
	type = new string[N];
	X = new double[N];
	Y = new double[N];
	Z = new double[N];
	Segid = new char[N][CHARMM_STR];
	Resid = new char[N][CHARMM_STR];
	Weight = new double[N];

	for (int atom =0; atom<N; atom++) {
		Atomno[atom] = Rhs.Atomno[atom];
		Resno[atom] = Rhs.Resno[atom];
		strcpy(Resname[atom], Rhs.Resname[atom]);
		Type=Rhs.Type;		// FIXME
		type[atom] = Rhs.type[atom];
		X[atom] = Rhs.X[atom];
		Y[atom] = Rhs.Y[atom];
		Z[atom] = Rhs.Z[atom];
		strcpy(Segid[atom], Rhs.Segid[atom]);
		strcpy(Resid[atom], Rhs.Resid[atom]);
		Weight[atom] = Rhs.Weight[atom];
	}
}
*/

// Copy Constructor with selection
// Sun May 27 07:18:34 BST 2001 - RB changed to renumber atoms and
//				  residues correctly!
coords::coords(const coords& Rhs, int * selection, int nselect,
		bool resrenum)
{
	N=nselect;
	Bound = Rhs.Bound;
	octcharmm = Rhs.octcharmm;
	boxlength = Rhs.boxlength;

	Atomno = new int[nselect];
	Resno = new int[nselect];
	Resname = new char[nselect][CHARMM_STR];
	Type = new AtomType[nselect];
	type = new string[nselect];
	X = new double[nselect];
	Y = new double[nselect];
	Z = new double[nselect];
	Segid = new char[nselect][CHARMM_STR];
	Resid = new char[nselect][CHARMM_STR];
	Weight = new double[nselect];

	for (int atom =0; atom<nselect; atom++) {
		Atomno[atom] = Rhs.Atomno[selection[atom]];
		Resno[atom] = Rhs.Resno[selection[atom]];
		strcpy(Resname[atom], Rhs.Resname[selection[atom]]);
		Type[atom]=Rhs.Type[selection[atom]];
		type[atom] = Rhs.type[selection[atom]];
		X[atom] = Rhs.X[selection[atom]];
		Y[atom] = Rhs.Y[selection[atom]];
		Z[atom] = Rhs.Z[selection[atom]];
		strcpy(Segid[atom], Rhs.Segid[selection[atom]]);
		strcpy(Resid[atom], Rhs.Resid[selection[atom]]);
		Weight[atom] = Rhs.Weight[selection[atom]];
	}
	// new renumbering bit
	int nures = 1; int lastres = Resno[0]; int curres = 1;
	for (int atom =0; atom<nselect; atom++) {
		Atomno[atom] = atom+1;
		if (resrenum) {
			if ( Resno[atom] != lastres ) {		// new residue
				curres++;
				lastres = Resno[atom];
			}
			Resno[atom] = curres;
			//strcpy(Resid[atom], itos(curres));
			sprintf(Resid[atom],"%4i",curres);
		}
	}
}

// die die die!
coords::~coords()
{
	delete [] Atomno;
	delete [] Resno;
	delete [] Resname;
	delete [] Type;
	delete [] type;
	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] Segid;
	delete [] Resid;
	delete [] Weight;
}

//******************* General Functions *******************//

int coords::Coorwrite(const char *outfile, const char *comment)
{
	ofstream outp(outfile);

	if (!outp.good()) {
		cerr << "\n!! Could not open output file " << outfile 
			<< " !!" << endl;
		return 1;
	}

	outp << "* " << comment << "\n* \n" << *this;
	outp.close();
	return 0;
}

double coords::xdiff(int i, int j)
{ 
	double xdif= X[j]-X[i]; 
	double ydif, zdif;

	switch (Bound) {
		case CUBIC:
			return xdif - boxlength * anint ( xdif / boxlength );
			break;
		case OCTAHEDRAL:
			double xd, yd, zd;
			xd= X[j]-X[i]; 
			yd = Y[j]-Y[i];
			zd = Z[j]-Z[i];
			if (octcharmm) {
				linalg::multiply_Ax( orot, xd, yd, zd, xdif, ydif,
						zdif); 
			} else {
				xdif = xd; ydif=yd; zdif=zd;
			}
			xdif -=  boxlength* anint ( xdif / boxlength );
			ydif -=  boxlength* anint ( ydif / boxlength );
			zdif -=  boxlength* anint ( zdif / boxlength );

			if((fabs(xdif)+fabs(ydif)+fabs(zdif))
					>=1.5*boxlength/2.0) {
				if (xdif>=0)
					xdif -= boxlength/2;
				else xdif += boxlength/2;
			}
			return xdif;
			break;
		case NONE:
			return xdif;
			break;
	}
}

double coords::ydiff(int i, int j)
{ 
	double ydif = Y[j]-Y[i]; 
	double xdif, zdif;

	switch (Bound)
	{
		case CUBIC:
			return ydif - boxlength * anint ( ydif / boxlength );
			break;
		case OCTAHEDRAL:
			double xd, yd, zd;
			xd= X[j]-X[i]; 
			yd = Y[j]-Y[i];
			zd = Z[j]-Z[i];
			if (octcharmm) {
				linalg::multiply_Ax( orot, xd, yd, zd, xdif, ydif,
						zdif); 
			} else {
				xdif = xd; ydif=yd; zdif=zd;
			}
			xdif -=  boxlength* anint ( xdif / boxlength );
			ydif -=  boxlength* anint ( ydif / boxlength );
			zdif -=  boxlength* anint ( zdif / boxlength );

			if((fabs(xdif)+fabs(ydif)+fabs(zdif))
					>=1.5*boxlength/2.0) {
				if (ydif>=0)
					ydif -= boxlength/2;
				else ydif += boxlength/2;
			}
			return ydif;
			break;
		case NONE:
			return ydif;
			break;
	}
}

double coords::zdiff(int i, int j)
{ 
	double zdif = Z[j]-Z[i]; 
	double xdif, ydif;

	switch (Bound) {
		case CUBIC:
			return zdif - boxlength * anint ( zdif / boxlength );
			break;
		case OCTAHEDRAL:
			double xd, yd, zd;
			xd= X[j]-X[i]; 
			yd = Y[j]-Y[i];
			zd = Z[j]-Z[i];
			if (octcharmm) {
				linalg::multiply_Ax( orot, xd, yd, zd, xdif, ydif,
						zdif); 
			} else {
				xdif = xd; ydif=yd; zdif=zd;
			}
			xdif -=  boxlength* anint ( xdif / boxlength );
			ydif -=  boxlength* anint ( ydif / boxlength );
			zdif -=  boxlength* anint ( zdif / boxlength );

			if((fabs(xdif)+fabs(ydif)+fabs(zdif))
					>=1.5*boxlength/2.0)
			{ 
				if (zdif>=0)
					zdif -= boxlength/2;
				else zdif += boxlength/2;
			}
			return zdif;
			break;
		case NONE:
			return zdif;
			break;
	}
}

void coords::diffvec(int i, int j, double & xdiff, 
		double & ydiff,	double & zdiff )
{ 
	xdiff = X[j]-X[i];
	ydiff = Y[j]-Y[i];
	zdiff = Z[j]-Z[i];

	double aint_x_l,aint_y_l,aint_z_l;

	double dis=0;
	switch(Bound) {
		case CUBIC:
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			//cubic periodic boundaries
			xdiff -=  boxlength*aint_x_l;
			ydiff -=  boxlength*aint_y_l;
			zdiff -=  boxlength*aint_z_l;
			break;

		case OCTAHEDRAL:
			double xd, yd, zd;
			if (octcharmm) {
				linalg::multiply_Ax( orot, xdiff, ydiff, zdiff, xd, 
						yd, zd); 
				xdiff=xd; ydiff=yd; zdiff=zd;
			} 
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff =  xdiff - boxlength*aint_x_l;
			ydiff =  ydiff - boxlength*aint_y_l;
			zdiff =  zdiff - boxlength*aint_z_l;

			if((fabs(xdiff)+fabs(ydiff)+fabs(zdiff))>=1.5*boxlength/2.0)
			{ 
				if (xdiff>=0)
					xdiff= xdiff - boxlength/2;
				else xdiff= xdiff + boxlength/2;

				if (ydiff>=0)
					ydiff= ydiff - boxlength/2;
				else ydiff= ydiff + boxlength/2;

				if (zdiff>=0)
					zdiff= zdiff - boxlength/2;
				else zdiff= zdiff + boxlength/2;
			}
			break;
		case NONE: 
			break;

	}

	return;
}

double coords::dist(int i, int j)
{ 
	double xdiff = X[j]-X[i];
	double ydiff = Y[j]-Y[i];
	double zdiff = Z[j]-Z[i];

	double aint_x_l,aint_y_l,aint_z_l;

	double dis=0;
	switch(Bound) {
		case CUBIC:
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff -=  boxlength*aint_x_l;
			ydiff -=  boxlength*aint_y_l;
			zdiff -=  boxlength*aint_z_l;
			break;

		case OCTAHEDRAL:
			double xd, yd, zd;
			if (octcharmm) {
				linalg::multiply_Ax( orot, xdiff, ydiff, zdiff, xd, 
						yd, zd); 
				xdiff=xd; ydiff=yd; zdiff=zd;
			} 
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff -= boxlength*aint_x_l;
			ydiff -= boxlength*aint_y_l;
			zdiff -= boxlength*aint_z_l;

			if((fabs(xdiff)+fabs(ydiff)+fabs(zdiff))
					>= 1.5*boxlength/2.0) {
				if (xdiff>=0)
					xdiff= xdiff - boxlength/2;
				else xdiff= xdiff + boxlength/2;

				if (ydiff>=0)
					ydiff= ydiff - boxlength/2;
				else ydiff= ydiff + boxlength/2;

				if (zdiff>=0)
					zdiff= zdiff - boxlength/2;
				else zdiff= zdiff + boxlength/2;

			}

			break;

	}

	dis = sqrt(xdiff*xdiff +  ydiff*ydiff + zdiff*zdiff); 
	return dis;
}

void coords::vector_from_point(double COM[3], int i, double & xdiff, 
		double & ydiff, double & zdiff) const
{ 
	xdiff = X[i]-COM[0];
	ydiff = Y[i]-COM[1];
	zdiff = Z[i]-COM[2];

	double aint_x_l,aint_y_l,aint_z_l;

	double dis=0;
	switch(Bound) {
		case CUBIC:

			//cubic periodic boundaries
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff -=  boxlength*aint_x_l;
			ydiff -=  boxlength*aint_y_l;
			zdiff -=  boxlength*aint_z_l;
			break;

		case OCTAHEDRAL:
			double xd, yd, zd;
			if (octcharmm) {
				linalg::multiply_Ax( orot, xdiff, ydiff, zdiff, xd, 
						yd, zd); 
				xdiff=xd; ydiff=yd; zdiff=zd;
			} 
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff =  xdiff - boxlength*aint_x_l;
			ydiff =  ydiff - boxlength*aint_y_l;
			zdiff =  zdiff - boxlength*aint_z_l;

			if((fabs(xdiff)+fabs(ydiff)+fabs(zdiff))
					>= 1.5*boxlength/2.0) {
				if (xdiff>=0)
					xdiff= xdiff - boxlength/2;
				else xdiff= xdiff + boxlength/2;

				if (ydiff>=0)
					ydiff= ydiff - boxlength/2;
				else ydiff= ydiff + boxlength/2;

				if (zdiff>=0)
					zdiff= zdiff - boxlength/2;
				else zdiff= zdiff + boxlength/2;
			}
			break;
	}

	return;
}

double coords::dist_from_point( double cm[3],int i)
{ 
	double xdiff = X[i]-cm[0];
	double ydiff = Y[i]-cm[1];
	double zdiff = Z[i]-cm[2];

	double aint_x_l,aint_y_l,aint_z_l;

	double dis=0;
	switch(Bound) {
		case CUBIC:
			//cubic periodic boundaries
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff =  xdiff - boxlength*aint_x_l;
			ydiff =  ydiff - boxlength*aint_y_l;
			zdiff =  zdiff - boxlength*aint_z_l;
			break;
		case OCTAHEDRAL:
			double xd, yd, zd;
			if (octcharmm) {
				linalg::multiply_Ax( orot, xdiff, ydiff, zdiff, xd, 
						yd, zd); 
				xdiff=xd; ydiff=yd; zdiff=zd;
			} 
			aint_x_l = anint( (xdiff)/boxlength );
			aint_y_l = anint( (ydiff)/boxlength );
			aint_z_l = anint( (zdiff)/boxlength );
			xdiff =  xdiff - boxlength*aint_x_l;
			ydiff =  ydiff - boxlength*aint_y_l;
			zdiff =  zdiff - boxlength*aint_z_l;
			if( (fabs(xdiff)+fabs(ydiff)+fabs(zdiff)) 
					>= 1.5*boxlength/2.0 ) {
				if (xdiff>=0)
					xdiff= xdiff - boxlength/2;
				else xdiff= xdiff + boxlength/2;

				if (ydiff>=0)
					ydiff= ydiff - boxlength/2;
				else ydiff= ydiff + boxlength/2;

				if (zdiff>=0)
					zdiff= zdiff - boxlength/2;
				else zdiff= zdiff + boxlength/2;

			}
			break;
	}

	dis = sqrt(xdiff*xdiff +  ydiff*ydiff + zdiff*zdiff); 
	return dis;
}

//m. kuttel - cribbed from KJN's calcdihe.f
// need to fix for actahedral periodic boundary conditions
double coords::dihedral_angle(int at1, int at2, int at3, int at4)
{
	const double RMIN =0.0001;
	const double RXMIN=0.005;

	//...Vector Aij
	double FX =  xdiff(at1,at2);
	double FY =  ydiff(at1,at2);
	double FZ =  zdiff(at1,at2);

	// ...Vector Ajk and Bjk
	double GX =  xdiff(at2,at3);
	double GY =  ydiff(at2,at3);
	double GZ =  zdiff(at2,at3);

	//...Vector Blk
	double HX =  xdiff(at4,at3);
	double HY =  ydiff(at4,at3);
	double HMMMM =  zdiff(at4,at3);

	//Calc. normal to planes
	//  normal_ijk : A = Aij cross Ajk
	//  normal_lkj : B = Blk cross Bjk

	//...Normal to plane A : plane_ijk
	double AX=FY*GZ-FZ*GY;
	double AY=FZ*GX-FX*GZ;
	double AZ=FX*GY-FY*GX;

	// ...Normal to plane B : plane_lkj
	double BX=HY*GZ-HMMMM*GY;
	double BY=HMMMM*GX-HX*GZ;
	double BZ=HX*GY-HY*GX;

	//---------- Calc. the dihedral angle -------------------

	double RA2=AX*AX+AY*AY+AZ*AZ;
	double RB2 =BX*BX+BY*BY+BZ*BZ;
	double RA=sqrt(RA2);
	double RB=sqrt(RB2);

	//---------- Check for Linear Angles ---------------------
	if(RA<=RXMIN) 
		RA = RMIN;

	double  RAR=1.0/RA;

	if(RB<=RXMIN) 
		RB = RMIN;

	double  RBR=1.0/RB;

	double AXR=AX*RAR;
	double AYR=AY*RAR;
	double AZR=AZ*RAR;
	double BXR=BX*RBR;
	double BYR=BY*RBR;
	double BZR=BZ*RBR;

	double CT=AXR*BXR+AYR*BYR+AZR*BZR;


	if (CT>1.00)  CT =  1.00;
	if (CT<-1.00) CT = -1.00;

	double AP = acos(CT);

	double  CX=AYR*BZR-AZR*BYR;
	double CY=AZR*BXR-AXR*BZR; 
	double CZ=AXR*BYR-AYR*BXR;  

	double  ST=sqrt(CX*CX + CY*CY + CZ*CZ);		
	double  S=GX*CX + GY*CY + GZ*CZ; 

	if (S< 0.0) {
		AP = -AP;
		ST = -ST;
	}

	double tors = (180.0/M_PI) * AP;
	return tors;

}

double coords::angle(int at1, int at2, int at3)
{
	double theta;

	const double RMIN =0.0001;
	const double RXMIN=0.005;

	//...Vector rij
	double AX =  xdiff(at2,at1);
	double AY =  ydiff(at2,at1);
	double AZ =  zdiff(at2,at1);

	// ...Vector rjk 
	double BX =  xdiff(at2,at3);
	double BY =  ydiff(at2,at3);
	double BZ =  zdiff(at2,at3);

	double dotp = AX*BX+AY*BY+AZ*BZ;
	double crossx = AY*BZ-AZ*BY;
	double crossy = AX*BZ-AZ*BX;
	double crossz = AX*BY-AY*BX;

	double normcross = sqrt(crossx*crossx+crossy*crossy+crossz*crossz);
	double norma = sqrt(AX*AX+AY*AY+AZ*AZ);
	double normb = sqrt(BX*BX+BY*BY+BZ*BZ);

	double costheta = dotp/(norma*normb);
	double sintheta = normcross/(norma*normb);

	if (costheta >= 0.) {
		theta = asin(sintheta);
	} else {
		theta = asin(sintheta);
		if (sintheta >=0.) {
			theta = M_PI - theta;
		} else {
			theta = -M_PI - theta;
		}
	}
	
	return theta*180./M_PI;

}

void coords::center_of_mass(int atom_list[], 
		int no_atoms, 
		double mass[], 
		double CM[3])
{
	CM[0] = 0;
	CM[1] = 0;
	CM[2] = 0;
	double total_mass=0;

	for(int i=0; i < no_atoms; i++) {
		CM[0] += mass[i]*X[atom_list[i]];
		CM[1] += mass[i]*Y[atom_list[i]];
		CM[2] += mass[i]*Z[atom_list[i]];
		total_mass += mass[i];
	}
	if (total_mass == 0) {
		print_err ( "coords::center_of_mass()", 
				"Error: division by zero in centre\
				of mass calculation. \
				Did you initialize the masses for your atoms?");
#ifndef SWIG
			throw DivZeroErr();
#else
			exit(1);
#endif
	}

	CM[0]=CM[0]/total_mass;
	CM[1]=CM[1]/total_mass;
	CM[2]=CM[2]/total_mass;  
}

int coords::Coorread(const char *infile, const char * type)
{
	if (N == 0) {
		*this = coords(infile);
		return 0;
	} else {
		ifstream inpf(infile);
		if (!inpf.good()) {
			print_err( string("coords::Coorread()"), 
					string("Could not open input file")
					+ string(infile)); 
			return 1;	// throw OpenErr();
		}
		inpf >> *this;
		inpf.close();
		return 0;
	} 
}

//******************* OPERATORS *************************//

// To read in coords from a CRDfile to an existing coords object.
// If default coords() operator used (i.e. no free store allocated), 
// this is done.

// FIXME: in the case that no crds defined, this should read 
// in all the info, not just coordinates.
istream& operator>> (istream& inpstr, coords& coors)
{
	char tmpstr[150];
	int new_N;

	while (inpstr.get() == '*') {
		inpstr.getline(tmpstr,150,'\n');
	}

	inpstr >> new_N;
	if (new_N != coors.N) {   // BEWARE: this is the only check performed!
		print_err( "operator>> (istream& inpstr, coords& coors)", 
				"New coordinate set has different number of atoms \
				from old !!");
#ifndef SWIG
				throw coords::WrongNumAtomsErr();
#else
				exit(1);
#endif
	}

	for (int atom = 0; atom < new_N; atom++) {
		inpstr >> tmpstr >> tmpstr >> tmpstr >> tmpstr
			>> coors.X[atom] >> coors.Y[atom] >> coors.Z[atom] 
			>> tmpstr
			>> tmpstr >> coors.Weight[atom];
	}
	return inpstr;
}

// Standard output stream operator.
ostream& operator<< (ostream& outstr, const coords& coors)
{
//#ifdef GCC295
//	long orig_flags = outstr.flags();   
//#else
	ios_base::fmtflags orig_flags = outstr.flags();   
//#endif
			// save original output format

	outstr.setf(ios::right);
	outstr << "* TITLE GOES HERE \n*" << endl;
	outstr << setw(5) << coors.N << endl;

	outstr.setf(ios::fixed | ios::right);
	outstr.precision(5);

	for (int atom=0; atom<coors.N; atom++) {
		outstr << setw(5) << coors.Atomno[atom] 
			<< setw(5) << coors.Resno[atom] << " ";
		outstr.setf(ios::left, ios::adjustfield);
		outstr << setw(4) << coors.Resname[atom] << " " 
			<< setw(4) << coors.type[atom].c_str(); 
		outstr.setf(ios::right, ios::adjustfield);
		outstr << setw(10) << coors.X[atom] << setw(10) << 
			coors.Y[atom] << setw(10)
			<< coors.Z[atom] << " ";
		outstr.setf(ios::left, ios::adjustfield);
		outstr << setw(4) << coors.Segid[atom] << " " 
			<< setw(4) << coors.Resid[atom];
		outstr.setf(ios::right, ios::adjustfield);
		outstr << setw(10) << coors.Weight[atom] << endl;    
	}

	outstr.flags(orig_flags);           // reset to original output format
	return outstr;                     
}

void operator>> (BaseITrajFile& dcdfile, coords& coors)
{
	//fprintf(stderr," in here now...1\n"); fflush(stderr);
	dcdfile.read_frame( coors.X, coors.Y, coors.Z, coors.N );
	//fprintf(stderr," in here now...2\n"); fflush(stderr);
	if ( dcdfile.crystal() ) {
		dcdfile.get_crystal_data( coors.Xtal_data);
		switch (coors.Bound) {
			case CUBIC:
				coors.boxlength = coors.Xtal_data[0];
				break;
			case OCTAHEDRAL:
				if (coors.octcharmm) {
					coors.boxlength = coors.Xtal_data[0] 
						* 1.2;
				} else {
					coors.boxlength = coors.Xtal_data[0];
				}
				break;
			default:
				coors.Bound = CUBIC;
				coors.boxlength = coors.Xtal_data[0];
				break;
				/*
			default:
				print_err( "operator>> \
(BaseITrajFile& dcdfile, coords& coors)", 
					"WARNING: you're \
					reading a crystal dcd file into a \
					coords with no boundary conditions \
					set!!!");
#ifndef SWIG
				throw coords::BoundaryErr();
#else
				exit(1);
#endif
*/
		}
	}
	//fprintf(stderr," in here now...3\n"); fflush(stderr);
}

// dcd output is a lot more low-level of necessity; beware!
void operator<< (coords& coors, BaseOTrajFile& dcdfile)
{
	dcdfile.write_frame( coors.X, coors.Y, coors.Z, coors.N );
}

ostream& operator<< (ostream & outstr , const AtomType & Tp)
{
	outstr << setw(4) << Tp.Type;
	return outstr;    
}

/* 
 * copies rhs to lhs, no questions asked!!
 * rather use more efficient copy function if you know
 * the coordinate sets are the same size
 */
coords & coords::operator= (const coords & rhs)
{
	N=rhs.N;
	Bound = rhs.Bound;
	boxlength = rhs.boxlength;

	Atomno = new int[N];
	Resno = new int[N];
	Resname = new char[N][CHARMM_STR];
	Type = new AtomType[N];
	type = new string[N];
	X = new double[N];
	Y = new double[N];
	Z = new double[N];
	Segid = new char[N][CHARMM_STR];
	Resid = new char[N][CHARMM_STR];
	Weight = new double[N];

	for (int atom =0; atom<N; atom++) {
		Atomno[atom] = rhs.Atomno[atom];
		Resno[atom] = rhs.Resno[atom];
		strcpy(Resname[atom], rhs.Resname[atom]);
		Type[atom]=rhs.Type[atom];		// FIXME
		type[atom] = rhs.type[atom];
		X[atom] = rhs.X[atom];
		Y[atom] = rhs.Y[atom];
		Z[atom] = rhs.Z[atom];
		strcpy(Segid[atom], rhs.Segid[atom]);
		strcpy(Resid[atom], rhs.Resid[atom]);
		Weight[atom] = rhs.Weight[atom];
	}
	return *this;
}

/* 
 * adds rhs to lhs, if coord set same size!
 */
coords & coords::operator+= (const coords & rhs)
{
	if (N != rhs.N) {
		print_err( "coords::coords()", 
				"Coords being added of different size!!!");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	for (int atom =0; atom<N; atom++) {
		X[atom] += rhs.X[atom];
		Y[atom] += rhs.Y[atom];
		Z[atom] += rhs.Z[atom];
	}
	return *this;
}

/* 
 * adds const rhs to lhs, if coord sets same size!
 */
coords & coords::operator*= (const double & rhs)
{
	for (int atom =0; atom<N; atom++) {
		X[atom] *= rhs;
		Y[atom] *= rhs; 
		Z[atom] *= rhs;
	}
	return *this;
}

// rms difference of entire coord sets
//
// does not account for periodic boundaries, since in all the
// cases I can think of this will not be necessary, and they are
// trivial to add
double rms ( const coords & set1, const coords & set2 )
{
	if ( set1.N != set2.N ) {
		print_err( "rms ( const coords & set1, const coords & set2 )", 
				"Coordinate sets being fitted have different\
				numbers of atoms");
#ifndef SWIG
		throw coords::WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	
	double RMS = 0.0;
	
	for (int i = 0; i<set1.N; i++) {
		double dx = set1.X[i] - set2.X[i];
		double dy = set1.Y[i] - set2.Y[i];
		double dz = set1.Z[i] - set2.Z[i];
		RMS += dx*dx + dy*dy + dz*dz;
	}
	
	RMS = sqrt(RMS/double(set1.N));
	
	return RMS;
	
}
	
// rms difference of a selection of atoms
double rms ( const coords & set1, const coords & set2,
		int * selection, int selected )
{
	if ( set1.N != set2.N ) {
		print_err( "rms ( const coords & set1, \
					const coords & set2, \
				       	int * selection, int selected )", 
				"Coordinate sets being fitted have different\
				numbers of atoms");
#ifndef SWIG
		throw coords::WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	
	double RMS;
	RMS = 0.0;
	
	for (int i = 0; i<selected; i++) {
		double dx = set1.X[selection[i]] - set2.X[selection[i]];
		double dy = set1.Y[selection[i]] - set2.Y[selection[i]];
		double dz = set1.Z[selection[i]] - set2.Z[selection[i]];
		RMS += dx*dx + dy*dy + dz*dz;
	}
	
	RMS = sqrt(RMS/double(selected));
	
	return RMS;
	
}

/*
 * calculates rms deviation of each atom in selection 
 */
double sel_rms2 ( const coords & set1, 
				const coords & set2,
				int * selection, int selected, 
				double * rmssel )
{

	if ( set1.N != set2.N ) {
		print_err( "sel_rms ( const coords & set1, \
					const coords & set2, \
				       	int * selection, int selected \
					double * rmssel)", 
				"Coordinate sets being fitted have different\
				numbers of atoms");
#ifndef SWIG
		throw coords::WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	
	double RMS2 = 0; double iRMS2;
	
	for (int i = 0; i<selected; i++) {
		double dx = set1.X[selection[i]] - set2.X[selection[i]];
		double dy = set1.Y[selection[i]] - set2.Y[selection[i]];
		double dz = set1.Z[selection[i]] - set2.Z[selection[i]];
		iRMS2 = dx*dx + dy*dy + dz*dz;
		RMS2 += iRMS2;
		rmssel[i] = iRMS2;
	}
	
	RMS2 = RMS2/double(selected);
	
	return RMS2;
}

// Assigns atomic masses of common biological elements based
// on first letter of atom name ! (works mostly, though)
void coords::assign_masses()
{
	for (int i=0; i< N; i++) {
		switch (type[i][0]) {
			case 'C':
			case 'c':
				Weight[i] = 12.011;
				break;
			case 'O':
			case 'o':
				Weight[i] = 15.9994;
				break;
			case 'H':
			case 'h':
				Weight[i] = 1.00794;
				break;
			case 'N':
			case 'n':
				Weight[i] = 14.007;
				break;
			case 'S':
			case 's':
				Weight[i] = 32.066;
				break;
			case 'P':
			case 'p':
				Weight[i] = 30.974;
				break;
			default:
				Weight[i] = 0.0;
		}
	}
}

void coords::center_of_mass(double CM[3])
{
	double cx, cy, cz, tot_mass, com, tmp_mass;
	tot_mass = 0;
	cx =0; cy = 0; cz = 0;
	com = 0;

	for (int i=0; i<N; i++) {
		tot_mass += tmp_mass = Weight[i];
		cx += X[i]*tmp_mass;
		cy += Y[i]*tmp_mass;
		cz += Z[i]*tmp_mass;
	}
	
	CM[0] = cx/tot_mass;
	CM[1] = cy/tot_mass;
	CM[2] = cz/tot_mass;
}

void coords::center_of_mass(double CM[3], int * selection, int n)
{
	double cx, cy, cz, tot_mass, com, tmp_mass;
	tot_mass = 0;
	cx =0; cy = 0; cz = 0;
	com = 0;

	for (int i=0; i<n; i++)
	{
		tot_mass += tmp_mass = Weight[selection[i]];
		cx += X[selection[i]]*tmp_mass;
		cy += Y[selection[i]]*tmp_mass;
		cz += Z[selection[i]]*tmp_mass;
	}
	
	CM[0] = cx/tot_mass;
	CM[1] = cy/tot_mass;
	CM[2] = cz/tot_mass;
}

void coords::assign_masses( int * massarray, int n )
{
	if ( n != N )
	{
		print_err( "coords::assign_masses()", 
				"Mass Array has the wrong size!");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	
	for (int i=0; i<n; i++)
	{
		Weight[i] = massarray[i];
	}
}

void coords::assign_masses ( int * massarray, int * selection,
				int nselect )
{
	if ( nselect >= N )
	{
		print_err( "coords::assign_masses()", 
				"Mass Array too big!");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}
	
	for (int i=0; i<nselect; i++)
	{
		Weight[selection[i]] = massarray[i];
	}
}

void coords::set_uniform_weights()
{
	for (int i=0; i<N; i++)
	{
		Weight[i] = 1.0;
	}
}

// Radius of gyration
double coords::rgyr()
{
	double gx, gy, gz;
	double diff;
	gx = 0; gy = 0; gz = 0;
	double tot_mass = 0;
	double tmp_mass;
	double CM[3];
	center_of_mass( CM );
	
	for ( int i = 0; i<N; i++)
	{
		tot_mass += tmp_mass = Weight[i];
		diff = X[i] - CM[0];
		gx += tmp_mass * diff * diff;
		diff = Y[i] - CM[1];
		gy += tmp_mass * diff * diff;
		diff = Z[i] - CM[2];
		gz += tmp_mass * diff * diff;
	}
	
	return sqrt( (gx + gy + gz)/tot_mass );
}
		
double coords::rgyr( int * selection, int n )
{
	double gx, gy, gz;
	double diff;
	gx = 0; gy = 0; gz = 0;
	double tot_mass = 0;
	double tmp_mass;
	double CM[3];
	center_of_mass( CM, selection, n );
	
	for ( int i = 0; i<n; i++)
	{
		tot_mass += tmp_mass = Weight[selection[i]];
		diff = X[selection[i]] - CM[0];
		gx += tmp_mass * diff * diff;
		diff = Y[selection[i]] - CM[1];
		gy += tmp_mass * diff * diff;
		diff = Z[selection[i]] - CM[2];
		gz += tmp_mass * diff * diff;
	}
	
	return sqrt( (gx + gy + gz)/tot_mass );
}

/*
double coords::dipole( int * selection, int n, double dip[3] )
{
	dip[0] = 0; dip[1] = 0; dip[2] = 0; 

}
*/

/*   woo!  */
void coords::translate( double x, double y, double z )
{
	for (int i = 0; i< N; i++)
	{	
		X[i] += x;
		Y[i] += y;
		Z[i] += z;
	}
}

void coords::translate( double x, double y, double z,
			 int * selection, int nselect )
{
	for (int i = 0; i< nselect; i++)
	{	
		X[selection[i]] += x;
		Y[selection[i]] += y;
		Z[selection[i]] += z;
	}
}

/*   WOO!  */
void coords::rotate( double rotmat[3][3], double &x, double &y, double &z )
{
	double tmpx, tmpy, tmpz;

	tmpx = rotmat[0][0] * x + rotmat[0][1] * y + rotmat[0][2] * z;
	tmpy = rotmat[1][0] * x + rotmat[1][1] * y + rotmat[1][2] * z;
	tmpz = rotmat[2][0] * x + rotmat[2][1] * y + rotmat[2][2] * z;
	
	x=tmpx; y=tmpy; z=tmpz;
}

void coords::rotate( double rotmat[3][3] )
{
	for (int i=0; i<N; i++)
	{
		rotate( rotmat, X[i], Y[i], Z[i] );
	}
}

void coords::rotate( double rotmat[3][3], int * selection, int nselect )
{
	for (int i=0; i<nselect; i++)
	{
		rotate( rotmat, X[selection[i]], Y[selection[i]], Z[selection[i]] );
	}
}

/* NOTE: 
   A lot of this stuff has been put in the linalg library which should
   preferably be used in future
   */
void normalize ( double vec[3] )
{
	double mod = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	for (int x=0; x<3; x++)
		vec[x] /= mod;
}

void normalize ( double & vecx, double & vecy, double & vecz )
{
	double mod = sqrt( vecx*vecx + vecy*vecy + vecz*vecz);
	vecx /= mod;
	vecy /= mod;
	vecz /= mod;
}

// y = Mx
void matrix_mult ( double M[3][3], double x[3], double y[3] )
{
	for ( int i=0; i<3; i++ )
	{
		y[i]=0;
		for (int j=0; j<3; j++)
			y[i] += M[i][j] * x[j];
	}
}

// A = MN
void matrix_mult ( double M[3][3], double N[3][3], double A[3][3] )
{
	for (int i=0; i<3; i++)
	{
		for ( int j=0; j<3; j++ )
		{
			A[i][j]=0;
			for (int k=0; k<3; k++)
				A[i][j] += M[i][k] * N[k][j];
		}
	}
}

void colswap( double M[3][3], int a, int b)
{
	double tmp;
	for (int i=0; i<3; i++)
	{
		tmp = M[i][a];
		M[i][a]=M[i][b];
		M[i][b]=tmp;
	}
}

void print_matrix( double M[3][3] )
{
	for (int i=0;i<3;i++)
	{
		cout << M[i][0] << '\t' << M[i][1] 
			<< '\t' << M[i][2] << endl;
	}
}

/************************************************************************
 *	kabsch_helper()							*
 * -------------------------------------------------------------------- *
 *	This function actually does the work for the kabsch fitting	*
 *	procedure after the functions below have made the appropriate	*
 *	atom selections and calculated the R matrix			*
 ************************************************************************/
double kabsch_helper( double R[3][3], double U[3][3] )
{
	double Rt[3][3], RtR[3][3];
	double low_thresh = pow(10.0,-5.0);

	// following kabsch, 
	// b) form RtR...
	for (int i=0; i<3; i++) 
	{
		for (int j=0; j<3; j++) 
		{
			Rt[i][j] = R[j][i];
		}
	}

	matrix_mult ( Rt, R, RtR );

#ifdef DEBUG
	cout << "R matrix is\n" << endl;
	print_matrix(R);
	cout << "det(R) = " << linalg::det(R) << endl;

	cout << "RtR matrix is\n" << endl;
	print_matrix(RtR);
#endif

	// catch 1-atom selection before we waste time
	// (i.e. all eigenvalues zero, tr(RtR = 0))
	double tr = RtR[0][0] + RtR[1][1] + RtR[2][2];
	if ( tr < low_thresh ) {
		linalg::ident_mat(U);
		return 0;
	}

	// ... determine e'vals and e'vecs
	// b) 2) find the eigenvalues and eigenvectors
	double evalue[3];
	double a[3][3];
	if(linalg::jacobi(RtR,evalue,a) != 0) return -1;
#ifdef DEBUG
	cout << "det(A) = " << linalg::det(a) << endl;
#endif
	// bubblesort ;-)
	if (evalue[0] < evalue[1]) 
		colswap(a, 0, 1);
	if (evalue[0] < evalue[2]) 
		colswap(a,0, 2);
	if (evalue[1] < evalue[2]) 
		colswap(a, 1, 2);

#ifdef DEBUG
	cout << "e'vals " << evalue[0] << " " << evalue[1] 
		<< " " << evalue[2] << endl;
#endif

	// ensure right-handed system a[-][2] = a[-][0] x a[-][1]
	a[0][2] = a[1][0]*a[2][1] - a[2][0]*a[1][1]; 
	a[1][2] = a[2][0]*a[0][1] - a[0][0]*a[2][1];
	a[2][2] = a[0][0]*a[1][1] - a[1][0]*a[0][1]; 

	// c) determine b(i) = R*a(i)
	double b[3][3];

	matrix_mult ( R, a, b );

	normalize ( b[0][0], b[1][0], b[2][0] );
	normalize ( b[0][1], b[1][1], b[2][1] );

	// ensure right-handed system b[-][2] = b[-][0] x b[-][1]
	b[0][2] = b[1][0]*b[2][1] - b[2][0]*b[1][1]; 
	b[1][2] = b[2][0]*b[0][1] - b[0][0]*b[2][1];
	b[2][2] = b[0][0]*b[1][1] - b[1][0]*b[0][1]; 

	// d) compute U = u(i,j) = sum(b(k,i) * a(k,j))
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			U[i][j]=0;
			for (int k =0; k<3; k++)
			{
				U[i][j] += b[i][k]*a[j][k];
			}
		}
	}
#ifdef DEBUG
	cout << "det(U) = " << linalg::det(U) << endl;
#endif

}

/************************************************************************
 *	coords::find_kabsch_fit()					*
 * -------------------------------------------------------------------- *
 *	DIFFERENT SELECTION FOR EACH SET				*
 * 	Routine to obtain a best fit rotation of two coordinate sets	*
 * 	about their respective centres of mass: i.e. both sets are	*
 * 	effectively 'translated' to the origin and the best fit		*
 *	rotation is found in this configuration. This is worth bearing	*
 *	in mind when applying the matrix. 				*
 *	The given matrix will rotate *this onto comp			*
 *	NOTE: both *this and comp must be assigned masses before	*
 *		this function is called.				*
 ************************************************************************/
double coords::find_kabsch_fit( coords & comp, int * selection, int n,
		double U[3][3] )
{
	double com_mine[3], com_comp[3];
	double R[3][3], Rt[3][3], RtR[3][3];
	double Eo = 0, xdiff, ydiff, zdiff;
	//assign_masses(); 
	//comp.assign_masses();
	center_of_mass( com_mine, selection, n);
	comp.center_of_mass( com_comp, selection, n);

	// following kabsch, 
	// a) calculate Eo=1/2 sum wn ( xn*xn + yn*yn )

	for ( int i=0; i<n ; i++)
	{
		double xi = X[selection[i]] - com_mine[0]; 
		double yi = Y[selection[i]] - com_mine[1];
		double zi = Z[selection[i]] - com_mine[2];
		double xci = comp.X[selection[i]] - com_comp[0]; 
		double yci = comp.Y[selection[i]] - com_comp[1];
		double zci = comp.Z[selection[i]] - com_comp[2];
		Eo += Weight[i] * ( xi*xi + xci*xci 
				+ yi*yi + yci*yci 
				+ zi*zi + zci*zci ) ;
	}
	Eo /= 0.5;

	// ...determine Rij = (sum wn * xni * ynj)
	for (int i=0; i<3; i++)
	{
		double * I;
		switch(i)
		{
			case 0:
				I=comp.X;
				break;
			case 1:
				I=comp.Y;
				break;
			case 2:
				I=comp.Z;
				break;
		}

		for (int j =0; j<3; j++)
		{
			double * J;
			switch(j)
			{
				case 0:
					J=X;
					break;
				case 1:
					J=Y;
					break;
				case 2:
					J=Z;
					break;
			}

			R[i][j] = 0;
			for (int x=0; x<n; x++)
			{
				int tmpx = selection[x];
				R[i][j] += Weight[tmpx]
					* (I[tmpx] - com_comp[i]) 
					* (J[tmpx] - com_mine[j]);
			}

		}
	}
	
	kabsch_helper( R, U );
	return 0;
}

/************************************************************************
 *	coords::find_kabsch_fit()					*
 * -------------------------------------------------------------------- *
 *	SAME SELECTION FOR BOTH SETS					*
 * 	Routine to obtain a best fit rotation of two coordinate sets	*
 * 	about their respective centres of mass: i.e. both sets are	*
 * 	effectively 'translated' to the origin and the best fit		*
 *	rotation is found in this configuration. This is worth bearing	*
 *	in mind when applying the matrix. 				*
 *	The given matrix will rotate *this onto comp			*
 *	NOTE: both *this and comp must be assigned massed before	*
 *		this function is called.				*
 ************************************************************************/
double coords::find_kabsch_fit( coords & comp, int * selection, 
		int * selection_comp, int n,
		double U[3][3] )
{
	double com_mine[3], com_comp[3];
	double R[3][3], Rt[3][3], RtR[3][3];
	double Eo = 0, xdiff, ydiff, zdiff;
	//assign_masses(); comp.assign_masses();
	center_of_mass( com_mine, selection, n);
	comp.center_of_mass( com_comp, selection_comp, n);

	// following kabsch, 
	// a) calculate Eo=1/2 sum wn ( xn*xn + yn*yn )

	for ( int i=0; i<n ; i++) {
		double xi = X[selection[i]] - com_mine[0]; 
		double yi = Y[selection[i]] - com_mine[1];
		double zi = Z[selection[i]] - com_mine[2];
		double xci = comp.X[selection_comp[i]] - com_comp[0]; 
		double yci = comp.Y[selection_comp[i]] - com_comp[1];
		double zci = comp.Z[selection_comp[i]] - com_comp[2];
		Eo += Weight[i] * ( xi*xi + xci*xci 
				+ yi*yi + yci*yci 
				+ zi*zi + zci*zci ) ;
	}
	Eo /= 0.5;

	// ...determine Rij = (sum wn * xni * ynj)
	for (int i=0; i<3; i++)
	{
		double * I;
		switch(i)
		{
			case 0:
				I=comp.X;
				break;
			case 1:
				I=comp.Y;
				break;
			case 2:
				I=comp.Z;
				break;
		}

		for (int j =0; j<3; j++)
		{
			double * J;
			switch(j)
			{
				case 0:
					J=X;
					break;
				case 1:
					J=Y;
					break;
				case 2:
					J=Z;
					break;
			}

			R[i][j] = 0;
			for (int x=0; x<n; x++)
			{
				int tmpx = selection[x];
				int tmpcx = selection_comp[x];
				R[i][j] += Weight[tmpx]
					* (I[tmpcx] - com_comp[i]) 
					* (J[tmpx] - com_mine[j]);
			}

		}
	}
	
#ifdef DEBUG
	cout << "R matrix is: " << endl;
	print_matrix(R);
#endif
	kabsch_helper( R, U );
	return 0;
}

/************************************************************************
 *	coords::Get_Volume()					*
 * -------------------------------------------------------------------- *
 * 	Returns volume of currently defined periodic box		*
 ************************************************************************/
double coords::Get_Volume() const
{
	if (Bound == CUBIC)
		return boxlength*boxlength*boxlength;
	else if (Bound == OCTAHEDRAL)
		return boxlength*boxlength*boxlength/2.0;
	else
	{
		print_err( "coords::Get_Volume()",
			      "Cannot calc volume without boundary cond!");
#ifndef SWIG
		throw BoundaryErr();
#else
		exit(1);
#endif
	}
}

void coords::Get_Subset( const coords & source )
{
	if ( source.Get_N() < N ) {
		print_err("coords::Get_Subset( const coords & source )",
				"Source has too few atoms");
#ifndef SWIG
		throw WrongNumAtomsErr();
#else
		exit(1);
#endif
	}

	for (int x=0; x<N; x++) {
		X[x] = source.X[x];
		Y[x] = source.Y[x];
		Z[x] = source.Z[x];
	}
}

int coords::Set_Atomtype(int i, string s)
{
	if (s.size() > 4) {
		print_err("int coords::Set_Atomtype(int i, string s);",
				"Word for atom type longer than 4 chars");
#ifndef SWIG
		throw StringLenErr();
#else
		exit(1);
#endif
	} else {
		type[i] = s;
		Type[i].Set(s.c_str());
	}
	return 0;
}

/*
   coords::Get_Selection()					
   -------------------------------------------------------------------- 
   Get coordinates for a selection	(float version)	
 */
void coords::Get_Selection( int * selection, int nselect, float * Xsel, 
		float * Ysel, float * Zsel) const
{
	for (int i=0; i<nselect; i++) {
		Xsel[i] = X[selection[i]];
		Ysel[i] = Y[selection[i]];
		Zsel[i] = Z[selection[i]];
	}
}

/*
   coords::Get_Selection()					
   -------------------------------------------------------------------- 
   Get coordinates for a selection	(double version)	
 */
void coords::Get_Selection( int * selection, int nselect, double * Xsel, 
		double * Ysel, double * Zsel) const
{
	for (int i=0; i<nselect; i++)
	{
		Xsel[i] = X[selection[i]];
		Ysel[i] = Y[selection[i]];
		Zsel[i] = Z[selection[i]];
	}
}

/*
   coords::Get_Selection()					
   -------------------------------------------------------------------- 
   Get coordinates for a selection	(double version)	
 */
void coords::Get_Selection( int * selection, int nselect, 
		vector<double> & XYZ) const
{
	for (int i=0; i<nselect; i++)
	{
		XYZ[i] = X[selection[i]];
		XYZ[nselect+i] = Y[selection[i]];
		XYZ[2*nselect+i] = Z[selection[i]];
	}
}

/*
   coords::Set_Selection()					
   -------------------------------------------------------------------- 
   Set coordinates for a selection	(double version)	
 */
void coords::Set_Selection( int * selection, int nselect, 
		vector<double> & XYZ) 
{
	for (int i=0; i<nselect; i++)
	{
		X[selection[i]] = XYZ[i];
		Y[selection[i]] = XYZ[nselect+i];
		Z[selection[i]] = XYZ[2*nselect+i];
	}
}


void coords::add( int i, double x, double y, double z )
{
	X[i] += x; Y[i] += y; Z[i] += z;
}

/*
   coords::cat ( const coords & set 2 )
   -------------------------------------------------------------------- 
   Adds set2's coordinates onto the end of *this
 */
void coords::cat ( const coords & set2)
{
	int old_N = N;
	N+=set2.N;
	if (Bound != set2.Bound) {
		print_err("void cat ( const coords & set2)",
				"Boundary types are different");
#ifndef SWIG
		throw BoundaryErr();
#else
		exit(1);
#endif
	}
	// for now, we simply use our own boxlength and ignore set2's
	int * nu_Atomno = new int[N];
	int * nu_Resno = new int[N];
	char (*nu_Resname)[5] = new char[N][CHARMM_STR];
	AtomType * nu_Type = new AtomType[N];
	string * nu_type = new string[N];
	double * nu_X = new double[N];
	double * nu_Y = new double[N];
	double * nu_Z = new double[N];
	char (*nu_Segid)[5] = new char[N][CHARMM_STR];
	char (*nu_Resid)[5] = new char[N][CHARMM_STR];
	double * nu_Weight = new double[N];

	for (int atom = 0; atom< old_N; atom++) {
		nu_Atomno[atom] = Atomno[atom];
		nu_Resno[atom] = Resno[atom];
		strcpy(nu_Resname[atom], Resname[atom]);
		nu_Type[atom]= Type[atom];		
		nu_type[atom] = type[atom];
		nu_X[atom] = X[atom];
		nu_Y[atom] = Y[atom];
		nu_Z[atom] = Z[atom];
		strcpy(nu_Segid[atom], Segid[atom]);
		strcpy(nu_Resid[atom], Resid[atom]);
		nu_Weight[atom] = Weight[atom];
	}
	int final_old_resno = Resno[old_N - 1];
	//cout << "copied old crds" << endl;
	for (int atom = 0; atom< set2.N; atom++) {
	//	cout << "copying new atom # " << atom << endl;
		int nu_atom = atom + old_N;
		nu_Atomno[nu_atom] = nu_atom + 1;
		nu_Resno[nu_atom] = set2.Resno[atom] + final_old_resno;
		strcpy(nu_Resname[nu_atom], set2.Resname[atom]);
		nu_Type[nu_atom]=set2.Type[atom];		// FIXME
		nu_type[nu_atom] = set2.type[atom];
		nu_X[nu_atom] = set2.X[atom];
		nu_Y[nu_atom] = set2.Y[atom];
		nu_Z[nu_atom] = set2.Z[atom];
		strcpy(nu_Segid[nu_atom], set2.Segid[atom]);
		strcpy(nu_Resid[nu_atom], set2.Resid[atom]);
		nu_Weight[nu_atom] = set2.Weight[atom];
	}
	
	//cout << "copied new crds" << endl;
	delete [] Atomno;
	delete [] Resno;
	delete [] Resname;
	delete [] Type;
	delete [] type;
	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] Segid;
	delete [] Resid;
	delete [] Weight;
	
	Atomno = &nu_Atomno[0];
	Resno = &nu_Resno[0];
	Resname = &nu_Resname[0];
	Type = &nu_Type[0];
	type = &nu_type[0];
	X = &nu_X[0];
	Y = &nu_Y[0];
	Z = &nu_Z[0];
	Segid = &nu_Segid[0];
	Resid = &nu_Resid[0];
	Weight = &nu_Weight[0];

}

/*
   coords::count_residues()					
   -------------------------------------------------------------------- 
   count residues in a segment!
 */
int coords::count_residues( const string & segid ) const {
	int atom = 0; int nres = 0;
	while ( strcmp( segid.c_str(), Segid[atom]) ) atom++;

	int curr_resno = -1; //invalid
	while ( ! strcmp( segid.c_str(), Segid[atom]) ) {
		if ( curr_resno != Resno[atom] ) {
			nres++;
			curr_resno = Resno[atom];
		}
		atom++;
	}

	return nres;
}

void coords::set_defaults( const coords & new_def )
{
	default_crd = new_def;
}

PDBOutput::PDBOutput(const string & pdbname, char *title)
{
	outfilename = pdbname;
	outfile = fopen(outfilename.c_str(),"w");
	if (outfile == NULL) {
		print_err( "PDBOutput::PDBOutput",
				"could not open pdb output file");
#ifndef SWIG
		throw FileErr();
#else
		exit(1);
#endif
	}
	//if (!outfile.good()) {
	//	cerr << "Could not open pdb output file " << outfilename
	//		<< endl;
	//	exit(1);
	//}
	natom=model=mdl_natom=0;
	//outfile << "REMARK  " << title << endl;
	fprintf(outfile, "REMARK  %s\n", title);
}

void PDBOutput::append_chain(const coords & coords, char chain)
{
	int oldnatom = natom;
	natom += coords.N;
	const char * format =
		"%-6s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n";
	//	"%2*%5i%*%4s%*%3s%2*%4i%4*%8.3f%8.3f%8.3f%6.2f%6.2f%6*%4s";
	for (int i=0; i<coords.N; i++) {
		//AtomType t = coords.AtomType[i];
		char tmpres[4];
		strncpy(tmpres, coords.Resname[i], 3);
		tmpres[3] = '\0';
		fprintf(outfile, format, "ATOM", oldnatom+i+1, 
				coords.type[i].c_str(), tmpres,
				chain, atoi(coords.Resid[i]), coords.X[i],
				coords.Y[i], coords.Z[i], 1.0, 
				coords.Weight[i], coords.Segid[i]);
	}
}

void PDBOutput::append_model(const coords & coords)
{
	int oldnatom = natom;
	natom += coords.N;
	fprintf(outfile, "MODEL     %4i\n", ++model);
	const char * format =
		"%-6s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n";
	//	"%2*%5i%*%4s%*%3s%2*%4i%4*%8.3f%8.3f%8.3f%6.2f%6.2f%6*%4s";
	for (int i=0; i<coords.N; i++) {
		//AtomType t = coords.AtomType[i];
		char tmpres[4];
		strncpy(tmpres, coords.Resname[i], 3);
		tmpres[3] = '\0';
		fprintf(outfile, format, "ATOM", i+1, 
				coords.type[i].c_str(), tmpres,
				' ', atoi(coords.Resid[i]), coords.X[i],
				coords.Y[i], coords.Z[i], 1.0, 
				coords.Weight[i], coords.Segid[i]);
	}
	fprintf(outfile, "ENDMDL\n");
}

PDBOutput::~PDBOutput()
{
	close();
}

void PDBOutput::close()
{
	if (outfile != NULL) {
		fprintf(outfile, "END\n");
		fclose(outfile);
		outfile = NULL;
	}
}

/*
void coords::add( const coords & Rhs )
{
	if (N != Rhs.N) {
		print_err( "coords::coords()", 
				"Coords being added of different size!!!");
		throw WrongNumAtomsErr();
	}

	for (int atom =0; atom<N; atom++) {
		X[atom] += Rhs.X[atom];
		Y[atom] += Rhs.Y[atom];
		Z[atom] += Rhs.Z[atom];
	}
}

void coords::multiply( double cnst )
{

	for (int atom =0; atom<N; atom++) {
		X[atom] *= cnst;
		Y[atom] *= cnst;
		Z[atom] *= cnst;
	}
}

*/


