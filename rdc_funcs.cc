
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctype.h>
#include "rdc_funcs.h"
#include "coords.h"
#include "traj.h"

double gyro(char atom)
{
	switch (atom) {
		case 'H':
			return 2.6752e8;
			break;
		case 'C':
			return 6.7287e7;
			break;
		case 'N':
			return 2.712e7;
			break;
		default:
			fprintf(stderr,"Unknown atom type in gyro %c\n",atom);
			exit(1);
	}
	return 0.0;
}

/*
double calc_Q_norm(gsl_vector *calc, gsl_vector *expt)
{
	int ndat;
	double sdev2,scalc2,tmpf;
	sdev2 = 0;
	scalc2 = 0;
	ndat = calc.size();
	for (int i=0; i<ndat; i++) {
		tmpf = calc[i]-exp[i];
		sdev2 += tmpf*tmpf;
		scalc2 += calc[i]*calc[i];
	}
	return sqrt(sdev2/scalc2);
}
*/

double calc_Q(gsl_vector *calc, gsl_vector *expt)
{
	int ndat,nexpt;
	double sdev2,scalc2,tmpf,expt_i,calc_i;
	sdev2 = 0;
	scalc2 = 0;
	ndat = calc->size;
	nexpt = expt->size;
	if (ndat != nexpt) {
		fprintf(stderr,"calc_Q: array size mismatch, %i:%i\n",ndat,nexpt);
		exit(1);
	}
	for (int i=0; i<ndat; i++) {
		calc_i = gsl_vector_get(calc,i);
		expt_i = gsl_vector_get(expt,i);
		//fprintf(stdout,"calc_i=%12.6e ; expt_i=%12.6e\n",calc_i,expt_i);
		//fprintf(stdout,"%12.5e %12.5e\n",calc_i,expt_i);
		tmpf = calc_i-expt_i;
		sdev2 += tmpf*tmpf;
		scalc2 += expt_i*expt_i;
	}
	//fprintf(stdout,"sqdev=%12.6e ; scalc2=%12.6e\n",sdev2,scalc2);
	fflush(stdout);
	return sqrt(sdev2/scalc2);
}

void strip(const string &s, string &st)
{
	st = "";
	for (int i=0;i<s.size(); i++) {
		if (isalnum(s[i]) || s[i]=='\'') {
			st += s[i];
		}
	}
}

double bondlen(const string &a, const string &b)
{
	char ati,atj;
	string at1, at2;
	string HN = "HN";
	string H = "H";
	string N = "N";
	string C = "C";
	string HA = "HA";
	string CA = "CA";
	string CB = "CB";
	if (a == H)
		at1 = HN;
	else
		at1 = a;
	if (b == H)
		at2 = HN;
	else
		at2 = b;
	//
	if ((at1==N && at2==HN) || (at1==HN && at2==N)) {
		return 1.04e-10;
	} else if ((at1==HN && at2==C) || (at1==C && at2==HN)) {
		return 2.04e-10;
	} else if ((at1==CA && at2==HA) || (at1==HA && at2==CA)) {
		return 1.12e-10;
	} else if ((at1==N && at2==C) || (at1==C && at2==N)) {
		return 1.33e-10;
	} else if ((at1==CA && at2==C) || (at1==C && at2==CA)) {
		return 1.53e-10;
	} else if ((at1==CA && at2==CB) || (at1==CB && at2==CA)) {
		return 1.53e-10;
	} else {
		ati = at1[0]; atj = at2[0];
		if (isdigit(ati)) 
			ati = at1[1];
		if (isdigit(atj)) 
			atj = at2[1];
		if ((ati=='C' && atj=='C')) {
			return 1.53e-10;
		} else if ((ati=='C' && atj=='H') || (ati=='H' && atj=='C')) {
			return 1.12e-10;
		} else if ((ati=='N' && atj=='H') || (ati=='H' && atj=='N')) {
			return 1.04e-10;
		} else if (ati=='H' && atj=='H') {
			return 1.633e-10; // need to check this value
		} else {
			fprintf(stderr, "Could not find bond length for atoms: %s - %s", 
					at1.c_str(), at2.c_str());
		}
	}
}


void parse_rdc(const string &rdc_file, const coords *crd, vector<rdc> &rdc_dat, bool s2)
{
	//double mu0 = 4*M_PI*1e-7;
	//double hcross = 1.05e-34;

	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
	int nrdc;
	int natom;
	bool foundi,foundj;
	char ati, atj;
	char *junk;
	string atom;
	FILE *f;
	//
	f = fopen(rdc_file.c_str(),"r");
	if (f == NULL) {
		fprintf(stderr,"Could not open file %s\n",rdc_file.c_str());
		exit(1);
	}
	nrdc = 0;
	while (feof(f) == 0) {
		nrdc++;
		junk = fgets(buf,buf_len,f);
	}
	nrdc --;
	rdc_dat.resize(nrdc);
	natom = crd->Get_N();
	f = fopen(rdc_file.c_str(),"r");
	for (int p=0; p<nrdc; p++) {
		junk = fgets(buf,buf_len,f);
		rdc_dat[p].resi = atoi(strtok(buf," \t"));
		rdc_dat[p].atomi = strtok(NULL," \t");
		rdc_dat[p].resj = atoi(strtok(NULL," \t"));
		rdc_dat[p].atomj = strtok(NULL," \t");
		if (!s2) {
			rdc_dat[p].Dij = atof(strtok(NULL," \t"));
			double blen = bondlen(rdc_dat[p].atomi,rdc_dat[p].atomj);
			ati  = rdc_dat[p].atomi[0];
			if ( isdigit (ati) )
				ati  = rdc_dat[p].atomi[1];
			atj  = rdc_dat[p].atomj[0];
			if ( isdigit (atj) )
				atj  = rdc_dat[p].atomj[1];
			rdc_dat[p].Dmax = (rdc_const::mu0*rdc_const::hcross)
				*gyro(ati)*gyro(atj)
				/ ( 4. *M_PI*M_PI*blen*blen*blen );
//				*gyro(rdc_dat[p].atomi[0])*gyro(rdc_dat[p].atomj[0])
		}
		foundi = false;
		foundj = false;
		for (int a=0;a<natom;a++) {
			strip(crd->Get_Atomtype(a),atom);
			if (atom == string("H"))
				atom = "HN";
			//atom.replace(atom.find(" "),1,"");
			int res = crd->Get_Resno(a);
			//fprintf(stderr,"[%s] %i\n",atom.c_str(),res);
			if (rdc_dat[p].atomi==atom && rdc_dat[p].resi==res) {
				rdc_dat[p].indi = a;
				foundi = true;
			} else if (rdc_dat[p].atomj==atom && rdc_dat[p].resj==res) {
				rdc_dat[p].indj = a;
				foundj = true;
			}
			if (foundi && foundj)
				break;
		}
		if ( !foundi ) {
			fprintf(stderr,"Couldn't find %i [%s]\n",rdc_dat[p].resi,
					rdc_dat[p].atomi.c_str());
			exit(1);
		} else if ( !foundj ) {
			fprintf(stderr,"Couldn't find %i [%s]\n",rdc_dat[p].resj,
					rdc_dat[p].atomj.c_str());
			exit(1);
		}
	}
	fclose(f);
}

void calc_coefs_trj(coords &a, const vector<string> &trjfiles,
		vector<rdc> &rdc_data, vector<gsl_matrix *> coef_mats, int *sel, int nsel)
{
	BaseITrajFile *trj;
	// these will be reference coordinates
	coords b(a);
	int ntrj,natom,frame,nframe,nrdc;
	double com[3],wsum;
	//
	b.set_uniform_weights();
	a.set_uniform_weights();
	if ( nsel != 0 ) {
		// centre reference
		b.center_of_mass(com, sel, nsel);
		//fprintf(stdout,"com: %12.5f %12.5f %12.5f\n",com[0],com[1],com[2]);
		b.translate( -com[0], -com[1], -com[2] );
	}
	//
	ntrj = trjfiles.size();
	nrdc = rdc_data.size();
	//
	// first count total number of frames ...
	frame = 0;
	for (int i=0; i<ntrj; i++) {
		int slen = trjfiles[i].size();
		if (slen-trjfiles[i].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			trj = new DCDITrajFile(trjfiles[i].c_str());
		} else if (slen-trjfiles[i].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			trj = new XTCITrajFile(trjfiles[i].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[i].c_str());
			exit(1);
		}
		natom = trj->num_atoms();
		fprintf(stdout,"Number of atoms = %i\n",natom);
		while (trj->frames_left()) {
			*trj >> a;
			frame++;
		}
		delete trj;
	}
	nframe = frame;
	// allocate a lot of  memory
	coef_mats.resize(nframe);
	for (int i = 0; i<nframe; i++) {
		coef_mats[i] = gsl_matrix_alloc(nrdc,5);
	}
	frame = 0;
	for (int i=0; i<ntrj; i++) {
		int slen = trjfiles[i].size();
		if (slen-trjfiles[i].rfind(string(".dcd"))==4) {
			trj = new DCDITrajFile(trjfiles[i].c_str());
		} else if (slen-trjfiles[i].rfind(string(".xtc"))==4) {
			trj = new XTCITrajFile(trjfiles[i].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[i].c_str());
			exit(1);
		}
		natom = trj->num_atoms();
		/*
		X = new float[natom];
		Y = new float[natom];
		Z = new float[natom];
		*/
		while (trj->frames_left()) {
			*trj >> a;
			calc_coef_single(a, b, rdc_data, coef_mats[frame], sel, nsel, 1.);
			frame++;
		}
		delete trj;
	}
}

void calc_coef_trj(coords &a, const vector<string> &trjfiles,
		vector<rdc> &rdc_data, gsl_matrix *coef_mat, int *sel, int nsel,
		vector<double> &weights)
{
	BaseITrajFile *trj;
	// these will be reference coordinates
	coords b(a);
	int ntrj,natom;
	double com[3];
	//
	fprintf(stderr,"calc_coefs_trj 0\n");
	fflush(stderr);
	b.set_uniform_weights();
	a.set_uniform_weights();
	fprintf(stderr,"calc_coefs_trj 1\n");
	fflush(stderr);
	if ( nsel != 0 ) {
		// centre reference
		b.center_of_mass(com, sel, nsel);
		//fprintf(stdout,"com: %12.5f %12.5f %12.5f\n",com[0],com[1],com[2]);
		b.translate( -com[0], -com[1], -com[2] );
	}
	fprintf(stderr,"calc_coefs_trj 2\n");
	fflush(stderr);
	//
	ntrj = trjfiles.size();
	//
	int frame = 0;
	double wsum = 0.0;
	for (int i=0; i<ntrj; i++) {
		int slen = trjfiles[i].size();
		if (slen-trjfiles[i].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			trj = new DCDITrajFile(trjfiles[i].c_str());
		} else if (slen-trjfiles[i].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			trj = new XTCITrajFile(trjfiles[i].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[i].c_str());
			exit(1);
		}
	fprintf(stderr,"opened trj %i\n",i);
	fflush(stderr);
		natom = trj->num_atoms();
		fprintf(stdout,"Number of atoms = %i\n",natom);
		fflush(stderr);
		/*
		X = new float[natom];
		Y = new float[natom];
		Z = new float[natom];
		*/
		while (trj->frames_left()) {
			(*trj) >> a;
	//fprintf(stderr,"trj %i\n",i);
	//fflush(stderr);
			/*
			fprintf(stdout,"X[1] = %8.3f\n",a.Get_X(1));
			fprintf(stdout,"Y[1] = %8.3f\n",a.Get_Y(1));
			fprintf(stdout,"Z[1] = %8.3f\n",a.Get_Z(1));
			*/
			//trj->read_frame(X,Y,Z,natom);
			//fprintf (stdout,"frame %5i: X[0] = %8.3f\n",frame++,X[0]);
			double w;
			if (weights.size() != 0) {
				w = weights[frame];
			} else {
				w = 1.0;
			}
			//fprintf(stdout,"%i %8.3f\n",frame,w);
			//calc_coef_single(a, b, rdc_data, coef_mat, NULL, 0, w);
			calc_coef_single(a, b, rdc_data, coef_mat, sel, nsel, w);
			wsum += w;
			frame++;
		}
		delete trj;
	}
	gsl_matrix_scale(coef_mat, 1./wsum);
}


void calc_coef_pdb(coords &a, vector<rdc> &rdc_data, gsl_matrix *coef_mat)
{
	coords b(a);
	calc_coef_single(a, b, rdc_data, coef_mat, NULL, 0);
//	fprintf(stdout,"coef_mat[0][0],%12.4e\n",gsl_matrix_get(coef_mat,0,0));
}

void calc_coef_single(coords &a, coords &b, vector<rdc> &rdc_data, gsl_matrix *coef_mat, int *sel, 
		int nsel, double w)
{
	int nrdc, i, j;
	double com[3], mat[3][3];
	nrdc = rdc_data.size();
	double mu_x, mu_y, mu_z,mu_r;
	//
	if ( nsel != 0 ) {
		// do alignment
		a.center_of_mass(com, sel, nsel);
		a.translate( -com[0], -com[1], -com[2] );
		//fprintf(stdout,"com: %12.5f %12.5f %12.5f\n",com[0],com[1],com[2]);
		a.find_kabsch_fit( b, sel, nsel, mat);
		a.rotate( mat );
	}
	//
	for (int d=0; d<nrdc; d++) {
		i = rdc_data[d].indi;
		j = rdc_data[d].indj;
		mu_x = a.Get_X(i)-a.Get_X(j);
		mu_y = a.Get_Y(i)-a.Get_Y(j);
		mu_z = a.Get_Z(i)-a.Get_Z(j);
		mu_r = sqrt(mu_x*mu_x+mu_y*mu_y+mu_z*mu_z);
		mu_x /= mu_r;
		mu_y /= mu_r;
		mu_z /= mu_r;
		gsl_matrix_set(coef_mat,d,0, 
				gsl_matrix_get(coef_mat,d,0)+w*(mu_x*mu_x-mu_z*mu_z));
		gsl_matrix_set(coef_mat,d,1, 
				gsl_matrix_get(coef_mat,d,1)+w*(mu_y*mu_y-mu_z*mu_z));
		gsl_matrix_set(coef_mat,d,2, 
				gsl_matrix_get(coef_mat,d,2)+w*(2.0*mu_x*mu_y));
		gsl_matrix_set(coef_mat,d,3, 
				gsl_matrix_get(coef_mat,d,3)+w*(2.0*mu_x*mu_z));
		gsl_matrix_set(coef_mat,d,4, 
				gsl_matrix_get(coef_mat,d,4)+w*(2.0*mu_y*mu_z));
		/*
		coef_mat[d][0] += mu_x*mu_x-mu_z*mu_z;
		coef_mat[d][1] += mu_y*mu_y-mu_z*mu_z;
		coef_mat[d][2] += 2.0*mu_x*mu_y;
		coef_mat[d][3] += 2.0*mu_x*mu_z;
		coef_mat[d][4] += 2.0*mu_y*mu_z;
		*/
	}
}

void fill_rdcvec(const vector<rdc> &R, gsl_vector *v)
{
	int nrdc = R.size();
	for (int i=0; i<nrdc; i++) {
		gsl_vector_set(v,i,R[i].Dij/R[i].Dmax);
	}
}

void unnorm_rdcvec(const vector<rdc> &R, gsl_vector *v, gsl_vector *w)
{
	int nrdc = R.size();
	double tmp;
	for (int i=0; i<nrdc; i++) {
		tmp = gsl_vector_get(v,i);
		gsl_vector_set(w,i,tmp*R[i].Dmax);
	}
}

void why_doesnt_gsl_have_matrix_multiplication(gsl_matrix *A, gsl_vector *x, gsl_vector *y)
{
	// y = A*x
	int A_m, A_n, x_n, y_m;
	double tmpf;
	//
	A_m = A->size1;
	A_n = A->size2;
	x_n = x->size;
	y_m = y->size;
	if ( A_n != x_n ) {
		fprintf(stderr,"matrix and vector of different dimensions!\n");
		exit(1);
	}
	if ( A_m != y_m ) {
		fprintf(stderr,"output vector of wrong dimension!\n");
		exit(1);
	}
	gsl_vector_set_zero(y);
	for (int i=0; i<A_m; i++) {
		tmpf = 0.0;
		for (int j=0; j<A_n; j++) {
			tmpf += gsl_matrix_get(A,i,j) * gsl_vector_get(x,j);
		}
		gsl_vector_set(y,i,tmpf);
	}
}

void read_weights(const string &weight_file, vector<double> &weights)
{
	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
	int n;
	char *junk;
	FILE *f;
	//
	f = fopen(weight_file.c_str(),"r");
	if (f == NULL) {
		fprintf(stderr,"Could not open file %s\n",weight_file.c_str());
		exit(1);
	}
	n = 0;
	while (feof(f) == 0) {
		n++;
		junk = fgets(buf,buf_len,f);
	}
	n--;
	fclose(f);
	f = fopen(weight_file.c_str(),"r");
	weights.resize(n);
	for (int p=0; p<n; p++) {
		junk = fgets(buf,buf_len,f);
		weights[p] = atof(strtok(buf," \t"));
		fprintf(stdout,"W[%i] = %8.3f\n",p,weights[p]);
	}
}

void get_align_res(vector<int> &align_res, const char *s)
{
	//fprintf(stdout,"GOT HERE 0\n");
	string ss = s;
	//fprintf(stdout,"GOT HERE 1\n");
	int i, j, state;
	int end = ss.size();
	int lo, hi;
	const char *sep1 = "-";
	const char *sep2 = ":";
	i=0;
	state = 0;
	align_res.resize(0);
	while (i>=0 && i<end) {
		if (state==0) {
			j = ss.find(sep1,i);
			//fprintf(stdout,"i,j : %i, %i\n",i,j);
			lo = atoi(ss.substr(i,j).c_str());
			state = 1;
			i = j+1;
		} else if (state == 1) {
			j = ss.find(sep2,i);
			if (j<0) {
				j=end;
			}
			//fprintf(stdout,"i,j : %i, %i\n",i,j);
			hi = atoi(ss.substr(i,j).c_str());
			int newres = hi-lo+1; //+1 for inclusive
			int cursize = align_res.size();
			align_res.resize(cursize+newres);
			for (int p=0; p<newres; p++) {
				align_res[cursize+p] = lo+p;
			}
			//fprintf(stdout,"lo -- hi : %i -- %i\n",lo,hi);
			i = j+1;
			state = 0;
		}
	}
	for (int i=0; i<align_res.size(); i++) 
		fprintf(stdout,"%i,",align_res[i]);
	fprintf(stdout,"\n");
}

void select_types( const coords & coor_main, int nselect, 
		int * lookup, AtomType * types, int ntype, vector<int> &align_res)
{
	int N = coor_main.Get_N();
	int l = 0;
	int res;

	for (int i=0; i < N; i++) {
		AtomType tmpat;
		coor_main.Get_Atomtype ( i, tmpat );
		res = coor_main.Get_Resno ( i );
		for (int t=0; t<ntype; t++) 
			if (tmpat == types[t]) {
				for (int p=0; p<nselect;p++) {
					if (res == align_res[p]) {
						lookup[l++] = i;
						break;
					}
				}
			}
	}
}

