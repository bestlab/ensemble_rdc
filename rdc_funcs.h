
#ifndef _RDC_FUNCS_H
#define _RDC_FUNCS_H

#include <vector>
#include <string>
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "coords.h"

using namespace std;

namespace rdc_const {
	const double mu0 = 4*M_PI*1e-7;
	const double hcross = 1.05e-34;
}

struct rdc {
	int resi; //residue
	int resj;
	int indi; //atom index
	int indj;
	string atomi; //atom name
	string atomj;
	double Dij;
	double Dmax; // max rdc
};

void parse_rdc(const string &rdc_file, const coords *crd, vector<rdc> &rdc_dat, bool s2=false);

//double calc_Q(vector<double> &calc, vector<double> &exp);
double calc_Q(gsl_vector *calc, gsl_vector *expt);

void read_weights(const string &weight_file, vector<double> &weights);

double gyro(char atom);

double bondlen(const string &a, const string &b);

void fill_rdcvec(const vector<rdc> &R, gsl_vector *v);

void unnorm_rdcvec(const vector<rdc> &R, gsl_vector *v, gsl_vector *w);

void calc_coef_pdb(coords &a, vector<rdc> &rdc_data, gsl_matrix *coef_mat);

void calc_coef_trj(coords &a, const vector<string> &trjfiles,
		vector<rdc> &rdc_data, gsl_matrix *coef_mat, int *sel, int nsel,
		vector<double> &weights);

void calc_coef_single(coords &a, coords &b, vector<rdc> &rdc_data, gsl_matrix *coef_mat, int *sel, 
		int nsel, double w=1.0);

void why_doesnt_gsl_have_matrix_multiplication(gsl_matrix *A, gsl_vector *x, gsl_vector *y);

void get_align_res(vector<int> &align_res, const char *s);

void select_types( const coords & coor_main, int nselect, 
		int * lookup, AtomType * types, int ntype, vector<int> &align_res);

#endif // _RDC_FUNCS_H

