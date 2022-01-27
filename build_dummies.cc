/*
 * build one dummy atom aligned with C=O
 *
 */

#include <unistd.h> // for getopt
#include <cstdio> 
#include <cstdlib> 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "coords.h"
#include "traj.h"
#include "rdc_funcs.h"

const char *usage = "\n\n"
"        Usage\n"
"             rdc -w w.dat -d rdcfit.dat -b rdcbackcalc.dat -o fit.dat -s xyz.pdb trj1 trj2 ...trjN\n\n"
"        ./rdc -d rdc.dat -o fit.dat -s xyz.pdb\n"
"             ...will fit the single structure in xyz.pdb\n"
"             to rdc.dat and print the alignment tensor (in comment line) and \n"
"             fitted data to output.dat\n"
"        ./rdc -d rdc.dat -o fit.dat -s xyz.pdb trj1 ... trjN\n"
"             ...will fit the ensemble of structures in trajectories trj1...trjN\n"
"             to rdc.dat and print the alignment tensor (in comment line) and \n"
"             fitted data to output.dat\n"
"         ./rdc -d rdc.dat -b bc.dat -o output.dat -s xyz.pdb trj1...trjN\n"
"             ... will do as above, but rdc.dat will be used to fit the alignment \n"
"             tensor and bc.dat will be used for backcalculating rdc's (useful if\n"
"             only a fraction of the rdc's (e.g. backbone rdcs in SS) are to be\n"
"             used to fit the alignment tensor).\n"
"         -w w.dat will weight the data with the weights in w.dat\n"
"         -a 11-20:33-45:56-72 will align frames based on the CA atoms of selected residues\n"
"\n\n";


int main(int argc, char **argv)
{
	int c,ntrj,natom,nrdc_fit,nrdc_bc;
	float *X,*Y,*Z;
	string itrj,rdc_fit,rdc_bc,output,pdbfile,weight_file;
	vector<string> trjfiles;
	vector<rdc> fit_data;
	vector<rdc> bc_data;
	vector<double> weight_vec;
	coords *crd;
	bool align;
	vector<int> align_res;
	gsl_vector *rdc_vec, *rdc_vec_bc, *S, *Stmp, *work, *bc;
	gsl_vector *rdc_vec_bc_unnorm, *bc_unnorm;
	gsl_matrix *coef_mat, *coef_mat_bc, *A, *V;
	double Sxx,Syy,Szz,Sxy,Sxz,Syz;
	double Q_norm, Q_raw;
	double Dmax_NH;
	int *sel, nsel;
	FILE *outp;

	rdc_bc = "NULL";
	rdc_fit = "NULL";
	output = "NULL";
	weight_file = "NULL";
	sel = NULL;
	nsel = 0;
	align = false;
	while (1) {
		c=getopt(argc,argv,"ha:d:b:o:s:w:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage);
				exit(0);
				break;
			case 'd':
				rdc_fit = optarg;
				break;
			case 'b':
				rdc_bc = optarg;
				break;
			case 'o':
				output = optarg;
				break;
			case 's':
				pdbfile = optarg;
				break;
			case 'w':
				weight_file = optarg;
				break;
			case 'a':
				align = true;
				//fprintf(stdout,"align: %s\n",optarg);
				get_align_res(align_res,optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}
	if (rdc_fit == string("NULL")) {
		fprintf(stderr,"Specify rdc_fit file with -d\n\n");
		fprintf(stderr,"%s\n",usage);
		exit(1);
	}
	if (pdbfile == string("NULL")) {
		fprintf(stderr,"Specify pdbfile file with -s\n\n");
		fprintf(stderr,"%s\n",usage);
		exit(1);
	}
	if (output == string("NULL")) {
		fprintf(stderr,"Specify output file with -s\n\n");
		fprintf(stderr,"%s\n",usage);
		exit(1);
	}
	if (rdc_bc == string("NULL")) {
		rdc_bc = rdc_fit;
	}
	ntrj = argc-optind;
	trjfiles.resize(ntrj);
	for (int i=0; i<ntrj; i++) {
		trjfiles[i] = argv[optind+i];
	}
	if (weight_file != string("NULL")) {
		read_weights(weight_file, weight_vec);
	}

	fprintf(stdout,"==============================================================\n");
	fprintf(stdout,"             Using PDB file: %s\n", pdbfile.c_str());
	fprintf(stdout,"Fitting alignment tensor to: %s\n", rdc_fit.c_str());
	fprintf(stdout,"Back-calculating RDC's from: %s\n", rdc_bc.c_str());
	fprintf(stdout,"              RDC output to: %s\n", output.c_str());
	if (ntrj > 0) {
		fprintf(stdout,"Using trajectories:\n");
		for (int i = 0; i<ntrj;i++) 
			fprintf(stdout,"\t%s\n",trjfiles[i].c_str());
	} else {
		fprintf(stdout,"Fitting to single structure in PDB file\n");
	}
	if (align) {
		fprintf(stdout,"Will align all frames to reference pdb\n");
	} else {
		fprintf(stdout,"Will not align frames to reference\n");
	}
	if (weight_file != string("NULL")) {
		fprintf(stdout,"Will weight structures using weights in file %s\n",
				weight_file.c_str());
	}
	fprintf(stdout,"==============================================================\n\n");

	// open for output:
	outp = fopen(output.c_str(),"w");
	crd = new coords(pdbfile.c_str(),"pdb");
	if (align) {
		nsel = align_res.size();
		sel = new int[nsel];

		AtomType tlst[1];
		tlst[0].Set(" CA ");
		select_types(*crd,nsel,sel, tlst, 1,align_res);
		/*
		for (int p=0;p<nsel;p++) {
			fprintf(stdout,"%i\n",sel[p]);
		}
		*/
	}

	parse_rdc(rdc_fit, crd, fit_data);
	nrdc_fit = fit_data.size();
	if (rdc_bc == rdc_fit) {
		parse_rdc(rdc_bc, crd, bc_data);
		nrdc_bc = bc_data.size();
	}
	fprintf(stdout,"%5s %4s %5s -- %5s %4s %5s : %8s (%12s)\n", "res_i", "at_i", "idx_i",
			"res_j", "at_j", "idx_j", "Dij_expt", "Dij_max");
	fprintf(stdout,"==============================================================\n");
	for (int i=0; i<nrdc_fit; i++) {
		fprintf(stdout,"%5i %4s %5i -- %5i %4s %5i : %8.3f (%12.3f)\n",
				fit_data[i].resi,fit_data[i].atomi.c_str(),fit_data[i].indi,
				fit_data[i].resj,fit_data[i].atomj.c_str(),fit_data[i].indj,
				fit_data[i].Dij,fit_data[i].Dmax);
	}
	rdc_vec = gsl_vector_alloc(nrdc_fit);
	rdc_vec_bc_unnorm = gsl_vector_alloc(nrdc_fit);
	bc = gsl_vector_alloc(nrdc_bc);
	bc_unnorm = gsl_vector_alloc(nrdc_fit);
	//tmp = gsl_vector_alloc(nrdc_bc);
	Stmp = gsl_vector_alloc(5);
	S = gsl_vector_alloc(5);
	work = gsl_vector_alloc(5);
	coef_mat = gsl_matrix_alloc(nrdc_fit,5);
	A = gsl_matrix_alloc(nrdc_fit,5);
	V = gsl_matrix_alloc(5,5);
	fill_rdcvec(fit_data, rdc_vec);
	gsl_matrix_set_zero(coef_mat);
	gsl_vector_set_zero(bc);

	if (ntrj==0) {
		// SINGLE STRUCTURE CASE:
		calc_coef_pdb(*crd,fit_data,coef_mat);
		if (rdc_bc == rdc_fit) {
			coef_mat_bc = coef_mat;
			rdc_vec_bc = rdc_vec;
		} else {
			rdc_vec_bc = gsl_vector_alloc(nrdc_fit);
			fill_rdcvec(bc_data, rdc_vec_bc);
			coef_mat_bc = gsl_matrix_alloc(nrdc_fit,5);
			gsl_matrix_set_zero(coef_mat_bc);
			calc_coef_pdb(*crd,bc_data,coef_mat_bc);
		}
	} else {
		// CALCULATE RDC'S FROM TRAJECTORIES
		calc_coef_trj(*crd, trjfiles, fit_data, coef_mat, sel, nsel, weight_vec);
		if (rdc_bc == rdc_fit) {
			coef_mat_bc = coef_mat;
			rdc_vec_bc = rdc_vec;
		} else {
			rdc_vec_bc = gsl_vector_alloc(nrdc_fit);
			fill_rdcvec(bc_data, rdc_vec_bc);
			coef_mat_bc = gsl_matrix_alloc(nrdc_fit,5);
			gsl_matrix_set_zero(coef_mat_bc);
		//	calc_coef_pdb(*crd,bc_data,coef_mat_bc);
		}
	}
	/*
	for (int i=0; i<nrdc_fit; i++) {
		fprintf(stdout,"%5i %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",i,
				gsl_matrix_get(coef_mat,i,0),
				gsl_matrix_get(coef_mat,i,1),
				gsl_matrix_get(coef_mat,i,2),
				gsl_matrix_get(coef_mat,i,3),
				gsl_matrix_get(coef_mat,i,4),
				gsl_vector_get(rdc_vec,i));
	}
	*/
	// let gsl do the work
	// ...first make a copy of coef_mat because it is modified by SVD functions
	gsl_matrix_memcpy(A,coef_mat);
	gsl_linalg_SV_decomp(A, V, Stmp, work);
	gsl_linalg_SV_solve(A, V, Stmp, rdc_vec, S);
	Sxx = gsl_vector_get(S,0);
	Syy = gsl_vector_get(S,1);
	Szz = -Sxx-Syy;
	Sxy = gsl_vector_get(S,2);
	Sxz = gsl_vector_get(S,3);
	Syz = gsl_vector_get(S,4);
	fprintf(outp,"# Alignment tensor\n");
	fprintf(outp,"# Sxx = %12.5e\n",Sxx);
	fprintf(outp,"# Syy = %12.5e\n",Syy);
	fprintf(outp,"# Szz = %12.5e\n",Szz);
	fprintf(outp,"# Sxy = %12.5e\n",Sxy);
	fprintf(outp,"# Sxz = %12.5e\n",Sxz);
	fprintf(outp,"# Syz = %12.5e\n",Syz);
	gsl_blas_dgemv (CblasNoTrans, 1.0, coef_mat_bc, S, 0., bc);
	//why_doesnt_gsl_have_matrix_multiplication(coef_mat_bc, S, bc);
	Q_norm = calc_Q(bc, rdc_vec_bc);
	fprintf(outp,"# Normalized Q = %12.6f\n",Q_norm);
	fprintf(stdout,"Normalized Q = %12.6f\n",Q_norm);
	unnorm_rdcvec(bc_data,bc,bc_unnorm);
	unnorm_rdcvec(bc_data,rdc_vec_bc,rdc_vec_bc_unnorm);
	Q_raw = calc_Q(bc_unnorm, rdc_vec_bc_unnorm);
	fprintf(outp,"# Raw Q = %12.6f\n",Q_raw);
	fprintf(stdout,"Raw Q = %12.6f\n",Q_raw);

	Dmax_NH = rdc_const::mu0*rdc_const::hcross*
		(gyro('N')*gyro('H'))/(4.*(M_PI*M_PI)*pow(bondlen(string("N"),string("H")),3.));

	fprintf(outp,"#%4s %4s %5s %4s %8s %8s %8s %8s\n",
		"res1","at1","res2","at2","D_calc","D_exp",
			"|D_calc|","|D_exp|");

	for (int i=0; i<nrdc_bc; i++) {
		fprintf(outp,"%5i %4s %5i %4s %8.3f %8.3f %8.3f %8.3f\n",
				bc_data[i].resi,bc_data[i].atomi.c_str(),
				bc_data[i].resj,bc_data[i].atomj.c_str(),
				gsl_vector_get(bc_unnorm,i),
				gsl_vector_get(rdc_vec_bc_unnorm,i),
				gsl_vector_get(bc,i)*Dmax_NH,
				gsl_vector_get(rdc_vec_bc,i)*Dmax_NH);
	}

	gsl_vector_free(rdc_vec);
	gsl_matrix_free(coef_mat);
	return 0;
}
