/*
 *
 * calculate order parameters from MD trajectory
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
"             s2 [-w w.dat] [-a selection] [-n nblock] -d pairs.dat -o s2out.dat -s xyz.pdb trj1 trj2 ...trjN\n\n"
"        e.g."
"        ./s2 -o s2out.dat -s xyz.pdb trj1 ... trjN\n"
"             will calculate order parameters from trajectories trj1...trjN\n"
"             for the atom pairs listed in pairs.dat in format: resi atomi resj atomj, e.g.\n"
"             3 N 3 HN\n"
"             after aligning structures with CA atoms in xyz.pdb\n"
"             writing output to s2out.dat\n"
"         -w w.dat will weight the data with the weights in w.dat\n"
"         -a 11-20:33-45:56-72 will align frames based on the CA atoms of selected residues\n"
"         -n nblock will divide the data into nblock blocks for error analysis (default: 10)\n"
"\n\n";


double calc_S2(const double &Sxx, const double &Syy,const double &Szz,
		const double &Sxy, const double &Sxz,const double &Syz,
		int ndat)
{
	double mu_xx, mu_yy, mu_zz, mu_xy, mu_xz, mu_yz, tmp;
	mu_xx = Sxx/float(ndat);
	mu_yy = Syy/float(ndat);
	mu_zz = Szz/float(ndat);
	mu_xy = Sxy/float(ndat);
	mu_xz = Sxz/float(ndat);
	mu_yz = Syz/float(ndat);
	tmp = 1.5*(mu_xx*mu_xx+mu_yy*mu_yy+mu_zz*mu_zz+2.*mu_xy*mu_xy+2.*mu_xz*mu_xz+2.*mu_yz*mu_yz) - 0.5;
	return tmp;
}

int main(int argc, char **argv)
{
	int c,ntrj,natom,nS2;
	float *X,*Y,*Z;
	string itrj,S2pairs,output,pdbfile,weight_file;
	vector<string> trjfiles;
	vector<rdc> pairs;
	vector<double> Sxx,Syy,Szz,Sxy,Sxz,Syz;
	vector<double> Sxx_b,Syy_b,Szz_b,Sxy_b,Sxz_b,Syz_b, S2_sum, S2_sum2;
	vector<double> weight_vec;
	coords *crd, *crd_ref;
	bool align;
	vector<int> align_res;
	int *sel, nsel,nblock,ndat,npair,blen,ii,jj;
	double x,y,z,xx,yy,zz,r2;
	FILE *outp;
	BaseITrajFile *trjx;

	S2pairs = "NULL";
	output = "NULL";
	nblock = 10;
	weight_file = "NULL";
	sel = NULL;
	nsel = 0;
	align = false;

	if (argc == 1) {
		fprintf(stdout,"%s\n",usage);
		exit(0);
	}
	while (1) {
		c=getopt(argc,argv,"ha:d:n:o:s:w:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage);
				exit(0);
				break;
			case 'd':
				S2pairs = optarg;
				break;
			case 'n':
				nblock = atoi(optarg);
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
	if (S2pairs == string("NULL")) {
		fprintf(stderr,"Specify pairs file with -d\n\n");
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
	ntrj = argc-optind;
	trjfiles.resize(ntrj);
	for (int i=0; i<ntrj; i++) {
		trjfiles[i] = argv[optind+i];
	}
	// count frames in all trajectories
	crd_ref = new coords(pdbfile.c_str(),"pdb");
	crd = new coords(pdbfile.c_str(),"pdb");
	crd->set_uniform_weights();
	crd_ref->set_uniform_weights();
	ndat = 0;
	for (int i=0; i<ntrj; i++) {
		int slen = trjfiles[i].size();
		if (slen-trjfiles[i].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			trjx = new DCDITrajFile(trjfiles[i].c_str());
			//ndat += trj->total_frames();
		} else if (slen-trjfiles[i].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			trjx = new XTCITrajFile(trjfiles[i].c_str());
			//ndat += trj->total_frames();
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[i].c_str());
			exit(1);
		}
		//natom = trj->num_atoms();
		//fprintf(stdout,"Number of atoms = %i\n",natom);
		while (trjx->frames_left()) {
			*trjx >> *crd;
			ndat++;
		}
		delete trjx;
	}
	blen = ndat/nblock;
	if (weight_file != string("NULL")) {
		read_weights(weight_file, weight_vec);
		if (weight_vec.size() != ndat) {
			fprintf(stderr,"Error!!");
			fprintf(stderr,"Weight vector of size: %i\n",weight_vec.size());
			fprintf(stderr,"Trajectory data of size: %i\n",ndat);
			exit(1);
		}
	}

	fprintf(stdout,"==============================================================\n");
	fprintf(stdout,"             Using PDB file: %s\n", pdbfile.c_str());
	fprintf(stdout,"  Pairs of nuclei from file: %s\n", S2pairs.c_str());
	fprintf(stdout,"              S2 output to: %s\n", output.c_str());
	if (ntrj > 0) {
		fprintf(stdout,"Using trajectories:\n");
		for (int i = 0; i<ntrj;i++) 
			fprintf(stdout,"\t%s\n",trjfiles[i].c_str());
		fprintf(stdout,"Total frames in all trajectories = %i\n",ndat);
		fprintf(stdout,"Number of blocks = %i\n",nblock);
		fprintf(stdout,"Block length = %i\n",blen);
	} else {
		fprintf(stderr,"No trajectories specified!!!\n");
		exit(1);
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

	parse_rdc(S2pairs, crd, pairs, true);
	npair = pairs.size();
	fprintf(stdout,"%5s %4s %5s -- %5s %4s %5s\n", "res_i", "at_i", "idx_i",
			"res_j", "at_j", "idx_j");
	fprintf(stdout,"==============================================================\n");
	for (int i=0; i<npair; i++) {
		fprintf(stdout,"%5i %4s %5i -- %5i %4s %5i\n",
				pairs[i].resi,pairs[i].atomi.c_str(),pairs[i].indi,
				pairs[i].resj,pairs[i].atomj.c_str(),pairs[i].indj);
	}

	Sxx.resize(npair); Sxx_b.resize(npair);
	Syy.resize(npair); Syy_b.resize(npair);
	Szz.resize(npair); Szz_b.resize(npair);
	Sxy.resize(npair); Sxy_b.resize(npair);
	Sxz.resize(npair); Sxz_b.resize(npair);
	Syz.resize(npair); Syz_b.resize(npair);
	S2_sum.resize(npair); S2_sum2.resize(npair);
	for (int p=0;p<npair;p++) {
		Sxx[p] = 0.; Syy[p] = 0.; Szz[p] = 0.; 
		Sxy[p] = 0.; Sxz[p] = 0.; Syz[p] = 0.; 
		S2_sum[p] = 0.; S2_sum2[p] = 0.;
	}
	int frame = 0;
	int pblock = -1;
	double com[3], mat[3][3];
	crd_ref->center_of_mass(com, sel, nsel);
	crd_ref->translate( -com[0], -com[1], -com[2] );
	int blockc = 0;
	for (int i=0; i<ntrj; i++) {
		BaseITrajFile *trj;
		int slen = trjfiles[i].size();
		if (slen-trjfiles[i].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			trj = new DCDITrajFile(trjfiles[i].c_str());
			//ndat += trj->total_frames();
		} else if (slen-trjfiles[i].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			trj = new XTCITrajFile(trjfiles[i].c_str());
			//ndat += trj->total_frames();
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[i].c_str());
			exit(1);
		}
		//natom = trj->num_atoms();
		//fprintf(stdout,"Number of atoms = %i\n",natom);
		while (trj->frames_left()) {
			//fprintf(stdout,"frame %i\n",frame);
			*trj >> *crd;
			if ( nsel != 0 ) {
				// do alignment
				crd->center_of_mass(com, sel, nsel);
				crd->translate( -com[0], -com[1], -com[2] );
				//fprintf(stdout,"com: %12.5f %12.5f %12.5f\n",com[0],com[1],com[2]);
				crd->find_kabsch_fit( *crd_ref, sel, nsel, mat);
				crd->rotate( mat );
			}
			int block = frame/blen;
			//fprintf(stdout,"frame %i block %i\n",frame,block);
			if (block>pblock) {
				pblock=block;
				if (block>0) {
					//fprintf(stdout,"block %i\n",block);
					for (int p=0;p<npair;p++) {
						double S2_i = calc_S2(Sxx_b[p],Syy_b[p],Szz_b[p],
								Sxy_b[p],Sxz_b[p],Syz_b[p],blen);
						//fprintf(stdout,"S2[%i] = %8.3f\n",p,S2_i);
						S2_sum[p] += S2_i;
						S2_sum2[p] += S2_i*S2_i;
						Sxx[p] += Sxx_b[p];
						Syy[p] += Syy_b[p];
						Szz[p] += Szz_b[p];
						Sxy[p] += Sxy_b[p];
						Sxz[p] += Sxz_b[p];
						Syz[p] += Syz_b[p];
					}
					blockc++;
				}
				for (int p=0;p<npair;p++) {
					Sxx_b[p] = 0.; Syy_b[p] = 0.; Szz_b[p] = 0.; 
					Sxy_b[p] = 0.; Sxz_b[p] = 0.; Syz_b[p] = 0.; 
				}
			}
			for (int p=0;p<npair;p++) {
				ii = pairs[p].indi;
				jj = pairs[p].indj;
				x = crd->Get_X(ii)-crd->Get_X(jj);
				y = crd->Get_Y(ii)-crd->Get_Y(jj);
				z = crd->Get_Z(ii)-crd->Get_Z(jj);
				xx = x*x;
				yy = y*y;
				zz = z*z;
				r2 = xx+yy+zz;
				Sxx_b[p] += xx/r2; 
				Syy_b[p] += yy/r2; 
				Szz_b[p] += zz/r2;
				Sxy_b[p] += x*y/r2; 
				Sxz_b[p] += x*z/r2; 
				Syz_b[p] += y*z/r2;
			}
			frame ++;
		}
		delete trj;
	}
	//fprintf(stdout,"frame = %i; nblock*blen = %i\n", frame, blen*nblock);
	if (frame==blen*nblock) {
		for (int p=0;p<npair;p++) {
			double S2_i = calc_S2(Sxx_b[p],Syy_b[p],Szz_b[p],
					Sxy_b[p],Sxz_b[p],Syz_b[p],blen);
			S2_sum[p] += S2_i;
			S2_sum2[p] += S2_i*S2_i;
			Sxx[p] += Sxx_b[p];
			Syy[p] += Syy_b[p];
			Szz[p] += Szz_b[p];
			Sxy[p] += Sxy_b[p];
			Sxz[p] += Sxz_b[p];
			Syz[p] += Syz_b[p];
		}
		blockc++;
	}

	for (int p=0; p<npair; p++) {
		double S2_block, S2_err, S2_glob;
		S2_glob = calc_S2(Sxx[p],Syy[p],Szz[p],
				Sxy[p],Sxz[p],Syz[p],nblock*blen);
		S2_block = S2_sum[p]/float(nblock);
		S2_err = S2_sum2[p]/float(nblock)-S2_block*S2_block;
		fprintf(stdout,"%4i %4s %4i %4s %8.3f %8.3f %8.3f\n",
				pairs[p].resi,pairs[p].atomi.c_str(),
				pairs[p].resj,pairs[p].atomj.c_str(),
				S2_glob,S2_err,S2_block);
		fprintf(outp,"%4i %4s %4i %4s %8.3f %8.3f %8.3f\n",
				pairs[p].resi,pairs[p].atomi.c_str(),
				pairs[p].resj,pairs[p].atomj.c_str(),
				S2_glob,S2_err,S2_block);
	}
	//fprintf(stdout,"blockc = %i\n",blockc);
	fclose(outp);

	return 0;
}
