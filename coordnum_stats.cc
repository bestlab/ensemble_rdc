/*
 * suitable comment here
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
//#include "rdc_funcs.h"

void strip(const string &s, string &st)
{
	st = "";
	for (int i=0;i<s.size(); i++) {
		if (isalnum(s[i]) || s[i]=='\'') {
			st += s[i];
		}
	}
}

void parse_pdb(const coords *crd, vector<int> &protein_heavy_atoms, vector<int> &water_oxygens)
{

	const int buf_len = 1024;
	char buf[buf_len];
	char tok[buf_len];
	int n_atom;
	bool foundi,foundj;
	char ati, atj;
	char *junk;
	int n_wat = 0;
	int n_protein_heavy = 0;
	string atom;
	n_atom = crd->Get_N();

	for (int a=0;a<n_atom;a++) {
		strip(crd->Get_Atomtype(a),atom);
		if (atom == string("OW")) { // water
			n_wat++;
		} else if (atom[0] == 'C' || atom[0] == 'N' || atom[0] == 'O' 
			|| atom[0] == 'S' ) { // heavy atom
			n_protein_heavy++;
		}
	}

	water_oxygens.resize(n_wat);
	protein_heavy_atoms.resize(n_protein_heavy);

	int w_ind = 0;
	int heavy_ind = 0;
	for (int a=0;a<n_atom;a++) {
		strip(crd->Get_Atomtype(a),atom);
		if (atom == string("OW")) { // water
			water_oxygens[w_ind++] = a;
		} else if (atom[0] == 'C' || atom[0] == 'N' || atom[0] == 'O' 
			|| atom[0] == 'S' ) { // heavy atom
			protein_heavy_atoms[heavy_ind++] = a;
		}
	}


	//fprintf(stdout,"Number of atoms = %i\n", n_atom);
	//fprintf(stdout,"Number of protein heavy atoms = %i\n", n_protein_heavy);
	//fprintf(stdout,"Number of water molecules = %i\n", n_wat);
	/*

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
	*/
}

const char *usage = "\n\n"
"        Usage\n"
"             coordnum -s pdb_file -o coordnum.dat trj1 trj2 ...trjN\n\n"
"\n\n";


int main(int argc, char **argv)
{
	string output,pdbfile;
	coords *crd;
	vector<int> protein_heavy_atoms, water_oxygens, neighbours, neighbour_stats;
	vector<string> trjfiles;
	int c,ntrj;
	double x_kl, y_kl, z_kl, r_kl2, r_cut2;
	double r_cut = 3.5;
	double xtal[6];

	while (1) {
		c=getopt(argc,argv,"ho:s:r:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage);
				exit(0);
				break;
			case 'o':
				output = optarg;
				break;
			case 'r':
				r_cut = atof(optarg);
				break;
			case 's':
				pdbfile = optarg;
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}
	r_cut2 = r_cut * r_cut;
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

	fprintf(stdout,"==============================================================\n");
	fprintf(stdout,"             Using PDB file: %s\n", pdbfile.c_str());
	fprintf(stdout,"                  output to: %s\n", output.c_str());
	if (ntrj > 0) {
		fprintf(stdout,"Using trajectories:\n");
		for (int i = 0; i<ntrj;i++) 
			fprintf(stdout,"\t%s\n",trjfiles[i].c_str());
	} else {
		fprintf(stdout,"Exiting, no trajectories specified\n");
	}
	fprintf(stdout,"==============================================================\n\n");

	// open for output:
	FILE *outp = fopen(output.c_str(),"w");
	crd = new coords(pdbfile.c_str(),"pdb");
	parse_pdb(crd,protein_heavy_atoms,water_oxygens);

	ntrj = trjfiles.size();
	int frame = 0;
	BaseITrajFile *trj;
	int nwat =  water_oxygens.size();
	int nprot =  protein_heavy_atoms.size();
	int maxn = 16;
	neighbours.resize(nwat);
	neighbour_stats.resize(maxn);
	
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
		int natom = trj->num_atoms();
		while (trj->frames_left()) {
			*trj >> *crd;
			trj-> get_crystal_data(xtal);
			double a = xtal[0]/sqrt(3); // assume tr oh for now!
			double L = 2.*a;

			//for (int q=0; q<6; q++) {
			//	fprintf(stdout,"%12.6f ",xtal[q]);
			//}
			//fprintf(stdout,"\n");
			//fflush(stdout);

			frame++;
			for (int k=0; k<nwat; k++) {
				neighbours[k] = 0;
			}


			for (int k=1; k<nwat; k++) {
				for (int l=0; l<k; l++) { 
					int wk = water_oxygens[k];
					int wl = water_oxygens[l];
					x_kl = crd->Get_X(wk)-crd->Get_X(wl);
					y_kl = crd->Get_Y(wk)-crd->Get_Y(wl);
					z_kl = crd->Get_Z(wk)-crd->Get_Z(wl);
					/*
					 // pbc - To be added properly
					x_kl -= L*round(x_kl/L);
					y_kl -= L*round(y_kl/L);
					z_kl -= L*round(z_kl/L);
					if((fabs(x_kl)+fabs(y_kl)+fabs(z_kl))
							>= 1.5*a) {
						if (x_kl>=0)
							x_kl= x_kl - a;
						else x_kl= x_kl + a;

						if (y_kl>=0)
							y_kl= y_kl - a;
						else y_kl= y_kl + a;

						if (z_kl>=0)
							z_kl= z_kl - a;
						else z_kl= z_kl + a;

					}
					*/
					r_kl2 = x_kl*x_kl + y_kl*y_kl + z_kl*z_kl; 
					if (r_kl2 < r_cut2) {
						neighbours[k]++;
						neighbours[l]++;
					}
				}
			}
			for (int k=0; k<maxn; k++) {
				neighbour_stats[k] = 0;
			}
			for (int k=0; k<nwat; k++) {
				bool isclose=false;
				for (int l=0; l<nprot; l++) {
					int wk = water_oxygens[k];
					int pl = protein_heavy_atoms[l];
					x_kl = crd->Get_X(wk)-crd->Get_X(pl);
					y_kl = crd->Get_Y(wk)-crd->Get_Y(pl);
					z_kl = crd->Get_Z(wk)-crd->Get_Z(pl);
					r_kl2 = x_kl*x_kl + y_kl*y_kl + z_kl*z_kl; 
					if (r_kl2 < r_cut2) {
						isclose=true;
						break;
					}
				}
				if (isclose) {
					if (neighbours[k] < maxn) {
						neighbour_stats[neighbours[k]]++;
					} else {
						fprintf(stderr,"Water %i has more than %i neighbours!!!!\n",
							k,maxn);
						fflush(stderr);
						exit(1);
					}
				}
			}
			for (int k=0; k<maxn; k++) {
				fprintf(outp,"%6i ",neighbour_stats[k]);
			}
			fprintf(outp,"\n");

			if (frame % 100==0) {
				fprintf(stdout,"frame %i \n",frame);
				fflush(stdout);
				fflush(outp);
			}
		}
		delete trj;
	}
	fprintf(stdout,"Total frames = %i\n",frame);

	delete crd;
	return 0;
}
