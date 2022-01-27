/*
 * calculate Q from MD trajectory
 *
 */

#define VERBOSE
#include <unistd.h> // for getopt
#include <cstdio> 
#include <cstdlib> 
#include <fstream> 
#include <cmath> 
#include "traj.h"

const char *usage = "\n\n"
"        Usage\n"
"             noe_corfun -n noelist.dat -s xyz.pdb -o output.dat trj1 trj2 ...trjN\n\n"
"\n\n";

void ReadNOE(const string &file, vector<int> &idx_i, vector<int > &idx_j,
		vector<int> &res_i, vector<int> &res_j,
		vector<string> &atom_i, vector<string> &atom_j  )
{
	ifstream inp(file.c_str());
	const int bufsz = 1024;
	char buf[bufsz];
	int p,q,nc,ni,nj;
	nc = 0;
	while (inp.good()) {
		inp.getline(buf,bufsz,'\n');
		nc++;
	}
	nc--;
	inp.close();
	idx_i.resize(nc);
	idx_j.resize(nc);
	res_i.resize(nc);
	res_j.resize(nc);
	atom_i.resize(nc);
	atom_j.resize(nc);
	ifstream inp2(file.c_str());
	for (int t=0; t<nc; t++) {
		inp2 >> res_i[t] >> atom_i[t] >> res_j[t] >> atom_j[t];
		inp2 >> idx_i[t] >> idx_j[t] ;
	}
	inp2.close();
	return;
}

int main(int argc, char **argv)
{
	const double hbar = 1.38e-23/2*M_PI; // J
	const double gamma_H = 2.6752e8; // C/kg
	int c,ntrj,natom,nnoe;
	float *X,*Y,*Z;
	double * rij_ave;
	int ii,jj;
	double dx,dy,dz,dr2;
	string itrj,output,noe_file;
	vector<string> trjfiles;
	vector<int> idx_i,idx_j;
	vector<int> res_i,res_j;
	vector<string> atom_i,atom_j;
	vector<double> rij_lo, rij_hi;
	FILE *outp;
	BaseITrajFile *trj;
	int L,nmon;
	bool nmer;

	nmer = false;
	output = "NULL";
	noe_file = "NULL";
	while (1) {
		c=getopt(argc,argv,"hn:q:o:");
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
			case 'n':
				noe_file = optarg;
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}
	if (output == string("NULL")) {
		fprintf(stderr,"Specify output file with -o\n\n");
		fprintf(stderr,"%s\n",usage);
		exit(1);
	}
	if (noe_file == string("NULL")) {
		fprintf(stderr,"Specify noe list file file with -n\n\n");
		fprintf(stderr,"%s\n",usage);
		exit(1);
	}
	ntrj = argc-optind;
	if (ntrj == 0) {
		fprintf(stderr,"No trajectories specified.\n\n");
		fprintf(stderr,"%s\n",usage);
		exit(1);
	}
	trjfiles.resize(ntrj);
	for (int i=0; i<ntrj; i++) {
		trjfiles[i] = argv[optind+i];
	}

	fprintf(stdout,"==============================================================\n");
	fprintf(stdout,"Reading atom pairs from: %s\n", noe_file.c_str());
	fprintf(stdout,"          NOE output to: %s\n", output.c_str());
	fprintf(stdout,"Using trajectories:\n");
	for (int i = 0; i<ntrj;i++) 
		fprintf(stdout,"\t%s\n",trjfiles[i].c_str());
	fprintf(stdout,"==============================================================\n\n");

	ReadNOE(noe_file,idx_i,idx_j,res_i,res_j,atom_i,atom_j);

	fprintf(stdout,"==============================================================\n");
	fprintf(stdout,"Reading atom pairs from: %s\n", noe_file.c_str());
	fprintf(stdout,"          NOE output to: %s\n", output.c_str());
	fprintf(stdout,"Using trajectories:\n");
	for (int i = 0; i<ntrj;i++) 
		fprintf(stdout,"\t%s\n",trjfiles[i].c_str());
	fprintf(stdout,"==============================================================\n\n");
	nnoe = idx_i.size();
#ifdef VERBOSE
	fprintf(stdout,"%5s %5s %5s %5s %6s %6s\n","Res_i","At_i","Res_j","At_j","idx_i","idx_j");
	for (int i=0; i<nnoe; i++) {
		fprintf(stdout,"%5i %5s %5i %5s %6i %6i\n",
				res_i[i],atom_i[i].c_str(),res_j[i],atom_j[i].c_str(),idx_i[i],idx_j[i]);
	}

#endif // VERBOSE
	
	// count the number of atoms
	int slen = trjfiles[0].size();
	if (slen-trjfiles[0].rfind(string(".dcd"))==4) {
		fprintf(stdout,"dcd file\n");
		trj = new DCDITrajFile(trjfiles[0].c_str());
	} else if (slen-trjfiles[0].rfind(string(".xtc"))==4) {
		fprintf(stdout,"xtc file\n");
		trj = new XTCITrajFile(trjfiles[0].c_str());
	} else {
		fprintf(stderr,"unknown trajectory file type:\n");
		fprintf(stderr,"%s\n",trjfiles[0].c_str());
		exit(1);
	}
	natom = trj->num_atoms();
	if (nmer)
		L = natom/nmon;
	delete trj;
	X = new float[natom];
	Y = new float[natom];
	Z = new float[natom];
	rij_ave = new double[nnoe];
	for (int i=0; i<nnoe; i++) 
		rij_ave[i] = 0.0;
	// open for output:
	outp = fopen(output.c_str(),"w");

	int total_fr = 0; // need this to allocate memory
	for (int T=0; T<ntrj; T++) {
		int slen = trjfiles[T].size();
		if (slen-trjfiles[T].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			trj = new DCDITrajFile(trjfiles[T].c_str());
		} else if (slen-trjfiles[T].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			trj = new XTCITrajFile(trjfiles[T].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[T].c_str());
			exit(1);
		}
		while (trj->frames_left()) {
			trj->read_frame(X, Y, Z, natom);
			total_fr++;
		}
		delete trj;
	}

	fprintf(stdout,"Total number of frames = %i\n",total_fr);
	fprintf(stdout,"MX_INT = %i\n",INT_MAX);

	/*
	int frame = 0;
	for (int T=0; T<ntrj; T++) {
		int slen = trjfiles[T].size();
		if (slen-trjfiles[T].rfind(string(".dcd"))==4) {
			fprintf(stdout,"dcd file\n");
			trj = new DCDITrajFile(trjfiles[T].c_str());
		} else if (slen-trjfiles[T].rfind(string(".xtc"))==4) {
			fprintf(stdout,"xtc file\n");
			trj = new XTCITrajFile(trjfiles[T].c_str());
		} else {
			fprintf(stderr,"unknown trajectory file type:\n");
			fprintf(stderr,"%s\n",trjfiles[T].c_str());
			exit(1);
		}
		natom = trj->num_atoms();
		fprintf(stdout,"Number of atoms = %i\n",natom);
		while (trj->frames_left()) {
			trj->read_frame(X, Y, Z, natom);
			frame++;
			for (int t=0; t<nnoe; t++) {
				double r2min = 10000.;
				for (int p=0; p<i[t].size(); p++) {
					ii = i[t][p] -1;
					for (int q=0; q<j[t].size(); q++) {
						jj = j[t][q] -1;
						if (nmer) {
							for (int a=0;a<nmon;a++) {
								int aoff = a*L;
								for (int b=0;b<nmon;b++) {
									int boff = b*L;
									dx = X[ii+aoff]-X[jj+boff];
									dy = Y[ii+aoff]-Y[jj+boff];
									dz = Z[ii+aoff]-Z[jj+boff];
									dr2 = dx*dx+dy*dy+dz*dz;
									if (dr2<r2min) 
										r2min = dr2;
								}
							}
								
						} else {
							dx = X[ii]-X[jj];
							dy = Y[ii]-Y[jj];
							dz = Z[ii]-Z[jj];
							dr2 = dx*dx+dy*dy+dz*dz;
							//fprintf(stdout,"i=%i ; j=%i; rij=%8.3f\n",ii,jj,sqrt(dr2));
							if (dr2<r2min) 
								r2min = dr2;
						}
					}
				}
				rij_ave[t] += 1./(r2min*r2min*r2min);
			}
		}
		fprintf(outp,"# %3s %3s %4s %3s %4s %8s %8s %8s\n","idx","ri","ati",
				"rj","atj","rij_ave","rij_lo","rij_hi");
		for (int p=0; p<nnoe; p++) {
			rij_ave[p] = pow(rij_ave[p]/double(frame),-1./6.);
			fprintf(outp,"%5i %3i %4s %3i %4s %8.3f %8.3f %8.3f\n",p,resi[p],atomi[p].c_str(),
					resj[p],atomj[p].c_str(),rij_ave[p],rij_lo[p],rij_hi[p]);
		}


		delete trj;
	}

	*/
	delete [] X;
	delete [] Y;
	delete [] Z;
	delete [] rij_ave;
	return 0;
}
