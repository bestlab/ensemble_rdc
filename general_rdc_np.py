#!/usr/bin/env python

# somewhat revamped script to fit an ensemble of structures (or a single
# structure) to a set of rdc's

#import PDBFile, sys,math,numarray,numarray.linear_algebra, getopt
import PDBFile, sys,math,numpy,numpy.linalg, getopt

gyro = { 'H': 2.6752e8, 'C': 6.7287e7, 'N': 2.712e7 }
mu0 = 4*math.pi*1e-7
hcross = 1.05e-34

def calc_Q(calc,exp):
	sum_dev2 = 0.0
	sum_calc2 = 0.0 
	for i in range(len(calc)):
		sum_dev2 += (calc[i]-exp[i])**2
		sum_calc2 += calc[i]**2
	return math.sqrt(sum_dev2/sum_calc2)

def gyromag(at1,at2):
	if at1 == "H":
		at1 = "HN"
	if at2 == "H":
		at2 = "HN"
	if (at1=="N" and at2=="HN") or (at1=="HN" and at2=="N"):
		return gyro[at1[0]]*gyro[at2[0]]
	elif (at1=="HN" and at2=="C") or (at1=="C" and at2=="HN"):
		return gyro[at1[0]]*gyro[at2[0]]
	elif (at1=="CA" and at2=="HA") or (at1=="HA" and at2=="CA"):
		return gyro[at1[0]]*gyro[at2[0]]
	elif (at1=="N" and at2=="C") or (at1=="C" and at2=="N"):
		return gyro[at1[0]]*gyro[at2[0]]
	elif (at1=="CA" and at2=="C") or (at1=="C" and at2=="CA"):
		return gyro[at1[0]]*gyro[at2[0]]
	elif (at1=="CA" and at2=="CB") or (at1=="CB" and at2=="CA"):
		return gyro[at1[0]]*gyro[at2[0]]
	elif (at1[0]=="C" and at2[0]=="C"):
		#print "GOT HERE", at1,at2
		return gyro["C"]*gyro["C"]
		#return gyro["C"]*gyro["H"]/(-3.17)
		#return gyro["C"]*gyro["H"]*(3.17)
		#return gyro[at1[0]]*gyro[at2[0]]
	elif (at1[0]=="C" and at2[0]=="H") or (at1[0]=="H" and at2[0]=="C"):
		return gyro[at1[0]]*gyro[at2[0]]
	else:
		print "could not find bond length for atoms: ", at1, at2
	
def bondlen(at1,at2):
	if at1 == "H":
		at1 = "HN"
	if at2 == "H":
		at2 = "HN"
	if (at1=="N" and at2=="HN") or (at1=="HN" and at2=="N"):
		return 1.04e-10
	elif (at1=="HN" and at2=="C") or (at1=="C" and at2=="HN"):
		return 2.04e-10
	elif (at1=="CA" and at2=="HA") or (at1=="HA" and at2=="CA"):
		return 1.12e-10
	elif (at1=="N" and at2=="C") or (at1=="C" and at2=="N"):
		return 1.33e-10
	elif (at1=="CA" and at2=="C") or (at1=="C" and at2=="CA"):
		return 1.53e-10
	elif (at1=="CA" and at2=="CB") or (at1=="CB" and at2=="CA"):
		return 1.53e-10
	elif (at1[0]=="C" and at2[0]=="C"):
		return 1.53e-10
		#return 1.12e-10
	elif (at1[0]=="C" and at2[0]=="H") or (at1[0]=="H" and at2[0]=="C"):
		return 1.12e-10
	else:
		print "could not find bond length for atoms: ", at1, at2
	
def parse_rdc(pdat,file):
	inp = open(file)
	rdc_list = []
	seq = pdat.get_seq()
	for line in inp.readlines():
		sl = line.split()
		a = sl[1][0]
		b = sl[3][0]
		r1 = int(sl[0])
		a1 = sl[1]
		r2 = int(sl[2])
		a2 = sl[3]
		if a1=='HN':
			a1 = 'H'
		if a2=='HN':
			a2 = 'H'
		#if seq[r1-1] == 'I' and a1 == "CD1":
		#	a1 = "CD"
		#if seq[r2-1] == 'I' and a2 == "CD1":
		#	a2 = "CD"
		i1 = pdat.get_index(r1,a1)
		i2 = pdat.get_index(r2,a2)
		if i1 < 0:
			print "could not find atom: ", r1, a1
			sys.exit(1)
		if i2 < 0:
			print "could not find atom: ", r2, a2
			sys.exit(1)
		blen = pdat.dist(i1,i2)*1e-10
		if a1[0] not in gyro.keys() or a2[0] not in gyro.keys():
			continue


		Dmax = mu0*hcross*(gyromag(a1,a2))/(4*(math.pi**2)*(bondlen(a1,a2))**3)
		#Dmax = mu0*hcross*(gyromag(a1,a2))/(4*(math.pi**2)*blen**3)


		#Dmax = mu0*hcross*(gyromag(a1,a2))/(4*(math.pi**2)*(blen**3))
		#Dmax = mu0*hcross*(gyro[a]*gyro[b])/(4*(math.pi**2)*(bond[a+b])**3)
		#Dmax = mu0*hcross*(gyro[a]*gyro[b])/(4*(math.pi**2)*(blen)**3)
		#print Dmax
		rdc_list.append((int(sl[0]),sl[1],int(sl[2]),sl[3],float(sl[4]),Dmax,i1,i2))
	return rdc_list

def calc_coef_nmr(nmr_file, rdc_list,refpdb):
	pfile = PDBFile.PDBFile(nmr_file)
	nmodel = pfile.get_nmodel()
	
	nrdc = len(rdc_list)
	coef_mat = numpy.zeros([nrdc,5],numpy.float64)
	rdc_vec = numpy.array(map(lambda x: x[4]/x[5],rdc_list),numpy.float64)
	
	refseq = refpdb.get_seq()
	plen = len(refseq)
	sel = range(1,plen+1)
	struc = 0
	for model in range(1,nmodel+1):
		print "...model " + str(model)
		pdat = PDBFile.PDBData(pfile,'ignore', 0, -1, model, 1)
		pdat.align(refpdb,sel,sel,"CA")
		for d in range(nrdc):
			i = rdc_list[d][6]
			j = rdc_list[d][7]
			mu_x = pdat.get_x(i)-pdat.get_x(j)
			mu_y = pdat.get_y(i)-pdat.get_y(j)
			mu_z = pdat.get_z(i)-pdat.get_z(j)
			mu_r = math.sqrt(mu_x**2+mu_y**2+mu_z**2)
			mu_x /= mu_r
			mu_y /= mu_r
			mu_z /= mu_r
			coef_mat[d][0] += mu_x**2-mu_z**2
			coef_mat[d][1] += mu_y**2-mu_z**2
			coef_mat[d][2] += 2.0*mu_x*mu_y
			coef_mat[d][3] += 2.0*mu_x*mu_z
			coef_mat[d][4] += 2.0*mu_y*mu_z
		struc+=1
	coef_mat = coef_mat / float(struc)
	return coef_mat, rdc_vec
	
# whether to fit structures to template (first structure)
# before computing intramolecular contrib to rdc
fit = 0

def calc_coef_ens(files, rdc_list,refpdb,convert_charmm):
	nrdc = len(rdc_list)
	coef_mat = numpy.zeros([nrdc,5],numpy.float64)
	rdc_vec = numpy.array(map(lambda x: x[4]/x[5],rdc_list),numpy.float64)
	
	refseq = refpdb.get_seq()
	plen = len(refseq)	
	sel = range(1,plen+1)
	struc = 0
	for file in files:
		print "...file " + file
		pfilew = PDBFile.PDBFile(file)
		pdat = PDBFile.PDBData(pfilew)
		if convert_charmm:
			pdat.pdbify()
		#pdat.align(refpdb,sel,sel,"CA")
		for d in range(nrdc):
			i = rdc_list[d][6]
			j = rdc_list[d][7]
			mu_x = pdat.get_x(i)-pdat.get_x(j)
			mu_y = pdat.get_y(i)-pdat.get_y(j)
			mu_z = pdat.get_z(i)-pdat.get_z(j)
			mu_r = math.sqrt(mu_x**2+mu_y**2+mu_z**2)
			mu_x /= mu_r
			mu_y /= mu_r
			mu_z /= mu_r
			coef_mat[d][0] += mu_x**2-mu_z**2
			coef_mat[d][1] += mu_y**2-mu_z**2
			coef_mat[d][2] += 2.0*mu_x*mu_y
			coef_mat[d][3] += 2.0*mu_x*mu_z
			coef_mat[d][4] += 2.0*mu_y*mu_z
		struc+=1

	coef_mat = coef_mat / float(struc)
	return coef_mat, rdc_vec

Usage = """

Usage:
	%s [-c] [-n nmr.pdb] [-d rdcfit.dat] [-b rdcbackcalc.dat] [-o output_fit.dat] [-p]

	examples:
		%s -n nmr.pdb -d rdc.dat -o fit.dat
			...will fit the NMR ensemble in nmr.pdb to rdc.dat and print
			the alignment tensor (in comment line) and fitted data to
			fit.dat
		%s -p -d rdc.dat -o fit.dat file1.pdb file2.pdb ... fileN.pdb
			...will fit the ensemble of structures given by file1.pdb...fileN.pdb
			to rdcfit.dat and print the alignment tensor (in comment line) and 
			fitted data to fit.dat
		%s -p -d rdc.dat -b bc.dat -o fit.dat file1.pdb file2.pdb ... fileN.pdb
		        ... will do as above, but rdc.dat will be used to fit the alignment 
			tensor and bc.dat will be used for backcalculating rdc's (useful if
			only a fraction of the rdc's (e.g. backbone rdcs in SS) are to be
			used to fit the alignment tensor).


""" % (sys.argv[0],sys.argv[0],sys.argv[0],sys.argv[0])

options, files = getopt.getopt(sys.argv[1:], 'chn:d:b:o:p')

rdc_fit_file = ""
rdc_bc_file = ""
output_file = "junk.pdb"
nmr = 0
pdb_ensemble = 0
charmm_names = 0
for opt, val in options:
	if opt in [ '-h' ] or len(sys.argv) == 1:
		print Usage
		sys.exit(0)
	elif opt in [ '-n' ]:
		nmr_file = val
		nmr = 1
	elif opt in [ '-d' ]:
		rdc_fit_file = val
	elif opt in [ '-b' ]:
		rdc_bc_file = val
	elif opt in [ '-o' ]:
		output_file = val
	elif opt in [ '-p' ]:
		pdb_ensemble = 1
	elif opt in [ '-c' ]:
		charmm_names = 1
	else:
		print "Unknown option %s" % opt
		print Usage
		sys.exit(1)

if rdc_fit_file == "":
	print Usage
	sys.exit(1)

if nmr and pdb_ensemble:
	print "can only specify one of -n or -p"
	sys.exit(1)

# default
if rdc_bc_file == "":
	rdc_bc_file = rdc_fit_file

sys.stdout.write("Will fit rdc's with the following input:\n")
sys.stdout.write("RDC fitting data: %s\n" % (rdc_fit_file))
sys.stdout.write("RDC back calculation data: %s\n" % (rdc_bc_file))
sys.stdout.write("Output file: %s\n" % (output_file))
if charmm_names:
	sys.stdout.write("Will convert atom names from CHARMM format\n")
if nmr:
	sys.stdout.write("Using nmr-style pdb file: %s\n" % (nmr_file))
else:
	sys.stdout.write("Using the following pdb files:\n")
	for file in files:
		sys.stdout.write("\t%s\n" % (file))

if nmr:
	pfile = PDBFile.PDBFile(nmr_file)
	refpdb = PDBFile.PDBData(pfile,'ignore', 0, -1, 1, 1)
else:
	pfile = PDBFile.PDBFile(files[0])
	refpdb = PDBFile.PDBData(pfile)
	if charmm_names:
		refpdb.pdbify()

print "GOT HERE"
rdc_fit = parse_rdc(refpdb,rdc_fit_file)
rdc_bc = parse_rdc(refpdb,rdc_bc_file)

if nmr:
	coef_mat_fit, rdc_vec_fit = calc_coef_nmr(nmr_file,rdc_fit,refpdb)
	if rdc_fit_file != rdc_bc_file:
		coef_mat_bc, rdc_vec_bc = calc_coef_nmr(nmr_file,rdc_bc,refpdb)
	else:
		coef_mat_bc, rdc_vec_bc = coef_mat_fit, rdc_vec_fit
else:
	coef_mat_fit, rdc_vec_fit = calc_coef_ens(files,rdc_fit,refpdb,charmm_names)
	if rdc_fit_file != rdc_bc_file:
		coef_mat_bc, rdc_vec_bc = calc_coef_ens(files,rdc_bc,refpdb,charmm_names)
	else:
		coef_mat_bc, rdc_vec_bc = coef_mat_fit, rdc_vec_fit

solution = numpy.linalg.lstsq(coef_mat_fit,rdc_vec_fit)

S = solution[0]

Sxx = S[0]
Syy = S[1]
Szz = -Sxx-Syy
Sxy = S[2]
Sxz = S[3]
Syz = S[4]

#bc = numpy.matrixmultiply(coef_mat_bc,S)
bc = numpy.dot(coef_mat_bc,S)

outp = open(output_file, "w")
outp.write("# Fitted rdc's with the following input:\n")
outp.write("# RDC fitting data: %s\n" % (rdc_fit_file))
outp.write("# RDC back calculation data: %s\n" % (rdc_bc_file))
if nmr:
	outp.write("# Using nmr-style pdb file: %s\n" % (nmr_file))
else:
	outp.write("# Using an ensemble of pdb files:\n")

outp.write("# Alignment tensor\n")
outp.write("# Sxx = %12.5e\n"%(Sxx))
outp.write("# Syy = %12.5e\n"%(Syy))
outp.write("# Szz = %12.5e\n"%(Szz))
outp.write("# Sxy = %12.5e\n"%(Sxy))
outp.write("# Sxz = %12.5e\n"%(Sxz))
outp.write("# Syz = %12.5e\n"%(Syz))

# calculate Q - unnormalized
exp_norm = []
exp_unnorm = []
calc_norm = []
calc_unnorm = []
for i in range(len(rdc_bc)):
	exp_norm.append(rdc_bc[i][4]/rdc_bc[i][5])
	exp_unnorm.append(rdc_bc[i][4])
	calc_unnorm.append(bc[i]*rdc_bc[i][5])
	calc_norm.append(bc[i])

Q_unnorm = calc_Q(calc_unnorm,exp_unnorm)
Q_norm = calc_Q(calc_norm,exp_norm)
outp.write("#\n# Q_unnorm = %12.5f\n"%(Q_unnorm))
outp.write("# Q_norm = %12.5f\n"%(Q_norm))

outp.write("#\n# |D| refers to the dipolar coupling rescaled such that\n")
outp.write("# D_max(A,X) = D_max(N,H)\n#\n")
outp.write("#%4s %4s %5s %4s %8s %8s %8s %8s\n"\
		% ("res1","at1","res2","at2","D_calc","D_exp",
			"|D_calc|","|D_exp|"))

# scale normalized RDC's by Dmax_NH to make them more friendly numbers
Dmax_NH = mu0*hcross*(gyromag("N","H"))/(4*(math.pi**2)*(bondlen("N","H"))**3)

for i in range(len(rdc_bc)):
	outp.write("%5i %4s %5i %4s %8.3f %8.3f %8.3f %8.3f\n"\
	%(rdc_bc[i][0],rdc_bc[i][1],rdc_bc[i][2],rdc_bc[i][3],\
	calc_unnorm[i],exp_unnorm[i],calc_norm[i]*Dmax_NH,exp_norm[i]*Dmax_NH))

outp.close()


