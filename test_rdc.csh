#!/bin/csh -f

# reference test - the python version
./general_rdc_np.py -p -d data/rdc_500_bb.dat -o testout/python_fit.dat data/1ubq_amber03_tip3p_prot.pdb 

# c++ version
# single structure
./rdc -s data/1ubq_amber03_tip3p_prot.pdb -d data/rdc_500_bb.dat -o testout/c_single_fit.dat 

exit
# trajectory
#./rdc -s 1ubq_amber03_tip3p_ions.pdb -d rdc_500_bb.dat -o c_single_fit.dat -a 3-20:30-45:50-60 \
#	1ubq_amber03_tip3p_unres.xtc

./rdc -s 1ubq_amber03_tip3p_ions.pdb -d rdc_500_bb.dat -o c_trj_fit_weight.dat  -a 3-65 \
	-w wrongweights.dat 1ubq_amber03_tip3p_unres.xtc 

./rdc -s 1ubq_amber03_tip3p_ions.pdb -d rdc_500_bb.dat -o c_trj_fit_oweight.dat  -a 3-65 \
	-w otherweights.dat 1ubq_amber03_tip3p_unres.xtc 

./rdc -s 1ubq_amber03_tip3p_ions.pdb -d rdc_500_bb.dat -o c_trj_fit_noweight.dat -a 3-65 \
	 1ubq_amber03_tip3p_unres.xtc 

./rdc -s 1ubq_amber03_tip3p_ions.pdb -d rdc_500_bb.dat -o c_trj_fit_noweight_noalign.dat \
	 1ubq_amber03_tip3p_unres.xtc 

./rdc -s test.pdb -d rdc_500_bb.dat -o testpdb_fit.dat 
