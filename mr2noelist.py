#!/usr/bin/env python

import sys,math

pdblines = filter(lambda x: x[0:4] == "ATOM" and x[17:20] not in ["SOL","Na+","Cl-"], \
		open(sys.argv[1]).readlines())

atom_map = {}
for l in pdblines:
	atom_no = int(l[6:11])
	atom_name = l[11:17].strip()
	res_no = int(l[23:26])
	if res_no not in atom_map.keys():
		atom_map[res_no] = {}
	atom_map[res_no][atom_name] = atom_no;

noelines = filter(lambda x: x[0:6] == "assign", open(sys.argv[2]).readlines())

#for x in map(lambda x: x.split(), noelines):
for l in noelines:
	bra = l.find('(')
	ket = l.find(')',bra+1)
	a = l[bra+1:ket].split()
	bra = l.find('(',ket)
	ket = l.find(')',bra+1)
	b = l[bra+1:ket].split()
	dc = l[ket+1:].split()
	resi = int(a[1])
	atomi = a[4].upper()
	if atomi == "HN":
		atomi = "H"
	resj = int(b[1])
	atomj = b[4].upper()
	if atomj == "HN":
		atomj = "H"
	d,dminus,dplus = map(float, dc)
	idx_i = []
	idx_j = []
	for k in atom_map[resi].keys():
		if k == atomi:
			idx_i.append(atom_map[resi][k])
		elif atomi.find("#") > 0:
			nc = atomi.find("#")
			if k[:nc] == atomi[:nc]:
				idx_i.append(atom_map[resi][k])
		elif atomi.find("*") > 0:
			nc = atomi.find("*")
			if k[:nc] == atomi[:nc]:
				idx_i.append(atom_map[resi][k])
	if len(idx_i) == 0:
		print "could not find %i %s" % (resi, atomi)
	
	for k in atom_map[resj].keys():
		if k == atomj:
			idx_j.append(atom_map[resj][k])
		elif atomj.find("#") > 0:
			nc = atomj.find("#")
			if k[:nc] == atomj[:nc]:
				idx_j.append(atom_map[resj][k])
		elif atomj.find("*") > 0:
			nc = atomj.find("*")
			if k[:nc] == atomj[:nc]:
				idx_j.append(atom_map[resj][k])
	if len(idx_j) == 0:
		print "could not find %i %s" % (resj, atomj)

	sys.stdout.write("%3i %4s %3i %4s" % (resi,atomi,resj,atomj))
	sys.stdout.write(" %2i" % len(idx_i))
	for xx in range(len(idx_i)):
		sys.stdout.write(" %5i"%(idx_i[xx]))
	sys.stdout.write(" %2i" % len(idx_j))
	for xx in range(len(idx_j)):
		sys.stdout.write(" %5i"%(idx_j[xx]))
	sys.stdout.write("%8.3f %8.3f\n"%(d-dminus,d+dplus))


