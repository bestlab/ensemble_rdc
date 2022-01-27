#!/usr/bin/env python2.1

#from Numeric import *
#from LinearAlgebra import *
#from numarray import *
from numpy import *
from numpy.linalg import *
#from numarray.linear_algebra import *
from string import *
import time, math, gzip

"""Utilities for accessing PDB files from Python
scripts in a sensible way"""

#class reference:
#	def __init__(self, authors, title, editor, journal, publisher):

# update this to include additional modified AA from pdb
# files
aalist = { 'ALA': 'A', 'GLY': 'G', 'THR': 'T', 'TYR': 'Y',
	   'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'TRP': 'W',
	   'GLU': 'E', 'ASP': 'D', 'SER': 'S', 'ASN': 'N',
	   'GLN': 'Q', 'PRO': 'P', 'PHE': 'F', 'ARG': 'R',
	   'CYS': 'C', 'HIS': 'H', 'LYS': 'K', 'MET': 'M',
	   'CGU': 'E', 'DBY': 'Y', 'GLZ': 'G', 'GLQ': 'E',
	   'HSD': 'H', 'HEM': 'X', 'ABA': 'B', 'CSO': 'C',
	   'ASPP': 'D'}
	   # ABA aminobutyric acid

# this does not need updating as protonation states will not
# be indicated
aamap = { 'A':'ALA', 'G':'GLY', 'T':'THR', 'Y':'TYR',
	  'V':'VAL', 'L':'LEU', 'I':'ILE', 'W':'TRP',
	  'E':'GLU', 'D':'ASP', 'S':'SER', 'N':'ASN',
	  'Q':'GLN', 'P':'PRO', 'F':'PHE', 'R':'ARG',
	  'C':'CYS', 'H':'HIS', 'K':'LYS', 'M':'MET' }

def onel (threel):
	"""Map a three letter amino acid code to one letter"""
	try:
		return aalist[ threel.upper() ]
	except KeyError:
		return 'X'
		#print "Unknown three-letter code: ", threel
		#print "Need to update aalist"
		#import sys
		#sys.exit(1)
		
def threel(onel):
	"""Map a one letter amino acid code to three letter"""
	try:
		return aamap[ onel.upper() ]
	except KeyError:
		print "Unknown one-letter code: ", onel
		import sys
		sys.exit(1)

def onel_seq(threel_seq):
	"""Map three letter AA sequence to one letter sequence"""
	os = ""
	for x in threel_seq:
		os += onel(x)
	return os

def translate( X, Y, Z, delta ):
	X += delta[0]
	Y += delta[1]
	Z += delta[2]

def transform( X, Y, Z, symop, box ):
	for i in range(len(X)):
		X[i] = X[i]*symop[0][0] + symop[0][1]*box[0]
		Y[i] = Y[i]*symop[1][0] + symop[1][1]*box[1]
		Z[i] = Z[i]*symop[2][0] + symop[2][1]*box[2]
		#X[i] = -X[i]+0.5*box[0]
		#Y[i] = -Y[i]+0.0*box[1]
		#Z[i] =  Z[i]+0.5*box[0]
	#xtr = 
	#print xtr
	#ytr = symop[1][1]*box[1]
	#print ytr
	#ztr = symop[2][1]*box[2]
	#print ztr
	#X += xtr
        #Y += ytr
        #Z += ztr

def calcrotmat( XS, YS, ZS, XR, YR, ZR, notrans = 0 ):
	R_trans = -array([ sum(XR)/len(XR), sum(YR)/len(YR), sum(ZR)/len(ZR) ], float64)
	S_trans = -array([ sum(XS)/len(XS), sum(YS)/len(YS), sum(ZS)/len(ZS) ], float64) 
	if not notrans:
		translate( XR, YR, ZR, R_trans )
		translate( XS, YS, ZS, S_trans )
	# Rij = sum( Rni * Snj )
	R = array( [[ dot(XR,XS), dot(XR,YS), dot(XR,ZS)],
			[ dot(YR,XS), dot(YR,YS), dot(YR,ZS)],
			[ dot(ZR,XS), dot(ZR,YS), dot(ZR,ZS)] ], float64 )
	RTR = dot( transpose(R), R )
	evalues, evectors = eig(RTR)
	evalues, evectors = esort(evalues, evectors)
	# force right-handed system
	evectors[2] = crossp( evectors[0], evectors[1] )
	b = transpose(dot( R, transpose(evectors) ))
	b[0] =  b[0] / sqrt(dot(b[0], b[0]))
	b[1] =  b[1] / sqrt(dot(b[1], b[1]))
	b[2] = crossp( b[0], b[1] )
	U = dot( transpose(b), evectors )
	return U

#class PDBFile(object):
class PDBFile:
	"""Represents a pdb file, gives header info, allows random
	access to coordinate records. Aimed at I/O rather than
	data manipulation"""
	def __init__( self, filename ):
		if filename[-2:] == 'gz':
			self.file = gzip.GzipFile( filename )
		else:
			self.file = open( filename )
		self.olist = []		
		self.title = ""
		self.expdata = ""
		self.seqmap = {}
		self.nmodel = 0 
		self.compnd = ""
		self.temperature = 0.0
		self.resolution = 0.0
		for line in self.file.readlines():
			self.parseline(line)
		if len(self.olist) != self.olist.count(""):
			print "Warning: pdb file ", self.idCode, "withdrawn"
			print "from database, replaced by "
		self.file.seek(0)

	def close(self):
		self.file.close()

	def get_nmodel(self):
		return self.nmodel

	def parseline(self, line):
		rkey = line[0:6]
		if rkey == "HEADER":
			self.depDate = line[50:59]
			self.classification = line[10:50]
			self.idCode = line[62:66]
		elif rkey == "OBSLTE":
			olist.append(line[31:35]).append(line[31:35])
			olist.append(line[36:40]).append(line[41:45])
			olist.append(line[46:50]).append(line[51:55])
			olist.append(line[56:60]).append(line[61:65])
			olist.append(line[66:70])
		elif rkey == "TITLE ":
			self.title += line[10:70].strip() + ' '
		elif rkey == "COMPND":
			self.compnd += line[10:70].strip() + ' '
		elif rkey == "EXPDTA":
			self.expdata += line[10:70].strip() + ' '
#		elif line[0:6] == "AUTHOR":
#			self.expdata += line[10:70].strip() + ' '
		elif rkey == "SEQRES":
			chainid = line[11]
			if chainid in self.seqmap.keys():
				self.seqmap[chainid] += line[19:70].split()
			else:
				self.seqmap[chainid] = line[19:70].split()
#		elif rkey == "HELIX ":
#		elif rkey == "SEQRES":
		elif rkey == "MODEL ":
			self.nmodel+=1
		elif rkey == "REMARK":
			if line[7:10] == '200' and line.find('TEMPERATURE') > 0:
					ls = line.split()
					T = ls[-1]
					if T != "NULL" and T[0].isdigit():
						self.temperature = float(T)
			if line[9] == '2' and line.find('RESOLUTION') > 0:
					ls = line.split()
					R = ls[3]
					if R != "NULL" and R[0].isdigit():
						self.resolution = float(R)

	def get_title(self):
		return self.title

	def get_compnd(self):
		return self.compnd

	def get_depdate(self):
		return self.depdate

	def get_resolution(self):
		return self.resolution

	def get_temperature(self):
		return self.temperature

	def get_expdata(self):
		return self.expdata

	def get_chains(self):
		return self.seqmap.keys()

	def get_seq(self, chainid):
		return self.seqmap[chainid]

	def find_subseq(self, subseq):
		"""Finds a single-letter code subsequence, by searching
		all chains. Returns id of first matching chain and offset
		within it"""
		# TODO: consider case where there are no SEQRES records
		if len(self.seqmap.keys()) > 0:
			for key in self.seqmap.keys():
				off = onel_seq(self.seqmap[key]).find(subseq)
				if (off >= 0):
					print "Found subseq in chain " + key + \
					" in SEQRES, looking in coordinates"
					# check coordinate records 
					# actually there!
					offset = \
					self.check_records( subseq, key )
					if offset >= 0:
						print "Found subseq in coords too!"
						return key, offset
					else:
						print "Could not find subseq in coords!"
						return "NULL", -1
		else:
			offset = self.check_records( subseq, "ignore" )
			if offset >= 0:
				return "ignore", offset
			else:
				return "NULL", -1

			
	def chop_model(self, model, lines):
		for line in lines:
			if line[0:6] == "MODEL " and model == int(line[10:14]):
				lines = lines[ lines.index(line)+1 : ]
				break
		for line in lines:
			if line[0:6] == "ENDMDL":
				lines = lines[ : lines.index(line) ]
				break
		return lines

	def check_records(self, subseq, chain, model = 0):
		"""Recursive, incremental search routine to find a subsequence
		in the actual coordinate records"""
		self.file.seek(0)
		lines = self.file.readlines()
		if model > 0:
			lines = self.chop_model(model, lines)
		sequence = ""
		resnum = -1
		if chain != "ignore":
			for x in lines:
				if x[0:6] == "ATOM  " and resnum != x[22:26] \
				and chain == x[21]:
					resnum = x[22:26]    
					sequence += onel(x[17:20])
		else: 
			for x in lines:
				if x[0:6] == "ATOM  " and resnum != x[22:26]:
					resnum = x[22:26]    
					sequence += onel(x[17:20])
		return sequence.find(subseq)

	def find_chain(self,chain):
		"""Is a certain chain actually in the ATOM records??!"""
		for line in self.file.readlines():
			if line[0:6] == "ATOM  ":
				c = line[21]
				if c == chain:
					return 1
		return 0

	def res_offset(self, chain, resn, model=0):
		"""What is the offset of residue number resn from
		the start of chain 'chain'?"""
		self.file.seek(0)
		lines = self.file.readlines()
		if model > 0:
			lines = self.chop_model(model, lines)
		resnum = -1
		offset = found = 0
		for x in lines:
			if x[0:6] == "ATOM  " and resnum != x[22:26] and chain == x[21]:
				if int(x[22:26]) == int(resn):
					found =1 
					break
				resnum = x[22:26]
				offset += 1
		if found:
			return offset
		else:
			return -1

		
#class PDBData(object):
class PDBData:
	"""Class to contain PDB atom records only; does not
	support multiple chains, models, etc!"""
	def __init__(self, pfile, chain='ignore', offset=0, nres=-1, model=0, resrenum=1,allowhet=0):
		"not very OO, but ..."
		if allowhet:
			atomrec = [ "ATOM  ", "HETATM" ]
		else:
			atomrec = [ "ATOM  " ]
		self.chain = chain
		self.nres = nres
		self.seq = ""
		# find right model
		pfile.file.seek(0)
		if model > 0:
			m=-1
			while m != model:
				input = pfile.file.readline()
				if input[0:6] == "MODEL ":
					m = int(input.split()[1]) # more robust
					#m = int(input[10:14])

		# now find right chain
		c = "ignore"

		if chain == "ignore":
			input = pfile.file.readline()
			while input[0:6] not in atomrec: #!= "ATOM  ":
				input = pfile.file.readline()
		else:
			while c != chain:
				input = pfile.file.readline()
				if input[0:6] in atomrec: # "ATOM  ":
					c = input[21]
		# and correct residue offset
		resseq = int(input[22:26])
		resnum=0
		while offset != resnum:
			input = pfile.file.readline()
			cres = int(input[22:26])
			if cres != resseq:
				resnum+=1       #(cres -  resseq)
				resseq = cres
		# rewind to start of line
		pfile.file.seek(pfile.file.tell() - len(input))
		# now read ze coordinates
		# for now, use fixed-size arrays; can change later
		self.X = zeros(20000, float64)
		self.Y = zeros(20000, float64)
		self.Z = zeros(20000, float64)
		self.wght = zeros(20000, float64)
		self.natom = 0
		self.atom = []
		self.resname = []
		self.chain = []
		self.resnum = []
		self.segid = []
		self.arec = []
		curres = 0
		curresseq = -1
		for x in pfile.file.readlines():
			if x[0:6] == "ANISOU":
				continue
			if x[0:6] == "TER   ":
				continue
			elif x[0:6] not in atomrec: #!= "ATOM  ":
				break

			if x[21] != chain and chain != 'ignore':
				break
			if curresseq != x[22:26]:
				oldresseq = curresseq
				curresseq = x[22:26]
				#if curres != 0:
				#	diff = int(curresseq) - int(oldresseq)
				#else:
				#	diff = 1
				diff=1
				curres+=diff
				if curres > nres and nres > 0:
					break
				else:
					if diff > 1:	## residues missing in pdb
						dummy = ""
						for i in range(diff-1):
							dummy += "X"
						self.seq += dummy
					self.seq += onel( x[17:20] ) 
			# skip lower of two partial occupancies
			occ = float(x[54:60])
			if model == 0:		## only consider occup. for x-ray structures
				if occ < 0.5 and occ != 0.0:
					continue
				elif occ == 0.5 and x[16] != 'A':
					continue 
			self.arec.append(x[0:6])
			self.atom.append(x[12:16])
			self.resname.append(x[17:20])
			self.chain.append(x[21])
			if resrenum:
				self.resnum.append(curres)
			else:
				self.resnum.append(int(x[22:26]))
			self.segid.append(x[72:76])
			self.X[self.natom] = float(x[30:38])
			self.Y[self.natom] = float(x[38:46])
			self.Z[self.natom] = float(x[46:54])
			self.wght[self.natom] = float(x[60:66])
			self.natom+=1
		if self.nres < 0:		# catch "select all" case
			self.nres = curres
	
	def zero(self):
		for i in range(self.natom):
			self.X[i] = 0.0
			self.Y[i] = 0.0
			self.Z[i] = 0.0
			self.wght[i] = 0.0
	
	def zero_weight(self):
		"""Zeros the entire weight array"""
		for i in range(self.natom):
			self.wght[i] = 0.0

	def assign_weight(self, val):
		"""Sets the entire weight array to val"""
		for i in range(self.natom):
			self.wght[i] = val

	def set_weight(self, wmap):
		"""Sets the weight array using the map wmap of the
		form weight[(residue,atomtype)]; residues not specified will be
		left unchanged"""
		wkeys = wmap.keys();
		print wkeys
		for i in range(self.natom):
			k = (str(self.resnum[i]),self.atom[i])
			if k in wkeys:
				self.wght[i] = wmap[k]

	def rms_from_bfact(self):
		rms = map(lambda x: math.sqrt(3.0*x/(8.0*(math.pi)**2)), self.wght)
		return rms

	def get_index(self,resn,at):
		for i in range(self.natom):
			if self.resnum[i] == resn and self.atom[i].strip() == at:
				return i
		#for i in range(self.natom):
		#	print i, self.resnum[i], self.atom[i]
		return -1

	def get_degenerate_indices(self,resn,at):
		indices = []
		if at.find('#') >= 0:
			at=at[:at.find('#')]
			for i in range(self.natom):
				if self.resnum[i] == resn and self.atom[i].find(at)>=0:
					indices.append(i)
		else:
			for i in range(self.natom):
				if self.resnum[i] == resn and self.atom[i].strip()==at:
					indices.append(i)

		#for i in range(self.natom):
		#	print i, self.resnum[i], self.atom[i]
		return indices

	def get_weight(self, key):
		for i in range(self.natom):
			k = (str(self.resnum[i]),self.atom[i])
			if k == key:
				return self.wght[i]
		return 0

	def angle(self,i,j,k):
		xji = self.X[j]-self.X[i]
		yji = self.Y[j]-self.Y[i]
		zji = self.Z[j]-self.Z[i]
		xjk = self.X[j]-self.X[k]
		yjk = self.Y[j]-self.Y[k]
		zjk = self.Z[j]-self.Z[k]
		rji = math.sqrt(xji**2+yji**2+zji**2)
		rjk = math.sqrt(xjk**2+yjk**2+zjk**2)
		xji /= rji; yji /= rji; zji /= rji; 
		xjk /= rjk; yjk /= rjk; zjk /= rjk; 
		dotp = xji*xjk + yji*yjk + zji*zjk;
		crossp = math.sqrt(1.0 - dotp**2)
		acos_angle = math.acos(dotp)
		asin_angle = math.asin(crossp)
		if asin_angle >= 0.0:
			angle = 180.0/math.pi*acos_angle
		else:
			angle = 360.0-180.0/math.pi*acos_angle
		return angle

	def out_of_plane_angle(self, i,j,k,l):
		"return angle by which jl is out of plane ijk"
		xij = self.X[i]-self.X[j]
		yij = self.Y[i]-self.Y[j]
		zij = self.Z[i]-self.Z[j]
		xjk = self.X[j]-self.X[k]
		yjk = self.Y[j]-self.Y[k]
		zjk = self.Z[j]-self.Z[k]
		xjl = self.X[j]-self.X[l]
		yjl = self.Y[j]-self.Y[l]
		zjl = self.Z[j]-self.Z[l]
		ax = yij*zjk-zij*yjk
		ay = zij*xjk-xij*zjk
		az = xij*yjk-yij*xjk
		#bx = ylk*zjk-zlk*yjk
		#by = zlk*xjk-xlk*zjk
		#bz = xlk*yjk-ylk*xjk
		ra2 = ax**2 + ay**2 + az**2
		rjl2 = xjl**2 + yjl**2 + zjl**2
		#ra = math.sqrt(ra2)
		#ax/=ra
		#ay/=ra
		#az/=ra
		cosphi = ( ax*xjl + ay*yjl + az*zjl ) / math.sqrt(ra2*rjl2)
		phi = math.acos(cosphi)
		return phi
		#rb2 = bx**2 + by**2 + bz**2
		#rjk2 = xjk**2 + yjk**2 + zjk**2
		#rjk = math.sqrt(rjk2)
		#rjkr = 1.0/rjk
		#ra2r = 1.0/ra2
		#rb2r = 1.0/rb2
		#rabr = math.sqrt(ra2r*rb2r)
		#cosphi = ( ax*bx + ay*by + az*bz ) * rabr
		#print 'cosphi = ', cosphi
		#sinphi = rjk*rabr*(ax*xlk + ay*ylk + az*zlk)
		#print 'sinphi = ', sinphi
		#acos = math.acos(cosphi)
		#asin = math.asin(sinphi)
		##print cosphi,sinphi,acos, asin
		#if acos < math.pi/2.0:
		#	return asin/math.pi*180.0
		#elif asin >= 0:
		#	return acos/math.pi*180.0	#(math.pi-asin)/math.pi*180.0
		#else:
		#	return -acos/math.pi*180.0	#(asin-math.pi)/math.pi*180.0

	def torsion(self, i,j,k,l):
		"return torsion angle between specified atoms"
		xij = self.X[i]-self.X[j]
		yij = self.Y[i]-self.Y[j]
		zij = self.Z[i]-self.Z[j]
		xjk = self.X[j]-self.X[k]
		yjk = self.Y[j]-self.Y[k]
		zjk = self.Z[j]-self.Z[k]
		xlk = self.X[l]-self.X[k]
		ylk = self.Y[l]-self.Y[k]
		zlk = self.Z[l]-self.Z[k]
		ax = yij*zjk-zij*yjk
		ay = zij*xjk-xij*zjk
		az = xij*yjk-yij*xjk
		bx = ylk*zjk-zlk*yjk
		by = zlk*xjk-xlk*zjk
		bz = xlk*yjk-ylk*xjk
		ra2 = ax**2 + ay**2 + az**2
		rb2 = bx**2 + by**2 + bz**2
		rjk2 = xjk**2 + yjk**2 + zjk**2
		rjk = math.sqrt(rjk2)
		rjkr = 1.0/rjk
		ra2r = 1.0/ra2
		rb2r = 1.0/rb2
		rabr = math.sqrt(ra2r*rb2r)
		cosphi = ( ax*bx + ay*by + az*bz ) * rabr
		#print 'cosphi = ', cosphi
		sinphi = rjk*rabr*(ax*xlk + ay*ylk + az*zlk)
		#print 'sinphi = ', sinphi
		acos = math.acos(cosphi)
		asin = math.asin(sinphi)
		#print cosphi,sinphi,acos, asin
		if acos < math.pi/2.0:
			return asin/math.pi*180.0
		elif asin >= 0:
			return acos/math.pi*180.0	#(math.pi-asin)/math.pi*180.0
		else:
			return -acos/math.pi*180.0	#(asin-math.pi)/math.pi*180.0

	def get_seq(self):
		return self.seq

	def get_wght(self, i):
		return self.wght[i]

	def set_wght(self, i, wght):
		self.wght[i] = float(wght)	# to be on the safe side

	def average(self, list):
		"""Makes self the average of structures in list,
		puts rms from avg in weight"""
		self.zero()
		nmodel = len(list)
		for x in list:
			if x.natom != self.natom:
				print "wrong number of atoms!!!"
				sys.exit(1)
			for i in range(self.natom):
				self.X[i] += x.X[i]
				self.Y[i] += x.Y[i]
				self.Z[i] += x.Z[i]
		for i in range(self.natom):
			self.X[i] /= nmodel
			self.Y[i] /= nmodel
			self.Z[i] /= nmodel
		for x in list:
			for i in range(self.natom):
				dx = x.X[i] - self.X[i]
				dy = x.Y[i] - self.Y[i]
				dz = x.Z[i] - self.Z[i]
				self.wght[i] += dx*dx + dy*dy + dz*dz
		for i in range(self.natom):
			self.wght[i] = sqrt(self.wght[i]/nmodel)

	def set_segid(self, s):
		for a in range(self.natom):
			self.segid[a] = s
		
	def map_segid(self, cur_s, new_s):
		for a in range(self.natom):
			if self.segid[a] == cur_s:
				self.segid[a] = new_s
	
	def charmm19ify(self):
		for a in range(self.natom):
			if self.resname[a] == 'ILE' and self.atom[a] == ' CD1':
				self.atom[a] = ' CD '
			if self.resname[a] == 'GLN' and self.atom[a] == '1HE2':
				self.atom[a] = 'HE21'
			if self.resname[a] == 'GLN' and self.atom[a] == '2HE2':
				self.atom[a] = 'HE22'
			if self.atom[a] == ' OXT':
				self.atom[a] = ' OT1'

	def charmm22ify(self):
		for a in range(self.natom):
			if self.resname[a] == 'ILE' and self.atom[a] == ' CD1':
				self.atom[a] = ' CD '
			if self.resname[a] == 'GLN' and self.atom[a] == '1HE2':
				self.atom[a] = 'HE21'
			if self.resname[a] == 'GLN' and self.atom[a] == '2HE2':
				self.atom[a] = 'HE22'
			if self.atom[a] == ' OXT':
				self.atom[a] = ' OT1'
			if self.resname[a] == 'HIS':
				self.resname[a] = 'HSD'
		
	def pdbify(self):
		for a in range(self.natom):
			if self.atom[a] == ' HN ':
				self.atom[a] = '  H '
			if self.resname[a] == 'ILE' and self.atom[a] == ' CD ':
				self.atom[a] = ' CD1'
			if self.resname[a] == 'GLN' and self.atom[a] == 'HE21':
				self.atom[a] = '1HE2'
			if self.resname[a] == 'GLN' and self.atom[a] == 'HE22':
				self.atom[a] = '2HE2'
			if self.atom[a] == ' OT1':
				self.atom[a] = ' OXT'
			if self.resname[a] == 'HSD':
				self.resname[a] = 'HIS'
		
	def mutate(self, res, target):
		"n.v. intelligent, but..."
		for a in range(self.natom):
			if self.resnum[a] == res:
				self.resname[a] = target #aamap[target]
		self.seq = self.seq[:res-1] + aalist[target] + self.seq[res:]

	def dist(self, a, b):
		"distance between atoms a and b"
		dx = self.X[a]-self.X[b]
		dy = self.Y[a]-self.Y[b]
		dz = self.Z[a]-self.Z[b]
		return sqrt( dx*dx + dy*dy + dz*dz )

	def selcontacts(self, sel, cutoff=6.0, excl=2, type="all"):
		"""calculate heavy atom contacts between a list of atoms
		sel and all other atoms (not in sel)"""
		contacts = 0;
		for i in range(self.natom):
			if ('H' not in self.atom[i]) and (i not in sel):
				for j in sel:
					if self.dist(i,j) < cutoff \
					and math.sqrt((self.resnum[i]-self.resnum[j])**2)\
					>= excl:
						contacts +=1
		return contacts

	def sin2contacts(self, p, q, cutoff=6.0, excl=2, type="all"):
		contacts = 0.0;
		xpq = self.X[p]-self.X[q]
		ypq = self.Y[p]-self.Y[q]
		zpq = self.Z[p]-self.Z[q]
		rpq = math.sqrt(xpq**2 + ypq**2 + zpq**2)
		xpq /= rpq
		ypq /= rpq
		zpq /= rpq
		for i in range(self.natom):
			if ('H' not in self.atom[i]) and i!=p and i!=q:
				if self.dist(i,p) < cutoff \
				and math.sqrt((self.resnum[i]-self.resnum[p])**2)\
				>= excl:
					xpi = self.X[p] - self.X[i]
					ypi = self.Y[p] - self.Y[i]
					zpi = self.Z[p] - self.Z[i]
					rpi = math.sqrt(xpi**2 + ypi**2 + zpi**2)
					xpi/=rpi
					ypi/=rpi
					zpi/=rpi
					dotp = xpi*xpq + ypi*ypq + zpi*zpq
					contacts += 1.0 - dotp**2
		return contacts

	def contacts(self, i, j, cutoff=6.0, type="all"):
		"""calculate all heavy atom contacts between residues i and j
		within a certain cutoff"""
		if i == j:
			print "can't calculate contacts within same residue"
			import sys
			sys.exit(1)
		bb = [ ' CA ', ' C  ', ' N  ', ' O  ' ]
		ilist, jlist = [], []
		clist = []
		for a in range(self.natom):
			if self.resnum[a] == i and self.atom[a][1] != 'H':
				if self.atom[a] in bb:
					if type == "bb" or type == "all":
						ilist.append(a)
				else:
					if type == "sc" or type == "all":
						ilist.append(a)
			elif self.resnum[a] == j and self.atom[a][1] != 'H':
				if self.atom[a] in bb:
					if type == "bb" or type == "all":
						jlist.append(a)
				else:
					if type == "sc" or type == "all":
						jlist.append(a)
		for ii in ilist:
			for jj in jlist:
				if ii in bb and jj in backbone and abs(ii-jj)==1:
					continue
				d = self.dist(ii,jj)
				if d < cutoff:
					clist.append((self.atom[ii],self.resname[ii],
					self.resnum[ii],self.atom[jj],self.resname[jj],
					self.resnum[jj], d))
		return clist

	def contacts2(self, i, cutoff=6.0, type="all"):
		"""calculate all heavy atom contacts between residue i and other
		atoms within a certain cutoff"""
		bb = [ ' CA ', ' C  ', ' N  ', ' O  ' ]
		ilist, jlist = [], []
		clist = []
		for a in range(self.natom):
			if self.resnum[a] == i and self.atom[a][1] != 'H':
				if self.atom[a] in bb:
					if type == "bb" or type == "all":
						ilist.append(a)
				else:
					if type == "sc" or type == "all":
						ilist.append(a)
			elif self.atom[a][1] != 'H':
				if self.atom[a] in bb:
					if type == "bb" or type == "all":
						jlist.append(a)
				else:
					if type == "sc" or type == "all":
						jlist.append(a)
		for ii in ilist:
			for jj in jlist:
				if ii in bb and jj in backbone and abs(ii-jj)==1:
					continue
				d = self.dist(ii,jj)
				if d < cutoff:
					clist.append((self.atom[ii],self.resname[ii],
					self.resnum[ii],self.atom[jj],self.resname[jj],
					self.resnum[jj], d))
		return clist

		
	def write(self):
		print "writing"

	def align(self, ref, slist, rlist, type="BB"):
		"""Align self to ref (rotate & translate) using residue selections
		slist and rlist for self and ref respectively. Only selected type
		of backbone atoms used for fitting. Type = "CA" or "BB".
		Works for any proteins (i.e. don't need matching residues)"""
		if type == "CA":
			stypes = [ ' CA ' ]
		elif type == "C1'":
			stypes = [ " C1'" ]
		elif type == "BB":
			stypes = [ ' CA ', ' C  ', ' N  ', ' O  ' ]
#		if type != "CA":
#			print "Only CA selection supported at present!"
#			import sys
#			sys.exit(1)
		if len(slist) != len(rlist):
			print "RESIDUE SELECTIONS OF DIFFERENT LENGTH, QUITTING!"
			import sys
			sys.exit(1)
		lsel = len(slist)
		ltype = len(stypes)
		latom = ltype * lsel
		XR = zeros(latom, float64)
		YR = zeros(latom, float64)
		ZR = zeros(latom, float64)
		XS = zeros(latom, float64)
		YS = zeros(latom, float64)
		ZS = zeros(latom, float64)
		# fill arrays with CA coordinates
		for i in range(lsel):
			for a in range(ltype):
				idx = i*ltype + a
				XS[idx],YS[idx],ZS[idx],u=self.get_crd(stypes[a],slist[i])
				XR[idx],YR[idx],ZR[idx],u=ref.get_crd(stypes[a],rlist[i])
				#print XS[idx],YS[idx],ZS[idx]
				#print XR[idx],YR[idx],ZR[idx]
			# coords missing for some reason ...
			if XR[idx] < -999 or XS[idx] < -999:
				print "GOT HERE"
				print stypes[a],rlist[i]
				rlist = rlist[:i] + rlist[i+1:]
				slist = slist[:i] + slist[i+1:]
				slist, rlist = self.align(ref, slist, rlist, type)
				return slist, rlist
		# now calculate translation and rotation operations
		#print sum(XS), len(XS)
		#print sum(XR), len(XR)
		S_trans = -array([ sum(XS)/len(XS), sum(YS)/len(YS), sum(ZS)/len(ZS) ], float64) 
		R_trans = -array([ sum(XR)/len(XR), sum(YR)/len(YR), sum(ZR)/len(ZR) ], float64)
		U = calcrotmat( XS, YS, ZS, XR, YR, ZR )
		# (I) translate self to centre of geometry
		translate( self.X, self.Y, self.Z, S_trans )
		# (II) rotate self
		self.rotate( U )
		# (III) translate self to centre of geometry of reference
		self.translate( -R_trans )
		return slist, rlist

	def homoalign(self, ref, slist, rlist,twofold=0):
		"""Align two structures of a protein using all available heavy
		atom coordinates. Can probably be used with any two proteins,
		but only makes sense for comparing identical sequences"""
		# get heavy atom lists from each protein
		atomlist = []
		satoms = []
		ratoms = []
		for i in range(self.natom):
			if self.resnum[i] in slist and self.atom[i].strip()[0] != 'H':
				atomlist.append( ( self.resnum[i], self.atom[i] ) )
				satoms.append(i)
		
		found = -1
		for j in atomlist:
			for x in range(ref.natom):
				if (ref.resnum[x], ref.atom[x]) == j:
					found = x
					break
			if found >= 0:
				ratoms.append(found)
			else:
				atomlist = atomlist[:j] + atomlist[j+1:]
				satoms = satoms[:j] + satoms[j+1:]

		if len(ratoms) != len(satoms):
			print ratoms
			print satoms
			print "woops!"
			import sys
			sys.exit(1)

		#lsel = len(slist)
		#ltype = len(stypes)
		latom = len(atomlist)
		XR = zeros(latom, float64)
		YR = zeros(latom, float64)
		ZR = zeros(latom, float64)
		XS = zeros(latom, float64)
		YS = zeros(latom, float64)
		ZS = zeros(latom, float64)
		# fill arrays with CA coordinates
		for i in range(latom):
			si, ri = satoms[i], ratoms[i]
			XS[i],YS[i],ZS[i]=self.X[si],self.Y[si],self.Z[si]
			XR[i],YR[i],ZR[i]=ref.X[ri],ref.Y[ri],ref.Z[ri]
		# now calculate translation and rotation operations
		S_trans = -array([ sum(XS)/len(XS), sum(YS)/len(YS), sum(ZS)/len(ZS) ], float64) 
		R_trans = -array([ sum(XR)/len(XR), sum(YR)/len(YR), sum(ZR)/len(ZR) ], float64)
		U = calcrotmat( XS, YS, ZS, XR, YR, ZR )
		# (I) translate self to centre of geometry
		translate( self.X, self.Y, self.Z, S_trans )
		# (II) rotate self
		if twofold == 1:
			nx = math.sqrt(0.5*(U[0][0]+1.0))
			ny = math.sqrt(0.5*(U[1][1]+1.0))
			nz = math.sqrt(0.5*(U[2][2]+1.0))
			print nx,ny,nz, math.sqrt(nx**2+ny**2+nz**2)
			if nx < ny:
				if nz < nx:
					nz = math.sqrt(1.0-nx**2-ny**2)
				else:
					nx = math.sqrt(1.0-ny**2-nz**2)
			else:
				if nz < ny:
					nz = math.sqrt(1.0-nx**2-ny**2)
				else:
					ny = math.sqrt(1.0-nx**2-nz**2)
			print nx,ny,nz, math.sqrt(nx**2+ny**2+nz**2)
			U[0][0] = 2.0*nx**2-1.0
			U[1][1] = 2.0*ny**2-1.0
			U[2][2] = 2.0*nz**2-1.0
			U[0][0] = 2.0*nx**2-1.0
			U[1][1] = 2.0*ny**2-1.0
			U[2][2] = 2.0*nz**2-1.0
			U[0][1] = 2.0*nx*ny
			U[0][2] = 2.0*nx*nz
			U[1][0] = 2.0*nx*ny
			U[1][2] = 2.0*ny*nz
			U[2][0] = 2.0*nx*nz
			U[2][1] = 2.0*ny*nz
			print U
			print dot(U,U)
		self.rotate( U )
		# (III) translate self to centre of geometry of reference
		self.translate( -R_trans )
		rms = 0.0
		for i in range(latom):
			si, ri = satoms[i], ratoms[i]
			dx = self.X[si] - ref.X[ri]
			dy = self.Y[si] - ref.Y[ri]
			dz = self.Z[si] - ref.Z[ri]
			rms += dx**2 + dy**2 + dz**2
		return math.sqrt(rms/float(latom))

	def symmetrize(self, ref, slist, rlist,twofold=0):
		"""Align two structures of a protein using all available heavy
		atom coordinates. Can probably be used with any two proteins,
		but only makes sense for comparing identical sequences"""
		# get heavy atom lists from each protein
		atomlist = []
		satoms = []
		ratoms = []
		for i in range(self.natom):
			if self.resnum[i] in slist and self.atom[i].strip()[0] != 'H':
				atomlist.append( ( self.resnum[i], self.atom[i] ) )
				satoms.append(i)
		
		found = -1
		for j in atomlist:
			for x in range(ref.natom):
				if (ref.resnum[x], ref.atom[x]) == j:
					found = x
					break
			if found >= 0:
				ratoms.append(found)
			else:
				atomlist = atomlist[:j] + atomlist[j+1:]
				satoms = satoms[:j] + satoms[j+1:]

		if len(ratoms) != len(satoms):
			print ratoms
			print satoms
			print "woops!"
			import sys
			sys.exit(1)

		#lsel = len(slist)
		#ltype = len(stypes)
		latom = len(atomlist)
		XR = zeros(latom, float64)
		YR = zeros(latom, float64)
		ZR = zeros(latom, float64)
		XS = zeros(latom, float64)
		YS = zeros(latom, float64)
		ZS = zeros(latom, float64)
		# fill arrays with CA coordinates
		for i in range(latom):
			si, ri = satoms[i], ratoms[i]
			XS[i],YS[i],ZS[i]=self.X[si],self.Y[si],self.Z[si]
			XR[i],YR[i],ZR[i]=ref.X[ri],ref.Y[ri],ref.Z[ri]
		# now calculate translation and rotation operations
		S_trans = -array([ sum(XS)/len(XS), sum(YS)/len(YS), sum(ZS)/len(ZS) ], float64) 
		R_trans = -array([ sum(XR)/len(XR), sum(YR)/len(YR), sum(ZR)/len(ZR) ], float64)
		print S_trans
		print R_trans
		trans_ave = (S_trans+R_trans)/2.0
		print trans_ave
		trans_ave = -array([(sum(XS)+sum(XR))/(len(XS)+len(XR)), \
				(sum(YS)+sum(YR))/(len(YS)+len(YR)), \
				(sum(ZS)+sum(ZR))/(len(ZS)+len(ZR)) ], float64)
		print trans_ave
		#self.translate( -R_trans )
		#self.translate( -S_trans )
		translate( XS, YS, ZS, trans_ave )
		translate( XR, YR, ZR, trans_ave )
		U = calcrotmat( XS, YS, ZS, XR, YR, ZR, 1 )
		print U
		# (II) rotate self
		if twofold == 1:
			nx = math.sqrt(0.5*(U[0][0]+1.0))
			ny = math.sqrt(0.5*(U[1][1]+1.0))
			nz = math.sqrt(0.5*(U[2][2]+1.0))
			print nx,ny,nz, math.sqrt(nx**2+ny**2+nz**2)
			if nx < ny:
				if nz < nx:
					nz = math.sqrt(1.0-nx**2-ny**2)
				else:
					nx = math.sqrt(1.0-ny**2-nz**2)
			else:
				if nz < ny:
					nz = math.sqrt(1.0-nx**2-ny**2)
				else:
					ny = math.sqrt(1.0-nx**2-nz**2)
			print nx,ny,nz, math.sqrt(nx**2+ny**2+nz**2)
			U[0][0] = 2.0*nx**2-1.0
			U[1][1] = 2.0*ny**2-1.0
			U[2][2] = 2.0*nz**2-1.0
			U[0][0] = 2.0*nx**2-1.0
			U[1][1] = 2.0*ny**2-1.0
			U[2][2] = 2.0*nz**2-1.0
			U[0][1] = 2.0*nx*ny
			U[0][2] = 2.0*nx*nz
			U[1][0] = 2.0*nx*ny
			U[1][2] = 2.0*ny*nz
			U[2][0] = 2.0*nx*nz
			U[2][1] = 2.0*ny*nz
			print U
			print dot(U,U)
		for i in range(len(self.X)):
			ref.X[i] = self.X[i]
			ref.Y[i] = self.Y[i]
			ref.Z[i] = self.Z[i]
		self.translate( trans_ave )
		ref.translate( trans_ave )
		self.rotate( U )
		print len(self.X), len(ref.X)

	def translate( self, transvec ):
		translate( self.X, self.Y, self.Z, transvec )

	def transform( self, symop, box ):
		transform( self.X, self.Y, self.Z, symop, box )

	def get_nres(self):
		return self.nres

	def get_atom(self, i):
		return self.atom[i]

	def get_resname(self, i):
		return self.resname[i]

	def get_reslett(self, i):
		return self.seq[i]

	def get_resnum(self, i):
		return self.resnum[i]

	def get_natom(self):
		return self.natom

	def calc_bb_rmsd( self, ref, self_align, ref_align, rs=0, rr=0 ):
		"""Calculates rmsd between corresponding backbone
		CA, C, O, N residues for two structures"""
		natom = ref.get_natom()
		nres = ref.get_nres()
		rmsarray = zeros( natom, float64 )
		countarray = zeros( natom, float64 )
		# need for calculating displacement correlations
		ms,mr = [],[]
		#rs=rr=0
		for x in range(len(self_align)):
			if letters.find(self_align[x]) >= 0:
				rs+=1
			if letters.find(ref_align[x]) >= 0:
				rr+=1
			if letters.find(self_align[x]) >= 0 and \
			letters.find(ref_align[x]) >= 0:
				ms.append(rs)
				mr.append(rr)
		for i in range( len(mr) ):
			for at in [ ' CA ', ' C  ', ' N  ', ' O  ' ]:
				xr, yr, zr, atr = ref.get_crd( at, mr[i] )
				xs, ys, zs, ats = self.get_crd( at, ms[i] )
				if xr < -999 or xs < -999:
					continue
				dx,dy,dz = xr-xs, yr-ys, zr-zs
				rms2 = dx*dx + dy*dy + dz*dz
				rmsarray[atr], countarray[atr] = rms2, 1
	        return rmsarray, countarray

	def inter_contacts( self, other, cut=4.0 ):
		"""Calculates heavy atom contacts between two structures"""
		ns = self.get_natom()
		no = other.get_natom()
		contact_list = []
		dr_min  = 1000.
		for i in range(ns):
			if self.atom[i].strip()[0] == "H":
				continue
			for j in range(no):
				if other.atom[j].strip()[0] == "H":
					continue
				dx = self.X[i]-other.X[j]
				dy = self.Y[i]-other.Y[j]
				dz = self.Z[i]-other.Z[j]
				dr = (dx**2+dy**2+dz**2)**0.5
				if dr_min > dr:
					dr_min = dr
				if dr < cut:
					contact_list.append((self.resname[i],self.resnum[i],
						self.atom[i],other.resname[j],other.resnum[j],
						other.atom[j],dr))

		#print dr_min
		return contact_list

	def calc_tot_bb_rmsd( self, ref, self_align, ref_align, rs=0, rr=0 ):
		"""Calculates rmsd between corresponding backbone
		CA, C, O, N residues for two structures"""
		natom = ref.get_natom()
		nres = ref.get_nres()
		# need for calculating displacement correlations
		ms,mr = [],[]
		#rs=rr=0
		for x in range(len(self_align)):
			if letters.find(self_align[x]) >= 0:
				rs+=1
			if letters.find(ref_align[x]) >= 0:
				rr+=1
			if letters.find(self_align[x]) >= 0 and \
			letters.find(ref_align[x]) >= 0:
				ms.append(rs)
				mr.append(rr)
		rms2 = 0.0
		nmatch = 0
		for i in range( len(mr) ):
			for at in [ ' CA ', ' C  ', ' N  ', ' O  ' ]:
				xr, yr, zr, atr = ref.get_crd( at, mr[i] )
				xs, ys, zs, ats = self.get_crd( at, ms[i] )
				if xr < -999 or xs < -999:
					continue
				dx,dy,dz = xr-xs, yr-ys, zr-zs
				rms2 += dx*dx + dy*dy + dz*dz
				nmatch += 1

		return math.sqrt(rms2/float(nmatch))

	def calc_ave_disp( self, ref, self_align, ref_align, rs=0, rr=0 ):
		"""Calculates average displacements between self and ref
		structures, using an average over all backbone atoms for
		each residue."""
		natom = ref.get_natom()
		nres = ref.get_nres()
		dr = zeros( [nres,3], float64 )
		dr_count = zeros( nres, float64 )
		ms,mr = [],[]
		#rs=rr=0
		for x in range(len(self_align)):
			if letters.find(self_align[x]) >= 0:
				rs+=1
			if letters.find(ref_align[x]) >= 0:
				rr+=1
			if letters.find(self_align[x]) >= 0 and \
			letters.find(ref_align[x]) >= 0:
				ms.append(rs)
				mr.append(rr)
		for i in range( len(mr) ):
			hits = 0
			r = mr[i] - 1
			for at in [ ' CA ', ' C  ', ' N  ', ' O  ' ]:
			#for at in [ ' CA ' ]:
				xr, yr, zr, atr = ref.get_crd( at, mr[i] )
				xs, ys, zs, ats = self.get_crd( at, ms[i] )
				if xr < -999 or xs < -999:
					continue
				dx,dy,dz = xr-xs, yr-ys, zr-zs
				dr[r][0] += dx
				dr[r][1] += dy
				dr[r][2] += dz
				hits += 1
			if hits > 0:
				dr[r][0] /= hits
				dr[r][1] /= hits
				dr[r][2] /= hits
				dr_count[r] = 1
		return dr, dr_count

	def get_crd(self, at, i):
		for x in range(self.natom):
			if self.atom[x] == at and self.resnum[x] == i:
				#print self.resnum[x], self.resname[x], self.atom[x]
				return [self.X[x], self.Y[x], self.Z[x], x]
		# No atom found?
		return -1000,-1000,-1000, -1000

	def get_x(self, i):
		return self.X[i];

	def get_y(self, i):
		return self.Y[i];

	def get_z(self, i):
		return self.Z[i];

	def get_CA_crd(self,i):
		for x in range(self.natom):
			if self.atom[x] == ' CA ' and self.resnum[x] == i:
				return [self.X[x], self.Y[x], self.Z[x]]
		# No CA found?
		return -1000,-1000,-1000

	def rotate(self, A):
		self.X, self.Y, self.Z = dot( A, array([self.X,self.Y,self.Z]) )

	def rrotate(self, A):
		for i in range(len(self.X)):
			self.X[i], self.Y[i], self.Z[i] = \
					dot( A, array([self.X[i],self.Y[i],self.Z[i]]))

	def centre(self):
		trans = -array([ sum(self.X)/len(self.X), sum(self.Y)/len(self.Y), sum(self.Z)/len(self.Z) ], float64) 
		self.translate(trans)

	def xrotate(self, phi):
		"""rotate around x axis by angle phi (in radians)"""
		sinphi = math.sin(phi)
		cosphi = math.cos(phi)
		U = array( [[ 1.0, 0.0, 0.0],
			[ 0.0, cosphi, -sinphi],
			[ 0.0, sinphi, cosphi] ], float64 )
		self.rotate(U)

	def yrotate(self, phi):
		"""rotate around y axis by angle phi (in radians)"""
		sinphi = math.sin(phi)
		cosphi = math.cos(phi)
		U = array( [[ cosphi, 0.0, -sinphi],
			[ 0.0, 1.0, 0.0],
			[ sinphi, 0.0, cosphi] ], float64 )
		self.rotate(U)

	def zrotate(self, phi):
		"""rotate around x axis by angle phi (in radians)"""
		sinphi = math.sin(phi)
		cosphi = math.cos(phi)
		U = array( [[ cosphi, -sinphi, 0.0],
			[ sinphi, cosphi, 0.0],
			[ 0.0, 0.0, 1.0] ], float64 )
		self.rotate(U)
	
def crossp( a, b ):
	"cross-product"
	return array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])
	
def ltswap( evals, ev, x, y ):
	if evals[x] < evals[y]:
		evals[x], evals[y] = evals[y], evals[x]
		ev[x,0],ev[x,1],ev[x,2], ev[y,0], ev[y,1], ev[y,2] = \
			ev[y,0],ev[y,1],ev[y,2], ev[x,0], ev[x,1], ev[x,2]

def esort( evals, evectors ):
	"Sort eigenvectors in descending order of eigenvalue"
	ltswap( evals, evectors, 0, 1 )
	ltswap( evals, evectors, 0, 2 )
	ltswap( evals, evectors, 1, 2 )
	return evals, evectors

	
#class PDBOutput(object):
class PDBOutput:
	"""The original PDBOutput class, now just a base for PDBModelOuput
	and PDBChainOutput which should be usind INSTEAD !!!"""
	def __init__(self, filename="junk.pdb"):
		self.ofile = open(filename, 'w')
		self.natom = 0
		self.nres = 0
		self.classification = "A PROTEIN"
		self.date = time.strftime("%d-%b-%y").upper()
		self.idcode = "1ROB"
		self.cryst1 = []
		self.sequence = {}
	
	def close(self):
		self.ofile.write("END\n")
		self.ofile.close()

	def write_header(self, classification="none", date="none"):
		if classification != "none":
			self.classification = classification.upper()
		if date != "none":
			self.date = date.upper()
		self.ofile.write("%-6s    %-40s%9s   %4s\n" 
		% ("HEADER", self.classification, self.date, self.idcode))

	def write_sequence(self):
		for k in sequence.keys():
			tmpstr = sequence[k]
			nr = len(tmpstr)
			lineno = 1
			while len(tmpstr) > 0:
				if len(tmpstr) < 13: 
					self.write_sequence_line(tmpstr)
					tmpstr = ""
				else:
					self.write_sequence_line(k,lineno,tmpstr[0:14],nr)
					tmpstr = tmpstr[14:]
				lineno += 1
		
	def write_sequence_line(self,chain,lineno,line,nres):
		self.ofile.write("%-6s  %2i %1s %4i " % ( "SEQRES", lineno, chain, nres ))
		nwrit = 0
		while len(line) > 0:
			aa = aamap[line[0].upper()]
			self.ofile.write(" %3s" % ( aa ))
			nwrit +=1
		for i in range(13 - nwrit):
			self.ofile.write(" %3s" % ( "   " ))
		self.ofile.write("\n")
		
#       deprecated, but kept for compatibility
	def append(self, pdata, chainid):
		oldnatom = self.natom
		self.natom += pdata.natom
		self.nres += pdata.nres
		for i in range(pdata.natom):
			format = "%-6s%5i %4s %3s %c%4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n"
			self.ofile.write( format
			% (pdata.arec[i], oldnatom + i +1, pdata.atom[i], pdata.resname[i],
			chainid, pdata.resnum[i], pdata.X[i], pdata.Y[i], pdata.Z[i], 
			1.0, pdata.wght[i], pdata.segid[i]))
		self.ofile.flush()
	
	def write(self,filename = "junk.pdb"):
		print "PDBOutput.write() not implemented yet, and "
		print "probably never will be..."


# Since the use of models and chains in PDB files is pretty much 
# mutually exclusive, the following two classes are naturally separate"""

class PDBModelOutput(PDBOutput):
	def __init__(self):
		PDBOutput.init(self, "none")
		self.nmodel = 0
	def append_model(self, pdata):
		if self.nmodel == 0:
			self.sequence[" "] = pdata.get_seq()
		else:
			if self.sequence[" "] != pdata.get_seq():
				print "New model sequence does not match old"
				sys.exit(1)

	def write(self, filename):
		self.write_header()
		print "PDBModelOutput.append_model() not implemented yet"

class PDBChainOutput(PDBOutput):
	def __init__(self):
		PDBOutput.init(self, "none")
	def append_chain(self, pdata):
		print "PDBChainOutput.append_chain() not implemented yet"

if __name__ == '__main__':
	"test case"
	#tpdb = PDBFile('pdb1bpv.ent')
	#tpdb = PDBFile('pdb1fak.ent')
	tpdb = PDBFile('pdb1ten.ent')
	print tpdb.depDate
	print tpdb.classification
	print tpdb.idCode
	print tpdb.title
	print tpdb.expdata
	print tpdb.seqmap.keys()
#	for key in tpdb.seqmap.keys():
#		print "SEQUENCE FOR CHAIN ", key
#		print tpdb.seqmap[key]
	print tpdb.nmodel
	tdata = PDBData(tpdb, ' ', 5, 80, 0)	#, 1)
	allc = tdata.contacts(18,20,6.0,"all")
	bbc = tdata.contacts(18,20,6.0,"bb")
	scc = tdata.contacts(18,20,6.0,"sc")
	print allc
	print bbc
	print scc
	pout = PDBOutput("testpdboutput.pdb")
	pout.append( tdata, 'I' )
	



