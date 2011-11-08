import re, xml, sys, string, os, optparse, math, copy
import xml.sax
import xml.sax.handler
from tempfile import mkstemp
from procrun import proc_run_exitOnError
from pdbr import protein, isPdbAAline, isPdbHydrogenLine, line2atomid, line2crd, line2resid, line2atomname, line2resn
from protinfo import getResidueCategory
from geom import vec_dist, vec_angle_deg, vec_angle_deg3, vec_dihed_deg, vec_add, vec_scale, cross_product, normalize, vec_diff, dot_product
from grid import makeGridZimmer

def pymol_line(p1, p2) :
	return (p1,p2)
	#return "VERTEX, %f,%f,%f, VERTEX, %f,%f,%f," % (p1[0], p1[1], p1[2], p2[0], p2[1], p2[2])
	
class GeomLim : ## geometric checks
	def __init__(self, tags, minimum, maximum) :
		self.tags, self.min, self.max = tags, copy.deepcopy(minimum), copy.deepcopy(maximum)
	def check(self, tag2points) :
		if len(self.tags) == 2 : x = vec_dist(tag2points[self.tags[0]], tag2points[self.tags[1]])
		elif len(self.tags) == 3 : x = vec_angle_deg3(tag2points[self.tags[0]], tag2points[self.tags[1]], tag2points[self.tags[2]])
		elif len(self.tags) == 4 : x = vec_dihed_deg(tag2points[self.tags[0]], tag2points[self.tags[1]], tag2points[self.tags[2]], tag2points[self.tags[3]])
		else : assert(0)
		ny = '?'
		for i in range(len(self.max)) :
			if self.min[i] <= x and x <= self.max[i] : ny = '!'
		print len(self.tags), "GeomCheck", self.min, x, self.max, ny
		if ny == '?' : return None
		else : return 1
	def info(self) :
		x = 'GeomLim'
		for i in range(len(self.min)) : x = x + " (%f %f)" % (self.min[i], self.max[i])
		for t in self.tags : x = x + " " + t
		print x

class TaggedAtom :
	def __init__(self, atomname, tag) :
		self.atn = atomname
		self.tag = tag
	def info(self) :
		print "atom " + self.atn + " " + self.tag

class TaggedGroup :
	def __init__(self, restype, tag, atoms) :
		assert tag == 'A' or tag == 'D'
		self.restype = restype
		self.ADtag = tag
		self.atoms = atoms
	def info(self) :
		print "group " + self.restype + " " + self.ADtag
		for a in self.atoms : a.info()

class Interaction :
	def __init__(self, name, type, groups, lims) :
		self.name, self.type = name, type
		self.groups, self.lims = groups, lims
	def acceptableResidue(self, resn, resADtag, atomnames) :
		ret = []
		for g in self.groups :
			retTagAtoms = []
			if g.ADtag != resADtag : continue
			carryon, rescat = 0, getResidueCategory(resn)
			if carryon==0 and g.restype == resn : carryon = 1
			elif carryon==0 and g.restype == 'any' : carryon = 1
			elif carryon==0 and g.restype == 'AA' and 'aa' in rescat : carryon = 1
			elif carryon==0 and g.restype == 'NU' and ('rna' in rescat or 'dna' in rescat) : carryon = 1
			if carryon==0 : continue
			for tat in g.atoms :
				for atn in atomnames :
					if tat.atn == atn : retTagAtoms.append(tat)
			if len(retTagAtoms) == len(g.atoms) : ## all atoms specified in xml are found in this residue
				ret.append(retTagAtoms)
		if len(ret) == 0 : return None
		return ret
	## justOnePerPair indicates whether multiple interactions of same name are permitted for a pair of residues
	def validPartners(self, reslines1, reslines2, atomnbrs, justOnePerPair=1) :
		a1, a2 = [], []
		res1, res2 = line2resid(reslines1[0]), line2resid(reslines2[0])
		resn1, resn2 = line2resn(reslines1[0]), line2resn(reslines2[0])
		print "Checking " + self.name + "-" + res1 + "-" + res2 + "-"
		for l in reslines1 :
			assert line2resid(l) == res1
			a1.append(line2atomname(l))
		for l in reslines2 :
			assert line2resid(l) == res2
			a2.append(line2atomname(l))
		## does res1 qualify to be D-partner in this interaction ?
		tags1 = self.acceptableResidue(resn1, 'D', a1)
		if not tags1 : print res1, 'doesnt qualify as D'
		## does res2 qualify to be A-partner in this interaction ? 
		tags2 = self.acceptableResidue(resn2, 'A', a2)
		if not tags2 : print res2, 'doesnt qualify as A'
		if not tags1 or not tags2 : return None
		lines = []
		## check if tags1, tags2 are voronoi nbrs
		for t1 in tags1 :
			atn1 = []
			for t in t1 : atn1.append(t.atn)
			for t2 in tags2 :
				atn2 = []
				for t in t2 : atn2.append(t.atn)
				vornbrs = None
				for an1,an2 in atomnbrs :
					if an1 in atn1 and an2 in atn2 :
						vornbrs = 1
						break
				if not vornbrs :
					print  "Reqd atoms are not voronoi nbrs"
					continue
				ret = None
				if self.type == 'DHAB' : ret = self.validGeometry_DHAB(t1, t2, reslines1, reslines2)
				elif self.type == 'DHR' : ret = self.validGeometry_DHR(t1, t2, reslines1, reslines2)
				elif self.type == 'RR' : ret = self.validGeometry_RR(t1, t2, reslines1, reslines2)
				elif self.type == 'contact' : ret = self.validGeometry_contact(t1, t2, reslines1, reslines2, atomnbrs)
				else : assert 1==0
				if ret :
					if justOnePerPair : return [ret]
					else : lines.append(ret)
		if len(lines) > 0 : return lines
		return None
	def Rcrd(self, reslines, tags) :
		Ratomnames = []
		for tat in tags :
			if tat.tag == 'R' : Ratomnames.append(tat.atn)
		Rcrd, AMcrd, ANcrd = [], [0.,0.,0.], []
		for rl in reslines :
			if line2atomname(rl) in Ratomnames : Rcrd.append(line2crd(rl))
		assert len(Rcrd) >= 3
		for rc in Rcrd : AMcrd = vec_add(AMcrd, rc)
		AMcrd = vec_scale(AMcrd, 1./len(Rcrd))
		AN = cross_product(vec_diff(Rcrd[0], Rcrd[1]), vec_diff(Rcrd[2], Rcrd[1]))
		ANcrd.append(vec_add(AMcrd, AN))
		ANcrd.append(vec_diff(AMcrd, AN))
		return AMcrd, ANcrd
	def DHcrd(self, reslines, tags) :
		Datomname = None
		for tat in tags :
			if tat.tag == 'D' : Datomname = tat.atn
		Hatomname = '.H' + Datomname[2:4]
		Dcrd, Hcrds = None, []
		for rl in reslines :
			if re.compile(Hatomname).search(line2atomname(rl)) : Hcrds.append(line2crd(rl))
			if Datomname == line2atomname(rl) : Dcrd = line2crd(rl)
		return Dcrd, Hcrds
	def validGeometry_contact(self, tags1, tags2, reslines1, reslines2, atomnbrs) : ## implement in subclasses
		contact = None
		for an1 in tags1 :
			if contact : break
			for an2 in tags2 :
				if contact : break
				if (an1.atn,an2.atn) in atomnbrs : contact = 1
		if not contact : return None
		c1,c2, H1crd, H2crd = 0,0, [0.,0.,0.,], [0.,0.,0.,]
		for rl in reslines1 :
			for an in tags1 :
				if line2atomname(rl) == an.atn : 
					H1crd = vec_add(H1crd, line2crd(rl))
					c1 = c1 + 1
		for rl in reslines2 :
			for an in tags2 :
				if line2atomname(rl) == an.atn : 
					H2crd = vec_add(H2crd, line2crd(rl))
					c2 = c2 + 1
		H1crd = vec_scale(H1crd, 1./c1)
		H2crd = vec_scale(H2crd, 1./c2)
		return pymol_line(H1crd, H2crd)
	def validGeometry_RR(self, tags1, tags2, reslines1, reslines2) : ## implement in subclasses
		## Dm and Dn from reslines1
		DMcrd, DNcrd = self.Rcrd(reslines1, tags1)
		## Am and An from reslines2
		AMcrd, ANcrd = self.Rcrd(reslines2, tags2)
		## limit checks
		tag2points = {}
		tag2points['Dm'] = DMcrd
		tag2points['Am'] = AMcrd
		for dncrd in DNcrd :
			for ancrd in ANcrd :
				tag2points['Dn'] = dncrd
				tag2points['An'] = ancrd
				all_lims_satisfied = 1
				for l in self.lims :
					if not l.check(tag2points) :
						all_lims_satisfied = 0
						break
				if all_lims_satisfied == 1 : return pymol_line(DMcrd, AMcrd)
		return None
#		all_lims_satisfied = 1
#		for l in self.lims :
#			if len(l.tags) == 2 and l.tags[0] in ["Dn","An"] and l.tags[1] in ["Dn","An"] : # angle betn normals
#				v1 = normalize(vec_diff(tag2points["Dn"],tag2points["Dm"]))
#				v2 = normalize(vec_diff(tag2points["An"],tag2points["Am"]))
#				dp = dot_product(v1,v2)
#				if dp > 0.9999 : dp = 0.9999
#				if dp < -0.9999 : dp = -0.9999
#				ang = math.acos(dp) * 180./math.pi
#				if ang > 90 : ang = 180 - ang # acute angle
#				if l.min > ang or ang > l.max : all_lims_satisfied = None
#			elif len(l.tags) == 2 :
#				if not l.check(tag2points) : all_lims_satisfied = None
#			elif len(l.tags) == 3 :
#				v1 = normalize(vec_diff(tag2points[l.tags[0]],tag2points[l.tags[1]]))
#				v2 = normalize(vec_diff(tag2points[l.tags[2]],tag2points[l.tags[1]]))
#				dp = dot_product(v1,v2)
#				if dp > 0.9999 : dp = 0.9999
#				if dp < -0.9999 : dp = -0.9999
#				ang = math.acos(dp) * 180./math.pi
#				if ang > 90 : ang = 180 - ang # acute angle
#				if l.min > ang or ang > l.max : all_lims_satisfied = None
#			else : assert 1==0
#		if all_lims_satisfied == 1 : return pymol_line(DMcrd, AMcrd)
#		return None
	def validGeometry_DHR(self, tags1, tags2, reslines1, reslines2) : ## implement in subclasses
		## D and H coordinates from reslines1
		Dcrd, Hcrds = self.DHcrd(reslines1, tags1)
		## Am and An from reslines2
		AMcrd, ANcrd = self.Rcrd(reslines2, tags2)
		## limit checks
		tag2points = {}
		tag2points['D'] = Dcrd
		tag2points['Am'] = AMcrd
		for hcrd in Hcrds :
			tag2points['H'] = hcrd
			for ancrd in ANcrd :
				tag2points['An'] = ancrd
				all_lims_satisfied = 1
				for l in self.lims :
					#print tag2points
					if not l.check(tag2points) :
						all_lims_satisfied = 0
						break
				if all_lims_satisfied == 1 : return pymol_line(Dcrd, AMcrd)
		return None
	def validGeometry_DHAB(self, tags1, tags2, reslines1, reslines2) : ## implement in subclasses
		## D and H coordinates from reslines1
		Dcrd, Hcrds = self.DHcrd(reslines1, tags1)
		## A and B coordinates from reslines2
		Acrd, Bcrd, Aatomname, Batomname = None, None, None, None
		for tat in tags2 :
			if tat.tag == 'A' : Aatomname = tat.atn
			if tat.tag == 'B' : Batomname = tat.atn
		for rl in reslines2 :
			if Aatomname == line2atomname(rl) : Acrd = line2crd(rl)
			if Batomname == line2atomname(rl) : Bcrd = line2crd(rl)
		## limit checks
		tag2points = {}
		tag2points['D'] = Dcrd
		tag2points['A'] = Acrd
		tag2points['B'] = Bcrd
		if len(Hcrds) == 0 : Hcrds = [(1e6,1e6,1e6)] ## probably hydrogens have no role in this intxn, thats why they were not discovered
		for hcrd in Hcrds :
			tag2points['H'] = hcrd
			all_lims_satisfied = 1
			for l in self.lims :
				#print tag2points
				if not l.check(tag2points) :
					all_lims_satisfied = 0
					break
			if all_lims_satisfied == 1 : return pymol_line(Dcrd, Acrd)
		return None
	def info(self) :
		print "-----------Interaction " + self.name + " " + self.type
		for l in self.lims : l.info()
		for g in self.groups : g.info()

class IntxHandler(xml.sax.handler.ContentHandler) :
	def __init__(self,verbose=None):
		self.reset()
		self.verbose = verbose
	def reset(self):
		self.intxn = []
		self.cur_intxn = None
	def startElement(self, name, attrs) :
		if self.verbose : print "START",name, attrs.getNames()
		if name == 'interaction_set' : pass
		elif name == 'interaction' :
			self.cur_intxn = Interaction(attrs.get("name"), attrs.get("type"), [], [])
		elif name == 'lim':
			ats, min, max = [], {},{}
			for at in attrs.getNames() :
				if at[0:3] == 'tag' : ats.append(attrs.get(at))
				elif re.compile('^max').search(at) :
					max[re.sub("max", "", at)+'_'] = string.atof(attrs.get(at))
				elif re.compile('^min').search(at) :
					min[re.sub("min", "", at)+'_'] = string.atof(attrs.get(at))
				#elif at == 'max' : max = string.atof(attrs.get(at))
				#elif at == 'min' : min = string.atof(attrs.get(at))
				else : assert 1==0
			min1,max1 = [],[]
			for k in min.keys() :
				assert k in max.keys()
				max1.append(max[k])
				min1.append(min[k])
			self.cur_intxn.lims.append( GeomLim(ats, min1, max1) )
		elif name == 'group' :
			self.cur_intxn.groups.append(TaggedGroup(attrs.get('res'), attrs.get('AD'), []))
		elif name == 'atom' :
			gi = len(self.cur_intxn.groups) - 1
			self.cur_intxn.groups[gi].atoms.append(TaggedAtom(attrs.get('atn'), attrs.get('tag')))
		else : assert 1==0
	def endElement(self, name) :
		if self.verbose : print "END",name
		if name == 'interaction_set' : pass
		elif name == 'interaction' : self.intxn.append(self.cur_intxn)
		elif name == 'lim' : pass
		elif name == 'group' : pass
		elif name == 'atom' : pass
		else : assert 1==0
	def endDocument(self) : pass
	def characters(self,content) :
		content = re.sub("\n", "", content)
		if content == '' : return
		content = content.encode("ascii")
		if self.verbose : print "-%s-" % content



if __name__=='__main__' :

	parser = optparse.OptionParser() 
	parser.add_option("--pdbin", action='store', type='string', dest="pdbin", help='pdb file')
	parser.add_option("--intxml", action='store', type='string', dest="intxml", help='xml file describing interactions to be detected')
	parser.add_option("--out", action='store', type='string', dest="out", help='output xml file containing interactions, in a format readable with pymol plugin')
	parser.add_option("--reprez", action='store', type='string', dest="reprez", help='pdb file')

	(options, args) = parser.parse_args()

	intxml = options.intxml
	outfn = options.out

	## solvate protein and decide voronoi adjacency in un-reduced protein
	prot = protein(options.pdbin, read_hydrogens=0, read_waters=1, read_hets=1) ## pdb file here
	HOHcrd = makeGridZimmer(prot.allcrds())
	vor_in = []
	vor_in.append("3")
	vor_in.append("%d" % (len(prot.atomlines)+len(HOHcrd)))
	for l in prot.atomlines : vor_in.append("%f %f %f" % line2crd(l))
	for crd in HOHcrd : vor_in.append("%f %f %f" % crd)
	dontcare_estatus,vorface_lines,dontcare_elines = proc_run_exitOnError("qvoronoi Fv", vor_in) # voronoi vertices
	vornbrs = {}
	for i in range(len(prot.reslines)) :
		if not isPdbAAline(prot.atomlines[ prot.reslines[i][0] ]) : continue
		resid1 = line2resid( prot.atomlines[prot.reslines[i][0]] )
		for j in range(len(prot.reslines)) :
			if not isPdbAAline(prot.atomlines[ prot.reslines[j][0] ]) : continue
			resid2 = line2resid( prot.atomlines[prot.reslines[j][0]] )
			vornbrs[(resid1,resid2)] = []
	for l in vorface_lines[1:] : #ignore 1st line of output
		flds = l.split()[1:] #ignore first field of each line
		a0, a1 = string.atoi(flds[0]), string.atoi(flds[1])
		if a0 >= len(prot.atomlines) or a1 >= len(prot.atomlines) : continue
		if not isPdbAAline(prot.atomlines[a0]) or not isPdbAAline(prot.atomlines[a1]) : continue
		resid1 = line2resid(prot.atomlines[a0])
		resid2 = line2resid(prot.atomlines[a1])
		an1 = line2atomname(prot.atomlines[a0])
		an2 = line2atomname(prot.atomlines[a1])
		vornbrs[(resid1,resid2)].append((an1,an2))
		vornbrs[(resid2,resid1)].append((an2,an1))

	## reduce pdb file
	pdbfile = mkstemp(suffix = ".pdb")[1]
	print pdbfile
	proc_run_exitOnError('reduce %s > %s' % (options.pdbin, pdbfile), [], [256])
	prot = protein(pdbfile, read_hydrogens=1, read_waters=0, read_hets=0) ## pdb file here
	os.unlink(pdbfile)

	# read interactions from xml file
	reader = xml.sax.make_parser()
	ih = IntxHandler()
	reader.setContentHandler(ih)
	reader.parse(intxml) ## xml file here

	for intx in ih.intxn : intx.info()

	#prot.info()
	reslines, reprez = [], []
	for ri1 in range(len(prot.reslines)) :
		s1,e1 = prot.reslines[ri1][0], prot.reslines[ri1][1]
		rl = []
		mcrd, cacrd = [0.,0.,0.], None
		for i in range(s1,e1) :
			crd = line2crd(prot.atomlines[i])
			if line2atomname(prot.atomlines[i]) == " CA " : cacrd = line2crd(prot.atomlines[i])
			for k in range(3) : mcrd[k] = mcrd[k] + crd[k]
			rl.append(prot.atomlines[i])
		reslines.append(rl)
		for k in range(3) : mcrd[k] = mcrd[k] / len(rl)
		if cacrd : reprez.append(cacrd)
		else : reprez.append(mcrd)

	f = open(outfn, 'w') ## xml out file here
	assert(f)
	f.write("<pync_interactions>\n")
	for ri1 in range(len(prot.reslines)) :
		resi1 = line2resid(reslines[ri1][0])
		for ri2 in range(len(prot.reslines)) :
			if ri2 == ri1 : continue
			resi2 = line2resid(reslines[ri2][0])
			if len(vornbrs[(resi1,resi2)]) == 0 : continue
			print "Residues", resi1, resi2, "are nbrs", vornbrs[(resi1,resi2)]
			#if re.compile('628').search(line2resid(reslines[ri1][0])) == None or re.compile('626').search(line2resid(reslines[ri2][0])) == None : continue
			for intx in ih.intxn :
				if intx.name == 'hydrophobic' and ri1 >= ri2 : continue ## undirectional interactions, count a pair just once
				ret = intx.validPartners(reslines[ri1], reslines[ri2], vornbrs[(resi1,resi2)])
				if not ret : continue
				else :
					if options.reprez == 'residue' :
						f.write('<intxn res1="%s" res2="%s" name="%s" p1="%f %f %f" p2="%f %f %f"></intxn>\n' %\
							(resi1, resi2, intx.name, reprez[ri1][0],reprez[ri1][1],reprez[ri1][2], reprez[ri2][0],reprez[ri2][1],reprez[ri2][2]))
					else :#print ret
						for p,q in ret :
							f.write('<intxn res1="%s" res2="%s" name="%s" p1="%f %f %f" p2="%f %f %f"></intxn>\n' %\
								(resi1, resi2, intx.name, p[0],p[1],p[2], q[0],q[1],q[2]))
		#		if not intx.validGeometry() : continue
				print "Possible", intx.name, "Interaction between %s and %s" % ( line2resid(reslines[ri1][0]), line2resid(reslines[ri2][0]))
	f.write("</pync_interactions>\n")
	f.close()
