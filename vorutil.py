#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################


import string, re, cPickle, os

from procrun import *
from pdbr import peptideCN, sameResidue, covBonded, line2resid, isPdbAAline, isPdbHOHline, isPdbHetLine, phosSugarO3P
from geom import *
from transform import calcTransform

def runVoronoi(scrds) :
	assert(len(scrds) > 0)
	#for si in range(len(scrds)) : print si, "SCRD", scrds[si]
	vor_in = []
	vor_in.append("3")
	vor_in.append("%d" % len(scrds))
	for s in scrds : vor_in.append("%f %f %f" % (s[0],s[1],s[2]))
	dontcare_estatus,vorcen_lines,dontcare_elines = proc_run_exitOnError("qvoronoi p", vor_in) # voronoi vertices
	vorcen = [ [0,0,0], ]
	for l in vorcen_lines[2:] :
		flds = l.split()
		vorcen.append( [ string.atof(flds[0]), string.atof(flds[1]), string.atof(flds[2]) ] )
	adjlist = []
	dontcare_estatus,vorface_lines,dontcare_elines = proc_run_exitOnError("qvoronoi Fv", vor_in) # sites and vor faces
	for l in vorface_lines[1:] :
		flds, a = l.split(), []
		for f in flds[1:] : a.append( string.atoi(f) )
		adjlist.append(a)
	return vorcen, adjlist

#predefine some standard colors
chem2rgb = {}
chem2rgb['mc'] = (0,1,0)
chem2rgb['sc'] = (0,1,1)
chem2rgb['HOH'] = (0,0,1)
chem2rgb['phos'] = (1,0.9,0.5)
chem2rgb['sugar'] = (0.2,0.35,0.35)
chem2rgb['base'] = (0.625,0.125,0.9375)
chem2rgb['het'] = (1,0,1)
chem2rgb['phob'] = (0,0,1)
chem2rgb['phil'] = (1,0,0)
chem2rgb['arom'] = (1,1,0)
# clr2rgb = {}
# clr2rgb['white'] = (1,1,1)
# clr2rgb['grey'] = (0.2,0.35,0.35)
# clr2rgb['black'] = (0,0,0)
# clr2rgb['red'] = (1,0,0)
# clr2rgb['green'] = (0,1,0)
# clr2rgb['blue'] = (0,0,1)
# clr2rgb['cyan'] = (0,1,1)
# clr2rgb['yellow'] = (1,1,0)
# clr2rgb['magenta'] = (1,0,1)
# clr2rgb['purple'] = (0.625,0.125,0.9375)
# clr2rgb['gold'] = (1,0.9,0.5)
# clr2rgb['pink'] = (1,0,0.5)

def areGroupsCovbonded(sn1, sn2, adjResids) :
	if len(adjResids) == 0 : return None
	for atmline1 in sn1 :
		if isPdbHOHline(atmline1) or isPdbHetLine(atmline1) : continue
		resid1 = line2resid(atmline1)
		for atmline2 in sn2 :
			if isPdbHOHline(atmline2) or isPdbHetLine(atmline2) : continue
			resid2 = line2resid(atmline2)
			if (resid1,resid2) in adjResids and (peptideCN(atmline1,atmline2) or phosSugarO3P(atmline1,atmline2)) : return 1
			if (resid2,resid1) in adjResids and (peptideCN(atmline2,atmline1) or phosSugarO3P(atmline2,atmline1)) : return 1
			if sameResidue(atmline1, atmline2) :
				if covBonded(atmline1,atmline2) : return 1
	return None
## if a metasite of removeChem is surrounded completely by other metasites of color removeChem,
## that metasite is no longer considered. contacts between 2 metasites of removeChem are neglected.
## removeChem is useful for neglecting unnecessary waters.
class Vorutil :
	def find_MSadj_faces(self,fadj, slave2master) :
		msadj = {}
		for fa in fadj :
			m0, m1, fi = slave2master[fa[0]], slave2master[fa[1]], fa[2]
			if m0 == m1 : continue
			if m0 > m1 :
				temp = m0
				m0 = m1
				m1 = temp
			try : msadj[(m0,m1)].append(fi)
			except KeyError :
			#if (m0,m1) in msadj.keys() : msadj[(m0,m1)].append(fi)
			#else :
				msadj[(m0,m1)] = [fi,]
		return msadj
	def adjlist2faceAdjlist(self,adjlist) :
		fadj, faces = [], []
		for a in adjlist :
			if 0 in a[2:] : continue
			fadj.append( (a[0],a[1],len(faces)) )
			faces.append(a[2:])
		return fadj, faces
	def slave_to_master(self,master2slave) :
		slave2master = []
		max = -1
		for msi in range(len(master2slave)) :
			for s in master2slave[msi] :
				if max < s : max = s
		for i in range(max+1) : slave2master.append(-1)
		for msi in range(len(master2slave)) :
			for s in master2slave[msi] : slave2master[s] = msi
		return slave2master
	def __init__(self, sites, msnames, master2slaves, mschems, slavenames, msframes, removeChem, calcVols=0):
		self.chem2clr = {}
		self.alpha = 0.5
		self.conwidth = 0.1
		slave2master = self.slave_to_master(master2slaves)
		## remove orphan sites
		newsites, newS2M = [],[]
		for si in range(len(slave2master)) :
			if slave2master[si] == -1 : continue
			newsites.append(sites[si])
			newS2M.append(slave2master[si])
		sites, slave2master = newsites, newS2M
		for msi in range(len(msnames)) :
			print "METASITE", msnames[msi], mschems[msi]
			for s in slavenames[msi] : print "SLAVE", s
			#for si in master2slaves[msi] : print sites[si]
		## voronoi
		vcen, adjlist = runVoronoi(sites)
		#for vc in vcen :
		#	print vc
		#	if math.fabs(vc[0]) > 1e5 : assert(1==0)
		#	if math.fabs(vc[1]) > 1e5 : assert(1==0)
		#	if math.fabs(vc[2]) > 1e5 : assert(1==0)
		print "DONE slave2master"
		fadj, faces = self.adjlist2faceAdjlist(adjlist)
		print "DONE adjlist2faceAdjlist"
		varea = []
		for f in faces : varea.append(-10000000000.)
		for fa in fadj :
			if re.compile('HOH').search(msnames[slave2master[fa[0]]]) and re.compile('HOH').search(msnames[slave2master[fa[1]]]) : continue
			vc = []
			for vi in faces[fa[2]] : vc.append(vcen[vi])
			varea[fa[2]] = poly_area(vc)
		print "DONE areas"
		self.faces, self.areas, self.vcen, self.slave2master, self.master2slave, self.sites = faces, varea, vcen, slave2master, master2slaves, sites
###################
		self.sitevols = []
		cells = {}
		if calcVols == 'no' :
			for si in range(len(sites)) : self.sitevols.append(-10000000000.)
		else :
			for a in adjlist :
				if 0 in a[2:] : continue
				cells[a[0]] = {}
				cells[a[1]] = {}
			for a in adjlist :
				if 0 in a[2:] : continue
				for i in a[2:] :
					cells[a[0]][i] = 0
					cells[a[1]][i] = 0
			for si in range(len(sites)) :
				if re.compile('HOH').search(slavenames[slave2master[si]][0]) : self.sitevols.append(-10000000000.)
				else :
					vertices = []
					for k in cells[si].keys() : vertices.append(vcen[k])
					self.sitevols.append(poly_vol_sa(vertices)[0])
		print "DONE volumes"
###################
		old_msadj = self.find_MSadj_faces(fadj, slave2master)
		print "DONE find_MSadj_faces"
		msadj = {}
		for (m0,m1) in old_msadj.keys() :
			if mschems[m0] == mschems[m1] and mschems[m0] == removeChem : continue
			if mschems[m0] != removeChem :
				if m0 not in msadj.keys() : msadj[m0] = {}
				msadj[m0][m1] = old_msadj[(m0,m1)]
			if mschems[m1] != removeChem :
				if m1 not in msadj.keys() : msadj[m1] = {}
				msadj[m1][m0] = old_msadj[(m0,m1)]
		self.msadj = msadj
		self.msnames = msnames
		self.msframes = msframes
		self.mschems = []
		for msc in mschems : self.mschems.append(msc)
		self.slavenames = slavenames

	#make faces into objects
	def makeFacePymolStr(self, facestyle) :
		facestr = []
		if facestyle == 'lines' :
			for fi in range(len(self.faces)) :
				vorinds = self.faces[fi]
				facestr.append("BEGIN,LINES,")
				for i in range(len(vorinds)) :
					p,q = i,i+1
					if q == len(vorinds) : q = 0
					p,q = vorinds[p],vorinds[q]
					facestr[fi] = facestr[fi] + "VERTEX,%f,%f,%f,VERTEX,%f,%f,%f," % \
						(self.vcen[p][0], self.vcen[p][1], self.vcen[p][2], self.vcen[q][0], self.vcen[q][1], self.vcen[q][2])
				facestr[fi] = facestr[fi] + "END"
		else :
			for fi in range(len(self.faces)) :
				vorinds = self.faces[fi]
				facestr.append("BEGIN,TRIANGLE_FAN,ALPHA,%f" % self.alpha)
				for i in range(len(vorinds)) :
					p = vorinds[i]
					facestr[fi] = facestr[fi] + "VERTEX,%f,%f,%f," % (self.vcen[p][0], self.vcen[p][1], self.vcen[p][2])
				facestr[fi] = facestr[fi] + "END"
		self.facestr = facestr

	def writeAreasNew(self, filename, style, MSI, MSI_1=None) :
		f = fileOpen(filename, 'w')
		if not f : return
		common_nbrlist = []
		for msi0 in MSI : # render msi0-th site
			avcrd = [0.,0.,0.]
			for si in self.master2slave[msi0] :
				for i in range(3) : avcrd[i] = avcrd[i] + self.sites[si][i]/len(self.master2slave[msi0])
			f.write("%s %f %f %f %f METASITE %s\n" % (self.mschems[msi0],self.sitevols[msi0], avcrd[0],avcrd[1],avcrd[2], self.msnames[msi0]))
			if style == 'exposed_surface' : ## surface of each metasite in MSI that isnt shared with a metasite in MSI
				area = 0
				for msi1 in self.msadj[msi0].keys() :
					if msi1 in MSI : continue
					for fi in self.msadj[msi0][msi1] : area = area + self.areas[fi]
				if area > 0 : f.write("%f\n" % area)
				else : f.write("0 buried\n")
			elif style == 'all_nbr' : ## nbrlist of each metasite in MSI # covbonded nbrs are filtered depending on size of self.adjResids
				nbrs = []
				areas, hohArea = [], 0.
				for msi1 in self.msadj[msi0].keys() :
					if areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) : continue
						#print "COV META", self.msnames[msi0], self.msnames[msi1]
					area = self.getCommonArea(msi0,msi1)
					nbrs.append(self.msnames[msi1])
					areas.append(area)
				totarea = 0
				for a in areas : totarea = totarea + a
				for ni in range(len(nbrs)) : f.write("%s NBR %f %f\n" % (nbrs[ni], areas[ni], areas[ni]/totarea))
				f.write('\n')
			elif style == 'all_nbr_orient' : ## nbrlist of each metasite in MSI # covbonded nbrs are filtered # orientations (angle,dihed) are calc in given frame
				print self.msframes[msi0], self.msnames[msi0]
				a,b,c = self.msframes[msi0]
				nbrs, areas, angles, diheds = [],[],[],[]
				for msi1 in self.msadj[msi0].keys() :
					if areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) : continue
						#print "COV META", self.msnames[msi0], self.msnames[msi1]
					area = self.getCommonArea(msi0,msi1)
					nbrs.append("%s %s" % (self.mschems[msi1],self.msnames[msi1]))
					areas.append(area)
					d = [0.,0.,0.] # average teh coordinates of nbr to get d
					if len(self.master2slave[msi1]) == 1 : d = tuple(self.sites[self.master2slave[msi1][0]])
					else :
						assert len(self.master2slave[msi1]) > 1
						for si in self.master2slave[msi1] : d[0],d[1],d[2] = d[0]+self.sites[si][0], d[1]+self.sites[si][1], d[2]+self.sites[si][2]
						d[0],d[1],d[2] = d[0]/len(self.master2slave[msi1]), d[1]/len(self.master2slave[msi1]), d[2]/len(self.master2slave[msi1])
					#print self.msnames[msi0], self.msnames[msi1],a,b,c,d
					angles.append(vec_angle_deg(b,c,d))
					if angles[len(angles)-1] < 1. or angles[len(angles)-1] > 178. : # dihed is ill-defined 
						diheds.append(0) # hence arbit set it to 0
					else : diheds.append(vec_dihed_deg(a,b,c,d))
					#print angles[len(angles)-1], diheds[len(diheds)-1]
				totarea = 0
				for ar in areas : totarea = totarea + ar
				for ni in range(len(nbrs)) : f.write("%s %f %f %f %f NBR %s\n" % (re.sub(" \[.*","",nbrs[ni]), areas[ni], 100.*areas[ni]/totarea, angles[ni], diheds[ni], re.sub(".* \[","[",nbrs[ni])))
				tr = calcTransform(a,b,c)
				f.write("TRANSFORM %f %f %f %f %f %f %f %f\n" % tr)
				f.write('\n')
			elif style == 'all_physchem' : ## nbrlist of all metasites in MSI, grouped on colour
				areas = {}
				for c in self.mschems : areas[c] = [0.,0] ## area-sum and number of neighbours
				for msi1 in self.msadj[msi0].keys() :
					area = self.getCommonArea(msi0,msi1)
					areas[self.mschems[msi1]][0] = areas[self.mschems[msi1]][0] + area
					areas[self.mschems[msi1]][1] = areas[self.mschems[msi1]][1] + 1
				totarea = 0
				for k,v in areas.items() : totarea = totarea + v[0]
				for k,v in areas.items() : f.write("%s %d %6.2f %6.2f\n" % (k, v[1],v[0],100*v[0]/totarea))
				f.write('\n')
			elif style == 'all_color' : ## nbrlist of all metasites in MSI, grouped on colour
				areas = {}
				for c in self.mschems : areas[c] = 0.
				for msi1 in self.msadj[msi0].keys() :
					if areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) : continue
					area = self.getCommonArea(msi0,msi1)
					areas[self.mschems[msi1]] = areas[self.mschems[msi1]] + area
				totarea = 0
				for k,v in areas.items() : totarea = totarea + v
				for k,v in areas.items() : f.write(" %6.2f(%6.2f)" % (v,100*v/totarea))
				f.write('\n')
			elif style == 'common_nbrlist' : ## list of non-HOH neighbours of MSI
				for msi1 in self.msadj[msi0].keys() :
					if re.compile('HOH').search(self.msnames[msi1]) : continue
					if msi1 not in common_nbrlist :
						common_nbrlist.append(msi1)
						f.write("%s\n" % self.msnames[msi1])
			elif style == 'inner_nbrs' : ## list of neighbours within MSI
				for msi1 in self.msadj[msi0].keys() :
					if msi1 < msi0 or msi1 not in MSI : continue
					if areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) : continue
					f.write("%s %f\n" % (self.msnames[msi1], self.getCommonArea(msi0, msi1)))
			elif style == 'ext_nbrs' : ## all nbr-pairs, one from MSI other from MSI_1
				assert(MSI_1)
				for msi1 in self.msadj[msi0].keys() :
					if msi1 not in MSI_1 : continue
					if areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) : continue
					f.write("%s %f\n" % (self.msnames[msi1], self.getCommonArea(msi0,msi1)))
			else : assert(1==0)
		f.close()


	def writeAreas(self, filename, style, removeChem, collapseWaters = 1) :
		assert(style == 'exposed_surface' or style == 'all_nbr' or style == 'all_color')
		f = fileOpen(filename, 'w')
		if not f : return
		for msi0 in self.msadj.keys() : # render msi0-th site
			if self.mschems[msi0] == removeChem : continue
			f.write("METASITE %s\n" % self.msnames[msi0])
			if style == 'exposed_surface' :
				area = 0
				for msi1 in self.msadj[msi0].keys() :
					if self.mschems[msi1] != removeChem : continue
					for fi in self.msadj[msi0][msi1] :
						area = area + self.areas[fi]
				if area > 0 : f.write("%f\n" % area)
			elif style == 'all_nbr' : # covbonded nbrs are filtered depending on size of self.adjResids
				nbrs, areas = [], []
				for msi1 in self.msadj[msi0].keys() :
					if len(self.adjResids) > 0 and areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) :
						#print "COV META", self.msnames[msi0], self.msnames[msi1]
						continue
					area = 0
					for fi in self.msadj[msi0][msi1] : area = area + self.areas[fi]
					nbrs.append(msi1)
					areas.append(area)
				f.write("NBRLIST")
				if collapseWaters == 1 :
					nbrs.append('HOH')
					areas.append(hohArea)
				for n in nbrs :
					f.write(" %s" % self.msnames[n])
				f.write('\n')
				totarea = 0
				for a in areas : totarea = totarea + a
				f.write("NBRAREA")
				for a in areas : f.write(" %6.2f(%6.2f)" % (a, 100*a/totarea))
				f.write('\n')
			elif style == 'all_color' :
				areas = {}
				for c in self.mschems : areas[c] = 0.
				for msi1 in self.msadj[msi0].keys() :
					if len(self.adjResids) > 0 and areGroupsCovbonded(self.slavenames[msi0],self.slavenames[msi1], self.adjResids) : continue
					area = 0
					for fi in self.msadj[msi0][msi1] : area = area + self.areas[fi]
					areas[self.mschems[msi1]] = areas[self.mschems[msi1]] + area
				totarea = 0
				for k,v in areas.items() : totarea = totarea + v
				for k,v in areas.items() : f.write(" %6.2f(%6.2f)" % (v,100*v/totarea))
				f.write('\n')
		else : assert(1==0)
		f.close()

	############# functions for running from within pymol
	def makeFacePymolObj(self, vcen, vorinds, faceclr, lineclr) :
		from pymol.cgo import BEGIN, TRIANGLE_FAN, END, VERTEX, COLOR, LINES, ALPHA
		faceobj = [BEGIN, TRIANGLE_FAN, ALPHA, self.alpha, COLOR, faceclr[0], faceclr[1], faceclr[2]]
		for vi in vorinds : faceobj = faceobj + [VERTEX, self.vcen[vi][0], self.vcen[vi][1], self.vcen[vi][2]]
		faceobj.append(END)
		faceobj = faceobj + [BEGIN, LINES, COLOR, lineclr[0], lineclr[1], lineclr[2]]
		for i in range(len(vorinds)) :
			vi1, vi2 = vorinds[i % len(vorinds)], vorinds[(i+1) % len(vorinds)]
			faceobj = faceobj + [VERTEX, self.vcen[vi1][0], self.vcen[vi1][1], self.vcen[vi1][2]]
			faceobj = faceobj + [VERTEX, self.vcen[vi2][0], self.vcen[vi2][1], self.vcen[vi2][2]]
		faceobj.append(END)
		return faceobj

	def showSurf(self,name,msis) :
		from pymol import cmd
		assert len(self.chem2clr.items()) > 0
		totarea, msobj, FIs, faceclr, lineclr = 0., [], [], [], []
		for msi0 in self.msadj.keys() :
			if re.compile('HOH').search(self.msnames[msi0]) : continue ## no waters
			if len(msis) > 0 and msi0 not in msis : continue
			for msi1 in self.msadj[msi0].keys() :
				#if len(msis) > 0 and msi1 not in msis : continue
				#if msi1 >= msi0 : continue
				#print self.msadj[msi0][msi1]
				for fi in self.msadj[msi0][msi1] :
					if fi in FIs : FIs.remove(fi)
					else : FIs.append(fi)
		for msi0 in self.msadj.keys() :
			if len(msis) > 0 and msi0 not in msis : continue
			if re.compile('HOH').search(self.msnames[msi0]) : continue ## no waters
			for msi1 in self.msadj[msi0].keys() :
				#if msi1 >= msi0 : continue
				#if len(msis) > 0 and msi1 not in msis : continue
				for fi in self.msadj[msi0][msi1] :
					if fi in FIs :
						totarea = totarea + self.areas[fi]
						faceclr, lineclr = self.chem2clr[self.mschems[msi0]], self.chem2clr[self.mschems[msi1]] ## own clr for surface
						msobj = msobj + self.makeFacePymolObj(self.vcen, self.faces[fi], faceclr, lineclr)
		cmd.load_cgo(msobj, name)
		return totarea
	def showvs(self, msi0, name) :
		assert len(self.chem2clr.items()) > 0
		from pymol import cmd
		msobj = []
		for msi1 in self.msadj[msi0].keys() :
			faceclr, lineclr, faceobj = self.chem2clr[self.mschems[msi1]], self.chem2clr[self.mschems[msi0]], []
			for fi in self.msadj[msi0][msi1] :
				faceobj = faceobj + self.makeFacePymolObj(self.vcen, self.faces[fi], faceclr, lineclr)
			msobj = msobj + faceobj
		cmd.load_cgo(msobj, name)
	def showLine(self, msi0, msi1, name='conn') :
		from pymol import cmd
		from pymol.cgo import BEGIN, END, LINES, VERTEX, COLOR, CYLINDER, SPHERE
		mean0, mean1, mid = [0.,0.,0.], [0., 0., 0.], [0.,0.,0.,]
		for si in self.master2slave[msi0] :
			mean0[0] = mean0[0] + self.sites[si][0]
			mean0[1] = mean0[1] + self.sites[si][1]
			mean0[2] = mean0[2] + self.sites[si][2]
		mean0[0] = mean0[0] / len(self.master2slave[msi0])
		mean0[1] = mean0[1] / len(self.master2slave[msi0])
		mean0[2] = mean0[2] / len(self.master2slave[msi0])
		for si in self.master2slave[msi1] :
			mean1[0] = mean1[0] + self.sites[si][0]
			mean1[1] = mean1[1] + self.sites[si][1]
			mean1[2] = mean1[2] + self.sites[si][2]
		mean1[0] = mean1[0] / len(self.master2slave[msi1])
		mean1[1] = mean1[1] / len(self.master2slave[msi1])
		mean1[2] = mean1[2] / len(self.master2slave[msi1])
		mid[0] = (mean0[0] + mean1[0]) / 2.
		mid[1] = (mean0[1] + mean1[1]) / 2.
		mid[2] = (mean0[2] + mean1[2]) / 2.
		clr0, clr1 = self.chem2clr[self.mschems[msi0]], self.chem2clr[self.mschems[msi1]]
		lobj = [CYLINDER, ] + list(mean0) + list(mean1) + [self.conwidth,] + list(clr0) + list(clr1)
		lobj = lobj + [COLOR,] + list(clr0) + [SPHERE,] + list(mean0) + [2*self.conwidth,] + [COLOR,] + list(clr1) + [SPHERE,] + list(mean1) + [2*self.conwidth,]
		print lobj
		#lobj = [BEGIN, LINES, COLOR,] + list(clr0) + [VERTEX,] + list(mean0) + [VERTEX,] + list(mid) + [COLOR,] + list(clr1) + [VERTEX,] + list(mid) + [VERTEX,] + list(mean1) + [END,]
		cmd.load_cgo(lobj, name)

	# show surface common to 2 sets, colored with colors of set1, ie area of set1 exposed to set2
	def showIntxn(self, nbrs, name='intxn') :
		from pymol import cmd
		faces1, faces2, msobj, totarea = [],[],[],0.
		for np in nbrs :
			msi0 = np[0]
			if re.compile('HOH').search(self.msnames[msi0]) : continue ## no waters
			for msi1 in self.msadj[msi0].keys() :
				for fi in self.msadj[msi0][msi1] :
					faces1.append(fi)
		for np in nbrs :
			msi0 = np[1]
			if re.compile('HOH').search(self.msnames[msi0]) : continue ## no waters
			for msi1 in self.msadj[msi0].keys() :
				for fi in self.msadj[msi0][msi1] :
					if fi in faces1 : faces2.append(fi)
		for np in nbrs :
			msi0, msi1 = np[0], np[1]
			for fi in self.msadj[msi0][msi1] :
				if not fi in faces2 : continue
				faceclr, lineclr = self.chem2clr[self.mschems[msi0]], self.chem2clr[self.mschems[msi1]] ## own clr for surface
				msobj = msobj + self.makeFacePymolObj(self.vcen, self.faces[fi], faceclr, lineclr)
				totarea = totarea + self.areas[fi]
		cmd.load_cgo(msobj, name)
		return totarea

	def getCommonArea(self, msi1, msi2) :
		totarea = 0.
		if msi1 not in self.msadj.keys() : return 0.
		if msi2 not in self.msadj[msi1].keys() : return 0.
		for fi in self.msadj[msi1][msi2] : totarea = totarea + self.areas[fi]
		return totarea

	def whichMSindices(self, resn='', chid='', resi='', ic='', atmn='', metasite='', mschem='') :
		print 'whichMS for', resn, ':', chid, ':', resi, ':', ic, ':', atmn, ':', metasite
		## '*' is wildcard
		if len(resn) == 0 : resn = '...'
		elif len(resn) == 1 : resn = '  ' + resn
		elif len(resn) == 2 : resn = ' ' + resn
		elif len(resn) > 3 : assert 1==0
		if chid == '' : chid = '.'
		if ic == '' : ic = '.'
		if resi == '' : resi = '....'
		elif len(resi) == 1 : resi = '   ' + resi
		elif len(resi) == 2 : resi = '  ' + resi
		elif len(resi) == 3 : resi = ' ' + resi
		elif len(resi) == 5 : ic = ''
		elif len(resi) > 5 : assert 0
		if atmn == '' : atmn = '....'
		elif len(atmn) == 1 : atmn = ' ' + atmn + '  '
		elif len(atmn) == 2 : atmn = ' ' + atmn + ' '
		elif len(atmn) == 3 : atmn = ' ' + atmn
		elif len(atmn) > 4 : assert 1==0
		restr = '\[' + resn + chid + resi + ic + ']_' + metasite
		print "LOOKING FOR", restr
		msis, rexp, at_rexp = [], re.compile(restr), re.compile('atmn')
		for msi in range(len(self.msnames)) :
			if re.compile('HOH').search(self.msnames[msi]) : continue
			if mschem != '' and self.mschems[msi] != mschem : continue
			if not rexp.search(self.msnames[msi]) : continue
			print self.msnames[msi]
			if atmn == '' : msis.append(msi)
			else :
				for slavename in self.slavenames[msi] :
					if re.compile('HOH').search(slavename) : continue
					if re.compile('............' + atmn).search(slavename) :
						msis.append(msi)
						break
		return msis
	def showMS(self, resn='', chid='', resi='', ic='', metasite='',atmn='', style='fused', name='voro') :
		msis = self.whichMSindices(resn, chid, resi, ic, metasite, atmn)
		for i in msis :
			if style == 'fused' : self.showvs(i, name)
			else : self.showvs(i,self.msnames[i])
