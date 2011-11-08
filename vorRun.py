#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################

import re, string, os, sys, cPickle, optparse, random
from protinfo import *
from grid import makeGrid, makeGridZimmer
from pdbr import *
from procrun import *
from geom import *
from vorutil import Vorutil
from myxml import myXmlParser
import xml, xml.sax

class metasite :
	def __init__(self, name, physchem, style, restype, slaveatoms, frA, frB, frC):
		self.name = name
		self.style = style
		self.physchem = physchem
		self.restype = restype
		self.slaveatoms = slaveatoms
		self.frA, self.frB, self.frC = frA, frB, frC
	def info(self):
		print "--------BEGIN METASITE INFO--------"
		print self.name, self.physchem, self.style
		print self.restype
		print self.slaveatoms
		print "---------END METASITE INFO---------"


class msHandler(xml.sax.handler.ContentHandler) :
#class msHandler :
	def __init__(self,verbose=None):
		self.metasites = []
		self.physchem2color = {}
		self.reset()
		self.verbose = verbose
	def reset(self):
		self.frA, self.frB, self.frC = None, None, None
		self.currEle = None
		self.slaveatoms = []
	def startElement(self, name, uri) :
		self.currEle = name
		if self.verbose : print "START",name
	def endElement(self, name) :
		if self.verbose : print "END",name
		if name == "metasite" :
			self.metasites.append( metasite(self.name, self.physchem, self.style, self.restype, self.slaveatoms, self.frA, self.frB, self.frC) )
			self.reset()
		elif name == "physchem2color" :
			self.physchem2color[self.physchem] = self.color
			self.reset()
		self.currEle = None
	def endDocument(self) :
		for ms in self.metasites : ms.info()
		for k,v in self.physchem2color.items() : print "PHYSCHEM2COLOR :",k,v
	def characters(self,content) :
		content = re.sub("\n", "", content)
		if content == '' : return
		content = content.encode("ascii")
		if self.verbose : print "-%s-" % content
		assert(self.currEle)
		if self.currEle == "metasite" : pass
		elif self.currEle == "name" : self.name = content
		elif self.currEle == "physchem" : self.physchem = content
		elif self.currEle == "restype" : self.restype = content
		elif self.currEle == "slaveatom" : self.slaveatoms.append(content)
		elif self.currEle == "physchem2color" : pass
		elif self.currEle == "color" : self.color = content
		elif self.currEle == "style" : self.style = content
		elif self.currEle == "frA" : self.frA = content
		elif self.currEle == "frB" : self.frB = content
		elif self.currEle == "frC" : self.frC = content
		else : assert(0)



## solvate and renumber
def solvate_prot(pdbfile, EMMDP, PRMDP, SPCGRO, water_thick, solv_md):
	assert(re.compile(".pdb$").search(pdbfile))
	origprot = protein(pdbfile, read_hydrogens=0, read_waters=0, read_hets=1)
	pdbf = re.sub(".pdb$", "", pdbfile)
	proc_run_exitOnError("pdb2gmx -ignh -f %s.pdb -o %s.gro -p %s.top" % (pdbf,pdbf,pdbf), ["4"])
	proc_run_exitOnError("editconf -f %s.gro -o %s.gro -bt tric -c -d %f" % (pdbf,pdbf,water_thick))
	if solv_md != 'yes' : ## no MD here
		proc_run_exitOnError("genbox -cp %s.gro -cs %s -o %s_b4em.pdb -p %s.top" % (pdbf,SPCGRO,pdbf,pdbf))
		proc_run_exitOnError("cat %s_b4em.pdb | sed 's/SOL/HOH/' >  %s_sol_displaced.pdb" % (pdbf,pdbf))
	else : ## MD
		proc_run_exitOnError("genbox -cp %s.gro -cs %s -o %s_b4em.gro -p %s.top" % (pdbf,SPCGRO,pdbf,pdbf))
		proc_run_exitOnError("grompp -f %s -c %s_b4em.gro -p %s.top -o %s_em.tpr" % (EMMDP,pdbf,pdbf,pdbf))
		proc_run_exitOnError("mdrun -s %s_em.tpr -o %s_em.trr -c %s_b4pr.gro" % (pdbf,pdbf,pdbf))
		proc_run_exitOnError("grompp -f %s -c %s_b4pr.gro -p %s.top -o %s_pr.tpr -r %s_b4pr.gro" % (PRMDP,pdbf,pdbf,pdbf,pdbf))
		proc_run_exitOnError("mdrun -s %s_pr.tpr -o %s_pr.trr -c %s_b4md.pdb" % (pdbf,pdbf,pdbf))
		proc_run_exitOnError("cat %s_b4md.pdb | sed 's/SOL/HOH/' >  %s_sol_displaced.pdb" % (pdbf,pdbf))
	proc_run_exitOnError("g_confrms -one -f1 %s.pdb -f2 %s_sol_displaced.pdb -o %s_sol.pdb" % (pdbf,pdbf,pdbf), ["3","3"])
	#proc_run_exitOnError("rm -Rf %s_sol_displaced.pdb" % (pdbf), [])
	os.unlink("%s_sol_displaced.pdb" % (pdbf))
	retpdbfn = "%s_sol.pdb" % pdbf
	proc_run_exitOnError("cat %s | grep HOH > %s_hoh.pdb" % (retpdbfn,pdbf))
	reslines = origprot.nonHOHsampleLines()
	prot = protein(retpdbfn, read_hydrogens=0, read_waters=0, read_hets=1)
	if not prot.changeNonHOHnumbering(reslines) : assert 1==0
	prot.write(retpdbfn)
	proc_run_exitOnError("cat %s_hoh.pdb >> %s" % (pdbf, retpdbfn))
	#proc_run_exitOnError("rm -Rf %s_hoh.pdb" % (pdbf))
	os.unlink("%s_hoh.pdb" % (pdbf))
	return retpdbfn

def prepareMDP(template_mdp, mdpout, solv_time) :
	lines = fileOpen(template_mdp, 'r').readlines()
	f = fileOpen(mdpout, 'w')
	for l in lines :
		l = re.sub('\n', '', l)
		if re.compile("dt *=").search(l) : l = "dt = 0.002"
		if re.compile("nsteps *=").search(l) : l = "nsteps = %f" % (solv_time/0.002)
		f.write(l + '\n')
	f.close()

def main(args) :
	parser = optparse.OptionParser()
	parser.add_option("--scratchdir", action='store', type='string', dest="scratchdir", help='all files are created here')
	parser.add_option("--pdbin", action='store', type='string', dest='pdbfile_in', help='input pdb file, used as read-only')
	parser.add_option("--recipe", action='store', type='string', dest="msXmlFile", help='xml file describing how to generate voronoi tesellations, see PROVATPATH/VorXmlRecipes for samples')
	parser.add_option("--checkMCcont", action='store', type='string', dest="checkMCcont", help='check whether consecutive CA atoms in same chain are closer than 5A. set to yes for checking (default), no for ignoring', default='yes')
	parser.add_option("--checkMCmissing", action='store', type='string', dest="checkMCmissing", help='check whether any of N,CA,C,O atoms are missing from amino acid residues. set to yes for checking (default) or no for ignoring', default='yes')
	parser.add_option("--checkSCmissing", action='store', type='string', dest="checkSCmissing", help='set to check, repair or ignore. check will make Provat halt if missing SC atoms, repair will rebuild bad sidechains with Scwrl (default), ignore will ignore', default='repair')
	parser.add_option("--read-het", action='store', type='string', dest="ReadHet", help='whether to read het atms from input pdb file, set to no to ignore hetatms, else to yes (default)', default='yes')
	parser.add_option("--read-hoh", action='store', type='string', dest="ReadHoh", help='whether to read waters from input pdb file, set to no to ignore, else to yes (default)', default='yes')
	parser.add_option("--useGMXsolv", action='store', type='string', dest="useGMXsolv", help='if set to yes, a water grid is put around/into protein, else without using gromacs (default)', default='no')
	parser.add_option("--emmdp", action='store', type='string', dest="EMMDP", help='path to em.mdp file for gromacs, treated read-only, sample is PROVATPATH/Gromacs/em.mdp')
	parser.add_option("--prmdp", action='store', type='string', dest="PRMDP", help='path to em.mdp file for gromacs, treated read-only, sample is PROVATPATH/Gromacs/pr.mdp')
	parser.add_option("--spcgro", action='store', type='string', dest="SPCGRO", help='path to water file for gromacs, treated read-only, sample is PROVATPATH/Gromacs/spc216.gro')
	parser.add_option("--water-thickness", action='store', type='float', dest="water_thick", help='thickness of water shell around protein in angstrom')
	parser.add_option("--em-time", action='store', type='float', dest="emTime", help='time for which energy minimization is carried out, in ps')
	parser.add_option("--pr-time", action='store', type='float', dest="prTime", help='time for which pos-restrained MD is carried out, in ps')
	parser.add_option("--gridDll", action='store', type='int', dest="gridDll", help='Zimmer grid parameter, min distance between any two solvent atoms, positive integer, defaults to 3', default=3)
	parser.add_option("--gridDpl", action='store', type='int', dest="gridDpl", help='Zimmer grid parameter, min distance between any solvent atom and any non-solvent atom, positive integer, defaults to 4', default=4)
	parser.add_option("--calcVols", action='store', type='string', dest="calcVols", help='decides whether to calculate volumes of individual polygons because that is a time consuming operation. no by default.', default='no')
	parser.add_option("--ignore-covbonded", action='store', type='string', dest="ignore_covbonded", help='set to yes if interfaces between metasites having covalent bonds are not to be reported. no by default.', default='no')
	parser.add_option("--datastyle", action='store', type='string', dest="datastyle", help='descibes the style in which interactions surface data is to be written : all_nbr (default), all_nbr_orient or all_color or exposed_surface', default='all_nbr')
	parser.add_option("--outChain", action='store', type='string', dest="outChain", help='use only this chain for writing out datafile. if all chains are desired, use "." for this argument (default).', default='.')

	(options, args) = parser.parse_args()

	if not os.path.isdir(options.scratchdir) : proc_run_exitOnError("mkdir %s" % options.scratchdir)

	pdbfile = "%s/%s" % (options.scratchdir,os.path.basename(options.pdbfile_in))
	assert(re.compile("\.pdb$").search(pdbfile))
	rhet = 0
	if options.ReadHet == 'yes' : rhet = 1
	else : assert options.ReadHet == 'no'
	rhoh = 0
	if options.ReadHoh == 'yes' : rhoh = 1
	else : assert options.ReadHoh == 'no'
	protein(options.pdbfile_in, read_hets = rhet, read_waters = rhoh, read_hydrogens=0, remove_duplicates=1).write(pdbfile, water=1, het=1, nonHOHnonHET=1, hydrogen=0)
	proc_run_exitOnError( "chmod 644 %s" % pdbfile )

	prot = protein(pdbfile, read_hets=1, read_waters=1, read_hydrogens=0)

	os.chdir(options.scratchdir)

	if options.checkMCcont == 'yes' :
		if checkMCcont(prot) == None :
			print "Mainchain discontinuous, exiting"
			sys.exit(256)
		else : print "Mainchain continuity ok"
	else : assert (options.checkMCcont == 'no')
	if options.checkMCmissing == 'yes' :
		if prot.findMissingAtoms('mc') :
			print "Mainchain atoms missing in %s, exiting" % pdbfile
			sys.exit(256)
		else : print "All expected mainchain atoms present"
	else : assert options.checkMCmissing == 'no'
	if options.checkSCmissing != 'ignore' :
		if options.checkSCmissing == 'check' :
			if prot.findMissingAtoms('sc') :
				print "Missing sidechain atoms, exiting"
				sys.exit(0)
			else : print "All expected sidechain atoms present"
		elif options.checkSCmissing == 'repair' :
			checkRepairMissingSCatoms(prot, pdbfile, options.scratchdir)
			prot = protein(pdbfile, read_hets=1, read_waters=1, read_hydrogens=0)
		else : assert 1==0
	else : assert options.checkSCmissing == 'ignore'

	#prot.info()
	if options.useGMXsolv == 'yes' :
		# prepare gromacs files from templates
		EMMDP = options.scratchdir + '/temp_' + os.path.basename(options.EMMDP)
		prepareMDP(options.EMMDP, EMMDP, options.emTime)
		PRMDP = options.scratchdir + '/temp_' + os.path.basename(options.PRMDP)
		prepareMDP(options.PRMDP, PRMDP, options.prTime)
		solv_md = 'yes'
		if options.emTime <= 0 and options.prTime <= 0 : solv_md = 'no'
		sol_pdbfile = solvate_prot(pdbfile, EMMDP, PRMDP, options.SPCGRO, options.water_thick/10., solv_md)
		proc_run_exitOnError("rm -Rf %s %s" % (EMMDP, PRMDP))
	else : ## use Zimmer et.al. cubic grid of water
		assert options.useGMXsolv == 'no'
		HOHcrd = makeGridZimmer(prot.allcrds(), dll=options.gridDll, dpl=options.gridDpl)
		sol_pdbfile = options.scratchdir + '/' + re.sub('\.pdb$', '_sol.pdb', os.path.basename(pdbfile))
		print pdbfile, sol_pdbfile
		#proc_run_exitOnError("cp %s %s" % (pdbfile, sol_pdbfile))
		prot.write(sol_pdbfile, water=1, het=rhet, nonHOHnonHET=1, hydrogen=0)
		f = fileOpen(sol_pdbfile, 'a')
		atmcnt = 0
		for crd in HOHcrd :
			## f.write("HETATM 2168  OW    C B 308      49.518  53.649  22.616  1.00 41.06      2DNJ2347"
			resnum = atmcnt % 10000
			atomnum = atmcnt % 100000
			f.write("HETATM%5d  OW  HOH  %4d    %8.3f%8.3f%8.3f                          \n" % (atomnum,resnum,crd[0],crd[1],crd[2]) )
			atmcnt = atmcnt + 1
		f.close()
	print "DONE solvation"

	#msreader = myXmlParser()
	msreader = xml.sax.make_parser()
	mshandler = msHandler()
	msreader.setContentHandler(mshandler)
	msreader.parse(options.msXmlFile)

	prot = protein(sol_pdbfile,read_waters=1,read_hydrogens=0,read_hets=1)
	#prot.info()
	msnames, msSlaves, mschems, scrds, msSlavenames, msframes = prot.makeSites(mshandler.metasites)
	assert(len(msnames) == len(mschems))
	assert(len(msnames) == len(msSlaves))

	assert('HOH' in mschems)
	assert options.calcVols == 'yes' or options.calcVols == 'no'
	vor = Vorutil(scrds, msnames, msSlaves, mschems, msSlavenames, msframes, 'HOH', options.calcVols)
	adjResids = []
	if options.ignore_covbonded == 'yes' : 
		adjacentResidues = prot.nonHOHsampleLines()
		for i in range(len(adjacentResidues)-1) :
			if not sameChain(adjacentResidues[i], adjacentResidues[i+1]) : continue
			if isPdbHOHline(adjacentResidues[i]) or isPdbHetLine(adjacentResidues[i]) : continue
			if isPdbHOHline(adjacentResidues[i+1]) or isPdbHetLine(adjacentResidues[i+1]) : continue
			if isPdbNUline(adjacentResidues[i]) and not isPdbNUline(adjacentResidues[i+1]) : continue
			if isPdbAAline(adjacentResidues[i]) and not isPdbAAline(adjacentResidues[i+1]) : continue
			adjResids.append( ( line2resid(adjacentResidues[i]),line2resid(adjacentResidues[i+1]) ) )
	else :  assert options.ignore_covbonded == 'no'
	vor.adjResids = adjResids
	print "DONE initiated vorutil"
	picklefilename = re.sub("\.pdb", "", pdbfile) + "_" + re.sub("\.xml","",os.path.basename(options.msXmlFile)) + '.pkl'
	#picklefilename = re.sub("\.pdb", ".pkl", pdbfile)
	cPickle.dump(vor,fileOpen(picklefilename,'w'))
	print "DONE wrote pickle", picklefilename
	#cgofilename = re.sub("\.pdb", "_"+re.sub("\.xml","",os.path.basename(options.msXmlFile))+".cgo", pdbfile)
	#print cgofilename
	#vor.writeCGO(cgofilename, all_poly_collapsed = collapsed, color=siteclr)
	#print "DONE cgo", cgofilename
	datafilename = re.sub("\.pdb", "", pdbfile) + "_" + re.sub("\.xml","",os.path.basename(options.msXmlFile)) + "_" + options.datastyle + '.data'
	MSI, pats = [], []
	hohPat = re.compile('HOH')
	chPat = re.compile("\[..."+options.outChain)
	for i in range(len(vor.msnames)) :
		if hohPat.search(vor.msnames[i]) : continue
		if not chPat.search(vor.msnames[i]) : continue
		if vor.mschems[i]=="het" : continue
		MSI.append(i)
	vor.writeAreasNew(datafilename, options.datastyle, MSI)
	#vor.writeAreas(datafilename, style=options.datastyle, removeChem = "HOH")
	print "DONE writing", datafilename

if __name__=='__main__' :
	main(sys.argv)
	sys.exit(0)
	msreader = myXmlParser()
	#msreader = xml.sax.make_parser()
	mshandler = msHandler()
	msreader.setContentHandler(mshandler)
	msreader.parse(sys.argv[3])
	sys.exit(0)
