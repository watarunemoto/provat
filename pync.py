

import cPickle, os, string, re, sys, random, math, tempfile, popen2, random
from tempfile import mktemp
from Tkinter import *
import Pmw
import tkSimpleDialog, tkFileDialog, tkColorChooser, tkFont
import tkMessageBox
import xml
import xml.sax
import xml.sax.handler


## preferred colors for interactions
## DHAB interactions get red shades, RR blue shades, DHR green shades
intxColors = {}
intxColors['hbond'] = [1., 0., 0.]
intxColors['CaHbond'] = [1., 0., 0.]

intxColors['NHpi'] = [0., 1., 0.]
intxColors['CaHpi'] = [0., 1., 0.1]
intxColors['cationPi'] = [0., 1., 0.2]
intxColors['thiolPi'] = [0., 1., 0.3]
intxColors['hydroPi'] = [0.1, 1., 0.]
intxColors['proCdPi'] = [0.2, 1., 0.]

intxColors['pipi'] = [0., 0., 1.]

intxColors['hydrophobic'] = [1., 1., 1.]

class Intxn :
	def __init__(self, res1, res2, name, p1, p2) :
		self.res1, self.res2, self.name, self.p1, self.p2 = res1, res2, name, list(p1), list(p2)
class IntxOutHandler(xml.sax.handler.ContentHandler) :
	def __init__(self,verbose=None):
		self.verbose = verbose
		self.intxns = []
	def startElement(self, name, attrs) :
		if self.verbose : print "START",name, attrs.getNames()
		if name == 'pync_interactions' : pass
		elif name == 'intxn' :
			self.intxns.append( Intxn(attrs.get("res1"), attrs.get("res2"), attrs.get("name"), attrs.get("p1"), attrs.get("p2")) )
		else : assert 1==0
	def endElement(self, name) :
		if self.verbose : print "END",name
		if name == 'pync_interactions' : pass
		elif name == 'intxn' : pass
		else : assert 1==0
	def endDocument(self) : pass
	def characters(self,content) :
		if self.verbose : print "-%s-" % content

class pyncData : pass

def get_interactions(outxmlfn) :
	interactions = []
## this is what shd ideally happen, make_parser shd work ! but it dsnt, so i have to write a dirty code of my own
#	reader = xml.sax.make_parser()
#	ih = IntxOutHandler()
#	reader.setContentHandler(ih)
#	reader.parse(outxmlfn)
#	interactions = ih.intxns
## my dirty code to substitute for missing make_parser, works only for a simple limited syntax
	for l in open(outxmlfn, 'r').readlines() :
		if not re.compile('^.intxn.*intxn.$').search(l) : continue
		flds = l.split('"')
		res1, res2, name, p1, p2 = None, None, None, None, None
		for i in range(len(flds)) :
			if i % 2 == 0 : continue
			if i==1 : res1 = flds[i]
			if i==3 : res2 = flds[i]
			if i==5 : name = flds[i]
			if i==7 : p1 = [string.atof(string.split(flds[i])[0]), string.atof(string.split(flds[i])[1]), string.atof(string.split(flds[i])[2])]
			if i==9 : p2 = [string.atof(string.split(flds[i])[0]), string.atof(string.split(flds[i])[1]), string.atof(string.split(flds[i])[2])]
		interactions.append(Intxn(res1, res2, name, p1, p2))
	return interactions

def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command', 'Pync', label='Pync', command = lambda s=self : initPyncWin(s))
	#initPyncWin(self)

def addFile() :
	fn = tkFileDialog.askopenfilename(title="Add file", parent=pyncData.pyncframe, filetypes=[("Pync output", "*.xml")])
	if not fn : return
	interactions = get_interactions(fn)
	colors = {}
	for int in interactions :
		if int.name not in colors.keys() :
			if int.name in intxColors.keys() : colors[int.name] = intxColors[int.name]
			else : colors[int.name] = [random.random(), random.random(), random.random()]
	pyncData.files[os.path.basename(fn)] = (interactions, colors)
	pyncData.filelist.delete(0,END)
	for k,v in pyncData.files.items() :
		pyncData.filelist.insert(END,k)

def removeFile() :
	for i in pyncData.filelist.curselection() :
		del pyncData.files[pyncData.filelist.get(i)]
	pyncData.filelist.delete(0,END)
	for k,v in pyncData.files.items() :
		pyncData.filelist.insert(END,k)

def refreshReslist() :
	assert len(pyncData.filelist.curselection()) == 1
	fn = pyncData.filelist.get(pyncData.filelist.curselection()[0])
	resnames = {}
	for int in pyncData.files[fn][0] :
		resnames[int.res1] = 1
		resnames[int.res2] = 1
	pyncData.reslist.delete(0,END)
	for rn in resnames : pyncData.reslist.insert(END,rn)
	pyncData.intlist.delete(0,END)
	for intname,clr in pyncData.files[fn][1].items() : pyncData.intlist.insert(END,intname)
	pyncData.curfile.set(fn)

def selRes() :
	pyncData.reslist.select_clear(0,END)
	resi = pyncData.sel_resi.get()
	if len(resi) == 0 : resi = '....'
	elif len(resi) == 1 : resi = '...' + resi
	elif len(resi) == 2 : resi = '..' + resi
	elif len(resi) == 3 : resi = '.' + resi
	elif len(resi) == 2 : pass
	else : assert 1==0
	resn = pyncData.sel_resn.get()
	if len(resn) == 0 : resn = '...'
	elif len(resn) == 1 : resn = '..' + resn
	elif len(resn) == 2 : resn = '.' + resn
	elif len(resn) == 3 : pass
	else : assert 1==0
	ic = pyncData.sel_ic.get()
	if len(ic) == 0 : ic = '.'
	elif len(ic) == 1 : pass
	else : assert 1==0
	ch = pyncData.sel_ch.get()
	if len(ch) == 0 : ch = '.'
	elif len(ch) == 1 : pass
	else : assert 1==0
	pat = re.compile(resn + ch + resi + ic)
	for i in range(pyncData.reslist.size()) :
		if pat.search(pyncData.reslist.get(i)) : pyncData.reslist.select_set(i)

def addSel1() :
	for k in pyncData.reslist.curselection() :
		pyncData.resSel1.insert(END,pyncData.reslist.get(k))
def addSel2() :
	for k in pyncData.reslist.curselection() :
		pyncData.resSel2.insert(END,pyncData.reslist.get(k))
def clrSel1() : pyncData.resSel1.delete(0,END)
def clrSel2() : pyncData.resSel2.delete(0,END)

def changeColor() :
	assert len(pyncData.intlist.curselection()) == 1
	intname = pyncData.intlist.get(pyncData.intlist.curselection()[0])
	curfile = pyncData.curfile.get()
	retclr = tkColorChooser.askcolor()
	if not retclr[0] : return
	pyncData.files[curfile][1][intname] = [retclr[0][0]/256., retclr[0][1]/256., retclr[0][2]/256.]

def posRepeat(pos, res) :
	for i in range(len(pos)) :
		for k in range(i+1, len(pos)) :
			x = pos[i][0]-pos[k][0]
			y = pos[i][1]-pos[k][1]
			z = pos[i][2]-pos[k][2]
			if math.sqrt(x*x+y*y+z*z) < 1 :
				print "DIST", math.sqrt(x*x+y*y+z*z), res[i], pos[i], res[k], pos[k]
				return 1
	return None


def getAllCycles(edgeCycleBase, edges) :
	edgeCycles = [[], edgeCycleBase[0]]
	for ci in range(1,len(edgeCycleBase)) :
		bc = edgeCycleBase[ci]
		newCycles = []
		## cycles bc, c have to have a common edge to be joinable. join using xor
		for c in edgeCycles :
			resultCycle = []
			for bce in bc :
				if not bce in c : resultCycle.append(bce)
			for ce in c :
				if not ce in bc : resultCycle.append(ce)
			if len(resultCycle) == 0 : continue
			if len(c) == 0 or len(resultCycle) < len(c) + len(bc) :
				newCycles.append(resultCycle)
		edgeCycles = edgeCycles + newCycles
	edgeCycles.pop(edgeCycles.index([]))
	## remove cycles with node-degree more than 2
	newEdgeCycles = []
	for c in edgeCycles :
		degree = {}
		for ec in c :
			try : degree[edges[ec][0]] = degree[edges[ec][0]] + 1
			except : degree[edges[ec][0]] = 1
			try : degree[edges[ec][1]] = degree[edges[ec][1]] + 1
			except : degree[edges[ec][1]] = 1
		goodedge = 1
		print degree
		for i in degree.values() :
			if i > 2 :
				goodedge = 0
				break
		if goodedge == 1 : newEdgeCycles.append(c)
	return newEdgeCycles

def showCycles() :
	resnames1, resnames2, intnames = [],[],[]
	for i in range(pyncData.resSel1.size()) : resnames1.append(pyncData.resSel1.get(i))
	for i in range(pyncData.resSel2.size()) : resnames2.append(pyncData.resSel2.get(i))
	for i in pyncData.intlist.curselection() : intnames.append(pyncData.intlist.get(i))
	res, pos, edges = [], [], []
	curfile = pyncData.curfile.get()
	for int in pyncData.files[curfile][0] :
		if not int.res1 in resnames1 : continue
		if not int.res2 in resnames2 : continue
		if len(intnames) > 0 and not int.name in intnames : continue
		if not int.res1 in res :
			res.append(int.res1)
			pos.append(list(int.p1))
		if not int.res2 in res :
			res.append(int.res2)
			pos.append(list(int.p2))
		edges.append((res.index(int.res1), res.index(int.res2)))
	if posRepeat(pos, res) : assert(0)
	fn = mktemp()[1]
	print fn
	f = open(fn, 'w')
	f.write("%d\n" % len(res))
	f.write("%d\n" % len(edges))
	for i in range(len(edges)) : f.write("%d %d\n" % edges[i])
	f.close()
	c = popen2.Popen3("/tmp/horton %s" % fn)
	edgeCycleBase = []
	for l in c.fromchild.readlines() :
		if not re.compile('Edges in cycle:').search(l) : continue
		cycle = []
		flds = string.split( re.sub('Edges in cycle:', '', l) )
		for fl in flds : cycle.append(string.atoi(fl))
		edgeCycleBase.append(cycle)
	if len(edgeCycleBase) == 0 :
		print "No cycles found !!"
		return
	edgeCycles = getAllCycles(edgeCycleBase, edges)
	from pymol import cmd
	from pymol.cgo import BEGIN,LINES,END,VERTEX,COLOR
	objname = pyncData.objname.get()
	if not objname or objname == '' : objname = 'cycles'
	for c in edgeCycleBase :
		obj = [BEGIN,LINES]
		for ce in c :
			obj = obj + [VERTEX,] + list(pos[edges[ce][0]]) + [VERTEX,] + list(pos[edges[ce][1]])
		obj = obj + [END]
		cmd.load_cgo(obj, objname)
	objname = pyncData.objname.get()
	if not objname or objname == '' : objname = 'cycles'
	objname = "all_" + objname
	for c in edgeCycles :
		obj = [BEGIN,LINES]
		for ce in c :
			obj = obj + [VERTEX,] + list(pos[edges[ce][0]]) + [VERTEX,] + list(pos[edges[ce][1]])
		obj = obj + [END]
		cmd.load_cgo(obj, objname)

def showCliques() :
	resnames1, resnames2, intnames = [],[],[]
	for i in range(pyncData.resSel1.size()) : resnames1.append(pyncData.resSel1.get(i))
	for i in range(pyncData.resSel2.size()) : resnames2.append(pyncData.resSel2.get(i))
	for i in pyncData.intlist.curselection() : intnames.append(pyncData.intlist.get(i))
	res, pos, edges = [], [], []
	curfile = pyncData.curfile.get()
	for int in pyncData.files[curfile][0] :
		if not int.res1 in resnames1 : continue
		if not int.res2 in resnames2 : continue
		if len(intnames) > 0 and not int.name in intnames : continue
		print int.res1, int.res2, int.p1, int.p2
		if not int.res1 in res :
			res.append(int.res1)
			pos.append(list(int.p1))
		if not int.res2 in res :
			res.append(int.res2)
			pos.append(list(int.p2))
		edges.append((res.index(int.res1), res.index(int.res2)))
	if posRepeat(pos, res) : assert(0)
	mat = []
	for i in range(len(res)) :
		x = []
		for i in range(len(res)) : x.append(0)
		mat.append(x)
	for e in edges :
		mat[e[0]][e[1]] = 1
		mat[e[1]][e[0]] = 1
	fn = mktemp()[1]
	print fn
	f = open(fn, 'w')
	f.write("%d\n" % len(mat[0]))
	for i in range(len(res)) :
		for k in range(i+1, len(res)) :
			f.write("%d " % mat[i][k])
		f.write("\n")
	f.close()
	c = popen2.Popen3("/tmp/maxclq %s 0.5" % fn)
	clqs = []
	for l in c.fromchild.readlines() :
		if not re.compile("^Maximal clique:").search(l) : continue
		l = re.sub("Maximal clique:", "", l)
		flds = l.split()
		if len(flds) < 3 : continue
		print "FLDS", flds
		clq = []
		for fl in flds : clq.append(string.atoi(fl)-1)
		clqs.append(clq)
	#os.unlink(fn)
	print "CLQS", clqs
	import random
	from pymol import cmd
	from pymol.cgo import BEGIN,LINES,END,VERTEX,COLOR
	objname = pyncData.objname.get()
	if not objname or objname == '' : objname = 'clq'
	for i in range(len(pos)) : print "POS", res[i], pos[i]
	for clq in clqs :
		print clq
		obj = [BEGIN,LINES,]
		for i in range(len(clq)) :
			for k in range(i+1,len(clq)) :
				#print "LINE BETN", clq, clq[i], clq[k],  pos[ clq[i] ], pos[ clq[k] ]
				obj = obj + [VERTEX,] + pos[ clq[i] ] + [VERTEX,] + pos[ clq[k] ]
		obj = obj + [END]
		cmd.load_cgo(obj, objname)
#############
#	curfile = pyncData.curfile.get()
#	obj = [BEGIN,LINES,]
#	for int in pyncData.files[curfile][0] :
#		if not int.res1 in resnames1 : continue
#		if not int.res2 in resnames2 : continue
#		if len(intnames) > 0 and not int.name in intnames : continue
#		obj = obj + [COLOR,] + pyncData.files[curfile][1][int.name] + [VERTEX,] + int.p1 + [VERTEX,] + int.p2

def showInt() :
	resnames1, resnames2, intnames = [],[],[]
	for i in range(pyncData.resSel1.size()) : resnames1.append(pyncData.resSel1.get(i))
	for i in range(pyncData.resSel2.size()) : resnames2.append(pyncData.resSel2.get(i))
	for i in pyncData.intlist.curselection() : intnames.append(pyncData.intlist.get(i))
	import random
	from pymol import cmd
	from pymol.cgo import BEGIN,LINES,END,VERTEX,COLOR,CYLINDER
	curfile = pyncData.curfile.get()
	obj = [] # [BEGIN,LINES,]
	for int in pyncData.files[curfile][0] :
		if not int.res1 in resnames1 : continue
		if not int.res2 in resnames2 : continue
		if len(intnames) > 0 and not int.name in intnames : continue
		obj = obj + [CYLINDER,] + int.p1 + int.p2 + [0.5] + pyncData.files[curfile][1][int.name] + pyncData.files[curfile][1][int.name]
		#obj = obj + [COLOR,] + pyncData.files[curfile][1][int.name] + [VERTEX,] + int.p1 + [VERTEX,] + int.p2
	#obj = obj + [END]
	objname = pyncData.objname.get()
	if not objname or objname == '' : objname = 'obj'
	cmd.load_cgo(obj, objname)

def initPyncWin(self) :
	pyncData.pyncframe = Toplevel(self.root)
	pyncData.files = {}
	
	fr = Frame(pyncData.pyncframe)
	vertscroll = Scrollbar(fr, orient=VERTICAL)
	horscroll = Scrollbar(fr, orient=HORIZONTAL)
	pyncData.filelist = Listbox(fr, selectmode=EXTENDED, xscrollcommand=horscroll.set, yscrollcommand=vertscroll.set, font=tkFont.Font(family='fixedsys',size=12))
	horscroll.config(command=pyncData.filelist.xview)
	vertscroll.config(command=pyncData.filelist.yview)
	vertscroll.grid(row=50,column=2,sticky=N+S+E)
	pyncData.filelist.grid(row=50,column=1,sticky=N+S+W)
	horscroll.grid(row=51,column=1,sticky=N+W+E)
	fr.grid(row=0,column=0)
	Button(pyncData.pyncframe, text="Add File", command=lambda:addFile()).grid(row=1,column=0, sticky=E+W)
	Button(pyncData.pyncframe, text="Remove File", command=lambda:removeFile()).grid(row=2,column=0, sticky=E+W)
	pyncData.curfile = StringVar()
	Label(pyncData.pyncframe, textvariable=pyncData.curfile, background='white').grid(row=3,column=0,sticky=E+W)

	fr = Frame(pyncData.pyncframe)
	vertscroll = Scrollbar(fr, orient=VERTICAL)
	horscroll = Scrollbar(fr, orient=HORIZONTAL)
	pyncData.reslist = Listbox(fr, selectmode=EXTENDED, xscrollcommand=horscroll.set, yscrollcommand=vertscroll.set, font=tkFont.Font(family='fixedsys',size=12))
	horscroll.config(command=pyncData.reslist.xview)
	vertscroll.config(command=pyncData.reslist.yview)
	vertscroll.grid(row=50,column=2,sticky=N+S+E)
	pyncData.reslist.grid(row=50,column=1,sticky=N+S+W)
	horscroll.grid(row=51,column=1,sticky=N+W+E)
	fr.grid(row=0,column=10)
	Button(pyncData.pyncframe, text="Refresh", command=lambda:refreshReslist()).grid(row=1,column=10, sticky=E+W)
	fr = Frame(pyncData.pyncframe)
	pyncData.sel_resn = Entry(fr,width=3)
	pyncData.sel_resn.grid(row=0, column=0)
	pyncData.sel_ch = Entry(fr,width=1)
	pyncData.sel_ch.grid(row=0, column=1)
	pyncData.sel_resi = Entry(fr,width=4)
	pyncData.sel_resi.grid(row=0, column=2)
	pyncData.sel_ic = Entry(fr,width=1)
	pyncData.sel_ic.grid(row=0, column=3)
	fr.grid(row=2,column=10, sticky=E+W)
	Button(pyncData.pyncframe, text="Select Residues", command=lambda:selRes()).grid(row=3, column=10, sticky=E+W)

	fr = Frame(pyncData.pyncframe)
	vertscroll = Scrollbar(fr, orient=VERTICAL)
	horscroll = Scrollbar(fr, orient=HORIZONTAL)
	pyncData.resSel1 = Listbox(fr, selectmode=EXTENDED, xscrollcommand=horscroll.set, yscrollcommand=vertscroll.set, font=tkFont.Font(family='fixedsys',size=12))
	horscroll.config(command=pyncData.resSel1.xview)
	vertscroll.config(command=pyncData.resSel1.yview)
	vertscroll.grid(row=50,column=2,sticky=N+S+E)
	pyncData.resSel1.grid(row=50,column=1,sticky=N+S+W)
	horscroll.grid(row=51,column=1,sticky=N+W+E)
	fr.grid(row=0,column=20)
	Button(pyncData.pyncframe, text="Add", command=lambda:addSel1()).grid(row=1, column=20, sticky=E+W)
	Button(pyncData.pyncframe, text="Clear", command=lambda:clrSel1()).grid(row=2, column=20, sticky=E+W)

	fr = Frame(pyncData.pyncframe)
	vertscroll = Scrollbar(fr, orient=VERTICAL)
	horscroll = Scrollbar(fr, orient=HORIZONTAL)
	pyncData.resSel2 = Listbox(fr, selectmode=EXTENDED, xscrollcommand=horscroll.set, yscrollcommand=vertscroll.set, font=tkFont.Font(family='fixedsys',size=12))
	horscroll.config(command=pyncData.resSel2.xview)
	vertscroll.config(command=pyncData.resSel2.yview)
	vertscroll.grid(row=50,column=2,sticky=N+S+E)
	pyncData.resSel2.grid(row=50,column=1,sticky=N+S+W)
	horscroll.grid(row=51,column=1,sticky=N+W+E)
	fr.grid(row=0,column=30)
	Button(pyncData.pyncframe, text="Add", command=lambda:addSel2()).grid(row=1, column=30, sticky=E+W)
	Button(pyncData.pyncframe, text="Clear", command=lambda:clrSel2()).grid(row=2, column=30, sticky=E+W)

	fr = Frame(pyncData.pyncframe)
	vertscroll = Scrollbar(fr, orient=VERTICAL)
	horscroll = Scrollbar(fr, orient=HORIZONTAL)
	pyncData.intlist = Listbox(fr, selectmode=EXTENDED, xscrollcommand=horscroll.set, yscrollcommand=vertscroll.set, font=tkFont.Font(family='fixedsys',size=12))
	horscroll.config(command=pyncData.intlist.xview)
	vertscroll.config(command=pyncData.intlist.yview)
	vertscroll.grid(row=50,column=2,sticky=N+S+E)
	pyncData.intlist.grid(row=50,column=1,sticky=N+S+W)
	horscroll.grid(row=51,column=1,sticky=N+W+E)
	fr.grid(row=0,column=40)
	Button(pyncData.pyncframe, text="Change Color", command=lambda:changeColor()).grid(row=1, column=40, sticky=E+W)
	fr = Frame(pyncData.pyncframe)
	Button(fr, text="Show", command=lambda:showInt()).grid(row=0, column=0, sticky=E+W)
	pyncData.objname = Entry(fr, width=5)
	pyncData.objname.grid(row=0, column=1, sticky=E+W)
	Label(fr, text="obj.name").grid(row=0, column=2, sticky=E+W)
	fr.grid(row=2, column=40, sticky=E+W)
	Button(pyncData.pyncframe, text="Show Cliques", command=lambda:showCliques()).grid(row=3, column=40, sticky=E+W)
	Button(pyncData.pyncframe, text="Show Cycles", command=lambda:showCycles()).grid(row=4, column=40, sticky=E+W)

