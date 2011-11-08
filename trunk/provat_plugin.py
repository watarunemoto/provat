#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################


import cPickle, os, string, re, sys, random, math
from pymol.cgo import *
from pymol import cmd
from Tkinter import *
import tkSimpleDialog, tkFileDialog, tkColorChooser, tkMessageBox

#predefine some standard colors
chem2rgb = {}
chem2rgb['mc'] = (0,0.99,0)
chem2rgb['peptide'] = (0,0.99,0)
chem2rgb['sc'] = (0,0.99,0.99)
chem2rgb['HOH'] = (0,0,0.99)
chem2rgb['phos'] = (0.99,0.9,0.5)
chem2rgb['sugar'] = (0.2,0.35,0.35)
chem2rgb['base'] = (0.625,0.5,0.9375) 
chem2rgb['het'] = (0.99,0,0.79)
chem2rgb['phob'] = (0.9,0.9,0.9)
chem2rgb['phil'] = (0.99,0,0) 
chem2rgb['arom'] = (0.99,0.99,0)

########### command line interface ##########

vorshows = {}

def newvs(filename, prname='') :
    if prname == '' : prname = re.sub('\..*$', '', os.path.basename(filename))
    assert(prname not in vorshows.keys())
    v = cPickle.load(open(filename,'r'))
    v.alpha = 0.5
    vorshows[prname] = v

def getPrname(prname) :
    if prname == '' :
        if len(vorshows.keys()) > 0 :
            prname = vorshows.keys()[0]
        else :
            print "protein name (prname) necessary"
            assert(0)
    assert(prname in vorshows.keys())
    return prname

def showvs(prname='', resn='', chid='', resi='', ic='', mets='', atmn='', style='fused',name='voro') :
    vorshows[getPrname(prname)].showMS(resn, chid, resi, ic, mets, atmn, style, name)

def makesurf(prname='', name='surf', MSIs = []) :
    vorshows[getPrname(prname)].showSurf(name, MSIs)

def setConwidth(cw, prname='') :
    vorshows[getPrname(prname)].conwidth = cw

def setAlpha(a, prname='') :
    vorshows[getPrname(prname)].alpha = a

def remvs(prname='') :
    if getPrname(prname) == '' : return
    del vorshows[getPrname(prname)]

#############################################

class provatData :
	def __init__() : pass

def __init__(self):
	print "initializing ProVAT viewer..."
	provatpath = ''
	if "PROVATPATH" not in os.environ.keys() :
		provatpath = tkFileDialog.askdirectory(title = 'Provide path to ProVAT distribution', parent=self.root)
		if not provatpath :
			print "Cannot initialize ProVAT :-("
			return
	else : provatpath = os.environ["PROVATPATH"]

	sys.path.append(provatpath)
	print 'SYS.PATH', sys.path
	from vorutil import Vorutil

	cmd.extend('remvs',remvs)
	cmd.extend('newvs',newvs)
	cmd.extend('setAlpha',setAlpha)
	cmd.extend('showvs',showvs)
	cmd.extend('makesurf', makesurf)

	print "initialized ProVAT successfully, enjoy !!"

	self.menuBar.addmenuitem('Plugin', 'command', 'ProVAT', label='ProVAT', command = lambda s=self : initProvatWin(s))
	provatData.provatframe = None
	#initProvatWin(self)

def addnewProt(self) :
	pklpath = tkFileDialog.askopenfilename(title = 'Choose pickle file', parent=self.root, filetypes=[("Python pickle files", "*.pkl")])
	if not pklpath : return
	protname = re.sub('\.[^\.]*', '', os.path.basename(pklpath))
	if protname in vorshows.keys() :
		protname = tkSimpleDialog.askstring('ProVAT', 'Choose new name for protein because %s is already used' % protname, parent=self.root)
	newvs(filename = pklpath, prname = protname)
	updateProtlist(self)
	assert len(vorshows[protname].msnames) == len(vorshows[protname].mschems)
	for msc in vorshows[protname].mschems :
		if msc in chem2rgb.keys() :
			vorshows[protname].chem2clr[msc] = chem2rgb[msc]
		else :
			vorshows[protname].chem2clr[msc] = ( random.random(), random.random(), random.random() )
	print vorshows[protname].chem2clr

def updateProtlist(self) :
	provatData.protlist.delete(0,END)
	for name in vorshows.keys() : provatData.protlist.insert(END,name)

def updateMslist(self) :
	if len(provatData.protlist.curselection()) != 1 : return
	provatData.curprname.set(provatData.protlist.get ( provatData.protlist.curselection()[0] ))
	provatData.mslist.delete(0,END)
	for msname in vorshows[ provatData.curprname.get() ].msnames :
		if re.compile('HOH').search(msname) : continue
		provatData.mslist.insert(END,msname)

def scrollListbox(master, selmode, label) :
	vertscroll = Scrollbar(master, orient=VERTICAL)
	vertscroll.grid(row=1, column=0, sticky=N+S)

	listlabel = Label(master, text=label)
	listlabel.grid(row=0,column=1)
	horscroll = Scrollbar(master, orient=HORIZONTAL)
	horscroll.grid(row=2, column=1, sticky=E+W)
	listbox = Listbox(master, selectmode=selmode, yscrollcommand=vertscroll.set, xscrollcommand=horscroll.set)
	listbox.grid(row=1,column=1,sticky=E+W)

	horscroll.config(command=listbox.xview)
	vertscroll.config(command=listbox.yview)

	return listbox, master

def showConn(self) :
	from vorutil import areGroupsCovbonded
	from pymol import cmd
	name = ''
	if provatData.objname.get() != '' : name = provatData.objname.get()
	else : name = 'conn'
	prn = getPrname(provatData.curprname.get())
	msis = getSelListMSIs(self)
	for mi in range(len(msis)) :
		msi0 = msis[mi]
		for mi1 in range(mi+1, len(msis)) :
			msi1 = msis[mi1]
			assert(len(vorshows[prn].adjResids) > 0)
			if provatData.ignCovbondedNbrs.get() == 0 and areGroupsCovbonded(vorshows[prn].slavenames[msi0],vorshows[prn].slavenames[msi1], vorshows[prn].adjResids) :
				print "COVBONDED"
				print vorshows[prn].slavenames[msi0]
				print vorshows[prn].slavenames[msi1]
				continue
			if msi1 not in vorshows[prn].msadj[msi0].keys() : continue
			vorshows[prn].showLine(msi0, msi1, name)
			cmd.set('all_states', 1, name)
	if provatData.writeObjToFile.get() != 0 : vorshows[prn].writeAreasNew(provatData.datafilename.get(), 'inner_nbrs', msis)

def showCommonSurf(self) :
	from vorutil import areGroupsCovbonded
	prn = getPrname(provatData.curprname.get())
	msis0 = getSelListMSIs(self)
	msis1 = getXSelListMSIs(self)
	surfname = ''
	if provatData.objname.get() == '' : surfname = 'Xsurf'
	else : surfname = provatData.objname.get()
	nbrs = []
	print "HERE", provatData.ignCovbondedNbrs.get()
	if provatData.ignCovbondedNbrs.get() == 0 : print 'doing covalency check, this may take time'
	for msi0 in msis0 :
		for msi1 in msis1 :
			if msi0 not in vorshows[prn].msadj.keys() : continue
			if msi1 not in vorshows[prn].msadj[msi0].keys() : continue
			if provatData.ignCovbondedNbrs.get() == 0 and areGroupsCovbonded(vorshows[prn].slavenames[msi0],vorshows[prn].slavenames[msi1], vorshows[prn].adjResids) : continue
			nbrs.append((msi0,msi1))
	totarea = vorshows[provatData.curprname.get()].showIntxn(nbrs, surfname)
	provatData.contAreaBox.set("%f" % totarea)
	if provatData.writeObjToFile.get() != 0 : vorshows[prn].writeAreasNew(provatData.datafilename.get(), 'ext_nbrs', msis0, msis1)

def showCrossConn(self) :
	from vorutil import areGroupsCovbonded
	name = ''
	if provatData.objname.get() != '' : name = provatData.objname.get()
	else : name = 'conn'
	prn, msis, msis1 = getPrname(provatData.curprname.get()), [], []
	msis = getSelListMSIs(self)
	msis1 = getXSelListMSIs(self)
	for msi0 in msis :
		for msi1 in msis1 :
			if provatData.ignCovbondedNbrs.get() == 0 and areGroupsCovbonded(vorshows[prn].slavenames[msi0],vorshows[prn].slavenames[msi1], vorshows[prn].adjResids) : continue
			#if msi0 not in vorshows[prn].msadj.keys() :
			#print 'here nums', msi0, msi1, len(msis), len(msis1), len(vorshows[prn].msnames)
			#print 'here', vorshows[prn].msnames[msi0], vorshows[prn].msnames[msi1]
			if msi1 not in vorshows[prn].msadj[msi0].keys() : continue
			vorshows[prn].showLine(msi0, msi1, name)
			cmd.set('all_states', 1, name)
	if provatData.writeObjToFile.get() != 0 : vorshows[prn].writeAreasNew(provatData.datafilename.get(), 'ext_nbrs', msis, msis1)

def getXSelListMSIs(self) :
	prn = getPrname(provatData.curprname.get())
	msis = []
	for i in provatData.xsellist.get(0,END) :
		if re.compile('HOH').search(i) : continue
		print i
		k = vorshows[prn].msnames.index(i)
		msis.append(k)
	return msis
def getSelListMSIs(self) :
	prn = getPrname(provatData.curprname.get())
	msis = []
	for i in provatData.sellist.get(0,END) :
		if re.compile('HOH').search(i) : continue
		print i
		k = vorshows[prn].msnames.index(i)
		msis.append(k)
	return msis

def showNbr(self) :
	prn = getPrname(provatData.curprname.get())
	style, name = 'fused', 'nbr'
	if provatData.isMSstyleFused.get() != 0 : #selected
		style = 'fused'
		name = provatData.objname.get()
		if name == '' :
			print "Using default name 'nbr' for the fused object because none is given"
			name = 'nbr'
	else : style = 'distinct'
	msis = getSelListMSIs(self)
	for k in msis :
		for ni in vorshows[prn].msadj[k].keys() :
			if re.compile('HOH').search(vorshows[prn].msnames[ni]) : continue
			if style == 'fused' : vorshows[prn].showvs(ni,name)
			else : vorshows[prn].showvs(ni, vorshows[prn].msnames[ni])
	if provatData.writeObjToFile.get() != 0 : vorshows[prn].writeAreasNew(provatData.datafilename.get(), 'common_nbrlist', msis)


def showSurface(self) :
	prn = getPrname(provatData.curprname.get())
	surfname = ''
	if provatData.objname.get() == '' : surfname = 'surf'
	else : surfname = provatData.objname.get()
	msis = getSelListMSIs(self)
	if len(msis) == 0 : return
	totarea = vorshows[prn].showSurf(surfname, msis)
	provatData.areabox.set("%f" % totarea)
	if provatData.writeObjToFile.get() != 0 : vorshows[prn].writeAreasNew(provatData.datafilename.get(), 'exposed_surface', msis)

def showPolyhedra(self) :
	prn = getPrname(provatData.curprname.get())
	style, name = 'fused', 'voro'
	if provatData.isMSstyleFused.get() != 0 : #selected
		style = 'fused'
		name = provatData.objname.get()
		if name == '' :
			print "Using default name 'voro' for the fused object because none is given"
			name = 'voro'
	else : style = 'distinct'
	totvol = 0.
	msis = getSelListMSIs(self)
	for k in msis :
		for si in vorshows[prn].master2slave[k] :
			totvol = totvol + vorshows[prn].sitevols[si]
		if style == 'fused' : vorshows[prn].showvs(k,name)
		else : vorshows[prn].showvs(k, vorshows[prn].msnames[k])
	provatData.volbox.set("%f" % totvol)
	if provatData.writeObjToFile.get() != 0 : vorshows[prn].writeAreasNew(provatData.datafilename.get(), 'all_nbr', msis)

def highlightMetasitesRE(self) :
	provatData.mslist.select_clear(0,END)
	print provatData.entry_resn.get(),':', provatData.entry_chid.get(),':', provatData.entry_resi.get(),':', provatData.entry_ic.get(),':', provatData.entry_atmn.get(),':', provatData.entry_mets.get()
	msis = vorshows[provatData.curprname.get()].whichMSindices(provatData.entry_resn.get(), provatData.entry_chid.get(), provatData.entry_resi.get(), provatData.entry_ic.get(), provatData.entry_atmn.get(),
				provatData.entry_mets.get(), provatData.entry_msc.get())
	for i in range(provatData.mslist.size()) :
		msi = vorshows[provatData.curprname.get()].msnames.index( provatData.mslist.get(i) )
		if msi in msis :
			print "Highlighting", provatData.mslist.get(i), vorshows[provatData.curprname.get()].msnames[msi]
			provatData.mslist.select_set(i)

def removeProtein(self) :
	selprot = provatData.protlist.curselection()
	assert(len(selprot) <= 1)
	print provatData.protlist.get(selprot[0])
	remvs(provatData.protlist.get(selprot[0]))
	updateProtlist(self)

def resetConwidth(self) :
	try :
		if provatData.conwid_val.get() == '' : cw = 0.1
		else : cw = string.atof(provatData.conwid_val.get())
		setConwidth(cw, provatData.curprname.get())
	except ValueError :
		print "Connection width invalid, enter a numeric value"
		tkMessageBox.showwarning(message="Connection width invalid, enter a numeric value")

def resetAlpha(self) :
	try :
		alpha = string.atof(provatData.alpha_val.get())
		if alpha < 0 or alpha > 1 : raise ValueError
		setAlpha(alpha, provatData.curprname.get())
	except ValueError :
		print "Alpha should be between 0 and 1"
		tkMessageBox.showwarning(message='Enter transparence between 0 and 1')

def chooseClr(button, protname, chem) :
	retclr = tkColorChooser.askcolor()
	if not retclr[0] : return
	clr = retclr[0]
	button.config(foreground=retclr[1])
	vorshows[protname].chem2clr[chem] = ( clr[0]/256., clr[1]/256., clr[2]/256. )
	print vorshows[protname].chem2clr

def adjustColors(self) :
	win = Toplevel(provatData.provatframe)
	curpr = getPrname(provatData.curprname.get())
	print vorshows[curpr].chem2clr.keys()
	for chem,clr in vorshows[curpr].chem2clr.items() :
		Rc = hex( int(math.floor(clr[0]*16)) )[2:]
		Gc = hex( int(math.floor(clr[1]*16)) )[2:]
		Bc = hex( int(math.floor(clr[2]*16)) )[2:]
		assert len(Rc) == 1 and len(Gc) == 1 and len(Bc) == 1
		b = Button(win, text='Color for ' + chem, foreground = '#' + Rc + Gc + Bc, background='white')
		b.config(command = lambda but=b, prn=curpr, c=chem : chooseClr(but, prn, c))
		b.grid(column=0, sticky=E+W)

# just append selection in metasite lisbox to selection listbox
def add2select(self) :
	for i in provatData.mslist.curselection() :
		provatData.sellist.insert(END, provatData.mslist.get(i))
def add2Xselect(self) :
	for i in provatData.mslist.curselection() :
		provatData.xsellist.insert(END, provatData.mslist.get(i))
# remove all entris from selection listbox
def clearSelect(self) :
	provatData.sellist.delete(0,END)
def clearXSelect(self) :
	provatData.xsellist.delete(0,END)

def initProvatWin(self) :
	print "Thanks for attempting to use ProVAT"
	if provatData.provatframe != None :
		print "ProVAT window is already initialized"
		tkMessageBox.showinfo(title="ProVAT", message="ProVAT window is already initialized", parent=self.root)
		return
	print "Making ProVAT window"
	provatData.provatframe = Toplevel(self.root)

	tfr = Frame(provatData.provatframe)
	provatData.protlist, prlistfr = scrollListbox(master=Frame(tfr), selmode=BROWSE, label='Protein List')
	prlistfr.grid(row=0,column=0)
	Button(tfr, text="Add New Protein", command= lambda s=self : addnewProt(s)).grid(column=0,sticky=E+W)
	Button(tfr, text="Remove Protein", command= lambda s=self : removeProtein(s)).grid(column=0,sticky=E+W)
	Button(tfr, text="Refresh Protein List", command= lambda s=self : updateProtlist(s)).grid(column=0,sticky=E+W)
	curPrFrame = Frame(tfr)
	Label(curPrFrame, text = 'Current Protein :').grid(row=0,column=0,sticky=E)
	provatData.curprname = StringVar()
	Label(curPrFrame, textvariable=provatData.curprname, background='white').grid(row=0,column=1,sticky=W)
	curPrFrame.grid(column=0,sticky=E+W)
	tfr.grid(row=0, column=0,sticky=E+W+N)
	
	tfr = Frame(provatData.provatframe)
	provatData.mslist, mslistfr = scrollListbox(master=Frame(tfr), selmode=EXTENDED, label='Metasite List')
	mslistfr.grid(row=0, column=1)
	Button(tfr, text="Refresh Metasites", command= lambda s=self : updateMslist(s)).grid(row=1,column=1,sticky=E+W)
	refr = Frame(tfr)
	Label(refr, text="resn").grid(row=0,column=0)#
	provatData.entry_resn = Entry(refr, width=3)
	provatData.entry_resn.grid(row=1,column=0)
	Label(refr, text="ch").grid(row=0,column=1)#
	provatData.entry_chid = Entry(refr, width=1)
	provatData.entry_chid.grid(row=1,column=1)
	Label(refr, text="resi").grid(row=0,column=2)#
	provatData.entry_resi = Entry(refr, width=4)
	provatData.entry_resi.grid(row=1,column=2)
	Label(refr, text="ic").grid(row=0,column=3)#
	provatData.entry_ic = Entry(refr, width=1)
	provatData.entry_ic.grid(row=1,column=3)
	Label(refr, text="atmn").grid(row=0,column=4)#
	provatData.entry_atmn = Entry(refr, width=4)
	provatData.entry_atmn.grid(row=1,column=4)
	Label(refr, text="mets").grid(row=0,column=5)#
	provatData.entry_mets = Entry(refr, width=6)
	provatData.entry_mets.grid(row=1,column=5)
	Label(refr, text="msc").grid(row=0,column=6)#
	provatData.entry_msc = Entry(refr, width=5)
	provatData.entry_msc.grid(row=1,column=6)
	refr.grid(row=2,column=1)
	Button(tfr, text="Highlight Metasites", command= lambda s=self : highlightMetasitesRE(s)).grid(row=3,column=1, sticky=E+W)
	tfr.grid(row=0, column=1,sticky=E+W+N)

	tfr = Frame(provatData.provatframe)
	provatData.sellist, sellistfr = scrollListbox(master=Frame(tfr), selmode=EXTENDED, label='Primary Selection')
	sellistfr.grid(row=0, column=0, sticky=E+W, columnspan=2)
	Button(tfr, text="Add to selection", command= lambda s=self : add2select(s)).grid(row=1,column=0, sticky=E+W)
	Button(tfr, text="Clear", command= lambda s=self : clearSelect(s)).grid(row=1,column=1, sticky=E+W)
	Button(tfr, text="Polyhedra", command= lambda s=self : showPolyhedra(s)).grid(row=2,column=0, sticky=E+W)
	provatData.volbox = StringVar()
	Label(tfr, textvariable=provatData.volbox, background='white').grid(row=2, column=1, sticky=E+W)
	Button(tfr, text="Surface", command= lambda s=self : showSurface(s)).grid(row=3,column=0, sticky=E+W)
	provatData.areabox = StringVar()
	Label(tfr, textvariable=provatData.areabox, background='white').grid(row=3, column=1, sticky=E+W)
	Button(tfr, text="Nbr Polyhedra", command= lambda s=self : showNbr(s)).grid(row=4,column=0, sticky=E+W)
	Button(tfr, text="Connections", command= lambda s=self : showConn(s)).grid(row=4,column=1, sticky=E+W)
	sellistfr = Frame(tfr)
	tfr.grid(row=0, column=2,sticky=E+W+N)

	tfr = Frame(provatData.provatframe)
	provatData.xsellist, xsellistfr = scrollListbox(master=Frame(tfr), selmode=EXTENDED, label='Additional Selection')
	xsellistfr.grid(row=0, column=0, sticky=E+W, columnspan=2)
	Button(tfr, text="Add to selection", command= lambda s=self : add2Xselect(s)).grid(row=1, column=0)
	Button(tfr, text="Clear", command= lambda s=self : clearXSelect(s)).grid(row=1, column=1)
	Button(tfr, text='Contact Surface', command= lambda s=self : showCommonSurf(s)).grid(row=2,column=0,sticky=E+W)
	provatData.contAreaBox = StringVar()
	Label(tfr, textvariable=provatData.contAreaBox, background='white', bd=1).grid(row=2, column=1, sticky=E+W)
	Button(tfr, text="Cross Connections", command= lambda s=self : showCrossConn(s)).grid(row=3,column=0,columnspan=2,sticky=E+W)
	provatData.ignCovbondedNbrs = IntVar()
	Checkbutton(tfr, text="ignore cov.bonded metasite nbrs?", variable=provatData.ignCovbondedNbrs).grid(row=4, column=0, columnspan=2, sticky=E+W+N)
	tfr.grid(row=0, column=3,sticky=E+W+N)



	tfr = Frame(provatData.provatframe)

	provatData.conwid_val = Entry(tfr, width=5)
	provatData.conwid_val.grid(row=0,column=0,sticky=E+W)
	Label(tfr, text = '(default .1)').grid(row=0,column=1,sticky=E+W)
	Button(tfr, text="Reset Connection Width", command= lambda s=self : resetConwidth(s)).grid(row=10, column=0, columnspan=2, sticky=E+W)

	Label(tfr).grid(row=11,column=0)

	provatData.alpha_val = Entry(tfr, width=5)
	provatData.alpha_val.grid(row=20, column=0,sticky=E+W)
	Label(tfr, text = '(within 0-1)').grid(row=20,column=1,sticky=E+W)
	Button(tfr, text="Reset Transparence", command= lambda s=self : resetAlpha(s)).grid(row=30, column=0, columnspan=2,sticky=E+W)

	Label(tfr).grid(row=31,column=0)

	provatData.isMSstyleFused = IntVar()
	Checkbutton(tfr, text="fused object?", variable=provatData.isMSstyleFused).grid(row=40,column=0, columnspan=2, sticky=E+W)
	Label(tfr, text='Object Name').grid(row=50,column=0,sticky=E+W)
	provatData.objname = Entry(tfr, width=10)
	provatData.objname.grid(row=50,column=1,sticky=E+W)

	Label(tfr).grid(row=51,column=0)

	provatData.writeObjToFile = IntVar()
	Checkbutton(tfr, text="write object to file?", variable=provatData.writeObjToFile).grid(row=60,column=0, columnspan=2, sticky=E+W)
	Button(tfr, text='Choose file', command= lambda s=self : chooseDatafile(s)).grid(row=61, column=0, columnspan=2, sticky=E+W)
	provatData.datafilename = StringVar()
	Label(tfr, textvariable=provatData.datafilename, background='white').grid(row=62,column=0,columnspan=2,sticky=E+W)

	Label(tfr).grid(row=69,column=0)
	Button(tfr, text='Adjust Colors', command= lambda s=self : adjustColors(s)).grid(row=70,column=0, columnspan=2,sticky=E+W)

	tfr.grid(row=0, column=4)

def chooseDatafile(self) :
	provatData.datafilename.set('')
	dfname = tkFileDialog.asksaveasfilename(title = 'Choose datafile', parent=provatData.root)
	print dfname
	if dfname : provatData.datafilename.set(dfname)
