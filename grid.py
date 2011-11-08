#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################


def surroundingGridPoints(gridmin, dgrid, pt) :
	base = (0,0,0)
	base[0] = math.floor( (pt[0] - gridmin[0]) / dgrid ) * dgrid + gridmin[0]
	base[1] = math.floor( (pt[1] - gridmin[1]) / dgrid ) * dgrid + gridmin[1]
	base[2] = math.floor( (pt[2] - gridmin[2]) / dgrid ) * dgrid + gridmin[2]
	surr = []
	surr.append( (base[0], base[1], base[2]) )
	surr.append( (base[0], base[1], base[2] + dgrid) )
	surr.append( (base[0], base[1] + dgrid, base[2]) )
	surr.append( (base[0], base[1] + dgrid, base[2] + dgrid) )
	surr.append( (base[0] + dgrid, base[1], base[2]) )
	surr.append( (base[0] + dgrid, base[1], base[2] + dgrid) )
	surr.append( (base[0] + dgrid, base[1] + dgrid, base[2]) )
	surr.append( (base[0] + dgrid, base[1] + dgrid, base[2] + dgrid) )
	return surr

def closestGridPoint(gridmin, dgrid, pt) :
	surr = surroundingGridPoints(gridmin, dgrid, pt)
	mindist = 100000000
	closest = (0,0,0)
	for s in surr :
		if vec_dist(s,pt) < mindist : closest = s
	return closest

def cubeSurf(dpl) :
	cs_dpl = {}
	for i in range(-1*dpl,dpl+1) :
		for j in range(-1*dpl,dpl+1) :
			cs_dpl[(i,j,dpl)] = 1
			cs_dpl[(i,j,-1*dpl)] = 1
			cs_dpl[(i,dpl,j)] = 1
			cs_dpl[(i,-1*dpl,j)] = 1
			cs_dpl[(dpl,i,j)] = 1
			cs_dpl[(-1*dpl,i,j)] = 1
	return list(cs_dpl.keys())
	cubesurf_dpl = []
	for i in range(-1*dpl,dpl+1) :
		for j in range(-1*dpl,dpl+1) :
			for k in range(-1*dpl,dpl+1) :
				if i*i + j+j + k*k >= dpl*dpl : cubesurf_dpl.append((i,j,k))

def makeGridZimmer(pcrds, dpl=4, dll=3) :
	# snap all pcrds to grid points and make a no-solvent zone around them for radius < dpl A
	Psites, noSsites = {}, {}
	rad_ltdpl_sph = [] ## which points around a grid-pt are less than 3A from it ?
	for i in range(-1*dpl,dpl+1) :
		for j in range(-1*dpl,dpl+1) :
			for k in range(-1*dpl,dpl+1) :
				if i*i + j+j + k*k < dpl*dpl : rad_ltdpl_sph.append((i,j,k))
	for crd in pcrds :
		newPsite = (round(crd[0]), round(crd[1]), round(crd[2]))
		Psites[ newPsite ] = 1
		for offset in rad_ltdpl_sph :
			noSsites[ (offset[0]+newPsite[0], offset[1]+newPsite[1], offset[2]+newPsite[2]) ] = 1
	cubesurf_dpl = cubeSurf(dpl)
	cubesurf_dll = cubeSurf(dll)
	rad_ltdll_sph = [] ## which points around a grid-pt are less than 4A from it ?
	for i in range(-1*dll,dll+1) :
		for j in range(-1*dll,dll+1) :
			for k in range(-1*dll,dll+1) :
				if i*i + j+j + k*k < dll*dll : rad_ltdll_sph.append((i,j,k))
	# for each Psite, find solvent sites s.t. each Ssite has no Psite < 3 from it and no Ssite < 4 from it
	Ssites = {}
	for site in Psites.keys() :
		for offset in cubesurf_dpl :
			Scandidate = (site[0]+offset[0], site[1]+offset[1], site[2]+offset[2])
			occ = None
			try : # does the candidate enter no-solvent zone defined by protein sites and existing Ssites ?
				if noSsites[Scandidate] == 1 : occ = 1
			except : pass
			if occ : continue
			Ssites[Scandidate] = 1
			for off in rad_ltdll_sph : # mark grid-points too close to new solvent site
				noSsites[ (Scandidate[0]+off[0], Scandidate[1]+off[1], Scandidate[2]+off[2]) ] = 1
	# add another layer of solvent around existing solvent
	sites = list(Ssites.keys())
	for site in sites :
		for offset in cubesurf_dll :
			Scandidate = (site[0]+offset[0], site[1]+offset[1], site[2]+offset[2])
			occ = None
			try : # does the candidate enter no-solvent zone defined by protein sites and existing Ssites ?
				if noSsites[Scandidate] == 1 : occ = 1
			except : pass
			if occ : continue
			Ssites[Scandidate] = 1
			for off in rad_ltdll_sph : # mark grid-points too close to new solvent site
				noSsites[ (Scandidate[0]+off[0], Scandidate[1]+off[1], Scandidate[2]+off[2]) ] = 1
	return list(Ssites.keys())

## make cubic grid around input coordinates and assign each input crd to a grid point
## make sure there are at least 3 layers of empty grid points on all sides of the grid
## return the empty grid points

def makeGrid(crd_in, dgrid = 3.) :
	maxcrd = [-100000000., -100000000., -100000000.]
	mincrd = [100000000., 100000000., 100000000.]
	for crd in crd_in :
		if crd[0] > maxcrd[0] : maxcrd[0] = crd[0]
		if crd[1] > maxcrd[1] : maxcrd[1] = crd[1]
		if crd[2] > maxcrd[2] : maxcrd[2] = crd[2]
		if crd[0] < mincrd[0] : mincrd[0] = crd[0]
		if crd[1] < mincrd[1] : mincrd[1] = crd[1]
		if crd[2] < mincrd[2] : mincrd[2] = crd[2]
	mincrd = (mincrd[0] - 3*dgrid, mincrd[1] - 3*dgrid, mincrd[2] - 3*dgrid)
	maxcrd = (maxcrd[0] + 3*dgrid, maxcrd[1] + 3*dgrid, maxcrd[2] + 3*dgrid)
	## for each point, find a grid index closest to it, from the eight arnd it
	occupiedIndices = []
	for crd in crd_in :
		occupiedIndices.append(
			(
			int( round((crd[0] - mincrd[0]) / dgrid) ),
			int( round((crd[1] - mincrd[1]) / dgrid) ),
			int( round((crd[2] - mincrd[2]) / dgrid) )
			)
		)
	#print occupiedIndices
	retcrd = []
	for i in range( int( (maxcrd[0]-mincrd[0]) / dgrid ) ) :
		for j in range( int( (maxcrd[1]-mincrd[1]) / dgrid ) ) :
			for k in range( int( (maxcrd[2]-mincrd[2]) / dgrid ) ) :
				if (i,j,k) in occupiedIndices : pass #print "OCCUPIED", i,j,k
				else :
					retcrd.append( (mincrd[0] + i*dgrid, mincrd[1] + j*dgrid, mincrd[2] + k*dgrid) )
	return retcrd
