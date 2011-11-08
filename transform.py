#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################

from geom import *
from random import random
import sys

# v3 shd become origin, v2 shd lie in -X axis, v1 shd lie in XY region with y>0
def calcTransform(u1,u2,u3, verbose=None) :
	if verbose : print "BEGIN-----------------------------"
	v1 = [u1[0],u1[1],u1[2]]
#list(u1)
	v2 = list(u2)
	v3 = list(u3)
	if verbose : print v1,v2,v3
	dx,dy,dz = -1*v3[0], -1*v3[1], -1*v3[2]
	v3[0],v3[1],v3[2] = v3[0]+dx, v3[1]+dy, v3[2]+dz
	v2[0],v2[1],v2[2] = v2[0]+dx, v2[1]+dy, v2[2]+dz
	v1[0],v1[1],v1[2] = v1[0]+dx, v1[1]+dy, v1[2]+dz
	if verbose : print v1,v2,v3
	alpha = vec_angle_deg( [-1.,0,0], [0,0,0], v2 )
	axis = normalize(cross_product(v2, [-1.,0,0]))
	if verbose : print alpha, axis
	v2 = rotate(v2, axis, alpha)
	v1 = rotate(v1, axis, alpha)
	if verbose : print v1,v2,v3
	theta = vec_dihed_deg(v1, [-1.,0.,0.], [0.,0.,0.], [0.,1.,0.])
	v1 = rotate(v1, [1.,0.,0.], theta)
	if verbose : print v1,v2,v3
	if verbose : checkTransformed(u1,u2,u3,v1,v2,v3)
	return dx,dy,dz,axis[0],axis[1],axis[2],alpha,theta

def checkTransformed(u1,u2,u3,v1,v2,v3) :
	acc = 0.15
	assert math.fabs(vec_dist(v1,v2) - vec_dist(u1,u2)) < acc
	assert math.fabs(vec_dist(v3,v2) - vec_dist(u3,u2)) < acc
	assert math.fabs(vec_angle_deg(v1,v2,v3) - vec_angle_deg(u1,u2,u3)) < acc
	assert math.fabs(v3[0]) < acc
	assert math.fabs(v3[1]) < acc
	assert math.fabs(v3[2]) < acc
	assert math.fabs(v2[1]) < acc
	assert math.fabs(v2[2]) < acc
	assert math.fabs(v1[2]) < acc

def doTransform(pt, dx,dy,dz, ax0,ax1,ax2, alpha,theta) :
	ret = [0.,0.,0.]
	ret[0] = pt[0] + dx
	ret[1] = pt[1] + dy
	ret[2] = pt[2] + dz
	ret = rotate(ret, [ax0,ax1,ax2], alpha)
	ret = rotate(ret, [1.,0.,0.], theta)
	return ret

def testTransform() :
	u1 = [random(), random(), random()]
	u2 = [random(), random(), random()]
	u3 = [random(), random(), random()]
	#u1 = [0.,0.,0.]
	#u2 = [0.,1.,0.]
	#u3 = [0.,0.,1.5]
	dx,dy,dz,ax1,ax2,ax3,alphaAxis,thetaX = calcTransform(u1,u2,u3)#(list(u1),list(u2),list(u3))
	v1,v2,v3 = None,None,None
	v1 = doTransform(u1, dx,dy,dz,ax1,ax2,ax3,alphaAxis,thetaX)
	v2 = doTransform(u2, dx,dy,dz,ax1,ax2,ax3,alphaAxis,thetaX)
	v3 = doTransform(u3, dx,dy,dz,ax1,ax2,ax3,alphaAxis,thetaX)
	print u1,u2,u3
	print v1,v2,v3
	checkTransformed(u1,u2,u3,v1,v2,v3)

if __name__ == "__main__" :
	for i in range(1000) : testTransform()

