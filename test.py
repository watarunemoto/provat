from geom import *

if __name__ == "__main__" :
	DCtoAlphaTheta = {}
	DCtoAlphaTheta[(1.,0.,0.)] = (180,0)
	DCtoAlphaTheta[(0.866,0.5,0.)] = (150,0)
	DCtoAlphaTheta[(0.5,0.866,0.)] = (120,0)
	DCtoAlphaTheta[(0.,1.,0.)] = (90,0)
	DCtoAlphaTheta[(-0.5,0.866,0.)] = (60,0)
	DCtoAlphaTheta[(-0.866,0.5,0.)] = (30,0)
	DCtoAlphaTheta[(-1.,0.,0.)] = (0,0)
	DCtoAlphaTheta[(-0.866,-0.5,0.)] = (30,180)
	DCtoAlphaTheta[(-0.5,-0.866,0.)] = (60,180)
	DCtoAlphaTheta[(0.,-1.,0.)] = (90,180)
	DCtoAlphaTheta[(0.5,-0.866,0.)] = (120,180)
	DCtoAlphaTheta[(0.866,-0.5,0.)] = (150,180)
	## now 45 degrees
	DCtoAlphaTheta[(0.707,0.,0.707)] = (180,0)
	DCtoAlphaTheta[(0.5,0.5,0.707)] = (180,0)
	DCtoAlphaTheta[(0.,0.707,0.707)] = (180,0)
	DCtoAlphaTheta[(-0.5,0.5,0.707)] = (180,0)
	DCtoAlphaTheta[(-0.707,0.,0.707)] = (180,0)
	DCtoAlphaTheta[(-0.5,-0.5,0.707)] = (180,0)
	DCtoAlphaTheta[(0.,-0.707,0.707)] = (180,0)
	DCtoAlphaTheta[(0.5,-0.5,0.707)] = (180,0)
	## now -45 degrees
	DCtoAlphaTheta[(0.707,0.,-0.707)] = (180,0)
	DCtoAlphaTheta[(0.5,0.5,-0.707)] = (180,0)
	DCtoAlphaTheta[(0.,0.707,-0.707)] = (180,0)
	DCtoAlphaTheta[(-0.5,0.5,-0.707)] = (180,0)
	DCtoAlphaTheta[(-0.707,0.,-0.707)] = (180,0)
	DCtoAlphaTheta[(-0.5,-0.5,-0.707)] = (180,0)
	DCtoAlphaTheta[(0.,-0.707,-0.707)] = (180,0)
	DCtoAlphaTheta[(0.5,-0.5,-0.707)] = (180,0)
	## now poles
	DCtoAlphaTheta[(0.,0.,1.)] = (180,0)
	DCtoAlphaTheta[(0.,0.,-1.)] = (180,0)

	a, b, c = (-2.,2.,0.), (-2.,0.,0.), (0.,0.,0.)
	for dc,at in DCtoAlphaTheta.items() :
		al, dh = vec_angle_deg(b,c,list(dc)), vec_dihed_deg(a,b,c,list(dc))
		print a, b, c, dc, al, dh
