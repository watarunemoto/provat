#############################################################################################
### THIS FILE IS PART OF PROVAT, A TOOL FOR TESELLATION AND VISUALIZATION OF MACROMOLECULES
### (c) Swanand Gore (swanand@cryst.bioc.cam.ac.uk)
### Structural Biology and Biocomputing Group, Dept of Biochemistry,
### University of Cambridge
### Tennis Court Road, Cambridge CB2 1GA, United Kingdom
#############################################################################################


########## some pdb format info
atmi = (6,11)
atmn = (12,16)
resn = (17,20)
chid = (21,22)
resi = (22,26)
ic = (26,27)
assert(atmn[0] < resn[0] and atmn[1] <= resn[0])
assert(resn[0] < chid[0] and resn[1] <= chid[0])
assert(chid[0] < resi[0] and chid[1] <= resi[0])
xcrd = (30,38)
ycrd = (38,46)
zcrd = (46,54)
