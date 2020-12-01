#**************************************************************************************
#	Sample Model Atlas for Nesting on a Truncated Icosahedral Star (MANTIS)
#**************************************************************************************
#
#	By Cameron D. Rae
#
#	Background:
#
#		Truncated Icosahedra are typically encountered as the 'Soccer ball' or 'Football' shape
#		possessing 20 regular hexagons sewn together by 12 regular pentagons. By calculating the
#		midpoints between each of the vertices and face centres, successively higher and higher
#		resolution grids may be constructed by projecting each of these midpoints back onto the
#		surface of the inscribing sphere (hence it can be interpreted as a 'star' shape). This
#		generating requirement limits, or quantizes, the available resolutions for this kind of
#		grid. 
#
#	Technical Implementation:
#
#		Let us define resolution 1 as the classic football shape, having 32 faces and 60 vertices
#		totaling 92 points spaced across the sphere, with an (approximate) equivalent length scale 
#		of sqrt(16/Npts) or about 0.417 radians (about 2650 km). Applying the midpoint algorithm 
#		described above, resolution 2 corresponds to 362 points on the sphere and a distance scale 
#		of 1340km. Continuing this algorithm a few times, this sample applies the grid with 1,474,562 
#		points and a distance scale of about 21km, or resolution 8. This can be interpreted as the 
#		distance traveled after walking for 4 hours, assuming an average walking pace of 5.25 km/hr,
#		and should represent a relatively familiar scale/resolution as it is 'within walking distance'. 
#
#	Motivation & Scope:
#
#		Projecting geospatial data onto this grid allows for modeling, analysis, and visualizations
#		with much less severe distortion as would be encountered in standard latitude/longitude grids
#		of a sphere. Each point is roughly the same size, and this is demonstrated in the png images
#		which represent each grid cell as a 15-pixel diameter circular texture, corresponding to a 
#		pre-calculated landcover type assigned to that point. These points can be clustered together
#		in what I've described as 'Provinces', which are simply customizable domains of the grid
#		which can be analyzed separately (the sample provinces which I have provided do NOT represent
#		political boundaries, are completely fictional and for demonstrative purposed only!!).
#
#	Data Storage:
#
#		To make the sample easier to interact with, I have provided the grid data as 5 text files, 
#		each of which is a simple list of values where the line number corresponds to the point ID 
#		number (0-59 are vertices of football, 60-91 are face centres, etc.). These files are:
#				
#			BRD8.txt is used to denote borders between defined Provinces, as in any point 
#				within a province which has at least one neighboring point in a different
#				province. Value is 1 if point is a border, or 0 if point is not a border.
#
#			CRDS8.txt is used to define the 3-D location of each grid cell on the unit sphere.
#				By convention, the North Pole occurs at (0,0,-1) and the South Pole occurs
#				at (0,0,1).
# 
#			FCE8.txt is used to identify the closest face center, and is therefore a number between
#				0 and 31 corresponding to the points on lines 60-91 in each file. Useful for
#				local point lookups, avoiding having to loop over all points in the global list.
#
#			PID8.txt is used to identify the Provincial affiliation of each point. Again, these do
#				not represent any political or physical boundaries, they're just to demonstrate
#				the capabilities of the grid.
#			
#			TYP8.txt is used to identify the IGBP land cover type corresponding to each point. The 
#				textures in the sample images provided do not necessarily reflect the IGBP 
#				types, but they will be labeled in PyPlot.
#
#	Providing the Sample & Script:
#
#		To demonstrate this within GitHub's data limits, I have provided tangent-plane projections
#		of each non-polar pentagonal face centre. The TangentPlane.png images can be viewed in most 
#		image-viewing software, and shouldn't require this script. The script is useful for showcasing
#		some of the algorithms needed for nesting models within the grid, like the calculus by 
#		regression which could be useful for solving differential equations, field operations, and 
#		nearest-neighbor identifications which link the grid together. Some basic algorithms are also
#		provided to extract data from more standard latitude/longitude grids. Proper re-gridding should
#		be done by equating open balls in each metric space & topology! 
#
#		For illustrative purposes, these sample images overlay output from a very crude sea ice & snow 
#		accumulation model which runs on MERRA2 daily average surface temperature and total precipitable
#		moisture (TQV). Data was randomly selected from the set of days in November occurring in years 
#		between 1980 and 2019. 
#	
#	References:
#
#		IGBP Data: Loveland, T., J. Brown, D. Ohlen, B. Reed, Z. Zhu, L. Yang, and S. Howard . 2009.ISLSCP 
#		II IGBP DISCover and SiB Land Cover, 1992-1993. In Hall, Forrest G., G. Collatz, B. Meeson, S. Los, 
#		E. Brown de Colstoun, and D. Landis (eds.). ISLSCP Initiative II Collection. Data set. Available on-
#		line [http://daac.ornl.gov/] from Oak Ridge National Laboratory Distributed Active Archive Center, 
#		Oak Ridge, Tennessee, U.S.A. doi:10.3334/ORNLDAAC/930
#
#		MERRA2 Data: The Modern-Era Retrospective Analysis for Research and Applications, Version 2 (MERRA-2),
#		Ronald Gelaro, et al., 2017, J. Clim., doi: 10.1175/JCLI-D-16-0758.1


#******************************************
#	Preliminary Imports & Setup
#******************************************
#	This sample version runs using python 2.7, the graphical algorithms will NOT work (as presented here) in python 3
#	Also, please run unzip.py or otherwise unzip data files, if you have not already done so!

import math
import numpy as np
import Point
import matplotlib.pyplot as plt
from PIL import Image

Point.LoadPoints('8')
Grid = Point.getPointList()
distancescale = math.sqrt(16.0/len(Grid))
imgdir = './imgcirc/IGBP/'
NNLIST = list(open('NNFULL.txt'))


#************************************
#	Basic Grid Algorithms
#************************************ 

#	Returns the dot product  between point1 and point2
def dot(XYZ1, XYZ2):
	return (XYZ1[0]*XYZ2[0]+XYZ1[1]*XYZ2[1]+XYZ1[2]*XYZ2[2])
#	Returns the cross product between point1 and point2 as a vector (list)
def cross(r1,r2):
	tot = 0.0
	rx = r1[1]*r2[2]-r1[2]*r2[1]
	ry = r1[2]*r2[0]-r1[0]*r2[2]
	rz = r1[0]*r2[1]-r1[1]*r2[0]
	R = (rx,ry,rz)
	return R

#	Returns the cross product between point1 and point2 as a string
def crossStr(r1,r2):
	tot = 0.0
	rx = r1[1]*r2[2]-r1[2]*r2[1]
	ry = r1[2]*r2[0]-r1[0]*r2[2]
	rz = r1[0]*r2[1]-r1[1]*r2[0]
	return str(rx)+' '+str(ry)+' '+str(rz)

#	Returns a vector from the corresponding string of XYZ coordinates
def toVector(strpoint):
	tmp = strpoint.split(' ')
	x0 = float(tmp[0])
	y0 = float(tmp[1])
	z0 = float(tmp[2])
	return (x0,y0,z0)	

#	Returns Euclidean distance between point1 and point2
def dEuclidStr(XYZ1s, XYZ2s):
	XYZ1 = toVector(XYZ1s)
	XYZ2 = toVector(XYZ2s)
	return math.sqrt((XYZ1[0]-XYZ2[0])*(XYZ1[0]-XYZ2[0])+(XYZ1[1]-XYZ2[1])*(XYZ1[1]-XYZ2[1])+(XYZ1[2]-XYZ2[2])*(XYZ1[2]-XYZ2[2]))

#	Returns Euclidean distance between point1 and point2
def dEuclid(XYZ1, XYZ2):
	return math.sqrt((XYZ1[0]-XYZ2[0])*(XYZ1[0]-XYZ2[0])+(XYZ1[1]-XYZ2[1])*(XYZ1[1]-XYZ2[1])+(XYZ1[2]-XYZ2[2])*(XYZ1[2]-XYZ2[2]))

#	Returns the distance on unit sphere between point1 and point 2
def dSurfStr(point1, point2):
	D0 = dEuclidStr(point1,point2)
	if(D0 < 2):
		return math.acos(1-D0*D0/2.0)
	else:
		return D0*math.pi/2.0

#	Returns the distance on unit sphere between point1 and point2
def dSurf(point1, point2):
	D0 = dEuclid(point1,point2)
	if(D0 < 2):
		return math.acos(1-D0*D0/2.0)
	else:
		return D0*math.pi/2.0

#	Returns the length of XYZ string vector
def modStr(point):
	XYZ = toVector(point)
	return math.sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2])

#	Returns the length of XYZ vector
def mod(XYZ):
	return math.sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2])

#	Returns the center of a collection of points
def getCentralPoint(points):
	xtot = 0.0
	ytot = 0.0
	ztot = 0.0
	for pts in points:
		ptXYZ = pts.getXYZPos()
		xtot += ptXYZ[0]
		ytot += ptXYZ[1]
		ztot += ptXYZ[2]
	xNew = xtot/len(points)
	yNew = ytot/len(points)
	zNew = ztot/len(points)
	size = math.sqrt(xNew*xNew+yNew*yNew+zNew*zNew)
	xNew = xNew/size
	yNew = yNew/size
	zNew = zNew/size
	return (xNew,yNew,zNew)

#Returns a list of all point IDs within radius of point
def NNID(point, radius):
	nearby = []
	for p in Grid:
		tmp = p.getXYZPos()
		x1 = tmp[0]
		y1 = tmp[1]
		z1 = tmp[2]
		if(dSurf(point,(x1,y1,z1)) <= 1.1*radius*distancescale):
			nearby.append(p.getID())
	return nearby

# Returns list of nearest neighbor points
def FindNN(ptID):
	NNpts = NNLIST[ptID]
	tmp = NNpts.split(' ')[0:-1]
	if(len(tmp) == 6):
		return (int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]))
	elif(len(tmp) == 7):
		return (int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]),int(tmp[6]))
	else:
		print('Wrong Number of Nearest Neighbors at point '+str(ptID))

#	Calculated and returns the meridional derivative (d/dy) of a field (list of values on grid). Calculates derivative using a regression of nearest neighbor points
def GridDY(field):
	DDYarr = []
	for pts in Grid:
		NXYZs = LocalCRDS(pts.getXYZPos(),FindNN(pts.getID()))
		dydy = 0.0
		dxdy = 0.0
		for nxyz in NXYZs:
			dxdy += (field[nyzy[0]]-field[pts.getID()])*(nxyz[2]*distancescale-0.0)
			dydy += (nxyz[2]*distancescale-0.0)*(nyxz[2]*distancescale-0.0)
		DDYarr.append(dxdy/dydy)
	return list(DDYarr)

#	Calculates and returns the zonal derivative (d/dx) of a field (list of values on grid). Calculates derivative using a regression of nearest neighbor points
def GridDX(field):
	DDYarr = []
	for pts in Grid:
		NXYZs = LocalCRDS(pts.getXYZPos(),FindNN(pts.getID()))
		dydy = 0.0
		dxdy = 0.0
		for nxyz in NXYZs:
			dxdy += (field[nyzy[0]]-field[pts.getID()])*(nxyz[1]*distancescale-0.0)
			dydy += (nxyz[1]*distancescale-0.0)*(nyxz[1]*distancescale-0.0)
		DDYarr.append(dxdy/dydy)
	return list(DDYarr)

#	Calculates and returns the average value of a field on global Grid
def GridAVG(field):
	tot = 0.0
	for pts in Grid:
		tot += field[pts.getID()]
	return tot/float(len(Grid))

#	Calculates and returns the average value of a field on a subset of Grid (domain)
def GridAreaAVG(field, domain):
	tot = 0.0
	for pts in domain:
		tot += field[pts.getID()]
	return tot/float(len(domain))

#	Converts the Z coordinate into a discrete latitude bin used for zonal averaging, residuals, or other conversion to a lat/lon grid
def ZtoLatGridCoordinate(zval,latgrid):
	if(zval == 1.0):
		return 0
	elif(zval == -1.0):
		return (len(latgrid)-1)
	else:
		return len(latgrid)/2-int(float(len(latgrid))*math.asin(zval)/math.pi)


#	Calculates and returns the zonal average of a field onto a specified latitude grid
def GridZonalAverage(field, latgrid):
	newlat = []
	for i in latgrid:
		newlat.append(0.0)
	for pts in field:
		xyz = pts.getXYZPos()
		newlat[ZtoLatGridCoordinate(xyz[2],len(latgrid))] += field[pts.getID()]
	return list(newlat)


#***********************************
#	Graphical Algorithms
#***********************************

#	Provides information about grid cell characteristics nearst to the mouse pointer in QuickVis pyplot method
class PointerData(object):
	def __init__(self,im,strdata):
		self.im = im
		self.data = strdata
	def __call__(self,x,y):
		outdat = self.data[int(x),int(y)]
		outint = outdat[0]*10000+outdat[1]*100+outdat[2]
		tmp0 = ''
		if(outint < 1474562):
			tmp0 = Grid[outint]
			return str(outint)+' '+tmp0.getProvince()+' '+tmp0.getXYZStr()+' '+tmp0.getXYStrEQR(2160)+' '+tmp0.getTerrain()
		else:
			return 'NullSpace'

#	Loads & Displays Pre-Existing Images/ Data to display in pyplot
#	Note: PyPlot will complain because image size exceeds the limit of 89478485 pixels, which is true for the provided sample data. This can be ignored, as pyplot will still open as long as image size is below twice that limit (which it is!)
def QuickVisPyPlot(plane):
	imgdat0 = Image.open('IMGDATA'+str(plane)+'.png')
	imgdat = imgdat0.load()
	img0 = Image.open('TangentPlane'+str(plane)+'.png')
	fig, ax = plt.subplots()
	im = ax.imshow(img0)
	ax.format_coord = PointerData(im,imgdat)
	plt.show()



