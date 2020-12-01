import math
import numpy as np
class Point:
	PointList = []
	PointCount = 0
	def __init__(self, IDnum, affiliation, xpos, ypos, zpos, gridface, isborder, terrain):#xypos, projxypos, NNlist ):
		self.ID = int(IDnum)					#Integer
		self.Province = affiliation				#String
		self.XYZPos = (float(xpos),float(ypos),float(zpos))	#Vector or list
		self.XYZStr = xpos+' '+ypos+' '+zpos			#String
		self.GridFace = int(gridface)				#Integer
		self.isBorder = int(isborder)				#Integer
		self.Terrain = terrain 					#String
		Point.PointList.append(self)				#Add new point to list of points
		Point.PointCount += 1					#Append the size of point list
		xpos = float(xpos)
		ypos = float(ypos)
		zpos = float(zpos)
		phi = math.asin(zpos/math.sqrt(xpos*xpos+ypos*ypos+zpos*zpos))
		lam = math.pi/2.0
		the = math.acos(-1*zpos)
		if(xpos	< 0 and ypos < 0):
			lam = math.atan(ypos/xpos)
		if(xpos == 0 and ypos < 0):
			lam = 1.0*math.pi/2.0
		if(xpos > 0 and ypos < 0):
			lam = math.atan(ypos/xpos)+math.pi
		if(xpos < 0 and ypos >= 0):
			lam = math.atan(ypos/xpos)
		if(xpos == 0 and ypos >= 0):
			lam = 3.0*math.pi/2.0
		if(xpos > 0 and ypos >= 0):
			lam = math.atan(ypos/xpos)+math.pi
		if(lam < 0):
			lam = lam+2*math.pi
		if(lam >= 2*math.pi):
			lam = lam-2*math.pi
		self.lat = phi
		self.lon = lam
		self.the = the
	def getID(self):
		return self.ID
	def getProvince(self):
		return self.Province
	def getXYZPos(self):
		return self.XYZPos
	def getXYZStr(self):
		return self.XYZStr
	def getGridFace(self):
		return self.GridFace
	def getBorder(self):
		return self.isBorder
	def isBorder(self):
		if(self.isBorder == 1):
			return True
		else:
			return False
	def getTerrain(self):
		return self.Terrain
	def getLat(self):
		return self.phi
	def getLatDeg(self):
		return 180.0*self.phi/math.pi
	def getLon(self):
		return self.lam
	def getLonDeg(self):
		return 180.0*self.lam/math.pi
	def getXYPosEQR(self,yframe):
		Y = int(self.the*2160/math.pi-0.5)
		X = int(self.lon*2160/math.pi-0.5)
		if(X < 0):
			X = 0
		if(X >= 2*yframe):
			X = 2*yframe-1
		if(Y < 0):
			Y = 0
		if(Y >= yframe):
			Y = yframe-1
		return (X,Y)
	def getXYStrEQR(self,yframe):
		Y = int(self.the*2160/math.pi-0.5)
		X = int(self.lon*2160/math.pi-0.5)
		if(X < 0):
			X = 0
		if(X >= 2*yframe):
			X = 2*yframe-1
		if(Y < 0):
			Y = 0
		if(Y >= yframe):
			Y = yframe-1
		return str(X)+' '+str(Y)

def getPoint(idnum):
	return Point.PointList[idnum]

def getPointList():
	return Point.PointList

#	Finds & returns the center of XYZ coordinates from Points list
def getCentralPoint(Points):
	xtot = 0.0
	ytot = 0.0
	ztot = 0.0
	for pt in Points:
		XYZ = pt.getXYZPos()
		xtot += XYZ[0]
		ytot += XYZ[1]
		ztot += XYZ[2]
	xNew = xtot/len(Points)
	yNew = ytot/len(Points)
	zNew = ztot/len(Points)
	size = math.sqrt(xNew*xNew+yNew*yNew+zNew*zNew)
	xNew = xNew/size
	yNew = yNew/size
	zNew = zNew/size
	return (xNew,yNew,zNew)

#	Loads list of points from text files at resolution number 'res'
def LoadPoints(res):
        PID = list(open('PID'+res+'.txt'))
        CRD = list(open('CRD'+res+'.txt'))
        FCE = list(open('FCE'+res+'.txt'))
        BRD = list(open('BRD'+res+'.txt'))
        TYP = list(open('TYP'+res+'.txt'))
	for ct in range(len(PID)):
		xyz = CRD[ct].split(' ')
                point = Point(str(ct),PID[ct][0:-1],xyz[0],xyz[1],xyz[2][0:-1],FCE[ct][0:-1],BRD[ct][0:-1],TYP[ct][0:-1])

