#
# nudraw
#
# A program for visualizing DNA Structures and Devices
# Copyright 2004
# Niles A. Pierce
# California Institute of Technology
# Version 1.0 coded 2004 (April 9,10,13,14)
#

import sys
from numpy import arange, array, cos, linspace, pi, sin, random, zeros, sqrt, ones, sqrt
try:
  import pylab
  from matplotlib import axes3d
except:
  pylab=None
from scipy import interpolate,optimize


from scipy.interpolate import splprep, splev, spalde, interp1d
import scipy
import pdb
try:
  import Coil
except:
  Coil=None
import pdb

#-------------------------------------------------------------------------------
def get_angle(radius, center, point):
  """ Returns the angle between the line point -> center and the line y==0
  """

  # cos TH=(point.x-center.x)
  #           point-center
  th = scipy.arccos((point[0] - center[0])/radius)
  #~ print "th1=%f"%th
  if th<0:
    th=scipy.pi+th
    #~ print "th2=%f"%th
  if center[1]>point[1]:
    th=-th
    #~ print "th3=%f"%th

  th=th*180/scipy.pi
  #~ print "get_point:(%f, %f) (%f, %f) -> %f"%(center[0], center[1],
    #~ point[0] - center[0], point[1] - center[1], th)
  return th

#-------------------------------------------------------------------------------
def get_angle_rad(radius, center, point):
  """ Returns the angle between the line point -> center and the line y==0
  """

  # cos TH=(point.x-center.x)
  #           point-center
  th = scipy.arccos((point[0] - center[0])/radius)
  #~ print "th1=%f"%th
  if th<0:
    th=scipy.pi+th
    #~ print "th2=%f"%th
  if center[1]>point[1]:
    th=-th
    #~ print "th3=%f"%th

  #~ print "get_point:(%f, %f) (%f, %f) -> %f"%(center[0], center[1],
    #~ point[0] - center[0], point[1] - center[1], th)
  return th


#-------------------------------------------------------------------------------
def spl3(x,y,x2,s=3):
  polycoeffs = array(scipy.polyfit(x, y, s)).flatten()

  y2 = scipy.polyval(polycoeffs, x2)
  dercoeffs=scipy.polyder(polycoeffs, 1)

  dy2=scipy.polyval(dercoeffs, x2)
  return y2,dy2

def spl2(x,y,x2,s=2.0):
  us = interpolate.UnivariateSpline(x, y, s=s)
  return us(x2)

def polyfit(yp,order=3,n=100):
  t=linspace(1,100,num=len(yp[0]))
  newt=linspace(1,100,num=n)
  newx,newdx=spl3(t,yp[0],newt,order)
  newy,newdy=spl3(t,yp[1],newt,order)
  newz,newdz=spl3(t,yp[2],newt,order)
  return newx,newy,newz,newdx,newdy,newdz

def polyfit2d(yp,order=3,n=100):
  t=linspace(1,100,num=len(yp))
  newt=linspace(1,100,num=n)
  newx,newdx=spl3(t,yp,newt)
  return newx,newdx


#-------------------------------------------------------------------------------
def dist(v1,v2):
  return sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + \
              (v2[2]-v1[2])**2)

#-------------------------------------------------------------------------------
def dist2d(v1,v2):
  return sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2)

#-------------------------------------------------------------------------------
def scalardist(v1,v2):
  return sqrt((v2-v1)**2)
  return d

#-------------------------------------------------------------------------------
# Search for index using bisection
def findpoint(x,p,si,ei,imin=False, imax=False, maxiter=50, debug=False, nofalse=False):
  if imin==False and imin!=0:
    imin=si
  if imin<0: imin=0
  if imax>=len(p): imax=len(p)-1
  if si<imin: si=imin
  if si>imax: si=imax

  if imax==False:
    imax=ei
  if ei>imax: ei=imax
  if ei<imin: ei=imin


  delta=(ei-si)/2.
  idelta=int(delta)
  ind=[si, si+idelta, ei]
  old_d=1e6
  d=0
  mycenters=[]
  eps=1e-3
  solvec=[False,False,False]
  count=0

  if x[2]==0.0:
    fdist=dist2d
  else:
    fdist=dist

  if debug:
    print "x=%s"%(x)
    print "\nsi3=%d ei=%d"%(si,ei)
    print "imin=%d imax=%d"%(imin,imax)

  while scalardist(old_d,d)>eps:
    mycenters.append(ind[1])
    if debug: print "\n",count, ind, solvec
    for i in range(3):
      if ind[i]>=len(p): ind[i]=len(p)-1
      if ind[i]<0: ind[i]=0
    #~ print "ind=%s, len(p)=%d, imin=%d, imax=%d"%(ind, len(p), imin, imax)
    for i in range(3):
      if not solvec[i]:
        solvec[i]=fdist(x,p[ind[i]])
      if debug: print "p[%d]=%s"%(ind[i],p[ind[i]])
    min_ind=0
    for i in [1,2]:
      if solvec[i]<solvec[min_ind]:
        min_ind=i

    if solvec[min_ind]<eps:
      return ind[min_ind],mycenters# We got lucky
    d_old=d
    d=solvec[min_ind]

    if debug: print solvec,"min_ind=",min_ind,d

    if min_ind==0: # Step left
      if debug: print "step left, idelta=%d"%idelta
      solvec[2]=solvec[1]
      ind[2]=ind[1]
      solvec[0]=False
      ind[0]=ind[0]-idelta
      if ind[0]<imin:
        ind[0]=imin
      solvec[0]=False
      idelta=(ind[2]-ind[0])/2
      delta=float(idelta)

      solvec[1]=False
      ind[1]=ind[0]+idelta
      if debug: print ind

    elif min_ind==1: # We bracketed the solution
      accel=1.;max(1,min(scalardist(solvec[0],solvec[1]),scalardist(solvec[1],solvec[2])))
      delta=min(0.5*delta,delta/(2.*accel))
      idelta=int(delta)
      if debug: print "bracket, idelta=%d"%idelta
      solvec[0]=False
      ind[0]=ind[1]-idelta
      solvec[2]=False
      ind[2]=ind[1]+idelta

    elif min_ind==2: #Step Right
      solvec[0]=solvec[1]
      ind[0]=ind[1]
      if debug: print "step right, idelta=%d"%idelta
      ind[2]=ind[2]+idelta
      if ind[2]>imax:
        ind[2]=imax-1
      solvec[2]=False

      idelta=(ind[2]-ind[0])/2
      delta=float(idelta)

      solvec[1]=False
      ind[1]=ind[0]+idelta
      if debug: print ind

    if idelta==0 or count>maxiter:
      if count>maxiter:
        print "Error: maxiter reached!"
      if min_ind==0 and ind[0]==False and nofalse:
        return 0,mycenters
      return ind[min_ind],mycenters
    count=count+1

# Find base points on interpolated set of points
def findPoints(points, bpoints):
  pdist=numpy.ones([len(points)+1])*1e6

  # First, calculate the distances between interpolated points
  map(lambda j: pdist[j].__setitem__(dist(points[j], points[j+1])), range(1,len(points)-1))

# spl("",sp,[xp, yp, zp, cp], ss ,xbc1, xbc2, 0.,0., order=4, stiff=20.0)
#-------------------------------------------------------------------------------
def spl(SPLINE,v,yp,xs,ybc1,ybc2,ksbc1,ksbc2,order=2, stiff=3.0):
  ys=[]
  dys=[]

  constant=False
  max_order = len(yp[0]) - 1

  if order > max_order:
    order = max_order
    #print("Setting order to "+str(max_order))

  #~ print "spline.spl:"
  #~ print "xp=\n%s"%yp[0]
  #~ print "yp=\n%s"%yp
  #~ print "zp=\n%s"%yp[2]

  #~ for y in yp:
    #~ t=(y==y[0])
    #~ if t.all():
      #~ constant=True
      #~ break
  if order==0:
    for y in yp:
      ys.append(y[0]*ones(xs.size))
      dys.append(zeros(xs.size))
  else:
    s=stiff
    nest=-1 # estimate of number of knots needed (-1 = maximal)
    k=order # spline order
    # find the knot points
    tckp,u = splprep(yp,u=v/v[-1],s=s,k=k,nest=len(yp[0]))

  # evaluate spline, including interpolated points

    ypoints=splev(xs,tckp)
    dypoints=splev(xs,tckp,1)

  return ypoints+dypoints


#-------------------------------------------------------------------------------
# Fit a circle through points (xp, yp), interpolate from z0 to zN
def circlefit(yp, center, radius, n, ax3d, colors):
  #if not len(yp)==3: raise ValueError("Points array needs to have 3xN dimensions")

  dz,z=polyfit2d(yp[2], order=min(len(yp)-1,2), n=n)
  points2d=array([yp[0],yp[1]]).T
  #~ print "points2d=\n%s"%points2d

  ang1=get_angle_rad(radius, center[:2], points2d[0])
  ang2=get_angle_rad(radius, center[:2], points2d[-1])

  #~ print "ang1=%s"%(ang1*180/pi)
  #~ print "ang2=%s"%(ang2*180/pi)

  th=linspace(ang1,ang2,num=n)
  #~ print "dth=%f\nth=\n%s"%((ang2-ang1)*180/pi,th)
  x=radius*cos(th)
  y=radius*sin(th)
  dx=-y.copy()
  dy=x.copy()
  x=x+center[0]
  y=y+center[1]

  #~ print "x=%s"%(x)
  #~ print "y=%s"%(y)

  return x, y, z, dx, dy, dz



#-------------------------------------------------------------------------------
def prep(x,y,z, n=100):
  sp=zeros(x.size)
  for i in range(1,x.size):
    sp[i] = sp[i-1] + sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + \
            (z[i]-z[i-1])**2)
  ss=linspace(0.,1.,n)
  x1,x2,x3, dx1,dx2,dx3 = spl(None,sp,[x,y,z],ss,0,0,0,0,order=1,stiff=10.)
  return x1,x2,x3, dx1,dx2,dx3

#-------------------------------------------------------------------------------
def test_sparse1():
  x=array([-11.53239257,  -7.71353944])
  y=array([-5.31375002,   -0.16157562])
  z=array([ 0.,           3.57089755])
  return  x,y,z,prep(x,y,z)

#-------------------------------------------------------------------------------
def test_sparse2():
  x=array([ 0.50598239,    13.21241212,    16.27147014,    17.29437447,  15.44431414])
  y=array([ 6.63842438,    10.70186058,    16.0003045,     22.0323034,   28.42015195])
  z=array([ 8.48492674e+00, 1.61800038e-15, 1.99261457e-15, 2.11788009e-15, -3.57089755e+00])
  return  x,y,z,prep(x,y,z)

#-------------------------------------------------------------------------------
def test_sparse3():
  x=array([ 3.76786639, -3.05905802, -8.79437447, -13.40138343])
  y=array([38.50775609, 38.88487358, 36.75473526,  31.95861174])
  z=array([-3.57089755, -3.57089755, -3.57089755,  -3.57089755])
  n=100
  newx,newy,newz,newdx,newdy,newdz=polyfit([x,y,z],order=3,n=n)
  points=array([newx,newy,newz]).T
  xmin=points[7]
  imin=findpoint(xmin,points,0,n-1)
  result=points[imin]

  return x,y,z,newx,newy,newz

#-------------------------------------------------------------------------------
def test_trig_gl(n=200):
  # make ascending spiral in 3-space
  t=linspace(0,1.75*2*pi,100)

  x = sin(t)
  y = cos(t)
  z = t

  # add noise
  #x+= random.normal(scale=0.1, size=x.shape)
  #y+= random.normal(scale=0.1, size=y.shape)
  #z+= random.normal(scale=0.1, size=z.shape)

  # spline parameters
  s=3 # smoothness parameter
  k=2 # spline order
  nest=-1 # estimate of number of knots needed (-1 = maximal)
  sp=zeros(x.size)
  sp[0]=0
  for i in range(1,x.size):
    sp[i] = sp[i-1] + sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + \
            (z[i]-z[i-1])**2)

  # find the knot points
  ss=linspace(0,1,n)
  return  spl(None,sp,[x,y,z],ss,0,0,0,0)

#-------------------------------------------------------------------------------
def test_trig(ax3d):
  # make ascending spiral in 3-space
  t=linspace(0,1.75*2*pi,100)

  x = sin(t)
  y = cos(t)
  z = t

  # add noise
  #x+= random.normal(scale=0.1, size=x.shape)
  #y+= random.normal(scale=0.1, size=y.shape)
  #z+= random.normal(scale=0.1, size=z.shape)

  # spline parameters
  s=3 # smoothness parameter
  k=2 # spline order
  nest=-1 # estimate of number of knots needed (-1 = maximal)
  sp=zeros(x.size)
  sp[0]=0
  for i in range(1,x.size):
    sp[i] = sp[i-1] + sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + \
            (z[i]-z[i-1])**2)

  # find the knot points
  tckp,u = splprep([x,y,z],u=sp/sp[-1],s=s,k=k,nest=-1)
  ss=linspace(0,1,400)
  # evaluate spline, including interpolated points
  xnew,ynew,znew = splev(ss,tckp)

  ax3d.plot3D(xnew, ynew, znew, 'b-')
    # evaluate spline, including interpolated points
  #ss=linspace(0,1,400)

  x1,x2,x3, dx1,dx2,dx3 = spl(None,sp,[x,y,z],ss,0,0,0,0)
  #print x1
  ax3d.plot3D(x1, x2, x3, 'r+')
  fig2d=pylab.figure()
  pylab.plot(x1,ss,'g-')
  pylab.plot(x2,ss,'r-')
  pylab.plot(x3,ss,'b-')
  pylab.plot(xnew,ss,'g+')
  pylab.plot(ynew,ss,'r+')
  pylab.plot(znew,ss,'b+')

#-------------------------------------------------------------------------------
if __name__=="__main__":
  print sys.argv
  if Coil and len(sys.argv)>1 and sys.argv[1]=='gl':
    mycoil=Coil.PolyCylinderMulti()
    x1,x2,x3, dx1,dx2,dx3 = test_trig_gl()
    mycoil.addPoly(array([x1,x2,x3]).T)
    mycoil.Display()
  else:
    fig3d=pylab.figure()
    ax3d = axes3d.Axes3D(fig3d)
    x,y,z, x1,x2,x3 = test_sparse3()
    ax3d.plot3D(x1, x2, x3, 'r+')
    #print x1,x2,x3
    ax3d.plot3D(x, y, z, 'g-')
    pylab.show()
