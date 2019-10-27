# 3D drawing utilities
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import numpy
import scipy
import Export
from spline import spl, polyfit, findpoint, circlefit
eye=[numpy.identity(0), numpy.identity(1),numpy.identity(2), numpy.identity(3)]
I3=eye[3]

#-------------------------------------------------------------------------------
# Stack base numbers of coil onto base array
def stackBases(i,coil, bindices, base, unit, cmd_options, gg):
  indices=numpy.arange(-1,len(coil[i].bases)+1)
  si=0
  ei=0
  st=unit[base[coil[i].bases[0]].iunit].type
  if st==-5:
    si=-1
  et=unit[base[coil[i].bases[-1]].iunit].type

  if et==-3:
    ei=+1
  indices=numpy.arange(si,len(coil[i].bases)+ei)
  #~ print indices, coil[i].bases

  blist=map(lambda x: coil[i].bases[0]+x, indices)
  barray=numpy.array(blist)

  if bindices==None: return barray
  return numpy.hstack([bindices, barray])

#-------------------------------------------------------------------------------
def findPointsInt(points, bindices, base, gg):
  bpoints=numpy.zeros([len(bindices),3])
  bpoints.__setitem__(0,numpy.array(base[bindices[0]].x3))
  map(lambda j: bpoints.__setitem__(j,numpy.array(base[bindices[j]].x3)), range(len(bindices)))

#-------------------------------------------------------------------------------
def findPoints(points2, ind_per_base, bindices, base, gg, debug=False):
  imin=0 # start index
  imax=ind_per_base*4 # end index
  si=sdelta=ind_per_base*0.8
  ei=edelta=ind_per_base*1.2
  pi=[] # indices of points closest to bases in bindices

  for bi in bindices:
    npi=findpoint(base[bi].x3, points2, si, ei, imin, imax, debug=debug, nofalse=True)[0]
    pi.append(npi)
    imax=min(npi+ind_per_base*4, len(points2)-1)
    imin=npi
    si=npi+sdelta
    ei=npi+edelta

  mi=[0] # middle indices, where color changes will occur
  for i in xrange(len(bindices)-1):
    mi.append(int(round((pi[i]+pi[i+1])/2)))
  mi.append(len(points2))

  texinds=[]

  lmi=len(mi)-1 # last middle index
  for i in xrange(lmi):
    #~ print mi[i],mi[i+1]
    for j in xrange(mi[i+1]-mi[i]):
      texinds.append(base[bindices[i]].domain)
    #~ print "len(texinds) =",len(texinds)
  #~ print "mi =",mi
  #~ print "len(bindices)=%d, lmi=%d, bindices=%s"%(len(bindices), lmi, bindices)
  #~ for j in xrange(mi[lmi]-mi[lmi-1]):
    #~ texinds.append(base[bindices[lmi]].domain)

    #~ cpoints=numpy.vstack([cpoints, numpy.ones([mi[i+1]-mi[i],3])*gg.domainCMap[base[bindices[i]].domain])
  #~ cpoints=numpy.vstack([cpoints, numpy.ones([mi[lmi]-mi[lmi-1],3])*numpy.array(gg.domainCMap[base[bindices[-1]].domain])])

  return gg.domainCMap, texinds

#-------------------------------------------------------------------------------
def sphere(N):
  # u and v are parametric variables.
  # u is an array from 0 to 2*pi, with N elements
  u=numpy.linspace(0, 2*scipy.pi, N)
  # v is an array from 0 to pi, with N elements
  v=numpy.linspace(0, scipy.pi, N)
  # x, y, and z are the coordinates of the points for plotting
  # each is arranged in a 100x100 array
  x=numpy.outer(scipy.cos(u),scipy.sin(v))
  y=numpy.outer(scipy.sin(u),scipy.sin(v))
  z=numpy.outer(numpy.ones(u.size),scipy.cos(v))
  return x,y,z

#-------------------------------------------------------------------------------
def dist(v1,v2):
  v1=numpy.array(v1).flatten()
  v2=numpy.array(v2).flatten()
  return numpy.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + \
              (v2[2]-v1[2])**2)

#-------------------------------------------------------------------------------
def remove_overlaps(v):
  newvec=v[0]
  i=0
  #print "%d: %s -> %s"%(i,0,v[i])
  for i in xrange(1,v.shape[0]):
    d=dist(newvec[-1],v[i])
    #print "%d: %s -> %s"%(i,d,v[i])
    if d>1e-5:
      newvec=numpy.vstack([newvec,v[i]])
  return newvec

#-------------------------------------------------------------------------------
def drawcoil(base, coilobj, nbase, xc1, nc1, thc1, xc2, nc2, xp, yp, zp, bctype,
             helical,random, shades, render, ccap, cdye, direction, seed, ax3d,
             loopcenter, loopradius, coiltype, gg, use_start=False, use_end=False,
             helpers=True):
  # nbase             # number of bases in helix
  # xc1(3)            # start coil position
  # nc1(3)            # start coil vector (into coil)
  # thc1              # start coil rotation around nc1
  # xc2(3)            # end coil position
  # nc2(3)            # end coil vector (into coil)
  # xp,yp,zp          # points to be splined between end points
  # bctype(2)         # 5' and 3' ends (1: first deriv bc, 2: second deriv bc)
  # helical           # make a helix around the spline
  # random            # make the helix angular progression and base orientations random
                      # in absence of helix, this still randomizes base
                      # orientations around spline
  # shades(2)         # map indices for coil, base pairs
  # render(4)         # render coil, first base, middle bases, last base
  # ccap(2)           # cap 5', 3' end of coil
  # cdye(2)           # dye 5', 3' end of coil
  # direction         # 1: 5' end at x1, 2: 5' end at x2
  # seed              # set seed if random


  # Put coordinates in base objects
  xc_p = numpy.array([xp, yp, zp])
  if len(xp):
    if render[1]>=0 and render[3]>=0:
      st=1
      base[coilobj.bases[0]-1].x3=xc_p[:,0]
    else:
      st=0

    if render[1]>=0 and render[3]<=-1:
      st=1

    for bi in coilobj.bases:
      base[bi].x3=xc_p[:,st]
      st=st+1

    if render[1]>=0 and render[3]>=0:
      base[coilobj.bases[-1]+1].x3=xc_p[:,-1]

  numpy.random.seed(seed)
  nhtot=[0,0]
  nhtot[0] = (gg.nh[0]-1)*(nbase-1) + 1  # total points along chain
  nhtot[1] = gg.nh[1]   # total points around chain

  SPLINE = 0

  ksbc1 = bctype[0]
  ksbc2 = bctype[1]
  xc_p = numpy.array([xp, yp, zp])

  if len(zp)<2: # Let the outer spline take care of this case
    print "Too few interpolation points, returning..."
    return 0., xc_p.T

  arcratio=0.0
  enough_points=True # just for now

  if enough_points:
    min=zp[0]
    max=zp[-1]

    if scipy.fabs(max-min)<1e-10:
      zp=numpy.ones(len(zp))*max
    else:
      n=len(zp)
      step=(max-min)/(n-1)
      zp=numpy.arange(zp[0],zp[-1]+step/2,step)

    sp=numpy.zeros(xp.size)
    sp[0] = 0.

    for i in range(1,xp.size):
        sp[i] = sp[i-1] + numpy.sqrt((xp[i]-xp[i-1])**2 + (yp[i]-yp[i-1])**2 + \
                (zp[i]-zp[i-1])**2)

    xbc1 = [nc1[0], nc1[1], nc1[2]]
    xbc2 = [-nc2[0], -nc2[1], -nc2[2]]

    ss = numpy.linspace(0,1,nhtot[0])

    try:
      xr,yr,zr, dxr, dyr, dzr=spl(SPLINE,sp,[xp, yp, zp],
        ss,xbc1,xbc2,ksbc1,ksbc2, order=3, stiff=4.0)
    except SystemError, v:
      print "Doing polyfit:"
      xr, yr, zr, dxr, dyr, dzr=polyfit([xp, yp, zp],order=2, n=nhtot[0])
    except TypeError, v:
      print "Doing polyfit:"
      xr, yr, zr, dxr, dyr, dzr=polyfit([xp, yp, zp],order=2, n=nhtot[0])

    if helpers and ax3d:
      ax3d.addPolyCylinder(numpy.array([xs,ys,zs]).T,colors=Export.colors["red"], radius=0.25)

    xc_r = numpy.array([xr, yr, zr])

  #-----------------------------------------------------------------------------
  # Now draw the bases

  divtotal=nbase-1
  line_int=len(xr)/(divtotal)
  if render[1]<-1:
    si=0
  else:
    si=1

  li=[]
  i=0
  for bi in coilobj.bases:
    ind=si*line_int
    if ind>0: ind-=1
    if ind>=len(xr): ind=len(xr)-1
    li.append(ind)
    base[bi].x3 = xc_r[:,ind]
    si+=1

  strand_points=xc_r.T
  strand_derivs=numpy.array([dxr,dyr,dzr]).T

  # Now calculate base angles
  thc_a=numpy.zeros(len(li))
  thc_a[0] = thc1*scipy.pi/180

  dthb_this=-gg.dthb

  for i in range(1,len(li)):
    if random:
      #~ dthb_rnd   = gg.dthb + numpy.random.randn(1)*numpy.sqrt(scipy.pi/8)
      dthb_rnd   = gg.dthb + numpy.random.randn(1)*numpy.sqrt(2*scipy.pi)
      #~ print "dthb_rnd=%s"%dthb_rnd
      thc_a[i]      = dthb_rnd
    else:
      thc_a[i]      = thc_a[i-1] + dthb_this

  xb = numpy.array([0.,0.3 * gg.rdh+gg.rhc])
  yb = numpy.zeros(len(xb))
  zb = yb.copy()


  for j in range(len(li)):
    i = li[j]
    l=j
    # 0) normalize chain axis vector
    ncj = numpy.array([strand_derivs[i,0], strand_derivs[i,1], strand_derivs[i,2]])
    ncj = ncj/numpy.linalg.norm(ncj)
    ncjm = numpy.matrix(ncj).T
    # 2) rotate helix axis to vector u_n \equiv ncj[2]
    # axis starts out as u_z
    # rotation is around vector u_rot = u_z x u_n
    sinth_rot = numpy.sqrt(ncj[0]**2 + ncj[1]**2)
    costh_rot = ncj[2]
    if sinth_rot>0:
      u_rot = numpy.array([-ncj[1], ncj[0], 0])
      u_rot = u_rot/numpy.linalg.norm(u_rot) # make unit vectors
      u_rotm=numpy.matrix(u_rot).T
      th_rot = scipy.arctan2(sinth_rot,costh_rot)
      #~ print "th_rot=%f"%(th_rot*180./scipy.pi)
      # rotate(h,u_rot,th_rot,[0 0 0]); # "rotate" won't work properly if xc1 \neq 0
      # (since intrinsic function won't do translation)

      # th_rot needs to be reversed compared to value for using
      # matlab intrinsic function "rotate"
      rmat    = scipy.cos(-th_rot)*I3 \
              + (1-scipy.cos(-th_rot))*u_rotm*u_rotm.T \
              +    scipy.sin(-th_rot) * \
              numpy.matrix([[ 0,         u_rot[2], -u_rot[1]], \
                            [-u_rot[2],  0,         u_rot[0]], \
                            [ u_rot[1], -u_rot[0],  0]])
    else:
      rmat  = I3

    rmat1     = scipy.cos(thc_a[l])*I3 + (1-scipy.cos(thc_a[l]))*ncjm*ncjm.T \
                + scipy.sin(thc_a[l])* \
                 numpy.matrix([[ 0,       ncj[2], -ncj[1]], \
                               [-ncj[2],  0,       ncj[0]], \
                               [ ncj[1], -ncj[0],  0]]) # rotate around chain tangent
    # rmat=I3
    xbtmp=numpy.zeros(xb.shape[0])
    ybtmp=numpy.zeros(xb.shape[0])
    zbtmp=numpy.zeros(xb.shape[0])

    xbtmp = numpy.array([0,ncj[0]]) + strand_points[i,0]
    ybtmp = numpy.array([0,ncj[1]]) + strand_points[i,1]
    zbtmp = numpy.array([0,ncj[2]]) + strand_points[i,2]

    if helpers:
      ax3d.addPolyCylinder(numpy.array([xbtmp.flatten(),ybtmp.flatten(), \
          zbtmp.flatten()]).T, colors=Export.colors["yellow"], radius=0.25)#,0*zbtmp+shades[-1])


    for k in range(xb.shape[0]):
        xtmp = rmat1*rmat*numpy.matrix([xb[k], yb[k], zb[k]]).T
        xbtmp[k] = xtmp[0,:] + strand_points[i,0]
        ybtmp[k] = xtmp[1,:] + strand_points[i,1]
        zbtmp[k] = xtmp[2,:] + strand_points[i,2]  #translate to base location

    ax3d.addPolyCylinder(numpy.array([xbtmp.flatten(),ybtmp.flatten(), \
        zbtmp.flatten()]).T, colors=Export.colors["light_gray"], radius=gg.rbc[0])#,0*zbtmp+shades[-1])
    ax3d.addSphere(numpy.array([xbtmp[1],ybtmp[1], zbtmp[1]]),
        color=Export.colors["light_gray"], radius=gg.rbc[0])

  coilobj.drawn_bases=True

  if render[1]>-1 and enough_points: si=gg.nh[0]/2
  else: si=0
  if render[3]>-1  and enough_points: ei=-gg.nh[0]/2
  else: ei=-1

  return arcratio, xc_r.T


#-------------------------------------------------------------------------------
def  drawcoil_bases(coil, nbase, thc1, shades, render, ax3d, strand_points,
                    strand_derivs, start, end, int_start, int_end,
                    gg, use_start=True, use_end=True, ipb=1, random=False,
                    helpers=True):

  # nbase             # number of bases in helix
  # thc1              # start coil rotation around nc1
  # points          # points to be splined between end points
  # shades(2)         # map indices for coil, base pairs
  # render(4)         # render coil, first base, middle bases, last base
  # ax3d              # Coil instance
  # strand points     # Points along interpolated strand
  # strand_derivs     # Gradient of strand
  # start             # Start position of coil
  # end               # End position of coil
  # int_start         # Start index of the search interval



  #~ global nh, nb, rhc, rbc, rhelix, rdh, rcap, ncap, thcap, dzb, bturn, brise, \
        #~ dthb, dsb, dz, dth, cmap, rdye, ddye, ndye, majorgroove, minorgroove, \
        #~ dthgroove, strutangle, strutlength, strutrise, proptwist, inclination, \
        #~ extend5prime, extend3prime, extendcoil, rcaphead, colors

  numpy.random.seed()
  npoints=strand_points.shape[0]

  int_end=min(npoints-1,int_end)
  int_start=max(0,int_start)
  di=3*(int_end-int_start)/nbase

  imin=dmin=0#max(0, int_start-di)
  dmax=npoints-1#min(npoints-1,int_start+di)
  centers=[]

  if use_start:
    imin,centers=findpoint(start, strand_points, int_start-di, int_start+di,
          imin=dmin,imax=dmax, debug=helpers)
  if helpers:
    ax3d.addSphere(strand_points[dmin], color=Export.colors["red"], radius=0.7)
    ax3d.addSphere(strand_points[dmax], color=Export.colors["brown"], radius=0.7)
    if use_start:
      for ci in centers:
        ax3d.addSphere(strand_points[ci], color=Export.colors["light_green"], radius=0.7)

  dmin=0#max(0, int_end-di)
  dmax=npoints-1#min(npoints-1,int_end+di)

  if use_end:
    imax,centers=findpoint(end, strand_points, int_end-di, int_end+di,
        imin=dmin, imax=dmax, debug=helpers)
  else:
      imax=min(npoints-1,imin+nbase*ipb)
      imax=npoints-1
  if not use_start and coil.bases[0]>0:
    imin=max(0,imax-nbase*ipb)
  #~ print "imin=%d, imax=%d, nbase=%s"%(imin, imax, nbase)


  if helpers:
    #~ print "render=%s"%render
    ax3d.addSphere(start, color=Export.colors["white"], radius=0.4)
    ax3d.addSphere(end, color=Export.colors["red"], radius=0.4)

    ax3d.addSphere(strand_points[dmin], color=Export.colors["dark_green"], radius=0.7)
    for ci in centers:
      ax3d.addSphere(strand_points[ci], color=Export.colors["turquoise"], radius=0.7)
    ax3d.addSphere(strand_points[dmax], color=Export.colors["dark_blue"], radius=0.7)


    ax3d.addSphere(strand_points[imin], color=Export.colors["yellow"], radius=2)
    ax3d.addSphere(strand_points[imax], color=Export.colors["orange"], radius=0.7)

  #~ print "nbase=",nbase

  #~ if imin>0:
    #~ imin=imin-ipb
    #~ imax=imax-ipb

  # Now calculate base angles
  thc_a=numpy.zeros(nbase)
  thc_a[0] = thc1*scipy.pi/180

  this_npoints=imax-imin+1
  this_ipb=this_npoints/nbase
  #~ print "npoints=%d imin=%d, imax=%d"%(npoints,imin,imax)
  #~ print "ipb=%d this_ipb=%d"%(ipb,this_ipb)
  dthb_this=-gg.dthb*this_ipb/ipb
  #~ print "dthb=%f dthb_this=%f"%(dthb, dthb_this)


  for i in range(1,nbase):
    if random:
      #~ dthb_rnd   = gg.dthb + numpy.random.randn(1)*numpy.sqrt(scipy.pi/8)
      dthb_rnd   = gg.dthb + numpy.random.randn(1)*numpy.sqrt(2*scipy.pi)
      #~ print "dthb_rnd=%s"%dthb_rnd
      thc_a[i]      = dthb_rnd
    else:
      thc_a[i]      = thc_a[i-1] + dthb_this

  #~ print "thc_a=%s"%thc_a
  # define base pair
  #
  xb = numpy.array([0.,0.3 * gg.rdh+gg.rhc])
  #~ print "xb=%s"%(xb)
  yb = numpy.zeros(len(xb))
  zb = yb.copy()

  # base pair struts
  #~ print "base pair struts"
  cpoints=imax-imin+1
  ind=numpy.linspace(imin,imax,nbase,endpoint=True)
  #~ print npoints,nbase,ind
  bi=coil.bases[0]
  #~ print coil.bases
  for j in range(nbase):
    if j==0 and nbase==2 and not use_end: continue
    if j==nbase-1 and nbase==2 and not use_start: continue
    i = int(ind[j])

    l=j%nbase
    #~ print i,j
    # 0) normalize chain axis vector
    ncj = numpy.array([strand_derivs[i,0], strand_derivs[i,1], strand_derivs[i,2]])
    ncj = ncj/numpy.linalg.norm(ncj)
    ncjm = numpy.matrix(ncj).T
    # 2) rotate helix axis to vector u_n \equiv ncj[2]
    # axis starts out as u_z
    # rotation is around vector u_rot = u_z x u_n
    sinth_rot = numpy.sqrt(ncj[0]**2 + ncj[1]**2)
    costh_rot = ncj[2]
    if sinth_rot>0:
      u_rot = numpy.array([-ncj[1], ncj[0], 0])
      u_rot = u_rot/numpy.linalg.norm(u_rot) # make unit vectors
      u_rotm=numpy.matrix(u_rot).T
      th_rot = scipy.arctan2(sinth_rot,costh_rot)
      #~ print "th_rot=%f"%(th_rot*180./scipy.pi)
      # rotate(h,u_rot,th_rot,[0 0 0]); # "rotate" won't work properly if xc1 \neq 0
      # (since intrinsic function won't do translation)

      # th_rot needs to be reversed compared to value for using
      # matlab intrinsic function "rotate"
      rmat    = scipy.cos(-th_rot)*I3 \
              + (1-scipy.cos(-th_rot))*u_rotm*u_rotm.T \
              +    scipy.sin(-th_rot) * \
              numpy.matrix([[ 0,         u_rot[2], -u_rot[1]], \
                            [-u_rot[2],  0,         u_rot[0]], \
                            [ u_rot[1], -u_rot[0],  0]])
    else:
      rmat  = I3
    #~ if j==0: print"rmat=\n%s"%(rmat)

    rmat1     = scipy.cos(thc_a[l])*I3 + (1-scipy.cos(thc_a[l]))*ncjm*ncjm.T \
                + scipy.sin(thc_a[l])* \
                 numpy.matrix([[ 0,       ncj[2], -ncj[1]], \
                               [-ncj[2],  0,       ncj[0]], \
                               [ ncj[1], -ncj[0],  0]]) # rotate around chain tangent
    #rmat=I3
    xbtmp=numpy.zeros(xb.shape[0])
    ybtmp=numpy.zeros(xb.shape[0])
    zbtmp=numpy.zeros(xb.shape[0])



    xbtmp = numpy.array([0,ncj[0]]) + strand_points[i,0]
    ybtmp = numpy.array([0,ncj[1]]) + strand_points[i,1]
    zbtmp = numpy.array([0,ncj[2]]) + strand_points[i,2]

    if helpers:
      ax3d.addPolyCylinder(numpy.array([xbtmp.flatten(),ybtmp.flatten(), \
          zbtmp.flatten()]).T, colors=Export.colors["yellow"], radius=0.25)#,0*zbtmp+shades[-1])


    for k in range(xb.shape[0]):
        xtmp = rmat1*rmat*numpy.matrix([xb[k], yb[k], zb[k]]).T
        xbtmp[k] = xtmp[0,:] + strand_points[i,0]
        ybtmp[k] = xtmp[1,:] + strand_points[i,1]
        zbtmp[k] = xtmp[2,:] + strand_points[i,2]  #translate to base location

    if helpers or (ax3d and ((j==0 and render[1]<-1) or (j==nbase-1 and render[3]<-1) or \
        (j > 0 and  j < nbase-1 and render[2]>0))):
      #~ print xbtmp,ybtmp,zbtmp,gg.rbc[0]
      ax3d.addPolyCylinder(numpy.array([xbtmp.flatten(),ybtmp.flatten(), \
          zbtmp.flatten()]).T, colors=Export.colors["light_gray"],
          radius=gg.rbc[0], drawall=True)
      ax3d.addSphere(numpy.array([xbtmp[1],ybtmp[1], zbtmp[1]]),
          color=Export.colors["light_gray"], radius=gg.rbc[0])
    else:
      pass
      #~ print "Not drawing base",j
    #~ if j==0:
      #~ ax3d.addSphere(numpy.array([xbtmp[0],ybtmp[0], zbtmp[0]]),
          #~ color=Export.colors["light_gray"], radius=rhc)


  return imax

#-------------------------------------------------------------------------------
def drawhelix(base, helixobj, nbase, xc1, nc1, thc1, dzc1, shades, render, acap, bcap,
              adye, bdye, ax3d, gg, helpers=True):
  # nbase             # number of bases in helix
  # xc1(3)            # start helix axis position
  # nc1(3)            # start helix axis vector
  # thc1              # start helix rotation around axis
  # dzc1              # translation along helix axis from start of helix axis
  # shades(3)         # colormap indices for chain a, chain b, base pairs
  # render(3)         # render chain a, chain b, base pair struts
  # acap(2)           # cap 5', 3' end of a chain
  # bcap(2)           # cap 5', 3' end of b chain
  # adye(2)           # dye 5', 3' end of a chain
  # bdye(2)           # dye 5', 3' end of b chain
  #~ print "base = ",base
  #~ print "helix = ",helixobj
  #~ print "nbase = ",nbase
  #~ print "xc1 = ",xc1
  #~ print "nc1 = ",nc1
  #~ print "thc1 = ",thc1
  #~ print "dzc1 = ",dzc1
  #~ print "shades = ",shades
  #~ print "render = ",render
  #~ print "acap = ",acap
  #~ print "bcap = ",bcap
  #~ print "adye = ",adye
  #~ print "bdye = ",bdye
  #~ print "ax3d = ",ax3d


  #~ print helixobj.bases
  nhtot=[0,0]
  nhtot[0]= (gg.nh[0]-1)*(nbase-1) + 1  # total points along chain
  if nhtot[0] < 2:
    nhtot[0] = 2
  nhtot[1] = gg.nh[1]   # total points around chain

  #~print "drawhelix "+"-"*50
  #pdb.set_trace()
  # 0) numpy.linalg.normalize target helix axis vector
  n=numpy.linalg.norm(nc1)
  if n>0:
    nc1 = nc1/n

  dz = gg.dzb/(gg.nh[0]-1)
  dth = 2*scipy.pi/(gg.nh[1]-1)
  x_a = numpy.zeros(nhtot)
  y_a = numpy.zeros(nhtot)
  z_a = numpy.zeros(nhtot)
  x_b = numpy.zeros(nhtot)
  y_b = numpy.zeros(nhtot)
  z_b = numpy.zeros(nhtot)

  th=numpy.arange(0, 2*scipy.pi+dth/2, dth)
  step=dz*gg.dthb/gg.dzb


  if nbase == 1:
    thc_a      = numpy.arange(0, 1.1*step, step)
    xc_a = numpy.zeros([3,thc_a.size])
    xc_b = numpy.zeros([3,thc_a.size])

    xc_a[0,:]  = gg.rdh*scipy.cos(thc_a)
    xc_a[1,:]  = gg.rdh*scipy.sin(thc_a)
                            # move start of a chain so rise due to inclination of
                            # base pair is centered on helix origin
    xc_a[2,:]  = numpy.arange(0, 1.1*dz, dz) - .5*gg.strutrise

    thc_b      = gg.dthgroove + numpy.arange(0, 1.1*step, step)
    xc_b[0,:]  = gg.rdh*scipy.cos(thc_b)
    xc_b[1,:]  = gg.rdh*scipy.sin(thc_b)
                          # move start of chain b so rise due to inclination of
                          # base pair is centered on the origin in the z direction
    xc_b[2,:]  = numpy.arange(0, 1.1*dz, dz) + .5*gg.strutrise
  else:
    thc_a      = numpy.arange(0, gg.dthb*(nbase-1)+step/2, step)
    xc_a = numpy.zeros([3,thc_a.size])
    xc_b = numpy.zeros([3,thc_a.size])

    xc_a[0,:]  = gg.rdh*scipy.cos(thc_a)
    xc_a[1,:]  = gg.rdh*scipy.sin(thc_a)
                            # move start of a chain so rise due to inclination of
                            # base pair is centered on helix origin
    xc_a[2,:]  = numpy.arange(0, gg.dzb*(nbase-1)+dz/2, dz) - .5*gg.strutrise

    thc_b      = gg.dthgroove + numpy.arange(0, gg.dthb*(nbase-1)+step/2, step)
    xc_b[0,:]  = gg.rdh*scipy.cos(thc_b)
    xc_b[1,:]  = gg.rdh*scipy.sin(thc_b)
                          # move start of chain b so rise due to inclination of
                          # base pair is centered on the origin in the z direction
    xc_b[2,:]  = numpy.arange(0, gg.dzb*(nbase-1)+dz/2, dz) + .5*gg.strutrise

  if helpers and ax3d:
    ax3d.addPolyCylinder(numpy.array([xc_a[0], xc_a[1], xc_a[2]]).T,
        colors=Export.colors["shady_blue"], radius=gg.rhc)
    ax3d.addPolyCylinder(numpy.array([xc_b[0], xc_b[1], xc_b[2]]).T,
        colors=Export.colors["shady_green"], radius=gg.rhc)
    minp=-20.
    maxp=20
    n=11
    d=maxp-minp
    step=d/(n-1)
    points3=numpy.zeros([n,3])
    points3[:,2]=points3[:,1]=numpy.zeros(n)
    z=numpy.arange(minp,maxp+step/2,step)
    #~ print z
    points3[:,0]=z
    #~ print points3
    ax3d.addPolyCylinder(points3, radius=1,colors=Export.colors["white"])

  #
  # define backbone
  #

  phi = scipy.pi/2 - scipy.arctan(gg.dzb/(gg.rdh*gg.dthb))
  x3_a=y3_a=z3_a=x5_b=y5_b=z5_b=None

  # convenient to rotate surface using matlab function rotate
  # however, still need to keep track of end positions and vectors using
  # rotation matrices, hence, might be more consistent just to explicitly
  # compute rotation matrix and do everything manually
  #
  # actually, would be useful reference check to keep moving surfaces using
  # "rotate" and move end info manually, unfortunately, since there is no
  # "translate" equivalent for translation, have to translate manually before
  # "rotate" and this makes it messy to rotate the end points since have
  # to change origin of rotation for them
  #
  # decided just to do everything manually in the end
  #


  # 1) first rotate around z axis amount thc1
  # rotate(h,[0 0 1],thc1,[0 0 0]);  # "rotate" won't work properly if xc1 \neq  (since "rotate" won't do translation)
  # (sign of sin terms seems reversed to me...????)
  rmat1   = numpy.matrix([[scipy.cos(thc1*scipy.pi/180), \
                          -scipy.sin(thc1*scipy.pi/180), 0], \
                         [scipy.sin(thc1*scipy.pi/180), \
                          scipy.cos(thc1*scipy.pi/180), 0], \
                         [0, 0, 1]])  # rotate around z axis

  # 2) rotate helix axis to vector u_n \equiv nc1(3)
  # axis starts out as u_z
  # rotation is around vector u_rot = u_z x u_n
  sinth_rot = numpy.sqrt(nc1[0]**2 + nc1[1]**2)
  costh_rot = nc1[2]
  if sinth_rot>0:
     u_rot = numpy.matrix([-nc1[1], nc1[0], 0]).T
     u_rot = u_rot/numpy.linalg.norm(u_rot) # make unit vectors
     th_rot = 180./scipy.pi * scipy.arctan2(sinth_rot,costh_rot)
     # rotate(h,u_rot,th_rot,[0 0 0]); # "rotate" won't work properly if xc1 \neq 0
     # (since intrinsic function won't do translation)

     # th_rot needs to be reversed compared to value for using
     # matlab intrinsic function "rotate"
     rmat2   = scipy.cos(-th_rot*scipy.pi/180)*I3  \
             + (1-scipy.cos(-th_rot*scipy.pi/180))*u_rot*u_rot.T \
             +    scipy.sin(-th_rot*scipy.pi/180) * \
             numpy.matrix([[0       ,  u_rot[2],  -u_rot[1]], \
                          [-u_rot[2],  0       ,   u_rot[0]], \
                          [u_rot[1] , -u_rot[0],         0]])
  elif costh_rot == -1:  # need special case for u_n = [0; 0; -1]
     u_rot = numpy.matrix([0, 1, 0]).T
     th_rot = 180.
     rmat2   = scipy.cos(-th_rot*scipy.pi/180)*I3  \
             + (1-scipy.cos(-th_rot*scipy.pi/180))*u_rot*u_rot.T \
             +    scipy.sin(-th_rot*scipy.pi/180) * \
             numpy.matrix([[0        , u_rot[2], -u_rot[1]], \
                           [-u_rot[2], 0       ,  u_rot[0]], \
                           [ u_rot[1], -u_rot[0], 0      ]])
  else: # special case for u_n = [0; 0; 1]
     rmat2 = numpy.matrix(I3)

  # 3) then translate the helix

  # chains
  for j in range(xc_a.shape[1]):
      xtmp = numpy.array(rmat2*rmat1*numpy.matrix([xc_a[0,j], xc_a[1,j], \
              xc_a[2,j]]).T).flatten() + xc1 + nc1*dzc1
      xc_a[0,j] = xtmp[0]
      xc_a[1,j] = xtmp[1]
      xc_a[2,j] = xtmp[2]
      xtmp = numpy.array(rmat2*rmat1*numpy.matrix([xc_b[0,j], xc_b[1,j], \
              xc_b[2,j]]).T).flatten() + xc1 + nc1*dzc1
      xc_b[0,j] = xtmp[0]
      xc_b[1,j] = xtmp[1]
      xc_b[2,j] = xtmp[2]


  # base pair struts
  #~ print "Calculating base positions"
  for j in range (1, nbase+1):
    i = (j-1)*(gg.nh[0]-1)
    bar=numpy.array([[xc_a[0,i], xc_a[1,i], xc_a[2,i]],
                               [xc_b[0,i], xc_b[1,i], xc_b[2,i]]])
    base[helixobj.bases[j-1][0]].x3=xc_a[:,i]
    base[helixobj.bases[j-1][1]].x3=xc_b[:,i]
    ax3d.addPolyCylinder(bar, radius=gg.rbc[0],colors=Export.colors["light_gray"])


  ### Draw cylinders
  if False and helpers and ax3d:
    ax3d.addPolyCylinder(numpy.array([xc_a[0], xc_a[1], xc_a[2]]).T,
        colors=Export.colors["yellow"], radius=gg.rhc)
    ax3d.addPolyCylinder(numpy.array([xc_b[0], xc_b[1], xc_b[2]]).T,
        colors=Export.colors["yellow"], radius=gg.rhc)

  offset=gg.nh[0]
  x1a  =  numpy.matrix(xc_a[:,0]).T

  n1a  = numpy.matrix([0, -gg.rdh*gg.dthb/numpy.sqrt(gg.dzb**2 + (gg.rdh*gg.dthb)**2), \
          -gg.dzb/numpy.sqrt(gg.dzb**2 + (gg.rdh*gg.dthb)**2)])

  x1b  = numpy.matrix(xc_b[:,0]).T
  ca   = scipy.cos(2*scipy.pi-gg.dthgroove)  # rotation matrix is cw but groove angle is ccw
  sa   = scipy.sin(2*scipy.pi-gg.dthgroove)  # rotate around z axis
  rmat = ca*I3 + (1-ca)*numpy.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]]) \
               + sa*numpy.matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
  n1b  = rmat*n1a.T

  # end chain information
  xc2     = numpy.matrix([0, 0, gg.dzb*(nbase-1)])  # end helix axis
  nc2     = nc1                                     # end helix axis vector

  thc2    = gg.dthb*(nbase-1)*180/scipy.pi          # end helix rotation around axis
  rmat    = numpy.matrix([[scipy.cos(thc2*scipy.pi/180), \
                          -scipy.sin(thc2*scipy.pi/180), 0], \
                         [scipy.sin(thc2*scipy.pi/180), \
                          scipy.cos(thc2*scipy.pi/180), 0], \
                         [0, 0, 1]])  # rotate around z axis
  x2a=numpy.matrix(xc_a[:,-1]).T

  n2a     = -rmat*n1a.T                  # end helix chain vector 5'->3' chain
  x2b=numpy.matrix(xc_b[:,-1]).T
  n2b     = -rmat*n1b                    # end helix chain vector 3'->5' chain

  n1a = numpy.array(n1a/numpy.linalg.norm(n1a)).flatten()  # normalize normal vectors
  n1b = numpy.array(n1b/numpy.linalg.norm(n1b)).flatten()
  n2a = numpy.array(n2a/numpy.linalg.norm(n2a)).flatten()
  n2b = numpy.array(n2b/numpy.linalg.norm(n2b)).flatten()

  thc2    = thc2 + thc1
  return xc2,nc2,thc2,x1a,n1a,x2a,n2a,x1b,n1b,x2b,n2b, nhtot, numpy.array([xc_a[0], xc_a[1], xc_a[2]]).T, numpy.array([xc_b[0], xc_b[1], xc_b[2]]).T

#-------------------------------------------------------------------------------
def  drawcap(xc_a, extend, rhc, dzb, ax3d, color, arrow=True):
  # xc_a(3)
  # extend            amount to extend chain at end of strand before cap
  #                   differs for 3' and 5' ends in helix and in coil
  # rhc               radius of polycylinder
  dvec=xc_a[-1]-xc_a[-2]
  end_vec=dvec/numpy.linalg.norm(dvec)
  endpt=xc_a[-1]+end_vec*extend
  lpts=numpy.vstack([xc_a[-4:],endpt])


  # Fit spline through points
  spoints=lpts.shape[0]
  sp=numpy.zeros(spoints)
  sp[0] = 0.

  xp=numpy.array(lpts[:,0]).flatten()
  yp=numpy.array(lpts[:,1]).flatten()
  zp=numpy.array(lpts[:,2]).flatten()

  for i in range(1,spoints):
      sp[i] = sp[i-1] + numpy.sqrt((xp[i]-xp[i-1])**2 + (yp[i]-yp[i-1])**2 + \
              (zp[i]-zp[i-1])**2)

  xbc1 = xbc2 = numpy.zeros(3)

  ss = numpy.linspace(0,1.,spoints*5)

  xs,ys,zs, dxs, dys, dzs=spl("",sp,[xp, yp, zp], ss ,xbc1, xbc2, 0.,0., order=2, stiff=10.0)
  fitpts=numpy.array([xs, ys, zs]).T

  ax3d.addPolyCylinder(fitpts, colors=color, radius=rhc, drawall=True)

  if arrow:
    r=0.3*rhc
    max=2.0
    a=r*0.66745721602838382 #r*scipy.arccos(scipy.py/4)

    th=numpy.linspace(8*scipy.pi/8,3*scipy.pi/8, 10)
    x=3*(r*a)*(scipy.cos(th)+1)
    y=r*scipy.sin(th)+max       # y is the cone radius

    rpts1=numpy.zeros([len(x),3])
    for i in [0,1,2]:
      rpts1[:,i]=end_vec[i]*x

    th2=numpy.linspace(3*scipy.pi/8, 0., 10)
    y2=1.25*r*scipy.sin(th2)

    x2=r*scipy.cos(th2)
    rpts2=numpy.zeros([len(x2),3])
    for i in [0,1,2]:
      rpts2[:,i]=end_vec[i]*(x2+rhc*1.5+x[-1])


    conepts=numpy.vstack([
      fitpts[-2],
      fitpts[-1]+rpts1,
      fitpts[-1]+rpts2,
      ])
    radii=numpy.hstack([numpy.array([0.0]),y,y2])
  else:
    r=0.3*rhc
    max=rhc-r
    a=r*0.66745721602838382 #r*scipy.arccos(scipy.py/4)

    th=numpy.linspace(4*scipy.pi/8,0, 10)
    x=r*(scipy.cos(th))
    y=r*scipy.sin(th)+max #  y is the cone radius

    rpts1=numpy.zeros([len(x),3])
    for i in [0,1,2]:
      rpts1[:,i]=end_vec[i]*x

    conepts=numpy.vstack([
      fitpts[-2],
      fitpts[-1]+rpts1,
      ])

    radii=numpy.hstack([numpy.array([rhc]),y,numpy.array([0.0])])

  ax3d.addCone(conepts, radii, colors=color, nsides=20)
