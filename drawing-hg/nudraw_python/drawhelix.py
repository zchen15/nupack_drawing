#
# nudraw::drawhelix
# 
# A program for visualizing DNA Structures and Devices
# Copyright 2004
# Niles A. Pierce
# California Institute of Technology
# Version 1.0 coded 2004 (April 9,10,13,14)
#
# Converted to Python/Numpy/SciPy/PyLab Jul 30, 2008
# Conrad Steenberg <conrad.steenberg@caltech.edu>

import numpy
import scipy
import pylab

from drawcap import *


#----------------------------------------------------------------------------------------
def drawhelix(nbase, xc1, nc1, thc1, dzc1, shades, render, acap, bcap, adye, bdye):
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

  global nh, nb, rhc, rbc, rhelix, rdh, rcap, ncap, thcap, dzb, bturn, brise, \
         dthb, dsb, dz, dth, cmap, rdye, ddye, ndye, majorgroove, minorgroove, \
         dthgroove, strutangle, strutlength, strutrise, proptwist, inclination, \
         extend5prime, extend3prime, extendcoil, rcaphead, I3

  nhtot=[0,0]
  nhtot[0] = (nh[0]-1)*(nbase-1) + 1  # total points along chain
  nhtot[1] = nh[1]   # total points around chain

  helpers = False

  # 0) pylab.normalize target helix axis vector
  nc1 = nc1/pylab.norm(nc1)

  dz = dzb/(nh[0]-1)
  dth = 2*scipy.pi/(nh[2]-1)
  x_a = numpy.zeros(nhtot)
  xc_a = numpy.zeros([3,nhtot[0]])
  xc_b = numpy.zeros([3,nhtot[0]])
  y_a = numpy.zeros(nhtot)
  z_a = numpy.zeros(nhtot)
  x_b = numpy.zeros(nhtot)
  y_b = numpy.zeros(nhtot)
  z_b = numpy.zeros(nhtot)

  x3_a = None
  x3_b = None
  x5_a = None
  x5_b = None
  

  th = numpy.array([0, dth, 2*scipy.pi])

  thc_a      = numpy.array([0, dz*dthb/dzb, dthb*(nbase-1)])
  xc_a[0,:]  = rdh*scipy.cos(thc_a)
  xc_a[1,:]  = rdh*scipy.sin(thc_a)
  xc_a[2,:]  = numpy.array([0, dz, dzb*(nbase-1)]) - .5*strutrise # move start of a chain so rise due to inclination of base pair is centered
                          # on helix origin

  thc_b      = dthgroove + numpy.array([0, dz*dthb/dzb, dthb*(nbase-1)])
  xc_b[0,:]  = rdh*scipy.cos(thc_b)
  xc_b[1,:]  = rdh*scipy.sin(thc_b)
  xc_b[2,:]  = numpy.array([0, dz, dzb*(nbase-1)]) + .5*strutrise  # move start of chain b so rise due to inclination of base pair is 
                          # centered on the origin in the z direction
  #
  # define backbone
  #

  phi = scipy.pi/2 - ascipy.arctan(dzb/(rdh*dthb))
  for i in range(nhtot[0]):
    xtmp[0,:] = rhc*scipy.cos(th)
    xtmp[1,:] = rhc*scipy.sin(th)
    xtmp[2,:] = 0  
    ca = scipy.cos(thc_a[i])
    sa = scipy.sin(thc_a[i])  # rotate circles so that cross section is circular normal to chain axis instead of helix axis
    # see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
    rmat_a = scipy.cos(phi)*I3 \
          + (1-scipy.cos(phi))*numpy.matrix([[ca**2, ca*sa, 0],[ca*sa, sa**2, 0], [0, 0,  0]]) \
          + scipy.sin(phi)*numpy.matrix([[0, 0, -sa], [0, 0, ca], [sa, -ca, 0]])
    xtmp = rmat_a*xtmp
    x_a[i,:] = xtmp[0,:] + xc_a[0,i] # translate
    y_a[i,:] = xtmp[1,:] + xc_a[1,i]
    z_a[i,:] = xtmp[2,:] + xc_a[2,i] 


    xtmp[0,:] = rhc*scipy.cos(th)
    xtmp[1,:] = rhc*scipy.sin(th)
    xtmp[2,:] = 0  

    cb = scipy.cos(thc_b[i])
    sb = scipy.sin(thc_b[i])

    # see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
    rmat_b = scipy.cos(phi)*I3 \
          + (1-scipy.cos(phi))*\
          numpy.matrix([[cb**2, cb*sb, 0], [cb*sb, sb**2, 0], [0, 0, 0]]) \
          + scipy.sin(phi)*\
          numpy.matrix([[0, 0, -sb], [0, 0, cb], [sb, -cb, 0]])

    xtmp = rmat_b*xtmp
    x_b[i,:] = xtmp[0,:] + xc_b[0,i] # translate
    y_b[i,:] = xtmp[1,:] + xc_b[1,i]
    z_b[i,:] = xtmp[2,:] + xc_b[2,i]

    if i == nhtot[0]: # must go before i==1 case so don't overwrite rmat_a and rmat_b
     endtype = 3
     x3_a,y3_a,z3_a=  drawcap(i,xc_a,rmat_a,endtype,extend3prime)
     endtype = 5
     x5_b,y5_b,z5_b =  drawcap(i,xc_b,rmat_b,endtype,extend5prime)

    if i == 0:
      # see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
      rmat_a = scipy.cos(phi+scipy.pi)*I3 \
             + (1-scipy.cos(phi+scipy.pi))*numpy.matrix([[ca**2, ca*sa, 0], [ca*sa, sa**2, 0], [0, 0, 0]]) \
             + scipy.sin(phi+scipy.pi)*numpy.matrix([[0, 0, -sa], [0, 0, ca], [sa, -ca, 0]])

      # see p. 58 of "Dynamics": crandall, karnopp, kurtz and pridmore-brown
      rmat_b = scipy.cos(phi+scipy.pi)*I3 \
             + (1-scipy.cos(phi+scipy.pi))* \
             numpy.matrix([[cb**2, cb*sb, 0], [cb*sb, sb**2, 0], [0, 0, 0]]) \
             + scipy.sin(phi+scipy.pi)* \
             numpy.matrix([[0, 0, -sb], [0, 0, cb], [sb, -cb, 0]])
      endtype = 5
      x5_a, y5_a, z5_a =  drawcap(i,xc_a,rmat_a,endtype,extend5prime)
      endtype = 3
      x3_b, y3_b, z3_b =  drawcap(i,xc_b,rmat_b,endtype,extend3prime)

  # if rotating with "rotate" intrinsic (but can't do translation)
  # h(1)     = surface(x_a,y_a,z_a, 0*z_a   +shades(1),'CDataMapping','direct'); hold on;
  # h(end+1) = surface(x5_a,y5_a,z5_a,0*z5_a+shades(1),'CDataMapping','direct');
  # h(end+1) = surface(x3_a,y3_a,z3_a,0*z3_a+shades(1),'CDataMapping','direct');
  # h(end+1) = surface(x_b,y_b,z_b, 0*z_b   +shades(2),'CDataMapping','direct'); 
  # h(end+1) = surface(x5_b,y5_b,z5_b,0*z5_a+shades(2),'CDataMapping','direct');
  # h(end+1) = surface(x3_b,y3_b,z3_b,0*z3_a+shades(2),'CDataMapping','direct');


  # define canonical base pair strut
  #
  xb = numpy.zeros(nb)
  yb = numpy.zeros(nb)
  zb = numpy.zeros(nb)




  dm1 = pylab.linspace(0, strutlength, nb[0])
  #dm1 = linspace(-rdh,rdh,nb(1))
  dm2 = pylab.linspace(0, 2*scipy.pi, nb[1]) 
  for i in range(nb[0]):
    xb[i,:] = dm1[i]
    yb[i,:] = rbc[0]*scipy.cos(dm2)
    zb[i,:] = rbc[1]*scipy.sin(dm2)

  #rotate propeller twist along the axis of the strut
  for k in range(nb[0]):
     twistang = (dm1(k) - strutlength/2)/strutlength * proptwist
     ca = scipy.cos(twistang)
     sa = scipy.sin(twistang)  #rotate around x axis in ccw direction
     rmat0 = ca*I3 + (1-ca)*numpy.matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]]) + \
              sa*numpy.matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
     for m in range(nb[1]):
        xtmp = rmat0*numpy.array([xb[k,m], [b[k,m]], zb[k,m]])
        xb[k,m] = xtmp[0]
        yb[k,m] = xtmp[1]
        zb[k,m] = xtmp[2]

  #rotate inclination around the y axis
  ca = scipy.cos(inclination)
  sa = scipy.sin(inclination)
  rmat0 = ca*I3 + (1-ca)*numpy.matrix([[0, 0, 0],[0, 1, 0], [0, 0, 0]]) + \
          sa**numpy.matrix([[0, 0, -1], [0, 0, 0], [1, 0, 0]])
  for k in range(nb[0]):
     for m in range(nb[1]):
        xtmp = rmat0*numpy.array([xb[k,m], yb[k,m], zb[k,m]])
        xb[k,m] = xtmp[0] 
        yb[k,m] = xtmp[1]
        zb[k,m] = xtmp[2]  

  #rotate and then translate canonical strut
  ca = scipy.cos(2*scipy.pi-strutangle)   # correct for fact that rotation is CW and strutangle is CCW
  sa = scipy.sin(2*scipy.pi-strutangle)   # rotate around z axis
  rmat0 = ca*I3 + (1-ca)*numpy.matrix([[0, 0, 0], [0, 0, 0] [0, 0, 1]]) + \
                        sa*numpy.matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
  for k in range(nb[0]):
     for m in range(nb[1]):
        xtmp = rmat0*numpy.array(xb[k,m], yb[k,m], zb[k,m])
        xb[k,m] = xtmp[0] + rdh # translate to canonical location  
        yb[k,m] = xtmp[1]
        zb[k,m] = xtmp[2]  #translate to proper base height


  # if rotating with "rotate" intrinsic, but can't do translation
  # for j=1:nbase
  #     i = 1 + (j-1)*(nh(1)-1);
  #    dbase = sqrt((xc_b(1,i)-xc_a(1,i))^2 + (xc_b(2,i)-xc_a(2,i))^2);
  #    nx = (xc_b(1,i)-xc_a(1,i))/dbase;
  #    ny = (xc_b(2,i)-xc_a(2,i))/dbase;
  #    nz = 0;
  #    xbtmp = xc_a(1,i) + dbase*nx*xb - ny*yb; 
  #    ybtmp = xc_a(2,i) + dbase*ny*xb + nx*yb;
  #    zbtmp = xc_a(3,i) + zb;
  #    h(end+1) = surface(xbtmp,ybtmp,zbtmp,0*zbtmp+58,'CDataMapping','direct'); hold on;
  # end

  # start chain information
  x1a  = xc_a[:,0]   # start helix chain position and vector 5'->3' chain
  n1a  = numpy.array([0, -rdh*dthb/scipy.sqrt(dzb**2 + (rdh*dthb)**2), -dzb/scipy.sqrt(dzb**2 + (rdh*dthb)^2)])
  x1b  = xc_b[:,0]   # start helix chain position and vector 3'->5' chain

  ca   = scipy.cos(2*scipy.pi-dthgroove) #rotation matrix is cw but groove angle is ccw
  sa   = scipy.sin(2*scipy.pi-dthgroove)  # rotate around z axis 
  rmat = ca*I3 + (1-ca)*numpy.matrix([[0, 0, 0], [0, 0, 0] [0, 0, 1]]) \
               + sa*numpy.matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
  n1b  = rmat*n1a

  #n1b  = [0;  rdh*dthb/sqrt(dzb^2 + (rdh*dthb)^2); -dzb/sqrt(dzb^2 + (rdh*dthb)^2)];

  # end chain information
  xc2     = numpy.array([0, 0, dzb*(nbase-1)])         # end helix axis
  nc2     = nc1                          # end helix axis vector
  thc2    = dthb*(nbase-1)*180/scipy.pi        # end helix rotation around axis
  rmat    = numpy.matrix([[scipy.cos(thc2*scipy.pi/180), -scipy.sin(thc2*scipy.pi/180), 0], \
                         [scipy.sin(thc2*scipy.pi/180), scipy.cos(thc2*scipy.pi/180), 0], \
                         [0, 0, 1]])  # rotate around z axis
  x2a     = xc_a[:,-1]                  # end helix chain position 5'->3' chain
  n2a     = -rmat*n1a                    # end helix chain vector 5'->3' chain
  x2b     = xc_b[:,-1]                  # end helix chain position 3'->5' chain
  n2b     = -rmat*n1b                    # end helix chain vector 3'->5' chain

  n1a = n1a/pylab.norm(n1a)  # normalize normal vectors
  n1b = n1b/pylab.norm(n1b)
  n2a = n2a/pylab.norm(n2a)
  n2b = n2b/pylab.norm(n2b)

  # convenient to rotate surface uscipy.sing matlab function rotate
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
  rmat1   = numpy.matrix([[scipy.cos(thc1*scipy.pi/180), -scipy.sin(thc1*scipy.pi/180), 0], \
                         [scipy.sin(thc1*scipy.pi/180), scipy.cos(thc1*scipy.pi/180), 0], \
                         [0, 0, 1]])  # rotate around z axis

  # 2) rotate helix axis to vector u_n \equiv nc1(3)
  # axis starts out as u_z
  # rotation is around vector u_rot = u_z x u_n
  sinth_rot = scipy.sqrt(nc1[0]**2 + nc1[1]**2)
  costh_rot = nc1[2] 
  if sinth_rot>0:
     u_rot = numpy.matrix([-nc1[1], nc1[0], 0]).T
     u_rot = u_rot/pylab.norm(u_rot) # make unit vectors
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
  for i in range(nhtot[0]):
     for j in range(nhtot[1]):
         xtmp = rmat2*rmat1*numpy.array([x_a[i,j], y_a[i,j], z_a[i,j]]) + \
                xc1 + nc1*dzc1
         x_a[i,j] = xtmp[0]
         y_a[i,j] = xtmp[1]
         z_a[i,j] = xtmp[2]
         xtmp = rmat2*rmat1*numpy.array([x_b[i,j], y_b[i,j], z_b[i,j]]) + \
                xc1 + nc1*dzc1
         x_b[i,j] = xtmp[0]
         y_b[i,j] = xtmp[1]
         z_b[i,j] = xtmp[2]
  
  # 5 caps
  for i in range(x5_a[0].size):
     for j in range(x5_a[1].size):
        xtmp = rmat2*numpy.array(rmat1*[x5_a[i,j], y5_a[i,j], z5_a[i,j]] )+ \
                xc1 + nc1*dzc1
        x5_a[i,j] = xtmp[0]
        y5_a[i,j] = xtmp[1] 
        z5_a[i,j] = xtmp[2]
        xtmp = rmat2*rmat2*numpy.array(rmat1*[x5_b[i,j], y5_b[i,j], z5_b[i,j]]) + \
                xc1 + nc1*dzc1
        x5_b[i,j] = xtmp[0]
        y5_b[i,j] = xtmp[1] 
        z5_b[i,j] = xtmp[2]

  # 3 caps
  for i in range(x3_a[0].size):
     for j in range(x3_a[1].size):
        xtmp = rmat2*rmat1*rmat2*numpy.array([x3_a[i,j], y3_a[i,j], z3_a[i,j]]) + \
                xc1 + nc1*dzc1
        x3_a[i,j] = xtmp[0]
        y3_a[i,j] = xtmp[1] 
        z3_a[i,j] = xtmp[2] 
        xtmp = rmat2*rmat1*rmat2*numpy.array([x3_b[i,j], y3_b[i,j], z3_b[i,j]]) + \
                xc1 + nc1*dzc1
        x3_b[i,j] = xtmp[0] 
        y3_b[i,j] = xtmp[1] 
        z3_b[i,j] = xtmp[2] 

  # base pair struts
  for j in range (nbase):
     i = 1 + (j-1)*(nh[0]-1)
     ca = scipy.cos(-thc_a(i))
     sa = scipy.sin(-thc_a(i))   # rotate around axis
     rmat0 = ca*I3 + (1-ca)*numpy.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]]) + \
            sa*numpy.matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
     for k in range(xb[0].size):
        for m in range(xb[1].size):
           xtmp = rmat0*numpy.array([xb[k,m], yb[k,m], zb[k,m]])
           xbtmp[k,m] = xtmp[0,:] 
           ybtmp[k,m] = xtmp[1,:]
           zbtmp[k,m] = xtmp[2,:] + xc_a[2,i]  #translate to proper base height

     for k in range(xbtmp[0].size):
         for m in range(xbtmp[1].size):
             xtmp = rmat2*rmat1*numpy.array([xbtmp[k,m], ybtmp[k,m], zbtmp[k,m]]) + \
                    xc1 + nc1*dzc1;
             xbtmp[k,m] = xtmp[0]
             ybtmp[k,m] = xtmp[1]
             zbtmp[k,m] = xtmp[2]

    #~ if render[2]
       #~ surface(xbtmp,ybtmp,zbtmp,0*zbtmp+shades[2],'CDataMapping','direct'); hold on;

  x1a     = rmat2*rmat1*x1a + xc1 + nc1*dzc1 # displace along axis a distance dzc1
  n1a     = rmat2*rmat1*n1a
  x1b     = rmat2*rmat1*x1b + xc1 + nc1*dzc1
  n1b     = rmat2*rmat1*n1b
  x2a     = rmat2*rmat1*x2a + xc1 + nc1*dzc1
  n2a     = rmat2*rmat1*n2a
  x2b     = rmat2*rmat1*x2b + xc1 + nc1*dzc1
  n2b     = rmat2*rmat1*n2b
  xc2     =       rmat2*xc2 + xc1 + nc1*dzc1
  nc2     =       nc2
  thc2    = thc2 + thc1


  #~ [xs, ys, zs] = sphere(ndye)


  #~ if render[0] 
     #~ if nhtot[0]>1
     #~ surface(x_a,y_a,z_a, 0*z_a   +shades[0],'CDataMapping','direct') #hold on
     #~ end
     #~ hold on
     #~ if acap[0]
       #~ surface(x5_a,y5_a,z5_a,0*z5_a+shades[0],'CDataMapping','direct')
     #~ end 
     #~ if acap[1]
        #~ surface(x3_a,y3_a,z3_a,0*z3_a+shades[0],'CDataMapping','direct')
     #~ end
     #~ if adye[0]
       #~ xtmp = x1a + n1a*(extend5prime + ddye) # need to add cap in addition to dye to get extension drawn out to dye location
       #~ surface(xtmp[0]+rdye*xs,xtmp[1]+rdye*ys,xtmp[2]+rdye*zs,0*zs+shades[0],'CDataMapping','direct')       
     #~ end
     #~ if adye[1]
       #~ xtmp = x2a + n2a*(extend3prime + ddye)
       #~ surface(xtmp[0]+rdye*xs,xtmp[1]+rdye*ys,xtmp[2]+rdye*zs,0*zs+shades[0],'CDataMapping','direct')              
     #~ end
     #~ if helpers
        #~ plot3([0 5*n1a[0]] + x1a[0],[0 5*n1a[1]] + x1a[1],[0 5*n1a[2]] + x1a[2],'b-','linewidth',2)
        #~ plot3([0 5*n2a[0]] + x2a[0],[0 5*n2a[1]] + x2a[1],[0 5*n2a[2]] + x2a[2],'b-','linewidth',2)
        #~ surf(2*xs+x1a[0],2*ys+x1a[1],2*zs+x1a[2])
        #~ surf(2*xs+x2a[0],2*ys+x2a[1],2*zs+x2a[2])
     #~ end
  #~ end

  #~ if render[1]
      #~ if nhtot[0] > 1
          #~ surface(x_b,y_b,z_b, 0*z_b   +shades[1],'CDataMapping','direct') #hold on
      #~ end
      #~ hold on
     #~ if bcap[0]
        #~ surface(x5_b,y5_b,z5_b,0*z5_b+shades[1],'CDataMapping','direct')
     #~ end
     #~ if bcap[1]
        #~ surface(x3_b,y3_b,z3_b,0*z3_b+shades[1],'CDataMapping','direct')
     #~ end
     #~ if bdye[0]
       #~ xtmp = x2b + n2b*(extend5prime + ddye)
       #~ surface(xtmp[0]+rdye*xs,xtmp[1]+rdye*ys,xtmp[2]+rdye*zs,0*zs+shades[1],'CDataMapping','direct')       
     #~ end
     #~ if bdye[1]
       #~ xtmp = x1b + n1b*(extend3prime + ddye)
       #~ surface(xtmp[0]+rdye*xs,xtmp[1]+rdye*ys,xtmp[2]+rdye*zs,0*zs+shades[1],'CDataMapping','direct')              
     #~ end
     #~ if helpers
        #~ plot3([0 5*n1b[0]] + x1b[0],[0 5*n1b[1]] + x1b[1],[0 5*n1b[2]] + x1b[2],'r-','linewidth',2)
        #~ plot3([0 5*n2b[0]] + x2b[0],[0 5*n2b[1]] + x2b[1],[0 5*n2b[2]] + x2b[2],'r-','linewidth',2)
        #~ surf(2*xs+x1b[0],2*ys+x1b[1],2*zs+x1b[2])
        #~ surf(2*xs+x2b[0],2*ys+x2b[1],2*zs+x2b[2])
    #~ end
  #~ end

  #~ #surf(2*xs+xc2[0],2*ys+xc2[1],2*zs+xc2[2])



  #~ #colorbar

  #~ xlabel('x')
  #~ ylabel('y')
  #~ zlabel('z')
  return xc2,nc2,thc2,x1a,n1a,x2a,n2a,x1b,n1b,x2b,n2b
