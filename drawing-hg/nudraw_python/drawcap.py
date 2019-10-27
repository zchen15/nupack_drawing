#
# nudraw
# 
# A program for visualizing DNA Structures and Devices
# Copyright 2004
# Niles A. Pierce
# California Institute of Technology
# Version 1.0 coded 2004 (April 9,10,13,14)

# Converted to Python/Numpy/SciPy/PyLab Jul 31, 2008
# Conrad Steenberg <conrad.steenberg@caltech.edu>


import numpy
import scipy
import pylab

def  drawcap(i, xc_a, rmat_a, endtype, extend):
  # i                 index along chain patches
  # xc_a(3)           axis of helix (z is same for both chains before rotation)
  # rcap              radius of end cap corner
  # rbc(2)            radii of base strut
  # th                angle array around the chain
  # ncap              number of points going over radius of cap
  # thcap              angle of conical cap at 3' end of chain
  # rhc               radius of helix chain
  # rmat_a            rotation matrix for cap to match angle of chain tube
  # endtype           end of chain
  # extend            amount to extend chain at end of strand before cap
                      # differs for 3' and 5' ends in helix and in coil

  global  nh, nb, rhc, rbc, rhelix, rdh, rcap, ncap, thcap, dzb, bturn, brise, \
          dthb, dsb, dz, dth, cmap , dye, ddye, ndye, majorgroove, minorgroove, \
          dthgroove, strutangle, strutlength, strutrise, proptwist, inclination, \
          extend5prime, extend3prime, extendcoil, rcaphead

  start=0.
  stop=2*scipy.pi
  steps=(stop-start)/dth
  th = pylab.linspace(start, stop, steps)


  if endtype ==  5:
  #if endtype ==  3 # quick fix for polarity issue
    xtmp[0,:] = rhc*scipy.cos(th)
    xtmp[1,:] = rhc*scipy.sin(th)
    xtmp[2,:] = 0
    xtmp = rmat_a*xtmp
    xcap_a[0,:] = xtmp[0,:] + xc_a[0,i]
    ycap_a[0,:] = xtmp[1,:] + xc_a[1,i]
    zcap_a[0,:] = xtmp[2,:] + xc_a[2,i]

    for j in range(2,ncap-1): # ncap -1 not included
      dcap = pi/2 * (j-2)/(ncap-3)
      xtmp[0,:] = ((rhc - rcap) + rcap*scipy.cos(dcap))*scipy.cos(th)
      xtmp[1,:] = ((rhc - rcap) + rcap*scipy.cos(dcap))*scipy.sin(th)
      xtmp[2,:] = rcap*scipy.sin(dcap) + extend # offset from last base by rbc(1) 
      xtmp = rmat_a*xtmp
      xcap_a[j,:] = xtmp[0,:] + xc_a[0,i]
      ycap_a[j,:] = xtmp[1,:] + xc_a[1,i]
      zcap_a[j,:] = xtmp[2,:] + xc_a[2,i]

    xtmp[0,:] = 0*scipy.cos(th)
    xtmp[1,:] = 0*scipy.sin(th)
    xtmp[2,:] = rcap + extend #*scipy.sin(phi)
    xtmp = rmat_a*xtmp
    xcap_a[ncap-1,:] = xtmp[0,:] + xc_a[0,i]
    ycap_a[ncap-1,:] = xtmp[1,:] + xc_a[1,i]
    zcap_a[ncap-1,:] = xtmp[2,:] + xc_a[2,i]   
  # elseif endtype == 3
  #       xtmp[0,:] = rhc*scipy.cos(th)
  #       xtmp[1,:] = rhc*scipy.sin(th)
  #       xtmp[2,:] = 0  
  #       xtmp = rmat_a*xtmp
  #       xcap_a[0,:] = xtmp[0,:] + xc_a[0,i]
  #       ycap_a[0,:] = xtmp[1,:] + xc_a[1,i]
  #       zcap_a[0,:] = xtmp[2,:] + xc_a[2,i]
  #    for j=2:ncap-1
  #       dcap = (pi/2 - thcap) * (j-2)/(ncap-3)
  #       xtmp[0,:] = ((rhc - rcap) + rcap*scipy.cos(dcap))*scipy.cos(th)
  #       xtmp[1,:] = ((rhc - rcap) + rcap*scipy.cos(dcap))*scipy.sin(th)
  #       xtmp[2,:] = rcap*scipy.sin(dcap) + extend #offset from last base by rbc(1) 
  #       xtmp = rmat_a*xtmp
  #       xcap_a[j,:] = xtmp[0,:] + xc_a[0,i]
  #       ycap_a[j,:] = xtmp[1,:] + xc_a[1,i]
  #       zcap_a[j,:] = xtmp[2,:] + xc_a[2,i]
  #    end
  #    for j = ncap:2*ncap
  #       dcap = thcap * ((2*ncap)-j)/ncap  # goes from thcap to zero
  #       xtmp[0,:] = rcap*scipy.sin(dcap)*scipy.cos(th)
  #       xtmp[1,:] = rcap*scipy.sin(dcap)*scipy.sin(th)
  #       xtmp[2,:] = (rhc - rcap)*scipy.tan(thcap) + rcap*scipy.cos(dcap) + extend 
  #       xtmp = rmat_a*xtmp
  #       xcap_a[j,:] = xtmp[0,:] + xc_a[0,i]
  #       ycap_a[j,:] = xtmp[1,:] + xc_a[1,i]
  #       zcap_a[j,:] = xtmp[2,:] + xc_a[2,i]
  #    end  

  elif endtype == 3:
  #elseif endtype == 5 # quick fix
    xtmp[0,:] = rhc*scipy.cos(th)
    xtmp[1,:] = rhc*scipy.sin(th)
    xtmp[2,:] = 0
    xtmp = rmat_a*xtmp
    xcap_a[0,:] = xtmp[0,:] + xc_a[0,i]
    ycap_a[0,:] = xtmp[1,:] + xc_a[1,i]
    zcap_a[0,:] = xtmp[2,:] + xc_a[2,i]

    xtmp[0,:] = rhc*scipy.cos(th)
    xtmp[1,:] = rhc*scipy.sin(th)
    xtmp[2,:] = extend
    xtmp = rmat_a*xtmp
    xcap_a[1,:] = xtmp[0,:] + xc_a[0,i]
    ycap_a[1,:] = xtmp[1,:] + xc_a[1,i]
    zcap_a[1,:] = xtmp[2,:] + xc_a[2,i]

    for j in range(3, ncap-1):
      dcap = (pi - thcap) * (j-3)/(ncap-4) - pi/2
      xtmp[0,:] = ((rcaphead - rcap) + rcap*scipy.cos(dcap))*scipy.cos(th)
      xtmp[1,:] = ((rcaphead - rcap) + rcap*scipy.cos(dcap))*scipy.sin(th)
      xtmp[2,:] = rcap + rcap*scipy.sin(dcap) + extend #offset from last base by rbc(1) 
      xtmp = rmat_a*xtmp
      xcap_a[j,:] = xtmp[0,:] + xc_a[0,i]
      ycap_a[j,:] = xtmp[1,:] + xc_a[1,i]
      zcap_a[j,:] = xtmp[2,:] + xc_a[2,i]

    for j in pylab.linspace(ncap,2*ncap):
      dcap = thcap * ((2*ncap)-j)/ncap  # goes from thcap to zero
      xtmp[0,:] = rcap*scipy.sin(dcap)*scipy.cos(th)
      xtmp[1,:] = rcap*scipy.sin(dcap)*scipy.sin(th)
      xtmp[2,:] = rcap + (rcaphead - rcap)*scipy.tan(thcap) + rcap*scipy.cos(dcap) + extend 
      xtmp = rmat_a*xtmp
      xcap_a[j,:] = xtmp[0,:] + xc_a[0,i]
      ycap_a[j,:] = xtmp[1,:] + xc_a[1,i]
      zcap_a[j,:] = xtmp[2,:] + xc_a[2,i]

  elif endtype == 0:  # cap for base pair
    rx = rbc[1]
    ry = rbc[2]
    dcap = pylab.linspace(pi/2,0,ncap)
    for j in range(1, ncap):
      xtmp[0,:] = rx*scipy.cos(dcap[j])*scipy.cos(th)
      xtmp[1,:] = ry*scipy.cos(dcap[j])*scipy.sin(th)
      xtmp[2,:] = rcap*scipy.sin(dcap[j]) # *scipy.sin(phi)  # radius starts at edge of last pair 
      xtmp = rmat_a*xtmp
      xcap_a[j,:] = xtmp[0,:] + xc_a
      ycap_a[j,:] = xtmp[1,:]
      zcap_a[j,:] = xtmp[2,:]

  return xcap_a,ycap_a,zcap_a
