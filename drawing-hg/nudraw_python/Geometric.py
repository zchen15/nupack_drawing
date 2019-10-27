# define loop geometries in terms of center point and direction of all
# stems on
# each loop
#
# convention: normal points out from current loop at each base pair
#

import numpy
import math
from Classes import *
import Newton

def calc(rdh, dzb, loop, base, stacknickedhelices, dsb, show3d=False, debug=False, circle=False):
  if debug:
    print "unpairedcircle=%s"%circle

  nside = 2
  sidelength=numpy.array([2*rdh, 2*rdh])
  ds  = dzb
  nds = 2
  # get helix parameters
  r, sideangle, ang_ds, f = Newton.loop(nside, sidelength, ds, nds)
  if debug:
    print "-------------- Geometric.calc:"
    print "sideangle=",sideangle
    print "ang_ds=",ang_ds
    print "f=",f


  stacked=stackedClass(r, sideangle, ang_ds)

  nicklength = 3*rdh
  safe_count = 0
  stem_safe=numpy.zeros([0,8])
  for iloop in range(len(loop)):
    loop[iloop].nickedhelix = False # put in zeros for all loops and overwrite for nicked helices below
    nside       = loop[iloop].nside # number of stems in loop

    is_stacked      = (nside == 2) and (loop[iloop].sidenbase[0] == 0) \
                      and (loop[iloop].sidenbase[1] == 0)
    is_smallhairpin = (nside == 1) and (loop[iloop].sidenbase[0] == 1)   # hairpin with 3 bases in loop
    is_hairpin = (nside == 1)    # hairpin
    is_blunthelix   = (nside == 2) and( loop[iloop].sidenbase[0] == -1) and \
                      (loop[iloop].sidenbase[1] == -1)  # exterior loop of size zero

    #~ print "loop[%d].nside=%d, sidenbase=%s"%(iloop, nside, loop[iloop].sidenbase)
    # find out if it is a nicked helix
    is_nickedhelix = False
    if nside == 3:
      # a nicked helix can be recognized by the nick between two of the
      # side ends which are 1 base apart
      total_size = 0
      has_nick = False

      for iside in range(0,nside):
        cur_nick = loop[iloop].nick[:,iside]
        cur_base = loop[iloop].sidebase[iside]
        cur_size = loop[iloop].sidenbase[iside]
        if cur_size >= 0:
          total_size += cur_size+1

        prev_side = (iside - 1) % nside
        last_nick = loop[iloop].nick[:,prev_side]
        last_base = loop[iloop].sidebase[prev_side]
        if cur_nick[0] == -5 and cur_nick[1] == 0 and \
           last_nick[0] == 0 and last_nick[1] == -3 and \
           abs(cur_base - last_base) == 1:
          has_nick=True

      if has_nick and total_size == 0:
        is_nickedhelix = True


    if debug:
      l=loop[iloop]
      print "L%d:"%iloop
      print "  nside: %s"%nside
      print "  is_stacked: %s"%is_stacked
      print "  is_smallhairpin: %s"%is_smallhairpin
      print "  is_hairpin: %s"%is_hairpin
      print "  is_blunthelix: %s"%is_blunthelix
      print "  is_nickedhelix: %s"%is_nickedhelix
      print "  sidenbase=%s\n  sideunit=%s"%(l.sidenbase, l.sideunit)
      print "  sidebase=%s\n  coilnum=%s"%(l.sidebase, l.coilnum)
      print "  toloop=%s\n  tohelix=%s"%(l.toloop, l.tohelix)
      print "  center=%s\n  radius=%s"%(l.center, l.radius)
      print "  strand=%s"%(l.strand)
      print "  loop[iloop].nick=%s"%(str(l.nick))


    sidenbase=numpy.zeros(nside)
    for iside in range(nside):
       sidenbase[iside] = loop[iloop].sidenbase[iside]

    if is_stacked:
        if debug: print "is_stacked:"
        r           = stacked.radius
        ang_ds      = stacked.ang_ds
        sideangle   = stacked.sideangle.copy()
        sidelength  = numpy.array([2*rdh, 2*rdh])
    elif is_smallhairpin:
        if debug: print "is_smallhairpin:"
        r = rdh
        ang_ds = math.pi/2 # for hairpin of length three array bases on half circle centered at center of paired base
        sideangle[0] = math.pi # angle between ends of helix
        sidelength[0] = 2*rdh
    elif is_blunthelix: # external loop of size zero
        if debug: print "is_blunthelix:"
        r = rdh
        ang_ds = 0    # no angle between bases along strand
        sideangle = [0,0]
        sideangle[0] = math.pi   # angle between strand ends
        sideangle[1] = math.pi
        sidelength   = numpy.array([2*rdh, 2*rdh])
    elif stacknickedhelices and is_nickedhelix: # Broken!
        r           = stacked.radius
        ang_ds      = stacked.ang_ds
        if loop[iloop].sidebase[0]==0:
          if loop[iloop].sidenbase[0] == -1:
            sideangle   = numpy.array([(math.pi-stacked.sideangle[1]), \
                          stacked.sideangle[0],  stacked.sideangle[1] ])
            sidelength  = numpy.array([2*rdh, 2*rdh, 2*rdh])
          elif loop[iloop].sidenbase[1] == -1:
            sideangle   = numpy.array([stacked.sideangle[0], \
                         (math.pi-stacked.sideangle[1]), stacked.sideangle[1] ])
            sidelength  = numpy.array([2*rdh, dzb, 2*rdh])
        else:
          if loop[iloop].sidenbase[0] == -1:
            sideangle   = numpy.array([(math.pi-stacked.sideangle[1]), \
                          stacked.sideangle[0],  stacked.sideangle[1] ])
            sidelength  = numpy.array([dzb, 2*rdh, 2*rdh])
          elif loop[iloop].sidenbase[1] == -1:
            sideangle   = numpy.array([stacked.sideangle[0], \
                         (math.pi-stacked.sideangle[1]), stacked.sideangle[1] ])
            sidelength  = numpy.array([2*rdh, dzb, 2*rdh])
        if debug:
          print "stacknickedhelices and is_nicked_helix:"
          print "  sidelength: %s"%sidelength
          print "  sideangle=",sideangle
          print "  ang_ds=",ang_ds
          print "r=",r


        loop[iloop].nickedhelix = True
    else:
        nds         = 0 # number of gaps
        if nside>2:
          sidelength=numpy.resize(sidelength,nside)

        for iside in range (nside):
            nds                 = nds + loop[iloop].sidenbase[iside] + 1
            if loop[iloop].nick[1, iside] == -3:    # space for nick
                loop[iloop].geo[iside].length = nicklength
            else:                            # space for stem
                loop[iloop].geo[iside].length = 2*rdh
            sidelength[iside]   = loop[iloop].geo[iside].length
        ds      = dsb #use arc length between bases in single-stranded regions

        r, sideangle, ang_ds, f = Newton.loop(nside, sidelength, ds, nds)
        if debug:
          print "After Newton.loop:"
          print "  sidelength: %s"%sidelength
          print "  sideangle=",sideangle
          print "  ang_ds=",ang_ds
          print "  r=",r
          print "show3d=%s, is_hairpin=%s"%(show3d, is_hairpin)

    if iloop == 0: # taking starting point from global origin
        isideref = 0
        loop[iloop].geo[isideref].xc = numpy.array([0, 0, 0]).T # starting position for secondary structure drawing
        loop[iloop].geo[isideref].nc = numpy.array([0, 1, 0]).T # starting direction for secondary structure drawing
        jsidestart = 1
        jsidestop = nside-1
    else: # take starting point from previous loop
        isideref = nside-1
        prevloop = loop[iloop].toloop[isideref]
        prevside = loop[iloop].toside[isideref]
        loop[iloop].geo[isideref].xc =  loop[prevloop].geo[prevside].xc.copy()
        loop[iloop].geo[isideref].nc = -loop[prevloop].geo[prevside].nc.copy()

        jsidestart = 0
        jsidestop = nside-2

    # assign stem position and direction for jside stem (by
    # convention this follows the jside single-stranded region)
    #
    xc = loop[iloop].geo[isideref].xc.copy()
    nc = loop[iloop].geo[isideref].nc.copy()
    if show3d and is_hairpin:
      loopcenter = xc - \
                 nc*math.sqrt(math.pow(r,2) - \
                 math.pow(sidelength[isideref]/3,2))
    else:
      loopcenter = xc - \
                 nc*math.sqrt(math.pow(r,2) - \
                 math.pow(sidelength[isideref]/2,2))
    loop[iloop].center = loopcenter
    if show3d and is_hairpin:
      fac=rdh*2.5/float(loop[iloop].sidenbase[0]+1)
      #~ print "fac=%s"%fac
      if fac<r:
        tempR = math.sqrt(math.pow(r,2) - \
                     math.pow(rdh*2.5/float(loop[iloop].sidenbase[0]+1),2))
        #~ print "tempR=%f, r=%f, nside=%d"%(tempR, r, nside)
      else:
        tempR=0.0
      tempR=max(tempR,rdh*0.5)
      #~ print "tempR=%f, r=%f"%(tempR, r)
      loop[iloop].radius=tempR
    else:
      loop[iloop].radius = r
    jprev = isideref
    stem_size=jsidestop-jsidestart+1

    for jside in range (jsidestart,jsidestop+1):
      jgap = loop[iloop].sidenbase[jside] + 1
      jang    = ang_ds*jgap + .5*(sideangle[jside] + sideangle[jprev])

      ca      = math.cos(jang)   # rotation matrix is CW around z axis
      sa      = math.sin(jang)
      rmat    = ca*numpy.eye(3) + (1-ca)*numpy.array([[0,0,0],[ 0,0,0],[0,0,1]]) + \
                            sa*numpy.array([[0,1,0],[-1,0,0],[0,0,0]])
      nc      = numpy.dot(rmat,nc)
      xc      = loopcenter + nc*math.sqrt(math.pow(r,2) - \
                math.pow(sidelength[jside]/2,2))

      loop[iloop].geo[jside].nc = nc.copy()
      loop[iloop].geo[jside].xc = xc.copy()

      jprev   = jside
      stem_safe=numpy.vstack([stem_safe, numpy.hstack([numpy.array([iloop, jside]),xc.T,nc.T])])
      safe_count = safe_count + 1

    #
    # assign base positions for jside single-stranded region
    #
    jprev = nside-1
    #~ pdb.set_trace()
    for jside in range(nside):
        jbasestart = loop[iloop].sidebase[jside]
        jbasestop = loop[iloop].sidebase[jside] + loop[iloop].sidenbase[jside]
        if loop[iloop].nick[1,jside] == -3 or (nside==1 and loop[iloop].nick[1,jside+1]==-3):
            jbasestop = jbasestop + 1 # handle base at nick that would not otherwise get assigned an x value

        nc = loop[iloop].geo[jprev].nc
        ang = .5*sideangle[jprev]
        if len(loop)>1 or circle:
          for jbase in range(jbasestart,jbasestop+1): #loop over paired base and following single-stranded region
              ca      = math.cos(ang)   # rotate CW by ang
              sa      = math.sin(ang)   # rotation matrix is CW around z axis
              rmat    = ca*numpy.eye(3) + (1-ca)*numpy.array([[0,0,0],[ 0,0,0],[0,0,1]]) + \
                                    sa*numpy.array([[0,1,0],[-1,0,0],[0,0,0]])
              ncstep = numpy.dot(rmat,nc)
              base[jbase].x = loop[iloop].center + ncstep*loop[iloop].radius
              ang = ang + ang_ds
        else:
          l=2*math.pi*loop[iloop].radius # Loop length
          inc=l/(jbasestop+1 - jbasestart)
          li=0
          for jbase in range(jbasestart,jbasestop+1):
            base[jbase].x = loop[iloop].center + li*inc*numpy.array([1,0,0])
            li+=1

        jprev = jside
