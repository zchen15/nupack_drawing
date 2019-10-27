# Update layout from a given JSON file
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 20, 2008

import json2 as json
import numpy
import math
import Export

#-------------------------------------------------------------------------------
def loadJson(fname):
  f=open(fname)
  data=f.read()
  json_data=json.read(data)
  if json_data.has_key('struc'):
    json_data=json_data['struc']

  #~ print json_data.keys()
  return json_data, \
         json_data["dotparens"], \
         json_data["sequence"]

#-------------------------------------------------------------------------------
def loadJson25d(fname):
  return loadJson(fname)

#-------------------------------------------------------------------------------
def trans(x, miny_pad, maxy_pad):
  if type(x)==list:
    retval=numpy.array(x)
  else:
    retval=x.copy()
  retval[1]=-(retval[1] - maxy_pad - miny_pad)
  return retval

#-------------------------------------------------------------------------------
def updateBases(json_data, base, unit, update_prob):
  json_base=json_data["bases"]
  miny_pad, maxy_pad=json_data["properties"]["ypad"]
  #~ print len(base), len(json_base)
  iu=0
  lu=len(unit)
  for i in range(len(json_base)):
    loc=json_base[i]["loc"]
    #~ bx=json_base[i]["loc"]
    #~ by=float(maxy_pad-u.baseInfo.x[1]+miny_pad)
    newx=trans([loc[0], loc[1], 0.], miny_pad, maxy_pad)
    #~ print "%3d: base.x=%s"%(i,base[i].x)
    #~ print "     newx  =%s"%newx
    base[i].x=newx
    while unit[iu].base<i: iu+=1
    while iu<lu and unit[iu].base==i:
      if update_prob:
        unit[iu].prob=json_base[i]["probability"]
      iu+=1
  #~ print json_base

#-------------------------------------------------------------------------------
def updateBases25d(json_data, base, unit, update_prob):
  json_base=json_data["bases"]
  miny_pad, maxy_pad=json_data["properties"]["ypad"]
  #~ print len(base), len(json_base)
  iu=0
  lu=len(unit)
  for i in range(len(json_base)):
    loc=json_base[i]["loc"]
    #~ bx=json_base[i]["loc"]
    #~ by=float(maxy_pad-u.baseInfo.x[1]+miny_pad)
    newx=trans(loc, miny_pad, maxy_pad)
    #~ print "%3d: base.x=%s"%(i,base[i].x)
    #~ print "     newx  =%s"%newx
    base[i].x3=newx
    base[i].x[:2]=newx[:2]
    while unit[iu].base<i: iu+=1
    while iu<lu and unit[iu].base==i:
      if update_prob:
        unit[iu].prob=json_base[i]["probability"]
      iu+=1


#-------------------------------------------------------------------------------
def avg_dist(loc1, loc2):
  "Calculates the point halfway between loc1,2 and half the separation distance"
  return (loc2+loc1)/2, numpy.sqrt((loc2[0]-loc1[0])**2 + (loc2[1]-loc1[1])**2)

#-------------------------------------------------------------------------------
def perp_vec(loc1, loc2):
  "Calculate a normalized vector perpendicular to two given points. I.e. rotated 90 deg"
  v1=numpy.array(loc2-loc1)
  nv1 = v1/numpy.linalg.norm(v1)
  return numpy.array([-nv1[1], nv1[0]])

#-------------------------------------------------------------------------------
def perp_vec3(loc1, loc2):
  "Calculate a normalized vector perpendicular to two given points, in the xy plane."
  v1=numpy.array(loc2-loc1)
  nv1 = v1/numpy.linalg.norm(v1)
  #~ print "v1=%s, nv1=%s"%(v1, nv1)
  return numpy.array([-nv1[1], nv1[0], nv1[2]])


#-------------------------------------------------------------------------------
def update_loop_from_jloop(loopObj, jloopObj, miny_pad, maxy_pad):
  loopObj.radius=jloopObj["radius"]
  loopObj.center[:2]=trans(jloopObj["center"], miny_pad, maxy_pad) # transform coords

#-------------------------------------------------------------------------------
def update_loop_from_helix(loop, jhelixObj, base, iloop, miny_pad, maxy_pad, loops_todo):
  hlen=len(jhelixObj["bases"])/2
  #~ print "\nupdate_loop_from_helix(%d), hlen=%d"%(iloop, hlen)
  if hlen==1:
    pass
    #~ for xci in [0,1]:
      #~ print "  old(%d,%d).xc=%s"%(iloop, xci, loop[iloop].geo[xci].xc)
  else:
    for j in range(1,hlen):
      j1=jhelixObj["bases"][j]
      j2=jhelixObj["bases"][j+hlen]
      #~ print "\nbase[%d, %d] -> loop[%d]"%(j1, j2, iloop)
      loc1=base[j1-1].x # A base pos, to use for radius
      csum=numpy.array([0.,0., 0])

      for jside in [0,1]:
        jl1=j1-1+jside
        jl2=j2+1-jside
        l1=base[jl1].x
        l2=base[jl2].x
        c,r=avg_dist(l1,l2)
        if jside==0:
          nc=perp_vec(l2,l1)
        else:
          nc=perp_vec(l1,l2)
        #~ print "loc(%d, %d)= %s -> %s"%(jl1, jl2, l1[:2],l2[:2])
        csum+=c
        xci=(jside+1)%2
        #~ print "  old(%d,%d).xc=%s nc=%s"%(iloop, xci, loop[iloop].geo[xci].xc,
                                          #~ loop[iloop].geo[xci].nc)
        loop[iloop].geo[xci].xc=c
        loop[iloop].geo[xci].nc=nc
        #~ print "  new(%d,%d).xc=%s nc=%s"%(iloop, xci, loop[iloop].geo[xci].xc,
                                          #~ loop[iloop].geo[xci].nc)

      #~ print "old c=%s, r=%f"%(loop[iloop].center[:2],loop[iloop].radius)
      loop[iloop].center=csum/2
      loc2=l2
      loop[iloop].radius=numpy.sqrt((loc2[0]-loc1[0])**2 + (loc2[1]-loc1[1])**2)/2
      #~ print "new c=%s, r=%f"%(loop[iloop].center[:2],loop[iloop].radius)

      iloop+=1

  return iloop

#-------------------------------------------------------------------------------
def updateLoops(json_data, loop, base):
  "Update loop radii and centers usin info from json file"
  # Get some info from json data structure
  json_base=json_data["bases"]
  json_loop=json_data["loops"]
  json_helix=json_data["helices"]
  miny_pad, maxy_pad=json_data["properties"]["ypad"]

  # Indices of helices left to do in json structure
  helices_todo=[]
  loops_todo=[0]
  loops_done={}

  # Iterate through loops in json file
  iloop=0 # counter for already laid out loops
  while len(loops_done) < len(json_loop) and loops_todo:
    #~ print "\nloops_todo=%s"%(loops_todo)
    jloop=loops_todo[-1]
    loops_todo=loops_todo[:-1]
    #~ print "jloop=%d"%jloop
    #~ print "updating jloop[%d] -> loop[%d]"%(jloop, iloop)
    update_loop_from_jloop(loop[iloop], json_loop[jloop], miny_pad, maxy_pad)

    # Update values of loop[].geo[side].xc
    iseg=0
    segments=json_loop[jloop]["segments"]
    #~ print "\nsegments=%s, sides=%s"%(segments,loop[iloop].nside)
    if segments:
      jsidemax=len(loop[iloop].geo)
      jside=0
      while jside<jsidemax:
        #~ print "segments[%d]=%s"%(iseg-1, segments[iseg-1])
        #~ print "segments[%d]=%s"%(iseg, segments[iseg])

        l1=base[segments[iseg-1][1]].x
        l2=base[segments[iseg][0]].x
        c,r=avg_dist(l1,l2)
        #~ print "old: xc[%d, %d] = %s"%(iseg, jside, loop[iloop].geo[jside].xc)
        #~ print "new: xc[%d, %d] = %s"%(iseg, jside, c)
        loop[iloop].geo[jside].xc=c
        iseg+=1
        while iseg<len(segments) and segments[iseg][0]==segments[iseg-1][1]: # look for next segment break
          iseg+=1
        if iseg==len(segments): break
        jside+=1

    iloop+=1
    helices_todo+=json_loop[jloop]["childHelices"][::-1]
    #~ print "helices_todo=%s"%(helices_todo)
    #~ print "helices_todo=%s"%(helices_todo)
    while len(helices_todo)>0:
      jhelix=helices_todo[-1]
      helices_todo=helices_todo[:-1]

      #~ print "updating jhelix[%d] -> loop[%d]"%(jhelix, iloop)
      # Iterate through loops in helix
      iloop=update_loop_from_helix(loop, json_helix[jhelix], base, iloop,
            miny_pad, maxy_pad, loops_todo)
      childLoop=json_helix[jhelix]["childLoop"]
      if childLoop>=0:
        loops_todo.append(childLoop)
        break #Do this loop first
      #~ print "helix update done, iloop=%d"%iloop


#-------------------------------------------------------------------------------
def updateHelix25d(json_data, loop, base, helix, gg):
  """Update helix .xc and .nc paramters using info from json file. Call after bases are updated"""
  #~ print "****updateHelix25d"
  # Get some info from json data structure
  json_base=json_data["bases"]
  json_loop=json_data["loops"]
  json_helix=json_data["helices"]
  miny_pad, maxy_pad=json_data["properties"]["ypad"]
  #~ print "miny_pad=%f, maxy_pad=%f"%(miny_pad, maxy_pad)
  for hi in range(len(json_helix)):
    bases=json_helix[hi]["bases"]
    pl=json_helix[hi]["parentLoop"]
    #~ print "parentLoop=%d"%pl
    nb=len(bases)

    l1=base[bases[0]].x
    l2=base[bases[nb/2]].x

    #~ print "h%d: b1=%d, b2=%d\nl1=%s, l2=%s angle=%f"%(hi, bases[0],bases[nb/2], l1, l2, json_helix[hi]["angle"])
    hAng=math.pi*2-json_helix[hi]["angle"] # This is more accurate
    # Calculate nc vector and xc position
    xc,r=avg_dist(l1,l2)
    #~ print "\nold: helix[%d].xc=%s"%(hi,helix[hi].xc)
    #~ print "new: helix[%d].xc=%s"%(hi,xc)
    helix[hi].xc=xc

    # Calculate vector nc
    nc=numpy.zeros(3)
    nc[0]=math.cos(hAng)
    nc[1]=math.sin(hAng)
    nc[2]=0.
    #~ print "old: helix[%d].nc=%s"%(hi,helix[hi].nc)
    #~ print "new1: helix[%d].nc=%s"%(hi,nc)

    helix[hi].nc=nc
    # Calculate thc
    helix[hi].thc    = 180./math.pi*(math.pi - gg.strutangle)
    adjust=90.

    # If parent was flipped, we need to flip helix too
    if pl>-1 and len(json_loop[pl]["seg_orient"])>0 and \
        json_loop[pl]["seg_orient"][0]==1:
      adjust=-90.

    helix[hi].thc = helix[hi].thc  + adjust + \
    180./math.pi*numpy.arctan2(nc[1],nc[0])

  return

#-------------------------------------------------------------------------------
class coilpts:
  def __init__(self, bases, nsplinepts, xsplinepts, center, radius):
    self.bases=bases
    self.nsplinepts=nsplinepts
    self.xsplinepts=xsplinepts
    self.center=center
    self.radius=radius

  #~ def __str__(self):
    #~ return "r=%s, c=%s, bases=%s"%(self.radius, self.center, self.bases)

#-------------------------------------------------------------------------------
def newCoil25d(newcoil, jloopObj, miny_pad, maxy_pad, base, jbase):
  #~ print "newCoil25d"
  segs=[]

  # Extract the full segment from the json loop object
  if len(jloopObj["segments"])>0:
    sb=jloopObj["segments"][0][0]
    eb_old=jloopObj["segments"][0][1]
    ss=jbase[sb]["strand"]
    seg=[sb]

    for i in range(len(jloopObj["segments"])):
      sbt=jloopObj["segments"][i][0]
      eb=jloopObj["segments"][i][1]
      es=jbase[eb]["strand"]
      if (i>0 and sbt!=eb_old) or es!=ss:
        segs.append(seg)
        seg=[sbt]
        sb=eb
        ss=jbase[sb]["strand"]
      seg.append(eb)
      eb_old=eb

    if seg:
      segs.append(seg)

  # Create temporary coil object and put it in an indexed structure
  for seg in segs:
    fseg=filter(lambda x: jbase[x]["structype"]=='.',seg)
    nbases=len(fseg)
    if nbases==0: continue

    splinepts=numpy.zeros([3, nbases])

    for i in range(len(fseg)):
      bi=fseg[i]
      #~ print bi,base[bi].x3
      splinepts[:,i]=base[bi].x3 # already transformed
    fcoil=coilpts(fseg, nbases, splinepts,
      trans(jloopObj["center"], miny_pad, maxy_pad),
      jloopObj["radius"])
    #~ print "bases2=",fseg
    #~ print "xsplinepts2=\n",splinepts
    newcoil[fseg[0]]=fcoil
    #~ print "newcoil[%s]=%s"%(fseg[0],fcoil)


#-------------------------------------------------------------------------------
def updateCoil25d(loop, oldcoil, newcoil):
  #~ print "old:",oldcoil.nsplinepts,oldcoil.xsplinepts.shape
  #~ print "new:",newcoil.nsplinepts,newcoil.xsplinepts.shape
  if newcoil.nsplinepts!=oldcoil.nsplinepts or \
    oldcoil.xsplinepts.shape!=newcoil.xsplinepts.shape:
    raise RuntimeError("Coil spline point mismatch")
  #~ print "loop[%d].bases=%s\n    - newbases=%s"%(oldcoil.loopnum, oldcoil.bases, newcoil.bases)
  #~ print "oldcenter=%s\nnewcenter=%s, coilobj=%s\n"%(loop[oldcoil.loopnum].center,newcoil.center,newcoil)
  oldcoil.nsplinepts=newcoil.nsplinepts
  oldcoil.xsplinepts=newcoil.xsplinepts
  oldcoil.bases=newcoil.bases
  oldcoil.center=newcoil.center
  loop[oldcoil.loopnum].radius=newcoil.radius
  #~ loop[oldcoil.loopnum].center=newcoil.center


#-------------------------------------------------------------------------------
def updateCoils25d(json_data, loop, base, coil):
  "Update coil radii centers and xsplinepts info from json file"
  # Get some info from json data structure
  json_base=json_data["bases"]
  json_loop=json_data["loops"]
  json_helix=json_data["helices"]
  miny_pad, maxy_pad=json_data["properties"]["ypad"]

  newcoil={}
  oldcoil={}
  #~ print ">>>>>>>>  Creating new coils"
  for i in range(len(json_loop)):
    newCoil25d(newcoil, json_loop[i], miny_pad, maxy_pad, base, json_base)
  #~ print "newcoil_bases=\n",newcoil.keys()

  for i in xrange(len(coil)):
    c=coil[i]
    #~ print "oldcoil[%d].bases=\n%s"%(i,c.bases)
    #~ print "oldcoil[%d].center=\n%s"%(i,c.center)
    if c.nsplinepts>0:
      oldcoil[c.bases[0]]=c
    #~ else:
      #~ print ">>>>>>>>  oldcoil[%d] NOT ADDED"%i

  for bi in oldcoil.keys():
    updateCoil25d(loop, oldcoil[bi], newcoil[bi])
    l=loop[oldcoil[bi].loopnum]
    #~ for li in filter(lambda x: x>-1, l.coilnum):
      #~ print "coil[%d].center=%s, bases=%s"%(li, coil[li].center,coil[li].bases)

#-------------------------------------------------------------------------------
def updateCoils25dSweep(json_data, loop, base, coil):
  "Update coil radii centers info from json file"
  # Get some info from json data structure
  json_base=json_data["bases"]
  json_loop=json_data["loops"]
  json_helix=json_data["helices"]
  miny_pad, maxy_pad=json_data["properties"]["ypad"]

  for i in range(len(json_loop)):
    jloop=json_loop[i]
    center=trans(jloop["center"], miny_pad, maxy_pad)
    radius=jloop["radius"]

    #~ print "jloop=%d:"%i
    for ci in jloop["coils"]:
      #~ print ci, coil[ci]
      coil[ci].center=numpy.array(center)
      coil[ci].radius=numpy.array(radius)

    li=jloop["loop"]
    #~ print "iloop=%d:"%li
    if loop[li].coilnum[0]<0:
      #~ print "Setting loop center..."
      loop[li].center=numpy.array(center)
      loop[li].radius=numpy.array(radius)


