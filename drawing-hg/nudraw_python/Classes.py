# nudraw class definitions
import numpy
#-------------------------------------------------------------------------------
class unitClass:
  type=1
  base=1
  strand=1
  pair=-1
  ppair=-1
  next=-1
  prevloop=-1
  prevside=-1

  def __init__(self, newtype=1, newbase=1, newstrand=1, baseInfo=None, prob=-1.0):
    self.type=newtype
    self.base=newbase
    self.strand=newstrand
    self.baseInfo=baseInfo
    self.prob=prob

  def __repr__(self):
    return "[%s, %02d, %02d, %02d, %02d]"%(self.type, self.base, self.strand,
                                           self.pair, self.next)

  def __str__(self):
    return self.__repr__()

#-------------------------------------------------------------------------------
class geoClass:
  length=0.0
  xc=numpy.zeros(3)
  nc=numpy.zeros(3)
  def __init__(self):
    pass

  def __str__(self):
    return "len=%s\nxc=%s\nnc=%s\n"%(self.length, self.xc, self.nc)

#-------------------------------------------------------------------------------
class refClass:
  def __init__(self, max):
    self.max=max+60
    self.angle  = numpy.zeros(max+60)
    self.origin = numpy.zeros([3,max+60]) # keep the origin the same as helix[1].xc(:,1)
    self.axis   = numpy.zeros([3,max+60]) # rotations work for z axis but not for arbitrary axes
    self.num    = 0         # number of reference rotations

#-------------------------------------------------------------------------------
class loopClass:
  sideunit=None
  sidebase=None
  sidenbase=None
  strand=None
  prevside=None
  toside=None
  nextside=None
  toloop=None
  strand=None
  nickedhelix=0
  coilnum=None
  radius=0.0
  nside=0
  globalmax=0

  geo=[]
  nick=[]

  def __init__(self):
    self.geo=[]
    self.nick=[]
    self.center=numpy.zeros(3)

    self.nside=0

  def setSize(self,in_max):
    if in_max>1:
      max=in_max
    else:
      max=2
    self.sideunit=-1*numpy.ones(max,numpy.int32)
    self.sidebase=-1*numpy.ones(max,numpy.int32)
    self.sidenbase=-1*numpy.ones(max,numpy.int32)
    self.prevside=-1*numpy.ones(max,numpy.int32)
    self.toside=-1*numpy.ones(max,numpy.int32)
    self.nextside=-1*numpy.ones(max,numpy.int32)
    self.toloop=-1*numpy.ones(max,numpy.int32)
    self.tohelix=-1*numpy.ones(max,numpy.int32)
    self.strand=-1*numpy.ones(max,numpy.int32)
    self.coilnum=-1*numpy.ones(max,numpy.int32)

    tempArray=[]
    for i in range(max):
      tempArray.append(-1.*numpy.ones(2))
      self.geo.append(geoClass())

    self.nick=numpy.array(tempArray).T
    self.nickedhelix=False

  def setRef(self,in_max=0):
    if in_max>1:
      max=in_max
    else:
      max=2

    self.ref=refClass(max)

#-------------------------------------------------------------------------------
class stackedClass:
  radius=1.0
  sideangle=None
  ang_ds=0.0

  def __init__(self, radius, sideangle, ang_ds):
    self.radius=radius
    self.sideangle=sideangle.copy()
    self.ang_ds=ang_ds

#-------------------------------------------------------------------------------
class baseClass:
  x=numpy.zeros(3)
  x3=numpy.zeros(3)
  color="#000000"
  domain=None
  domainid=None
  def __init__(self, iunit, type=""):
    self.type=type
    x=numpy.zeros(3)
    x3=numpy.zeros(3)
    self.iunit=iunit

  def __str__(self):
    return "%s  at %s and %s"%(self.type, self.x, self.x3)

#-------------------------------------------------------------------------------
class helixClass:
  prevloop=0
  prevside=0
  npair=0
  thc=0.
  xc=[]
  nc=[]

  def __init__(self):
    self.xc=numpy.zeros(3).T
    self.nc=numpy.zeros(3).T
    self.cap = numpy.array([0, 0, 0, 0]).T
    self.tocoil = numpy.array([4,2],numpy.int32)
    self.strand = numpy.array(2,numpy.int32)
    self.points = None
    self.added=[False, False]
    self.bases=[]

#-------------------------------------------------------------------------------
class coilClass:
  nbase=0
  loopnum=0
  loopside=0
  nsplinepts=0

  def __init__(self, max):
    self.nbase=0
    self.loopnum=0
    self.loopside=0
    self.tohelix=-1*numpy.ones([2,2] ,numpy.int32)
    self.strand=0

    tempArray=[]
    for i in range(max):
      tempArray.append(numpy.zeros(2))

    self.nick=numpy.array(tempArray).T
    self.strand=0
    self.nsplinepts=0
    self.center=numpy.zeros(3)
    self.x=numpy.zeros([3,2])
    self.n=numpy.zeros([3,2])
    self.bctype=numpy.ones(2)
    self.th=numpy.zeros([1,1])
    self.points=None
    self.added=[False,False]
    self.start_base=-1
    self.end_index=-1
    self.bases=[]
    self.drawn_bases=False

  def setSplineSize(self,max):
    self.xsplinepts=numpy.zeros([3,max])

