#!/usr/bin/python
#
# Coil scene renderer for OpenGL and producer for Tachyon ray tracer
#
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# 12 Aug 2008

# This is statement is required by the build system to query build info
if __name__ == '__build__':
  raise Exception

import math
import numpy
try:
  from OpenGL.GL import *
  from OpenGL.GLU import *
  from OpenGL.GLUT import *
  import OpenGLdisplay
except:
  OpenGL=None
  class OpenGLdisplay:
    class Base:
      def __init__(self):
        pass

  GL_POINTS=1
  GL_LINES=1

try:
  if OpenGL:
    from gle import *
except:
  print "importing OpenGL.GLE"
  from OpenGL.GLE import *

import json2 as json
mesh=None

# Constants
white=     [1.0,      1.0,      1.0      ]
orange=    [1.,       0.2857,    0.      ]
green=     [109./255,  216./255, 45./255 ]
red=       [164./255,  0.0,      0       ]     # red 3
dark_green=[0./255,    100./255, 0./255  ]     # dark green 4
yellow=    [1.,        0.9375,   0       ]     # yellow 5
turquoise= [0.,        0.9286,   1       ]     # turquoise 6
dark_blue= [0.,        0.,       0.8571  ]     # dark blue 7
brown=     [142./255,   73./255,  5./255 ]     # brown 8

# Tachyon scene information
tachyon_scene_header="""
BEGIN_SCENE
  RESOLUTION 512 512

TEXDEF txt001 AMBIENT 0.2  DIFFUSE 0.8  SPECULAR 0  OPACITY 1
 COLOR 1 0.75 0.33
 TEXFUNC 0

TEXDEF txt002 AMBIENT 0.6 DIFFUSE 0.9 SPECULAR 0.0 OPACITY 1.0
  COLOR 1.0 1.0 1.0 TEXFUNC 0

TEXDEF txt003 AMBIENT 0.3 DIFFUSE 0.9 SPECULAR 0.5 OPACITY 0.5
  COLOR 0.24313725490196078 0.44313725490196076 0.5490196078431373 TEXFUNC 0


Background 1.000 1.000 1.000

"""

tachyon_scene_camera="""

CAMERA
  ZOOM 0.9
  ASPECTRATIO 1.0
  ANTIALIASING 1
  RAYDEPTH 12
  CENTER %f %f %f
  VIEWDIR 0 0 1
  UPDIR 0 1 0

END_CAMERA

"""

tachyon_scene_light="""


LIGHT CENTER %f %f %f RAD 10 COLOR 1.0 1.0 1.0

"""


tachyon_scene_footer="""
END_SCENE
"""

# POVRay scene information
povray_scene_header="""
#include "colors.inc"

"""

povray_scene_light="""
light_source { <%f, %f, %f> color White }

"""

povray_scene_camera="""
camera {
  location <%f, %f, %f>
  look_at <0, 0, 1>
}

"""

povray_scene_footer=""

def findNext(x1, x2):
  """Finds the next point in parametric space after points x1 -> x2
  """
  # Set t1=1, t2=-1, t0=3, t3=-3
  t3=-3
  c=(x2+x1)/2
  m=x1 -c
  xn=m*t3+c
  return xn

class Object:
  tachyon_textures={}
  json_scene={ "camera"    : {},
                "lights"    : [],
                "objects"   : {
                  "cones"   : [],
                  "spheres" : [],
                  "polycylinders": []
                },
                "textures"  : []
              }
  texture_map={}
  texcolors=[]
  tachcolors=[]

  """Generic Base Object"""
  def get_texture_tachyon(self, color, f):
    texkey="texture_%f%f%f"%(color[0], color[1], color[2])
    texval=self.tachyon_textures.get(texkey,False)
    if texval: return texkey # tex is a texture spec string
    else:
      texval="\nTEXDEF %s AMBIENT 0.3  DIFFUSE 0.8  SPECULAR 0.15  OPACITY 1.0\nPHONG PLASTIC 0 PHONG_SIZE 100000\n"%texkey
      texval=texval+"  COLOR %f %f %f TEXFUNC 0\n\n"%(color[0], color[1], color[2])
      f.write(texval)
      texind=len(self.texture_map)
      self.tachyon_textures[texkey]=texval
      self.tachcolors.append(texkey)
    return texkey

  def get_texture_json(self, color):
    texkey="texture_%f%f%f"%(color[0], color[1], color[2]) # TODO: proper hash value
    texind=self.texture_map.get(texkey,-1)
    if texind>=0: return texind # tex is a texture spec string
    else:
      texval={"color"     : [float(color[0]), float(color[1]), float(color[2])],
              "ambient"   : 0.3,
              "diffuse"   : 0.8,
              "specular"  : 0.15,
              "opacity"   : 1.0
              }
      texind=len(self.texture_map)
      self.texcolors.append([float(color[0]), float(color[1]), float(color[2])])
      self.json_scene["textures"].append(texval)
      self.texture_map[texkey]=texind

    return texind


  def export_tachyon(self, f):
    pass

  def export_povray(self, f):
    pass

  def export_json(self):
    pass

  def get_maxpoints(self):
    return numpy.array([self.points.T[0].max(), self.points.T[1].max(), self.points.T[2].max()])

  def get_minpoints(self):
    return numpy.array([self.points.T[0].min(), self.points.T[1].min(), self.points.T[2].min()])

class ConeObject(Object):
  """Cone Object"""
  obj_type=GL_POINTS
  def __init__(self, points, radii, colors=None, sides=20):
    self.points=points
    self.radii=radii
    self.sides=sides

    if colors!=None:
      if len(colors)==3:
        self.colors=None
        self.color=colors
      else:
        self.colors=colors
    else:
      self.colors=colors
      self.color=white
  def renderObject(self):
    if self.color!=None:
      glColor3f(self.color[0], self.color[1], self.color[2])
    gleSetNumSides(self.sides)
    glePolyCone(self.points, self.colors, self.radii)

  def export_tachyon_triangle(self, f, pts, texkey):
    f.write("TRI\n")
    for i in range(3):
      f.write("  V%d %f %f %f\n"%(i, pts[i][0], pts[i][1], pts[i][2]))
    f.write("  %s\n"%texkey)

  def export_tachyon(self,f):
    texkey=self.get_texture_tachyon(self.color, f)
    texkey2=self.get_texture_tachyon(white, f)
    points=numpy.zeros([self.sides+1,3])
    f.write("# Cone object\n")

    for j in range(1, len(self.points)-1):
      f.write("#   Circle %d\n"%j)
      for i in range(self.sides+1):
        step=float(i%self.sides)
        points[i,0]=self.radii[j]*math.sin( 2*math.pi*(step/self.sides))
        points[i,1]=self.radii[j]*math.cos( 2*math.pi*(step/self.sides))
        points[i,2]=0

      v1=self.points[j]
      v2=self.points[j-1]
      up=numpy.array([0,0,1])
      rotMatrix=numpy.ones([4,4])
      uviewpoint(rotMatrix, v1, v2, up)

      rotMatrix=rotMatrix[:3,:3]
      points2=numpy.zeros([self.sides+1,3])
      for i in range(self.sides+1):
        points2[i]=numpy.dot(rotMatrix,points[i])+v1

      if j>1:
        for i in range(1, self.sides+1):
          self.export_tachyon_triangle(f,[prevpoints2[i-1], prevpoints2[i], points2[i]], texkey)
          if self.radii[j]>0.0: # points2[i] == points2[i-1]
            self.export_tachyon_triangle(f,[prevpoints2[i-1], points2[i], points2[i-1]], texkey)
      f.write("\n")
      prevpoints2=points2

  def export_json(self):
    cones=self.json_scene["objects"]["cones"]
    texind=self.get_texture_json(self.color)

    # Map polyline coords to floats
    pts=[]
    for p in self.points:
      pts.append(map(float,p))

    conedict= { "points"  : pts,
                "radii"   : map(float,self.radii),
                "texture" : texind
              }

    cones.append(conedict)



class PointsObject(Object):
  """Points Object"""
  obj_type=GL_POINTS
  def __init__(self, points, width=1., colors=None):
    self.points=points
    self.width=width

    if colors!=None:
      if len(colors)==3:
        self.colors=None
        self.color=colors
      else:
        self.colors=colors
    else:
      self.colors=colors
      self.color=white

  def renderObject(self):
    if self.color:
      glColor3f(self.color[0], self.color[1], self.color[2])
    glBegin(self.obj_type)
    for i in range(self.points.shape[0]):
      glVertex3fv(self.points[i])
    glEnd()

class LineObject(PointsObject):
  """Line Object"""
  obj_type=GL_LINES
  def __init__(self, points, width=1., colors=None):
    self.width=width
    self.points=points

    if colors!=None:
      if len(colors)==3:
        self.colors=None
        self.color=colors
      else:
        self.colors=colors
    else:
      self.colors=colors
      self.color=white

  def renderObject(self):
    if self.color:
      glColor3f(self.color[0], self.color[1], self.color[2])
    if len(self.points.shape)==2:
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointerf(self.points)
      glDrawElementsui(GL_LINES, numpy.arange(0, self.points.shape[0],1))
      glVertexPointerf(self.points[1:-1])
      glDrawElementsui(GL_LINES, numpy.arange(0, self.points.shape[0]-2,1))
      glDisableClientState(GL_VERTEX_ARRAY)
    elif len(self.points.shape)==3:
      glEnableClientState(GL_VERTEX_ARRAY);
      for i in xrange(self.points.shape[0]):
        glVertexPointerf(self.points[i])
        glDrawElementsui(GL_LINES, numpy.arange(0, self.points[i].shape[0],1))

        glVertexPointerf(self.points[i][1:-1])
        glDrawElementsui(GL_LINES, numpy.arange(0, self.points[i].shape[0]-2,1))

      glDisableClientState(GL_VERTEX_ARRAY)

  def get_maxpoints(self):
    if len(self.points.shape)==2:
      return numpy.array([self.points[0].max(), self.points[1].max(), self.points[2].max()])
    elif len(self.points.shape)==3:
      for i in xrange(self.points.shape[0]):
        return numpy.array([self.points[i][0].max(), self.points[i][1].max(), self.points[i][2].max()])

  def get_minpoints(self):
    if len(self.points.shape)==2:
      return numpy.array([self.points[0].min(), self.points[1].min(), self.points[2].min()])
    elif len(self.points.shape)==3:
      for i in xrange(self.points.shape[0]):
        return numpy.array([self.points[i][0].min(), self.points[i][1].min(), self.points[i][2].min()])

class PolyCylinderObject(Object):
  """Polycylinder Object"""
  def __init__(self, points, radius=0.5, colors=None, drawall=True, nsides=20, texinds=None):
    self.radius=radius
    self.drawall=drawall
    self.nsides=nsides
    self.texinds=[]
    self.tachinds=[]
    texlookup=[]
    try:
      f=float(colors[0])
      single_color=True
    except:
      single_color=False

    if colors!=None:
      if single_color:
        self.colors=None
        self.color=colors
        self.texind=self.get_texture_json(colors)
      else:
        self.color=None
        self.colors=[]
        self.texind=-1
        map(texlookup.append, map(self.get_texture_json, colors)) # Fill texture map
        map(self.texinds.append, map(lambda i:texlookup[i], texinds))
        map(self.colors.append, map(lambda i: colors[i], texinds))

    else:
      self.colors=None
      self.color=white
      self.texind=self.get_texture_json(green)

    #~ print "texlookup =", texlookup
    #~ print "texinds =",texinds
    #~ print "self.texinds =",self.texinds
    #~ if self.colors:
      #~ print "len(self.colors) =",len(self.colors)

    if drawall:
      x1=points[0]
      x2=points[1]

      # Set t1=1, t2=-1, t0=3, t3=-3
      t0=3
      t3=-3
      c1=(x2+x1)/2
      m1=x1 -c1
      x0=m1*t0+c1

      xm1=points[-2]
      xm2=points[-1]

      c2=(xm2+xm1)/2
      m2=xm1 -c2
      xl=m2*t3+c2

      points=numpy.vstack([x0, points, xl])
      if self.colors:
        self.texinds=[self.texinds[0]]+self.texinds+[self.texinds[-1]]
        col=[colors[0]]
        col+=self.colors
        col.append(colors[-1])
        self.colors=col
      #~ print "colors=%s/%s, radius=%s, points=\n%s"%(self.colors, self.color, self.radius,self.points)

    self.points=points
    if type(self.points)==numpy.matrix:
      self.points=numpy.array(self.points)
    #~ print "len(points)=%d"%len(self.points)
    #~ if self.colors:
      #~ print "len(self.colors) =",len(self.colors)

  def write_tachyon_segment(self, col, points, f):
    texkey=self.get_texture_tachyon(col,f)

    f.write("\nPOLYCYLINDER\n  POINTS %d\n"%(len(points)))
    for p in points:
      f.write("    %f %f %f\n"%(p[0],p[1],-p[2]))
    f.write("  RAD %f\n"%self.radius)
    f.write("  %s\n\n"%texkey)

  def export_tachyon(self, f):
    if self.texind>=0:
      col=self.color
      points=self.points[1:-1]
      self.write_tachyon_segment(col, points, f)
    else:
      prevind=-1
      startind=1
      lpmi=len(self.points)-1
      segment=0
      f.write("\n\n# Polycylinder object "+"-"*40+"\n")
      for i in xrange(1,lpmi):
        ind=self.texinds[i]
        #~ print "%d: ind=%d, prevind=%d, color=%s"%(i, ind, prevind, self.colors[i])
        if (ind!=prevind and i>1) or i==lpmi-1:
          col=self.colors[i-1]
          #~ print "col=%s"%col

          mag=0
          si=startind
          while si>1 and si<i and mag<self.radius/5:
            x=self.points[si]-self.points[startind-1]
            mag=numpy.sqrt(numpy.vdot(x,x))
            si+=1
            #~ print "rad/100=%f, mag=%f"%(self.radius/5, mag)

          startind=si-1
          points=self.points[startind:i]
          #~ print "writing segment %d -> %d, len=%f"%(startind,i, mag)
          f.write("# Polycylinder segment %d\n"%segment)
          self.write_tachyon_segment(col, points, f)
          startind=i
          segment+=1
        prevind=ind

  def export_json(self):
    polys=self.json_scene["objects"]["polycylinders"]
    #~ print "self.color=",self.color
    #~ print "self.colors=",self.colors

    # Map polyline coords to floats
    pts=[]

    for p in self.points:
      try:
        pts.append(map(float,p))
      except:
        print p.shape,type(p),p
        raise

    #~ if self.texind<0:
      #~ print "exporting polycylinder len(pts)=%d, len(texinds)=%d"%(len(pts), len(self.texinds))

    polydict= { "points"  : pts,
                "radius"  : float(self.radius),
                "texture" : self.texind,
                "textures": self.texinds
              }

    polys.append(polydict)


  def export_povray(self, f):
    if len(self.points)>3:
      if len(self.points)>4:f.write("merge {\n")
      for i in range(1,len(self.points)-2):
        q=self.points[i]
        r=self.points[i+1]
        s=\
"""  cylinder {
              <%f, %f, %f>, <%f, %f, %f>, %f
              pigment { color Yellow }
     }

"""%(q[0],q[1],q[2], r[0],r[1],r[2], self.radius)

        f.write(s)
      if len(self.points)>4: f.write("}\n")
      f.write("\n")

  def renderObject(self):
    if self.color!=None:
      glColor3f(self.color[0], self.color[1], self.color[2])
      gleSetNumSides(self.nsides);
    glePolyCylinder(self.points, self.colors, self.radius)

class SphereObject(Object):
  """Sphere Object"""
  def __init__(self, quadric, origin, radius=0.5, color=None, slices=20,
               stacks=20):
    self.radius=radius
    self.origin=origin
    self.slices=slices
    self.stacks=stacks
    self.quadric=quadric

    if color!=None:
      if len(color)==3:
        self.color=None
        self.color=color
      else:
        self.color=color
    else:
      self.color=white

  def export_tachyon(self, f):
    #~ print "exporting SphereObject(%s)"%id(self)
    texkey=self.get_texture_tachyon(self.color, f)
    f.write("SPHERE CENTER ")
    f.write(" %f %f %f"%(self.origin[0], self.origin[1], -self.origin[2]))
    f.write(" RAD %f\n"%self.radius)
    f.write("  %s\n"%texkey)

  def export_povray(self,f):
    s=\
"""  sphere {
      <%f, %f, %f>, %f
      pigment { color Blue }
    }\n"""%(self.origin[0], self.origin[1], self.origin[2], self.radius)

    f.write(s)

  def export_json(self):
    spheres=self.json_scene["objects"]["spheres"]
    texind=self.get_texture_json(self.color)
    spheredict={"center"  : map(float, self.origin),
                "radius"  : float(self.radius),
                "texture" : texind
                }
    spheres.append(spheredict)

  def renderObject(self):
    if self.color!=None:
      glColor3f(self.color[0], self.color[1], self.color[2])

    glTranslatef(self.origin[0], self.origin[1], self.origin[2])
    gluSphere(self.quadric, self.radius, self.slices, self.stacks)

  def get_maxpoints(self):
    return numpy.array([self.origin[0]+self.radius, self.origin[1]+self.radius, self.origin[2]+self.radius])

  def get_minpoints(self):
    return numpy.array([self.origin[0]-self.radius, self.origin[1]-self.radius, self.origin[2]-self.radius])


class DnaCoil(OpenGLdisplay.Base):
  points=None
  radius=None

  def __init__(self):
    self.objects=[]
    self.quadric=None
    if OpenGL:
      OpenGLdisplay.Base.__init__(self)
      print "Doing GL init..."
      self.quadric = gluNewQuadric()
      gluQuadricNormals(self.quadric, GLU_SMOOTH)         # Create Smooth Normals (NEW)
      gluQuadricTexture(self.quadric, GL_TRUE)                    # Create Texture Coords (NEW)

      glEnable(GL_TEXTURE_2D)
      glClearColor(0.0, 0.0, 0.0, 0.0)                    # This Will Clear The Background Color To Black
      #glClearDepth(1.0)                                  # Enables Clearing Of The Depth Buffer
      #glDepthFunc(GL_LESS)                               # The Type Of Depth Test To Do
      #glEnable(GL_DEPTH_TEST)                            # Enables Depth Testing
      glShadeModel(GL_SMOOTH)                             # Enables Smooth Color Shading

      glMatrixMode(GL_PROJECTION)
      glLoadIdentity()                                    # Reset The Projection Matrix

      #glutReshapeFunc(self.reshape)

  def addPolyCylinder(self,points ,colors=None, radius=0.5, drawall=True, nsides=20, texinds=None):
    self.objects.append(PolyCylinderObject(points, radius=radius,
        colors=colors, drawall=drawall, nsides=nsides, texinds=texinds))

  def addPoints(self,points,colors=None, width=1):
    self.objects.append(PointsObject(points, width=width, colors=colors))

  def addPolyLine(self,points, colors=None, width=1):
    self.objects.append(LineObject(points, width=width, colors=colors))

  def addSphere(self, origin, radius=0.5, color=None, slices=20,
               stacks=20):
    self.objects.append(SphereObject(self.quadric, origin, radius=radius, \
                        color=color, slices=slices,stacks=stacks))
  def addCone(self,points, radii, colors=None, nsides=20):
      self.objects.append(ConeObject(points, radii, colors=colors, sides=nsides))

  def get_extremepoints(self):
    maxpoints=-1000*numpy.zeros(3)
    minpoints=1000*numpy.zeros(3)
    for obj in self.objects:
      #~ print obj.__doc__
      mp=obj.get_maxpoints()
      #~ print "maxpoints=%s"%repr(mp)
      maxpoints=numpy.vstack([maxpoints, mp])
      mp=obj.get_minpoints()
      #~ print "minpoints=%s"%repr(mp)
      minpoints=numpy.vstack([minpoints, mp])


    maxp=numpy.zeros(3)
    minp=numpy.zeros(3)
    for i in range(3):
      maxp[i]=maxpoints.T[i].max()
      minp[i]=minpoints.T[i].min()

    return maxp,minp


  def export_tachyon(self, fname="nudraw.dat"):
    #~ print "Exporting to Tachyon"
    f=open(fname,"w")
    f.write(tachyon_scene_header)
    maxp,minp=self.get_extremepoints()

    width=maxp[0]-minp[0]
    height=maxp[1]-minp[1]

    mx=max([width,height])
    mn=min([width,height])
    fac=0.65+math.e**(-0.65*numpy.log10(math.sqrt(mn)))

    camera=numpy.array([(maxp[0]+minp[0])/2, (maxp[1]+minp[1])/2, -fac*mx])
    light=numpy.array([(maxp[0]+minp[0])/2, 150+(maxp[1]+minp[1])/2, -max([width,height])/2])


    #~ print "camera=%s"%camera

    f.write(tachyon_scene_camera%(camera[0], camera[1], camera[2]))


    f.write(tachyon_scene_light%(light[0], light[1], light[2]))

    for obj in self.objects:
      obj.export_tachyon(f)
    f.write(tachyon_scene_footer)
    f.close()

  def export_json(self, fname="nudraw.json"):
    if len(self.objects)==0:
      raise ValueError("No objects defined")
    f=open(fname,"w")
    maxp,minp=self.get_extremepoints()

    width=maxp[0]-minp[0]
    height=maxp[1]-minp[1]

    mx=max([width,height])
    mn=min([width,height])
    fac=0.65+math.e**(-0.65*numpy.log10(math.sqrt(mn)))
    #~ print "width=%f, height=%f, mn=%f, fac=%f"%(width,height, mn, fac)

    camera=numpy.array([(maxp[0]+minp[0])/2, (maxp[1]+minp[1])/2, fac*mx])
    light=numpy.array([(maxp[0]+minp[0])/2, 150+(maxp[1]+minp[1])/2, -max([width,height])/2])

    #~ print "camera=%s"%camera

    scene=self.objects[0].json_scene
    lookat=map(float, camera)
    lookat[2]=0.0
    scene["camera"]={ "position": map(float, camera),
                      "lookat"  : lookat,
                      "updir"   : [0., 1., 0.]
                    }

    scene["lights"].append({"position": map(float, camera),
                            "color"   : [0., 1., 0.]})


    for obj in self.objects:
      obj.export_json()
    #~ print json
    #~ print dir(json)
    f.write(json.write(scene))
    f.close()


  def export_povray(self, fname="nudraw.pov"):
    #~ print "Exporting to Tachyon"
    f=open(fname,"w")
    f.write(povray_scene_header)
    maxp,minp=self.get_extremepoints()

    width=maxp[0]-minp[0]
    height=maxp[1]-minp[1]
    #~ print "width=%f, height=%f"%(width,height)
    camera=numpy.array([(maxp[0]+minp[0])/2, (maxp[1]+minp[1])/2, -max([width,height])])
    light=numpy.array([(maxp[0]+minp[0])/2, 150+(maxp[1]+minp[1])/2, -max([width,height])/2])

    #~ print "camera=%s"%camera

    f.write(povray_scene_camera%(camera[0], camera[1], camera[2]))


    f.write(povray_scene_light%(light[0], light[1], light[2]))

    for obj in self.objects:
      obj.export_povray(f)
    f.write(povray_scene_footer)
    f.close()


  # draw the polys
  def Draw(self):
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    # set up some matrices so that the object spins with the mouse
    gleSetJoinStyle (TUBE_NORM_EDGE | TUBE_JN_ANGLE | TUBE_JN_CAP)
    glPushMatrix ()
    glTranslatef (0.0, 0.0, 30.0)
    glRotatef (self.lastx, 0.0, 1.0, 0.0)
    glRotatef (self.lasty, 1.0, 0.0, 0.0)

    for obj in self.objects:
      glPushMatrix()
      obj.renderObject()
      glPopMatrix()

    glPopMatrix ()
    glFlush()
    glutSwapBuffers ()

  def reshape(self, w, h):
    glViewport(0, 0, w, h)
    glMatrixMode (GL_PROJECTION)
    glLoadIdentity()
    if w <= h:
      glOrtho(-2.5, 2.5, -2.5*h/w,
               2.5*h/w, -10.0, 10.0)
    else:
      glOrtho(-2.5*w/h,
               2.5*w/h, -2.5, 2.5, -10.0, 10.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

def drawCone(coil):
  #~ cone_points=numpy.array([[0,-1,0.0], [0,2.5,0.0], [0,5,0.0], [0,7.5,0.0], [0,10.0,0.0]])
  cone_points=numpy.array([[0,0.0,-1], [0,0.0,5], [0,0.0,7.5], [0,0.0,10.0]])
  w1=5.0
  w2=1.0
  cone_radii=numpy.array([0.0, w1, 0.0, 0.0])
  #~ print "cone_points=\n%s"%cone_points
  coil.addCone(cone_points, cone_radii, colors=red, nsides=10)


def drawCircle(coil):
  NPOINTS=10
  points=numpy.zeros([NPOINTS,3])
  for i in range(NPOINTS):
      points[i,0]=5*math.sin( 2*math.pi*(float(i)/NPOINTS))
      points[i,1]=5*math.cos( 2*math.pi*(float(i)/NPOINTS))
      points[i,2]=0
      coil.addSphere(points[i], radius=0.12)


  #~ coil.addSphere(numpy.array([0,0,10]), color=red, radius=0.12)
  #~ coil.addSphere(numpy.array([0,0,-5]), color=yellow, radius=0.12)

def drawConeReal(coil):
  points=numpy.array([[ 10.61856106, -17.64698602,  -0.53877544],
          [ 10.57411878, -17.69977316,  -0.56916634],
          [ 10.56841763, -17.70655436,  -0.5728937 ],
          [ 10.55158445, -17.72657646,  -0.58389907],
          [ 10.52441728, -17.75889026,  -0.60166071],
          [ 10.48820405, -17.80196383,  -0.62533658],
          [ 10.44466155, -17.85375514,  -0.65380425],
          [ 10.39585406, -17.91180886,  -0.68571412],
          [ 10.34409544, -17.97337278,  -0.71955342],
          [ 10.29183945, -18.03552829,  -0.75371789],
          [ 10.24156345, -18.09532871,  -0.78658786],
          [  8.80480819, -19.80426691,  -1.72592483],
          [  8.77349385, -19.84151353,  -1.74639785],
          [  8.74446306, -19.876044,    -1.76537791],
          [  8.71821255, -19.90726748,  -1.78254025],
          [  8.69519146, -19.93464974,  -1.79759121],
          [  8.67579371, -19.95772226,  -1.81027328],
          [  8.66035118, -19.97609026,  -1.82036946],
          [  8.64912811, -19.98943946,  -1.82770699],
          [  8.64231653, -19.99754145,  -1.83216034],
          [  8.64003298, -20.00025761,  -1.83365331]])
  radii=numpy.array([ 0.,          2.,          2.09739783,  2.19017822,  2.27394264,  2.34472,
          2.39915487,  2.43466662,  2.4495717,   2.44316349,  2.41574579,  0.51968224,
          0.48713929,  0.44626125,  0.39774756,  0.3424283,   0.28125,     0.21525943,
          0.14558571,  0.07342098,  0.        ])

  coil.addPolyCylinder(points, radius=0.05, colors=white, nsides=5)
  coil.addCone(points, radii, colors=red, nsides=20)

if __name__ == '__main__':
  NPOINTS=400
  points=numpy.zeros([NPOINTS,3])
  for i in range(NPOINTS):
      points[i,0]=5*math.sin( 6*math.pi*(float(i)/NPOINTS))
      points[i,1]=3*math.cos( 6*math.pi*(float(i)/NPOINTS))
      points[i,2]=float(i)/NPOINTS*5. -2.5


  coil=DnaCoil()
  drawConeReal(coil)
  drawCircle(coil)
  #~ coil.addPolyCylinder(points, colors=[109./255, 216./255,  45./255], radius=0.5, drawall=False)
#~
  #~ points2=numpy.zeros([NPOINTS,3])
  #~ for i in range(NPOINTS):
    #~ points2[i,0]=1.25*math.cos( 6*math.pi*(float(i)/NPOINTS))
    #~ points2[i,1]=1.25*math.sin( 6*math.pi*(float(i)/NPOINTS))
    #~ points2[i,2]=float(i)/NPOINTS*20. -10
  #~ coil.addPolyCylinder(points2, colors=[1, .2857, 0], radius=0.5, drawall=False)
#~
  #~ points2=numpy.zeros([NPOINTS+1,3])
  #~ for i in range(NPOINTS+1):
    #~ points2[i,0]=8*math.cos( 2*math.pi*(float(i)/NPOINTS))
    #~ points2[i,1]=8*math.sin( 2*math.pi*(float(i)/NPOINTS))
    #~ points2[i,2]=0.;
  #~ coil.addPolyLine(points2, colors=red)
#~
  #~ minp=-20.
  #~ maxp=20
  #~ n=11
  #~ d=maxp-minp
  #~ step=d/(n-1)
  #~ points3=numpy.zeros([n,3])
  #~ points3[:,0]=points3[:,1]=numpy.zeros(n)
  #~ z=numpy.arange(minp,maxp+step/2,step)
  #~ print z
  #~ points3[:,2]=z
  #~ print points3
  #~ coil.addPolyCylinder(points3, radius=1,colors=white)
#~
  #~ points4=numpy.zeros([NPOINTS,3])
  #~ for i in range(NPOINTS):
    #~ points4[i,0]=7*math.cos( 7*math.pi*(float(i)/NPOINTS))
    #~ points4[i,1]=7*math.sin( 7*math.pi*(float(i)/NPOINTS))
    #~ points4[i,2]=float(i)/NPOINTS*20. -10
  #~ coil.addPolyLine(points4, colors=orange, width=0.5)
  #~ coil.addPolyLine(numpy.vstack([[points4+numpy.array([0,0,2])], [points4+numpy.array([2,0,0])]]), colors=white, width=0.5)
#~
  #~ coil.addSphere([0,0,20],2,color=green)
  #~ coil.addSphere([0,0,-20],2,color=white)
  #~ cone_points=numpy.array([[0,0,20], [5,0,20], [10,0,20], [15,0,20]])
  #~ w1=0.5
  #~ w2=1.0
  #~ cone_radii=numpy.array([w1,w2,0.0,0.0])
  #~ print "cone_points=%s"%cone_points
  #~ coil.addCone(cone_points, cone_radii, colors=turquoise, sides=50)
  print("exporting %s.dat"%__name__)
  coil.export_tachyon(fname="%s.dat"%__name__)
  #~ coil.export_json(fname="%s.json"%__name__)
#~
  if OpenGL:
    coil.Display()
