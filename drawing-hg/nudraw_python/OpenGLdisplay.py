#!/usr/bin/python

# This is statement is required by the build system to query build info
if __name__ == '__build__':
  raise Exception


import sys
import math
from OpenGL.GL import *
from OpenGL.GLE import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

class Base:
  zoom_default=1.0
  def __init__(self):
    self.lastx=self.pressx=0
    self.lasty=self.pressy=0
    self.button=None
    self.state=None
    self.zoom=self.zoom_default
    self.xoffset=0.0
    self.yoffset=0.0

  def MouseEvent(self, button, state, x, y):
    self.button=button
    self.state=state

    self.mods=mods=glutGetModifiers()
    if button==0 and state==0:
      self.pressx = x
      self.pressy = y

    if button==1 and state==0:
      if mods == GLUT_ACTIVE_CTRL:
        self.xoffset=self.yoffset=0.0
      if mods == GLUT_ACTIVE_SHIFT:
        pass
      else:
        self.lastx=self.lasty=x=y=0.0
      self.pressx = x
      self.pressy = y

    if button==2 and state==0:
      if mods == GLUT_ACTIVE_CTRL:
        self.zoom=self.zoom_default
      self.pressx = x
      self.pressy = y

    if button==3 and state==0 and self.zoom>1.:
      if mods == GLUT_ACTIVE_CTRL:
        self.zoom=self.zoom-1.

    if button==4 and state==0:
      if mods == GLUT_ACTIVE_CTRL:
        self.zoom=self.zoom+1.

    x,y,w,h = glGetDoublev(GL_VIEWPORT)
    self.Reshape(int(w),int(h))

    glutPostRedisplay ()

  # get notified of mouse motions
  def MouseMotion (self, x, y):
    reshape=False
    if self.button==0 and self.state==0:
      if self.mods == GLUT_ACTIVE_CTRL:
        dx=x-self.pressx
        dy=y-self.pressy

        adx=math.fabs(dx)
        ady=math.fabs(dy)

        try:
          sgnx=dx/adx
          dx=sgnx*math.log(adx)
        except ZeroDivisionError:
          dx=0

        try:
          sgny=dy/ady
          dy=sgny*math.log(ady)
        except ZeroDivisionError:
          dy=0


        self.xoffset = self.xoffset + dx
        self.yoffset = self.yoffset + dy
        reshape=True
      elif self.mods == GLUT_ACTIVE_SHIFT:
        reshape=True
        self.zoom = self.zoom + (y-self.pressy)
        if self.zoom<4.0:
          self.zoom=4.0
      else:
        self.lastx = self.lastx + (x-self.pressx)
        self.lasty = self.lasty + (y-self.pressy)
      self.pressx=x
      self.pressy=y

    if self.button==2 and self.state==0:
#      dist=math.sqrt((x-self.pressx)**2 + (y-self.pressy)**2)
      self.zoom = self.zoom + (y-self.pressy)
      if self.zoom<4.0:
        self.zoom=4.0
      #self.lasty = self.lasty + (y-self.pressy)
      self.pressx=x
      self.pressy=y
      reshape=True

    if reshape:
      x,y,w,h = glGetDoublev(GL_VIEWPORT)
      self.Reshape(int(w),int(h))

    glutPostRedisplay ()

  def Reshape(self, width, height):
    h = float(height) / float(width);

    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    glFrustum (-self.zoom-self.xoffset, self.zoom-self.xoffset,
    -h*self.zoom+self.yoffset, h*self.zoom+self.yoffset, 20, 1000.0)
    x,y,width,height = glGetDoublev(GL_VIEWPORT)

    glMatrixMode(GL_MODELVIEW)

    glLoadIdentity()
    glTranslatef(0.0, 0.0, -200.0)
    glFlush()

  def JoinStyle (msg):
    sys.exit(0)

  # Overload this function
  def Draw(self):
    pass

  def Display(self):
    global glutDisplayFunc, glutMotionFunc
    # initialize glut
    glutInit(sys.argv)
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
    glutCreateWindow("Nudraw Display")
    glutDisplayFunc(self.Draw)
    glutMotionFunc(self.MouseMotion)
    glutMouseFunc(self.MouseEvent)
    glutReshapeFunc(self.Reshape)


    # initialize GL */
    glClearDepth (1.0)
    glEnable (GL_DEPTH_TEST)
    glClearColor (0.0, 0.0, 0.0, 0.0)
    glShadeModel (GL_SMOOTH)

    glMatrixMode(GL_MODELVIEW)

    # set up a light
    lightOnePosition = (40.0, 40, 100.0, 0.0)
    lightOneColor = (0.99, 0.99, 0.99, 1.0)

    lightTwoPosition = (-40.0, 40, 100.0, 0.0)
    lightTwoColor = (0.99, 0.99, 0.99, 1.0)

    # initialize lighting */
    glLightfv (GL_LIGHT0, GL_POSITION, lightOnePosition)
    glLightfv (GL_LIGHT0, GL_DIFFUSE, lightOneColor)
    glEnable (GL_LIGHT0)
    glLightfv (GL_LIGHT1, GL_POSITION, lightTwoPosition)
    glLightfv (GL_LIGHT1, GL_DIFFUSE, lightTwoColor)
    glEnable (GL_LIGHT1)
    glEnable (GL_LIGHTING)
    glColorMaterial (GL_FRONT_AND_BACK, GL_DIFFUSE)
    glEnable (GL_COLOR_MATERIAL)

    glutMainLoop ()
