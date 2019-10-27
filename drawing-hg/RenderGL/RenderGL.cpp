// C++ OpenGL Renderer
// Conrad Steenberg <conrad.steenberg@caltech.edu>
// Oct 28, 2008

#include <iostream>
#include <GL/glut.h>
#include <GL/gle.h>
#include <vector>
#include "Objects.hpp"

float zoom = 80.0f;
float rotx = 0;
float roty = 0.001f;
float tx = 0;
float ty = 0;
int lastx=0;
int lasty=0;
unsigned char Buttons[5] = {0};

using namespace std;
vector<Object *> objects;
int mods;

//------------------------------------------------------------------------------
/// \brief  Initialises the openGL scene
///
void Init()
{
  GLfloat LightDiffuse[]= { 1.0f, 1.0f, 1.0f, 1.0f };
  GLfloat LightDiffuseDim[]= { 0.5f, 0.5f, 0.5f, 1.0f };
  GLfloat LightPosition1[]= { 40.0, 40, 100.0, 0.0 };
  GLfloat LightPosition2[]= { -40.0, 40, 100.0, 0.0 };
  glShadeModel(GL_SMOOTH);
//  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glClearColor(0.0, 0.0, 0.0, 1.0);                    // This Will Clear The Background Color To Black
  glClearDepth(1.0);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_TEXTURE_2D);

  glDisable( GL_CULL_FACE );
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuseDim);
  glLightfv(GL_LIGHT0, GL_POSITION,LightPosition1);
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuseDim);
  glLightfv(GL_LIGHT1, GL_POSITION,LightPosition2);
  glEnable(GL_LIGHT1);

  glEnable(GL_LIGHTING);
  glColorMaterial (GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable (GL_COLOR_MATERIAL);

  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glFlush();
}

//-------------------------------------------------------------------------------
/// \brief  Draws the scene
///
void display()
{

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity();

  glTranslatef(0,0,-zoom);
  glTranslatef(tx,ty,0);
  glRotatef(rotx,1,0,0);
  glRotatef(roty,0,1,0);

  for (int i=0; i<objects.size();i++){
    //~ cout<<"Drawing Object "<<i<<endl;
    objects[i]->Draw();
  }

  glutSwapBuffers();
}

//-------------------------------------------------------------------------------
/// \brief  Called when the screen gets resized
/// \param  w - the new width
/// \param  h - the new height
///
void reshape(int w, int h)
{
  // prevent divide by 0 error when minimised
  if(w==0)
    h = 1;

  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45,(float)w/h,0.1,400);
  //~ gluPerspective(
            //~ #~ 45, # field of view in degrees
            //~ #~ width/float(height or 1), # aspect ratio
            //~ #~ 100, # near clipping plane
            //~ #~ 200, # far clipping plane
    //~ #~ )
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}


//-------------------------------------------------------------------------------
//
void Motion(int x,int y)
{
  int diffx=x-lastx;
  int diffy=y-lasty;

  lastx=x;
  lasty=y;

  if( Buttons[0] && mods == GLUT_ACTIVE_CTRL )
  {
    zoom -= (float) 0.2f * diffx;
  }
  else if( Buttons[0] && mods == GLUT_ACTIVE_SHIFT)
  {
    tx += (float) 0.2f * diffx;
    ty -= (float) 0.2f * diffy;
  }  
  else if( Buttons[0] )
  {
    rotx += (float) 0.5f * diffy;
    roty += (float) 0.5f * diffx;
  }

  glutPostRedisplay();
}

//-------------------------------------------------------------------------------
//
void Mouse(int b,int s,int x,int y)
{
  lastx=x;
  lasty=y;
  mods=glutGetModifiers();

  switch(b)
  {
  case GLUT_LEFT_BUTTON:
    Buttons[0] = ((GLUT_DOWN==s)?1:0);
    break;
  case GLUT_MIDDLE_BUTTON:
    Buttons[1] = ((GLUT_DOWN==s)?1:0);
    break;
  case GLUT_RIGHT_BUTTON:
    Buttons[2] = ((GLUT_DOWN==s)?1:0);
    break;
  case 3:
    zoom -= (float) 0.2f;
    Buttons[3] = ((GLUT_DOWN==s)?1:0);
    break;
  case 4:
    zoom += (float) 0.2f;
    Buttons[4] = ((GLUT_DOWN==s)?1:0);
    break;
  default:
    break;
  }
  glutPostRedisplay();
}

//-------------------------------------------------------------------------------
///
int InitViewer(int argc,char** argv)
{
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_ALPHA);
  glutInitWindowSize(640,480);
  glutInitWindowPosition(100,100);
  glutCreateWindow("Nudraw Interactive Viewer");


  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(Mouse);
  glutMotionFunc(Motion);

  Init();
}
