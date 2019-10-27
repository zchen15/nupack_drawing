// Json scene file viewer
// Conrad Steenberg <conrad.steenberg@caltech.edu>
// Nov 6 2008

#include <iostream>
#include <cstdio>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/gle.h>
#include <GL/osmesa.h>
#include <stdlib.h>
#include <vector>
#include <gd.h>
#include "Objects.hpp"
#include "RenderGL.hpp"
#include "JsonScene.hpp"

using namespace std;
JsonScene Scene; // Don't you love it when callbacks force you to use global
                 // variables ;-)

extern float zoom;
extern float rotx;
extern float roty;
extern float tx;
extern float ty;
int live;
int printhelp;

void InitJsonGL()
{
  GLfloat LightDiffuse[]= { 1.0f, 1.0f, 1.0f, 1.0f };
  GLfloat LightDiffuseDim[]= { 0.4f, 0.4f, 0.4f, 1.0f };
  GLfloat LightDiffuseDim2[]= { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat LightPosition1[]= { -40.0, -40.0, 50.0, 0.0 };
  GLfloat LightPosition2[]= { 0.0, 0.0, 50.0, 0.0 };
  
  GLfloat LightPosition3[]= { -20.0, -20.0, 50.0, 0.0 };
  GLfloat LightPosition4[]= { -20.0,  20.0, 50.0, 0.0 };

//  glShadeModel(GL_SMOOTH);
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClearDepth(1.0);
  
  glEnable(GL_LIGHTING);
  //~ glEnable(GL_TEXTURE_2D);

  glDisable( GL_CULL_FACE );
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuseDim2);
  glLightfv(GL_LIGHT0, GL_SPECULAR, LightDiffuseDim);
  glLightfv(GL_LIGHT0, GL_AMBIENT, LightDiffuseDim2);
  glLightfv(GL_LIGHT0, GL_POSITION, LightPosition2);
  glEnable(GL_LIGHT0);

  glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuseDim2);
  glLightfv(GL_LIGHT1, GL_SPECULAR, LightDiffuseDim);
  //~ glLightfv(GL_LIGHT1, GL_AMBIENT, LightDiffuseDim2);  
  glLightfv(GL_LIGHT1, GL_POSITION, LightPosition1);
  glEnable(GL_LIGHT1);

  glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuseDim2);
  glLightfv(GL_LIGHT2, GL_SPECULAR, LightDiffuseDim);
  glLightfv(GL_LIGHT2, GL_POSITION, LightPosition3);
  glEnable(GL_LIGHT2);

  glLightfv(GL_LIGHT3, GL_DIFFUSE, LightDiffuseDim2);
  glLightfv(GL_LIGHT3, GL_SPECULAR, LightDiffuseDim);
  glLightfv(GL_LIGHT3, GL_POSITION, LightPosition4);
  glEnable(GL_LIGHT3);



  glEnable (GL_COLOR_MATERIAL);
  glColorMaterial (GL_FRONT, GL_SPECULAR);
  glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
  glMateriali(GL_FRONT, GL_SHININESS, 127);

  glEnable(GL_LIGHTING);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, LightDiffuseDim);
  glShadeModel(GL_SMOOTH);


  //~ glEnable (GL_LINE_SMOOTH);
  //~ glEnable (GL_POINT_SMOOTH);
  //~ glEnable(GL_POLYGON_SMOOTH);

//  glEnable (GL_BLEND);
//  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  //~ glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_NORMALIZE);

  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glFlush();
}


//------------------------------------------------------------------------------
/// \brief  Draws the scene
///
void displayJsonScene()
{
  int numObjects=Scene.objects.size();
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glClearColor (1.0, 1.0, 1.0, 1.0);
  glLoadIdentity(); 

  glPushMatrix();

  gluLookAt (Scene.eyeCam[0], Scene.eyeCam[1], Scene.eyeCam[2],
               Scene.centerCam[0], Scene.centerCam[1], Scene.centerCam[2],
               Scene.upCam[0], Scene.upCam[1], Scene.upCam[2]);

  glMatrixMode(GL_MODELVIEW);
  if (!live)
    glScalef(1.0, 1.0, -1.0);

  glTranslatef(0,0,-zoom);
  glTranslatef(tx,ty,0);
  glRotatef(rotx,1,0,0);
  glRotatef(roty,0,1,0);

  vector<Object *>::iterator it=Scene.objects.begin();
  for (int i=0; i<numObjects;i++, it++){
    //~ cout<<"Drawing Object "<<i<<endl;
    glPushMatrix();
    (*it)->Draw();
    glPopMatrix();
  }

  glPopMatrix();
  //~ glFlush();
  if (live)
    glutSwapBuffers();
}

//------------------------------------------------------------------------------
void reshapeJsonScene(int w, int h)
{
  // prevent divide by 0 error when minimised
  if(w==0)
    h = 1;

  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(55,(float)w/h,0.1,Scene.eyeCam[2]+100.0);
  //~ gluPerspective(
            //~ #~ 45, # field of view in degrees
            //~ #~ width/float(height or 1), # aspect ratio
            //~ #~ 100, # near clipping plane
            //~ #~ 200, # far clipping plane
    //~ #~ )
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

//------------------------------------------------------------------------------
int initJsonViewer(int width, int height)
{
  zoom=30.0;
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_ALPHA);
  glutInitWindowSize(width, height);
  glutCreateWindow("Nudraw Interactive Viewer");
  glClearColor (1.0, 1.0, 1.0, 0.5);

  glutDisplayFunc(displayJsonScene);
  glutReshapeFunc(reshapeJsonScene);
  glutMouseFunc(Mouse);
  glutMotionFunc(Motion);

  InitJsonGL();
}

//------------------------------------------------------------------------------
GLubyte *initJsonRenderer(int width, int height)
{
  zoom=15.0;
  GLubyte *buffer;

#ifndef USE_OSMESA
  cout<<"OSMESA support disabled"<<endl;
  return 0;
#else
  /* Create an RGBA-mode context */
#if OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
  /* specify Z, stencil, accum sizes */
  OSMesaContext ctx = OSMesaCreateContextExt( GL_RGBA, 24, 0, 0, NULL );
#else
  OSMesaContext ctx = OSMesaCreateContext( GL_RGBA, NULL );
#endif

  if (!ctx) {
    printf("OSMesaCreateContext failed!\n");
    return 0;
  }
  

  /* Allocate the image buffer */
  buffer = (GLubyte *)malloc( width * height * 4 * sizeof(GLubyte));
  if (!buffer) {
    printf("Alloc image buffer failed!\n");
    return 0;
  }

  /* Bind the buffer to the context and make it current */
  if (!OSMesaMakeCurrent( ctx, buffer, GL_UNSIGNED_BYTE, width, height )) {
    printf("OSMesaMakeCurrent failed!\n");
    return 0;
  }

//  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_ALPHA);
  //~ glutInitWindowSize(512,512);

  glClearColor (1.0, 1.0, 1.0, 1.0);

  reshapeJsonScene(width, height);
  InitJsonGL();
  return buffer;
#endif
}


//------------------------------------------------------------------------------
// Code from mesa/tests/osdemo.c, marked as public domain in their distribution
static void
write_targa_float(const char *filename, const GLfloat *buffer, int width, int height)
{
   FILE *f = fopen( filename, "w" );
   if (f) {
      int i, x, y;
      const GLfloat *ptr = buffer;
      fputc (0x00, f);  /* ID Length, 0 => No ID  */
      fputc (0x00, f);  /* Color Map Type, 0 => No color map included */
      fputc (0x02, f);  /* Image Type, 2 => Uncompressed, True-color Image */
      fputc (0x00, f);  /* Next five bytes are about the color map entries */
      fputc (0x00, f);  /* 2 bytes Index, 2 bytes length, 1 byte size */
      fputc (0x00, f);
      fputc (0x00, f);
      fputc (0x00, f);
      fputc (0x00, f);  /* X-origin of Image  */
      fputc (0x00, f);
      fputc (0x00, f);  /* Y-origin of Image  */
      fputc (0x00, f);
      fputc (width & 0xff, f);      /* Image Width  */
      fputc ((width>>8) & 0xff, f);
      fputc (height & 0xff, f);     /* Image Height */
      fputc ((height>>8) & 0xff, f);
      fputc (0x18, f);    /* Pixel Depth, 0x18 => 24 Bits */
      fputc (0x20, f);    /* Image Descriptor */
      fclose(f);
      f = fopen( filename, "ab" );  /* reopen in binary append mode */
      for (y=height-1; y>=0; y--) {
         for (x=0; x<width; x++) {
            int r, g, b;
            i = (y*width + x) * 4;
            r = (int) (ptr[i+0] * 255.0);
            g = (int) (ptr[i+1] * 255.0);
            b = (int) (ptr[i+2] * 255.0);
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;
            fputc(b, f); /* write blue */
            fputc(g, f); /* write green */
            fputc(r, f); /* write red */
         }
      }
   }
}

static void
write_targa_byte(const char *filename, const GLubyte *buffer, int width, int height)
{
   FILE *f = fopen( filename, "w" );
   if (f) {
      int i, x, y;
      const GLubyte *ptr = buffer;
      fputc (0x00, f);  /* ID Length, 0 => No ID  */
      fputc (0x00, f);  /* Color Map Type, 0 => No color map included */
      fputc (0x02, f);  /* Image Type, 2 => Uncompressed, True-color Image */
      fputc (0x00, f);  /* Next five bytes are about the color map entries */
      fputc (0x00, f);  /* 2 bytes Index, 2 bytes length, 1 byte size */
      fputc (0x00, f);
      fputc (0x00, f);
      fputc (0x00, f);
      fputc (0x00, f);  /* X-origin of Image  */
      fputc (0x00, f);
      fputc (0x00, f);  /* Y-origin of Image  */
      fputc (0x00, f);
      fputc (width & 0xff, f);      /* Image Width  */
      fputc ((width>>8) & 0xff, f);
      fputc (height & 0xff, f);     /* Image Height */
      fputc ((height>>8) & 0xff, f);
      fputc (0x18, f);    /* Pixel Depth, 0x18 => 24 Bits */
      fputc (0x20, f);    /* Image Descriptor */
      fclose(f);
      f = fopen( filename, "ab" );  /* reopen in binary append mode */
      for (y=height-1; y>=0; y--) {
         for (x=0; x<width; x++) {
            i = (y*width + x) * 4;
            fputc(ptr[i+2], f); /* write blue */
            fputc(ptr[i+1], f); /* write green */
            fputc(ptr[i], f);   /* write red */
         }
      }
   }
}

static void
write_png_byte(const char *filename, const GLubyte *buffer, int width, int height,
               int pngAA)
{
  int color=0xffffffff;
  int x,y,i;
  const GLubyte *ptr = buffer;

  gdImagePtr im = gdImageCreateTrueColor(width, height);
  gdImagePtr imscaled;

  if (pngAA>1)
    imscaled = gdImageCreateTrueColor(width/pngAA, height/pngAA);
  gdImageSaveAlpha (im, 1);

  for (y=height-1; y>=0; y--) {
    for (x=0; x<width; x++) {
      i = (y*width + x) * 4;
      color=gdTrueColor(ptr[i], ptr[i+1], ptr[i+2]);
      gdImageSetPixel (im, x, y, color);

    }
  }

  if (pngAA>1)
    gdImageCopyResampled (imscaled, im, 0, 0,
                          0, 0, width/pngAA, height/pngAA, width, height);


  FILE *out = fopen( filename, "w" );

  if (pngAA==1)
    gdImagePng (im, out);
  else{
    gdImagePng (imscaled, out);
    gdFree(imscaled);
  }

  gdFree(im);
  fclose(out);
}

//------------------------------------------------------------------------------
void printHelp(int argc, char** argv){
  cout<<"Usage: "<<argv[0]<<" [OPT] [OPT2=VAL]...\n\n";
  cout<<"Options:\n";
  cout<<
"\t--help\t\t\tPrint this message\n"
//"\t--show\t\t\tDisplay JSON scene using interactive viewer\n"
"\t--jsonfile=[FILENAME]\tJSON scene file to use\n"
#ifdef USE_OSMESA
"\t--tgafile=[FILENAME]\tName of TGA file to create, in non-interactive mode\n"
"\t--pngfile=[FILENAME]\tName of PNG file to create, in non-interactive mode\n"
"\t--pngaa=[INTEGER]\tOversampling number\n"
#endif
"\t--width=[INTEGER]\tWidth of window or image in pixels\n"
"\t--height=[INTEGER]\tHeight of window or image in pixels\n";
}

//------------------------------------------------------------------------------
int parseArguments(int argc, char** argv, string &jsonFileName,
                    string & tgaFileName, int *width, int *height,
                    string & pngFileName, int *pngaa){

  for (int i=1; i<argc; i++){
    string item(argv[i]);
    size_t l=item.length();
    bool toggle=false;
    string option;
    string value;
    size_t posmin=item.find("--",0)+2;
    size_t pos=item.find('=',posmin);
    if (posmin>l) posmin=0;
    if (pos>l) pos=l;


    //~ cout<<"\nl="<<l<<" pos="<<pos<<" posmin="<<posmin<<endl;
    if (pos==l){ // A switch found
      toggle=true;
      option=item.substr(posmin,pos);
      value="";
      //~ cout<<"toggle found"<<endl;
    }
    else{
      option=item.substr(posmin,pos-posmin);
      value=item.substr(pos+1,l-pos);
    }
    
    //~ cout<<"option: "<<option<<endl;
    //~ cout<<"value: "<<value<<endl;

    // Assign values to variables
    if (option=="help" && toggle){
      printHelp(argc, argv);
      //~ cout<<"Interactive viewer"<<endl;
      return 1;
    }
    //~ else if (option=="show" && toggle){
      //~ *live=true;
      //~ cout<<"Interactive viewer"<<endl;
    //~ }
    else if(option=="jsonfile"){
      if (value.length()==0){
        cerr<<"Empty JSON filename supplied"<<endl;
        return 1;
      }
      jsonFileName=value;
    }
    else if(option=="tgafile"){
      if (value.length()==0){
        cerr<<"Empty TGA filename supplied"<<endl;
        return 1;
      }
      tgaFileName=value;
    }
    else if(option=="pngfile"){
      if (value.length()==0){
        cerr<<"Empty PNG filename supplied"<<endl;
        return 1;
      }
      pngFileName=value;
    }
    else if(option=="width"){
      if (value.length()==0){
        cerr<<"Empty width supplied"<<endl;
        return 1;
      }
      if(EOF == sscanf(value.c_str(), "%d", width)){
        cerr<<"Invalid width '"<<value<<"'"<<endl;
        return 1;
      }

    }
    else if(option=="height"){
      if (value.length()==0){
        cerr<<"Empty height supplied"<<endl;
        return 1;
      }
      if(EOF == sscanf(value.c_str(), "%d", height)){
        cerr<<"Invalid height '"<<value<<"'"<<endl;
        return 1;
      }
    }
    else if(option=="pngaa"){
      bool error=false;
      if (value.length()==0){
        cerr<<"Empty pngaa value supplied"<<endl;
        return 1;
      }
      if(EOF == sscanf(value.c_str(), "%d", pngaa)) error=true;
      if (!error and (*pngaa<1 or *pngaa>10)) error = true;

      if (error)
      {
        cerr<<"Invalid pngaa value"<<value<<"'. Must be 1,2,3..10."<<endl;
        return 1;
      }      
    }
    else
    cerr<<"Ignoring option '"<<item<<"'"<<endl;
  }
  return 0;
}

//------------------------------------------------------------------------------
int main(int argc, char** argv){
  int width=512;
  int height=512;
  string jsonFileName("../../python/nudraw3d.json");
  string tgaFileName("");
  string pngFileName("nudraw3d.png");
  int pngAA=1;
#ifdef USE_OSMESA
  live=0;
#else
  live=1;
#endif

  GLubyte *buffer=0;
  if (live)
    glutInit(&argc, argv);
  
  if(parseArguments(argc, argv, jsonFileName, tgaFileName, &width, &height,
                                              pngFileName, &pngAA)){
    exit(1);
  }

  if (tgaFileName.length()>0 and pngAA!=1) {
    cerr<<"Ignoring option '--pngaa'"<<endl;
    pngAA=1;
  }

  Scene.parseFile(jsonFileName);

  if (live)
    initJsonViewer(width, height);
  else{
    Scene.setupCamera();
    buffer = initJsonRenderer(width*pngAA, height*pngAA);
    if (!buffer) exit(1);
  }

  Scene.buildAll();

  if (live){
    Scene.setupCamera();
    glutMainLoop();
  }
  else if (buffer){
    Scene.eyeCam[2]=-Scene.eyeCam[2]; // Coordinate system is different on vs off-screen
    Scene.upCam[1]=-Scene.upCam[1];
    displayJsonScene();
    if (tgaFileName.length()>0)
      write_targa_byte(tgaFileName.c_str(), buffer, width*pngAA, height*pngAA);
    else if (pngFileName.length()>0)
      write_png_byte(pngFileName.c_str(), buffer, width*pngAA, height*pngAA, pngAA);
  }

  return 0;
}
