#include <iostream>
#include <GL/glut.h>
#include <GL/gle.h>
#include <stdlib.h>
#include <vector>
#include "Objects.hpp"
#include "RenderGL.hpp"
#include "JsonScene.hpp"

using namespace std;
extern vector<Object *> objects;


int DemoScene(){
  // Draw a Sphere

  GLUquadric* quadric = gluNewQuadric();
  gluQuadricNormals(quadric, GLU_SMOOTH);         // Create Smooth Normals (NEW)
  gluQuadricTexture(quadric, GL_TRUE);            // Create Texture Coords (NEW)
  gluQuadricOrientation(quadric, GLU_OUTSIDE);
  
  triplet cOrigin, cColor;
  setTriplet(cOrigin, 0.0, 0.0, 0.0);
  setTriplet(cColor, 0, 0.0, 1.0);
  Sphere *sp=new Sphere(cOrigin, 5.0, cColor, 40, 40, quadric);
  //~ sp->setQuadric(quadric);
  //~ objects.push_back(sp);

  // Draw a PolyCone
  gleDouble coneData[][3]={{0,0,11}, {0,0,12}, {0,0,16}, {0,0,17}};
  gleDouble coneRadii[]={0.0, 2.0, 0.0, 0.0};
  
  PolyCone *cone=new PolyCone(4, (gleDouble **)coneData, (gleDouble *)coneRadii,
                              cColor, 40, false);
  objects.push_back(cone);  

  // Draw a PolyCylinder
  gleDouble cylData[][3]={{0,0,-1}, {0,0,0}, {0,0,10}, {0,0,12}, {0,0,13}};
  PolyCylinder *pc=new PolyCylinder(5, (gleDouble **)cylData, 2.0, cColor,
                                    false, 40, quadric);
  //~ sp->setQuadric(quadric);
  objects.push_back(pc);


  setTriplet(cOrigin, 5.0, 0.0, 0.0);
  setTriplet(cColor, 1.0, 0.0, 0);
  Sphere *sp2=new Sphere(cOrigin, 5.0, cColor, 40, 40, quadric);
  //~ sp->setQuadric(quadric);
  //~ cout<<sp->color[0]<<endl;
  objects.push_back(sp2);

  gleDouble cylData2[][3]={{0,7,-1}, {0,7,5}, {0,7,10}, {0,7,12}, {0,7,13}};
  PolyCylinder *pc2=new PolyCylinder(5, (gleDouble **)cylData2, 2.0, cColor,
                                      false, 40, quadric);
  //~ sp->setQuadric(quadric);
  objects.push_back(pc2);

  int npoints3=10;
  gleDouble **cylData3=(gleDouble **)new gleDouble[npoints3][3];

  gleDouble x=0.0, y=0.0, z=0.0;
  for (int i=0;i<npoints3;i++){
    z=i*4;
    y=z*2;
    ((gleDouble (*)[3])cylData3)[i][0]=x;
    ((gleDouble (*)[3])cylData3)[i][1]=y;
    ((gleDouble (*)[3])cylData3)[i][2]=z;
  }
  setTriplet(cColor, 0.0, 1.0, 0);
  PolyCylinder *pc3=new PolyCylinder(npoints3, (gleDouble **)cylData3, 2.0,
                                      cColor, true, 40, quadric);
  objects.push_back(pc3);

  glutMainLoop();
  delete sp;
  delete sp2;
  delete cone;
  delete pc;
  delete pc2;
  delete pc3;
  free(quadric);
}

int main(int argc,char** argv){
  InitViewer(argc, argv);
  DemoScene();
  return 0;
}
