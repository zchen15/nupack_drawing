// Object class definitions

#include <vector>
#include <iostream>
#include <GL/glut.h>
#include <GL/gle.h>
#include <stdlib.h>
using namespace std;

#ifndef OBJECTS_HPP
#define OBJECTS_HPP

//------------------------------------------------------------------------------
typedef vector<float> triplet;
typedef vector<GLdouble> tripletDouble;
void  setTriplet (triplet &t, float x, float y, float z);
#define SPHERE 1
#define POLYCYLINDER 2

//------------------------------------------------------------------------------
class Object{
public:
  Object(){
  }

  virtual int Draw(){
    //~ cout<<"Base class Draw() called"<<endl;
  }

  ~Object(){
    //~ cout<<"~Object Destructor"<<endl;
  }
};

//------------------------------------------------------------------------------
class QuadricObject:public Object{
protected:
  GLUquadric* quadric;

public:
  QuadricObject(GLUquadric* in_quadric){
    quadric=in_quadric;
    //~ cout<<"QuadricObject constructor called"<<endl;
  }

  void setQuadric(GLUquadric* in_quadric){
    quadric=in_quadric;
  }

  virtual ~QuadricObject(){
    cout<<"~QuadricObject Destructor"<<endl;
  }
};

//------------------------------------------------------------------------------
class Sphere:public QuadricObject{
protected:
  triplet origin;
  triplet color;
  float radius;
  int slices;
  int stacks;

public:
  Sphere(triplet in_origin, float in_radius, triplet in_color,
          int in_slices, int in_stacks, GLUquadric* quadric);

  int Draw();

  virtual ~Sphere(){
    //~ cout<<"~Sphere Destructor"<<endl;
  }
};

//------------------------------------------------------------------------------
class PolyCone:public Object{
protected:
  triplet color;
  int nsides;
  int npoints;
  gleDouble **points;
  gleDouble *radii;
  bool ownsdata;

public:
  PolyCone(int in_npoints, gleDouble **in_points, gleDouble *in_radii,
    triplet in_color, int in_sides, bool owner);

  int Draw();

  virtual ~PolyCone(){
    if (ownsdata){
      delete[] points;
      delete[] radii;
    }
  }
};

//------------------------------------------------------------------------------
class PolyCylinder:public QuadricObject{
protected:
  triplet color;
  int nsides;
  int npoints;
  gleDouble **points;
  gleDouble radius;
  bool ownsdata;

public:
  PolyCylinder(int in_npoints, gleDouble **in_points, gleDouble in_radius,
    triplet in_color, bool owner, int sides, GLUquadric* quadric);
  int Draw();

  virtual ~PolyCylinder(){
    if (ownsdata){
      //~ cout<<"Deleting points..."<<endl;
      delete[] points;
    }
  }
};

#endif
