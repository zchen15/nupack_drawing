// Objects.cpp
#include "Objects.hpp"

float no_mat[] = {0.0f, 0.0f, 0.0f, 1.0f};
float mat_ambient[] = {0.7f, 0.7f, 0.7f, 1.0f};
float mat_ambient_color[] = {0.8f, 0.8f, 0.2f, 1.0f};
float mat_diffuse[] = {0.1f, 0.5f, 0.8f, 1.0f};
float mat_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
float no_shininess = 0.0f;
float low_shininess = 5.0f;
float high_shininess = 100.0f;
float mat_emission[] = {0.3f, 0.2f, 0.2f, 0.0f};

//------------------------------------------------------------------------------
void setTriplet (triplet &t, float x, float y, float z){

  if(!t.empty()) t.clear();
  t.push_back(x);
  t.push_back(y);
  t.push_back(z);
  //~ cout<<"color[0]="<<t[0]<<endl;
  //~ cout <<"size="<<t.size()<<endl;
}


//------------------------------------------------------------------------------
Sphere::Sphere(triplet in_origin, float in_radius, triplet in_color,
                int in_slices, int in_stacks, GLUquadric *in_quadric):
                slices(in_slices),
                stacks(in_stacks),
                color(in_color),
                origin(in_origin),
                radius(in_radius),
                QuadricObject(in_quadric){

}

int Sphere::Draw(){
  //~ cout<<"\nslices="<<slices<<endl;
//~
  //~ cout<<"colors=";
  //~ for(int i=0;i<color.size();i++) cout<<color[i]<<" ";
  //~ cout<<endl;
//~
  //~ cout<<"origin=";
  //~ for(int i=0;i<origin.size();i++) cout<<origin[i]<<" ";
  //~ cout<<endl;
//~
  //~ cout<<"radius="<<radius<<endl;

  glColor3f(color[0], color[1], color[2]);
  glTranslatef(origin[0], origin[1], origin[2]);
  gluSphere(quadric, radius, slices, stacks);

  return 0;
}

//------------------------------------------------------------------------------
PolyCone::PolyCone(int in_npoints, gleDouble **in_points, gleDouble *in_radii,
                    triplet in_color, int in_sides, bool owner=false):
                          npoints(in_npoints),
                          points(in_points),
                          radii(in_radii),
                          color(in_color),
                          nsides(in_sides),
                          ownsdata(owner){
}

int PolyCone::Draw(){
  //~ cout<<"\nPolyCone"<<endl;
  //~ cout<<"\nsides="<<nsides<<endl;
//~
  //~ cout<<"colors=";
  //~ for(int i=0;i<color.size();i++) cout<<color[i]<<" ";
  //~ cout<<endl;

  //~ cout<<"points="<<endl;
  //~ for(int i=0;i<npoints;i++)
    //~ cout<<((gleDouble (*)[3])points)[i][0]<<", "<<((gleDouble (*)[3])points)[i][1]<<", "<<((gleDouble (*)[3])points)[i][2]<<endl;
  //~ cout<<endl;
//~
  //~ cout<<"radius="<<endl;
  //~ for(int i=0;i<npoints;i++) cout<<radii[i]<<" ";
  //~ cout<<endl;

  glColor3f(color[0], color[1], color[2]);
  if (nsides!=20) gleSetNumSides(nsides);
  glePolyCone(npoints, (gleDouble (*)[3])points, NULL, radii);

  return 0;
}


//------------------------------------------------------------------------------
PolyCylinder::PolyCylinder(int in_npoints, gleDouble **in_points,
                                gleDouble in_radius, triplet in_color,
                                bool owner, int in_sides, GLUquadric *in_quadric):
                                points(in_points),
                                npoints(in_npoints),
                                radius(in_radius),
                                color(in_color),
                                nsides(in_sides),
                                ownsdata(owner),
                                QuadricObject(in_quadric) {
}

int PolyCylinder::Draw(){
  //~ cout<<"\nPolyCylinder"<<endl;
  //~ cout<<"\nsides="<<nsides<<endl;
//~
  //~ cout<<"colors=";
  //~ for(int i=0;i<color.size();i++) cout<<color[i]<<" ";
  //~ cout<<endl;
  //~
  //~ cout<<"points="<<endl;
  //~ for(int i=0;i<npoints;i++)
    //~ cout<<((gleDouble (*)[3])points)[i][0]<<", "<<((gleDouble (*)[3])points)[i][1]<<", "<<((gleDouble (*)[3])points)[i][2]<<endl;
  //~ cout<<endl;
//~
  //~ cout<<"radius="<<radius<<endl;

  glColor3f(color[0], color[1], color[2]);
  if (nsides!=20) gleSetNumSides(nsides);
    glePolyCylinder(npoints, (gleDouble (*)[3])points, NULL, radius);

  return 0;
}
