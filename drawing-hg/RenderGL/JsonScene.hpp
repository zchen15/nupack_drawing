// Json file parser/scene builder
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <GL/glut.h>
#include <GL/gle.h>
#include "Objects.hpp"
#include "json/json.h"
using namespace std;

class JsonScene{
protected:
  Json::Reader jsonReader;
  fstream fileStream;
  Json::Value jsonRoot;
  
  void buildSpheres();
  void buildCones();
  void buildPolys();

  GLUquadric *quadric;

public:
  vector<Object *> objects;
  triplet eyeCam, centerCam, upCam;

  JsonScene();
  bool parseFile(string fname);
  int setupCamera();
  int buildAll();
};

