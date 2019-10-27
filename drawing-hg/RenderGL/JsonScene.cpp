#include "JsonScene.hpp"

//------------------------------------------------------------------------------
JsonScene::JsonScene(){
  quadric = gluNewQuadric();
  gluQuadricNormals(quadric, GLU_SMOOTH);         // Create Smooth Normals
  gluQuadricTexture(quadric, GL_TRUE);            // Create Texture Coords
  gluQuadricOrientation(quadric, GLU_OUTSIDE);
}

bool JsonScene::parseFile(string fname){
  fileStream.open(fname.c_str(), fstream::in);
  bool success=jsonReader.parse(fileStream, jsonRoot, true);
  if ( !success ){
    // report to the user the failure and their locations in the document.
    cerr  << "Failed to parse scene file "<<fname<<"\n"
          << jsonReader.getFormatedErrorMessages();
  }

  return success;
}

//------------------------------------------------------------------------------
void setTripletJsonArray(triplet &t, Json::Value v){
  setTriplet(t, (gleDouble)v[(unsigned int)0].asDouble(),
                (gleDouble)v[(unsigned int)1].asDouble(),
                (gleDouble)v[(unsigned int)2].asDouble());
}

//------------------------------------------------------------------------------
void setPointsJsonArrayN3(gleDouble points[][3], Json::Value v, int npoints, int Offset=0){
  for (int i=0; i<npoints;i++){
    Json::Value row=v[i+Offset];
    points[i][0]=row[(unsigned int)0].asDouble();
    points[i][1]=row[(unsigned int)1].asDouble();
    points[i][2]=row[(unsigned int)2].asDouble();
  }
}

//------------------------------------------------------------------------------
void setPointsJsonArrayN(gleDouble points[], Json::Value v, int npoints){

  for (int i=0; i<npoints;i++){
    points[i]=v[(unsigned int)i].asDouble();
  }
}


//------------------------------------------------------------------------------
void JsonScene::buildSpheres(){
  const Json::Value spheres = jsonRoot["objects"]["spheres"];
  const Json::Value textures = jsonRoot["textures"];
  int numSpheres=spheres.size();
  int i,texind;
  Sphere *sp;
  triplet color;
  triplet origin;

  for (i=0; i<numSpheres;i++){
    texind=spheres[i]["texture"].asInt();
    setTripletJsonArray(color, textures[texind]["color"]);
    setTripletJsonArray(origin, spheres[i]["center"]);

    sp = new Sphere(origin, spheres[i]["radius"].asDouble(), color, 40, 40, quadric);
    objects.push_back(sp);
  }
}

//------------------------------------------------------------------------------
void JsonScene::buildCones(){
  const Json::Value cones = jsonRoot["objects"]["cones"];
  const Json::Value textures = jsonRoot["textures"];
  int numCones=cones.size();
  int numPoints;
  int i,texind;
  PolyCone *pc;
  triplet color;
  gleDouble **points;
  gleDouble *radii;

  for (i=0; i<numCones;i++){
    numPoints=cones[i]["points"].size();
    Json::Value dataPoints=cones[i]["points"];
    Json::Value dataRadii=cones[i]["radii"];

    // Get color
    texind=cones[i]["texture"].asInt();
    setTripletJsonArray(color, textures[texind]["color"]);

    // Get points
    points=(gleDouble **) new gleDouble[numPoints][3];
    setPointsJsonArrayN3((gleDouble (*)[3])points, dataPoints, numPoints);

    //Get radii
    radii=(gleDouble *) new gleDouble[numPoints];
    setPointsJsonArrayN((gleDouble *)radii, dataRadii, numPoints);

    pc=new PolyCone(numPoints, points, radii, color, 160, true);

    objects.push_back(pc);
  }
}

//------------------------------------------------------------------------------
void JsonScene::buildPolys(){
  const Json::Value polys = jsonRoot["objects"]["polycylinders"];
  const Json::Value textures = jsonRoot["textures"];
  int numPolys=polys.size();
  int numPoints;
  int i,texind;
  PolyCylinder *pc;
  triplet color;
  gleDouble **points;

  for (i=0; i<numPolys;i++){
    numPoints=polys[i]["points"].size();
    Json::Value dataPoints=polys[i]["points"];

    texind=polys[i]["texture"].asInt();
    if (texind>=0){
      setTripletJsonArray(color, textures[texind]["color"]);
      points=(gleDouble **) new gleDouble[numPoints][3];
      setPointsJsonArrayN3((gleDouble (*)[3])points, dataPoints, numPoints);

      pc=new PolyCylinder(numPoints, points,
                      polys[i]["radius"].asDouble(), color, true, 160, quadric);

      objects.push_back(pc);
    }
    else{
      Json::Value texInds=polys[i]["textures"];
      int numTextures=texInds.size();
      cout << "numTextures="<<numTextures<<endl;
      cout << "numPoints="<<numPoints<<endl;
      unsigned int j=1, segStart=0, // j = segment End
        prevInd=texInds[(unsigned int)0].asInt(), curInd;
      while (j<numPoints){
        for (;j<numPoints;j++){
          curInd=texInds[j].asInt();
          if (prevInd != curInd) // new segment
            break;
          prevInd=curInd;
        }
        unsigned int realSegstart=segStart, realNumPoints=j-segStart;
        if (realSegstart>1) {
          realSegstart--;
          realNumPoints++;
        }
        if (j<numPoints-1)
          realNumPoints+=2;
        cout << "realSegstart="<<realSegstart<<", j="<<j<<", realNumPoints="<<realNumPoints<<endl;
        setTripletJsonArray(color, textures[prevInd]["color"]);
        points=(gleDouble **) new gleDouble[realNumPoints][3];
        setPointsJsonArrayN3((gleDouble (*)[3])points, dataPoints, realNumPoints, realSegstart);

        pc=new PolyCylinder(realNumPoints, points,
                        polys[i]["radius"].asDouble(), color, true, 160, quadric);

        objects.push_back(pc);
        prevInd=curInd;
        segStart=j;
      }
    }

  }
}

//------------------------------------------------------------------------------
int JsonScene::setupCamera(){
  const Json::Value camera = jsonRoot["camera"];
  const Json::Value pos=camera["position"];
  const Json::Value lookat=camera["lookat"];
  const Json::Value updir=camera["updir"];

  setTripletJsonArray(eyeCam, pos);
  setTripletJsonArray(centerCam, lookat);
  setTripletJsonArray(upCam, updir);
}

//------------------------------------------------------------------------------
int JsonScene::buildAll(){
  buildSpheres();
  buildCones();
  buildPolys();
}
