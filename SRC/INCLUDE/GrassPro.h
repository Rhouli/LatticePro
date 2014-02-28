#ifndef GrassPro_H
#define GrassPro_H

#include "LatUtil.h"

class GrassPro{
 private:
  double* dim;
  double radius;
  double atomDensity;
  double scaleFactor;
  double cutOff;
  double sphereSep;
  double width;
  double length;
  double height;

  string dataFName;
  string sphFName;

  int atomNum;
  int objNum;
  int objAtoms;
  int atomGroup;
  int atomType;
  int bondType;
  int angleType;
  int radiusCnt[RADIUS_CNT];

  FILE *dataFile, *infoFile;

  // Flags
  int FLAG_TRANSFORM;
  int FLAG_SPH;
  int FLAG_RAND;

 public: 
  GrassPro(){};
  GrassPro(double*, double, string, string,
	   double, double, double, double,
	   int, int, int);
  ~GrassPro(){};
  
  // functions
  int createRegion();
  int gatherInfo();
  void printInfo();
  void printStderr();
  int maxStalks();
  int atomPerStalk();
};

#endif
