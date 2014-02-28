#ifndef GranPro_H
#define GranPro_H

#include "eulerparameters.h"
#include "mat3x3.h"
#include "colmatrix.h"
#include "LatUtil.h"

class GranPro{
 private:
  double *dim;
  double *radius;
  double **sphObj;
  double *boundSphere;
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
  double atomDensity;
  FILE *dataFile, *infoFile;

  // Flags
  int FLAG_TRANSFORM;
  int FLAG_SPH;
  int FLAG_VAR;
  int FLAG_WIGGLE;

 public: 
  GranPro(){};
  GranPro(double*, double*, string, string, double, int, int, int, int);
  ~GranPro();
  
  // functions
  int createRegion();
  int gatherInfo();
  void printInfo();
  void printStderr();
  void wiggle(double*, double*, double*, double, double);
};

#endif
