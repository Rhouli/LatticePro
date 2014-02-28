#include "INCLUDE/LatUtil.h"

// read the coordinates of the objects atoms into a double**
int sphObjRead(string sphFName, double** sphObj){
  int x=0; 
  string atomInput;
  ifstream sphFile(sphFName.c_str());

  if(sphFile.good()){
    while(getline(sphFile, atomInput)){
      istringstream iss(atomInput);
      string s;
      int y = 0;
      while(getline(iss, s, ' ')){
	if(y != 4){
	if(s[0] == '-' || isdigit(s[0]))
	  sphObj[x][y] = atof(s.c_str());
	}
	y++;
      }
      x++;
    }
    
    sphFile.close();
    return 0;
  }
  sphFile.close();
  return -1;
}

// find the maximum spheres that can be produced in the given region 
int maxObj(double width, double length, double height, double radius){
  double maxXAtoms, maxYAtoms, maxZAtoms;
  
  maxXAtoms = static_cast<int>((width)/(2.0*radius));
  maxYAtoms = static_cast<int>((length)/(2.0*radius));
  maxZAtoms = static_cast<int>((height)/(2.0*radius));

  return static_cast<int>(maxXAtoms*maxYAtoms*maxZAtoms);
}

// find the number of atoms per object
int atomPerObj(string sphFName){
  int atomNum = 0; 
  string temp;
  ifstream sphFile(sphFName.c_str());

  if(sphFile.good()){
    while(getline(sphFile,temp)){
      if(temp[0] == '-' || isdigit(temp[0]))
	atomNum++;
    }
    
    sphFile.close();

    return atomNum;
  }
  sphFile.close();
  return -1;
}

// confirm spheres radius will fit in defined region
int valRadius(double *dim, double radius){
  if(((abs(dim[0])+abs(dim[1]) < 2*radius) || ((abs(dim[2])+abs(dim[3])) < 2*radius) || (abs(dim[4])+abs(dim[5])) < 2*radius))
    return 1;
  return 0;
}

// find best fitting sphere of the inputed object 
double* boundingSphere(double **sphObj, int atomNum){
  int i;
  double* boundsph = (double*) malloc(sizeof(double)*4);
  double maxX=0, maxY=0, maxZ=0, maxR=0;
  
  // find center of the atoms
  for(i=0; i < atomNum; i++){
    boundsph[0]=sphObj[i][0];
    boundsph[1]=sphObj[i][1];
    boundsph[2]=sphObj[i][2];
  }
  boundsph[0]/=atomNum;
  boundsph[1]/=atomNum;
  boundsph[2]/=atomNum;
 
  // find max X, Y, and Z distance from center
  for(i=0; i < atomNum; i++){
    if((abs(sphObj[i][0])-abs(boundsph[0])) > maxX)
      maxX = abs(sphObj[i][0]);
    if((abs(sphObj[i][1])-abs(boundsph[1])) > maxY)
      maxY = abs(sphObj[i][1]);
    if((abs(sphObj[i][2])-abs(boundsph[2])) > maxZ)
      maxZ = abs(sphObj[i][2]);
    if(abs(sphObj[i][3]) > maxR)
      maxR = abs(sphObj[i][3]);
  }
  
  // find max of the max X, Y, and Z
  if(maxX > maxY && maxX > maxZ)
    boundsph[3]= maxX + maxR;
  else if(maxY > maxZ)
    boundsph[3]=maxY + maxR;
  else
    boundsph[3]= maxZ + maxR;

  return boundsph;
}

// Scale bounding sphere and sphObj so radius is as requested by user
void scaleRad(double* boundSphere, double** objSph, double radius, int objAtoms){
  double scaleFactor = radius/boundSphere[3];
  int x, y;
  for(x = 0; x < 4; x++)
    boundSphere[x] *= scaleFactor;

  for(x = 0; x < objAtoms; x++){
    for(y = 0; y < 4; y++)
      objSph[x][y] *= scaleFactor;
  }
}

bool areSame(double a, double b)
{
  double EPSILON = 0.00000001;
  return (fabs(a - b) < EPSILON);
}

double estimateZHeight(int* radiusCnt, double* radius, double* dim){
  double granVol=0.0, pi = atan(1)*4;
  
  int x;
  for(x=0; x < RADIUS_CNT; x++)
    granVol += (4/3)*pi*(((double)radiusCnt[x])*pow(radius[x], 3));
  
  return granVol/((abs(dim[0])+abs(dim[1]))*(abs(dim[2])+abs(dim[3])));
}

