#include "INCLUDE/GranPro.h"

GranPro::GranPro(double *dim_in, double *radius_in, string dataFName_in, 
		 string sphFName_in, double atomDensity_in, int FLAG_TRANSFORM_in, 
		 int FLAG_SPH_in, int FLAG_VAR_in, int FLAG_WIGGLE_in){
  dim = dim_in;
  radius = radius_in;
  dataFName = dataFName_in;
  sphFName = sphFName_in;
  FLAG_TRANSFORM = FLAG_TRANSFORM_in;
  FLAG_SPH = FLAG_SPH_in;
  FLAG_VAR = FLAG_VAR_in;
  FLAG_WIGGLE = FLAG_WIGGLE_in;
 
 int j;
  for(j=0;j<RADIUS_CNT;j++)
    radiusCnt[j] = 0;
  
  atomGroup = 1;
  atomType = 2;
  bondType =1;
  angleType = 1;

  atomDensity = atomDensity_in;
  // open the data file for writing
  dataFile = fopen(dataFName.c_str(), "w");
}

GranPro::~GranPro(){
  fclose(dataFile);
}

int GranPro::createRegion(){
  int j, k;
  double x, z, y;
  int atomStep = 0, moleculeId=1;

  gatherInfo();

  fprintf(stderr, "Creating Lattice, please wait...\n");

  // Print atom bond and angle info
  fprintf(dataFile,"\n%d atoms\n", atomNum);
  fprintf(dataFile,"\n%d atom types\n", atomType);

  if(FLAG_SPH){
    fprintf(dataFile,"%d bonds\n", objNum*(objAtoms-1));
    fprintf(dataFile,"%d bond types\n", bondType);
  }
  // Print box bounds
  fprintf(dataFile, "\n%g %g xlo xhi\n", dim[0], dim[1]);
  fprintf(dataFile, "%g %g ylo yhi\n", dim[2], dim[3]);
  fprintf(dataFile, "%g %g zlo zhi\n", dim[4], dim[5]);

  // Print Atoms
  fprintf(dataFile, "\nAtoms\n\n");

  // Setup for Matrix FLAG_TRANSFORM

  Mat3x3 M;  
  double point[4];
  double rad = radius[0], random; 
  double off_x, off_y, off_z;
  srand(time(0));

  if(FLAG_TRANSFORM){
    M.Zeros();
    M(1,1) = 99;
  }
  for(z=dim[4]+radius[0]; z<(dim[5]-radius[0]) || areSame(z, (dim[5]-radius[0])); z+=(2.0*radius[0])){
    for(y=dim[2]+radius[0]; y<(dim[3]-radius[0]) || areSame(y, (dim[3]-radius[0])); y+=(2.0*radius[0])){
      for(x=dim[0]+radius[0]; x<(dim[1]-radius[0]) || areSame(x,(dim[1]-radius[0])); x+=(2.0*radius[0])){
	if(FLAG_SPH){
	  if(FLAG_TRANSFORM){
	    // Object translation
	    ColMatrix A(4);
	    A.Zeros();
	    A(1)=(double) rand();
	    A(2)=(double) rand();
	    A(3)=(double) rand();
	    A(4)=(double) rand();
	    
	    //EP_FLAG_TRANSFORMation(ColMatrix& q, Mat3x3& C);
	    EP_Transformation(A, M);
	  }
	  for(k=0; k < objAtoms; k++){
	    if(FLAG_TRANSFORM){
	      // Add transformeded points to point
	      point[0] = M(1,1)*(sphObj[k][0]-boundSphere[0])
		+ M(1,2)*(sphObj[k][1]-boundSphere[1])
		+ M(1,3)*(sphObj[k][2]-boundSphere[2]);
	      point[1] = M(2,1)*(sphObj[k][0]-boundSphere[0])
		+ M(2,2)*(sphObj[k][1]-boundSphere[1])
		+ M(2,3)*(sphObj[k][2]-boundSphere[2]);
	      point[2] = M(3,1)*(sphObj[k][0]-boundSphere[0])
		+ M(3,2)*(sphObj[k][1]-boundSphere[1])
		+ M(3,3)*(sphObj[k][2]-boundSphere[2]);
	      point[0] += x;
	      point[1] += y;
	      point[2] += z;
	      point[3] = sphObj[k][3];
	    } else {
	      // Add normal points to point
	      point[0] = x+sphObj[k][0]-boundSphere[0];
	      point[1] = y+sphObj[k][1]-boundSphere[1];
	      point[2] = z+sphObj[k][2]-boundSphere[2];
	      point[3] = sphObj[k][3];
	    }
	    atomStep++;
	    fprintf(dataFile, "%d %d %g %g %g %g %g %d\n", atomStep, 1,
		    point[0], point[1], point[2], 2.0*point[3], 
		    atomDensity, moleculeId);
	  }
	  moleculeId++;
	} else{
	  if(FLAG_VAR){
	    off_x = 0.0;
	    off_y = 0.0;
	    off_z = 0.0;
	    random = (double) rand()/RAND_MAX;
	      if(random < (0.15)){
		rad = radius[0];
		radiusCnt[0]++;
	      } else if((random > (0.15)) && (random <(0.25))){
		rad = radius[1];
		radiusCnt[1]++;
		if(FLAG_WIGGLE)
		  wiggle(&off_x, &off_y,  &off_z, radius[0], radius[1]);
	      } else if((random > 0.25) && (random < 0.5)){
		rad = radius[2];
		radiusCnt[2]++;
		if(FLAG_WIGGLE)
		  wiggle(&off_x, &off_y, &off_z, radius[0], radius[2]);
	      } else{
		rad = radius[3];
		radiusCnt[3]++;
		if(FLAG_WIGGLE)
		  wiggle(&off_x, &off_y, &off_z, radius[0], radius[3]);
	      }
	  }
	  atomStep++;
	  if(areSame(y+off_y, 0.0) && areSame(x+off_x, 0.0) && areSame(z+off_z, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, 0.0, 0.0);
	  else if(areSame(y+off_y, 0.0) && areSame(x+off_x,0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, 0.0, z+off_z);
	  else if(areSame(y+off_y, 0.0) && areSame(z+off_z,0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, x+off_x, 0.0, 0.0);
	  else if(areSame(x+off_x, 0.0) && areSame(z+off_z, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, y+off_y, 0.0);
	  else if(areSame(x+off_x, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, y+off_y, z+off_z);
	  else if(areSame(y+off_y, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, x+off_x, 0.0, z+off_z);
	  else if(areSame(z+off_z, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, x+off_x, y+off_y, 0.0);
	  else
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, x+off_x, y+off_y, z+off_z);
	}
      }
    }
  }
  if(FLAG_SPH){
    // Print Bonds
    fprintf(dataFile, "\nBonds\n\n");
    for(j=0; j < objNum; j++){
      for(k=1; k < objAtoms; k++){
	fprintf(dataFile, "%d %d %d %d\n", j*objAtoms+(k-j), 1, j*objAtoms+k, j*objAtoms+k+1); 
      }
    }
  }
  printInfo();
  printStderr();
}

int GranPro::gatherInfo(){
  int k;

  if(FLAG_SPH){
    // Find the number of atoms per object
    objAtoms = atomPerObj(sphFName);
    if(objAtoms < 0){
      fprintf(stderr, "ERROR: Cannot open SPH file, exiting program.\n");
      return 0;
    }
    // Read in the object individual atoms coordinates
    sphObj = (double**) malloc(objAtoms*sizeof(double*));
    for(k = 0; k < objAtoms; k++)
      sphObj[k] = (double*) malloc(4*sizeof(double));
    
    if(sphObjRead(sphFName, sphObj)){
      fprintf(stderr, "ERROR: Cannot open SPH file, exiting program.\n");
      return 0;
    }
    // Find the best fit bounding sphere to the object
    boundSphere = boundingSphere(sphObj, objAtoms);
    scaleRad(boundSphere, sphObj, radius[0], objAtoms);
    radius[0] = boundSphere[3];
  }
  // error if radius for the given dimensions is to small
  if(valRadius(dim, radius[0])){    
      fprintf(stderr, "ERROR: Radius value of %g is to large. Exiting program now...\n", radius[0]);
      return 1; 
  }
  if(FLAG_SPH){
    // Number of objects in simulation domain
    objNum = maxObj(abs(dim[0])+abs(dim[1]), abs(dim[2])+abs(dim[3]), abs(dim[4])+abs(dim[5]), radius[0]);
    atomNum = objNum*objAtoms;
    radiusCnt[0]=objNum;

  }else{
    // Number of objects in simulation domain
    atomNum = maxObj(abs(dim[0])+abs(dim[1]), abs(dim[2])+abs(dim[3]), abs(dim[4])+abs(dim[5]), radius[0]);
  }
  return 0;
}

void GranPro::printInfo(){
  string temp = dataFName + "_info";
  int k;
  infoFile = fopen(temp.c_str(), "w");
  // Print to info file
  fprintf(infoFile, "Data file: %s\n", dataFName.c_str());
  fprintf(infoFile, "Atoms produced: %d\n", atomNum);
  fprintf(infoFile, "Atom Density: %g\n", atomDensity);
  if(FLAG_SPH){
    fprintf(infoFile, "Objects produced: %d\n", objNum);
    fprintf(infoFile, "Bonds produced: %d\n", (objAtoms-1)*objNum);
  }
  if(FLAG_VAR)
    fprintf(infoFile, "Radii: %d of size %g\n       %d of size %g\n       %d of size %g\n       %d of size %g\n", 
	    radiusCnt[0], radius[0], radiusCnt[1], radius[1], radiusCnt[2], radius[2], radiusCnt[3], radius[3]);
  else
    fprintf(infoFile, "Radius: %g\n", radius[0]);
  fprintf(infoFile, "Box Boundarys:\n");
  for(k = 0; k < 6; k+=2)
    fprintf(infoFile, "%g %g\n", dim[k], dim[k+1]);
  if(FLAG_WIGGLE)
    fprintf(infoFile, "Granules translated inside bounding sphere\n");
  if(FLAG_TRANSFORM)
    fprintf(infoFile, "Granules were randomly transformed\n");
  fprintf(infoFile, "Estimated z height at zero porosity: %g m\n", estimateZHeight(radiusCnt, radius, dim));
// Close info File
  fclose(infoFile);
}

// Print to screen
void GranPro::printStderr(){
  int k; 
  fprintf(stderr, "#######*********************#######\n");
  fprintf(stderr, "#######   LATTICE COMPLETE  #######\n");
  fprintf(stderr, "#######*********************#######\n");
  fprintf(stderr, "Data file: %s\n", dataFName.c_str());
  fprintf(stderr, "Atoms produced: %d\n", atomNum);
  fprintf(infoFile, "Atom Density: %g\n", atomDensity);
  if(FLAG_SPH){
    fprintf(stderr, "Objects produced: %d\n", objNum);
    fprintf(stderr, "Bonds produced: %d\n", (objAtoms-1)*objNum);
  }
  if(FLAG_VAR)
    fprintf(stderr, "Radii: %d of size %g\n       %d of size %g\n       %d of size %g\n       %d of size %g\n", 
	    radiusCnt[0], radius[0], radiusCnt[1], radius[1], radiusCnt[2], radius[2], radiusCnt[3], radius[3]);
  else
    fprintf(stderr, "Radius: %g\n", radius[0]);
  fprintf(stderr, "Box Boundarys:\n");
  for(k = 0; k < 6; k+=2)
    fprintf(stderr, "%g %g\n", dim[k], dim[k+1]);
  if(FLAG_WIGGLE)
    fprintf(stderr, "Granules translated inside bounding sphere\n");
  if(FLAG_TRANSFORM)
    fprintf(stderr, "Granules were randomly transformed\n");
  fprintf(stderr, "Estimated z height at zero porosity: %g m\n", estimateZHeight(radiusCnt, radius, dim));
}

void GranPro::wiggle(double* off_x, double* off_y, double* off_z, double rad1, double rad2){
  double x=0, y=0, z=0;
  x = (rad1 - rad2) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
  do{
    y = (rad1 - rad2) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);\
  }while(sqrt(pow(x, 2)+pow(y,2)) > (rad1-rad2));
  do{    
    z = (rad1 - rad2) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
  }while(sqrt(pow(z, 2)+pow(y,2))> (rad1-rad2) || sqrt(pow(x, 2)+pow(z,2))> (rad1-rad2));
  
  *off_x = x;
  *off_y = y;
  *off_z = z;
}
