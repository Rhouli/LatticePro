#include "INCLUDE/GrassPro.h"
#include "INCLUDE/LatUtil.h"

GrassPro::GrassPro(double *dim_in, double radius_in, string dataFName_in, 
		   string sphFName_in, double atomDensity_in, double scaleFactor_in, 
		   double cutOff_in, double sphereSep_in, int FLAG_TRANSFORM_in, 
		   int FLAG_SPH_in,int FLAG_RAND_in){
  dim = dim_in;
  radius = radius_in;
  dataFName = dataFName_in;
  sphFName = sphFName_in;
  atomDensity = atomDensity_in;
  scaleFactor = scaleFactor_in;
  cutOff = cutOff_in;
  sphereSep = sphereSep_in;
  FLAG_TRANSFORM = FLAG_TRANSFORM_in;
  FLAG_SPH = FLAG_SPH_in;
  FLAG_RAND = FLAG_RAND_in;
 
 int j;
 for(j=0;j<RADIUS_CNT;j++)
   radiusCnt[j] = 0;
  
 atomGroup = 1;
 atomType = 3;
 bondType = 1;
 angleType = 1;
 atomNum = 0;
 objNum = 0;
 objAtoms = 0;
 width = abs(dim[0])+abs(dim[1]);
 length = abs(dim[2])+abs(dim[3]);
 height = abs(dim[4])+abs(dim[5]);
}

int GrassPro::createRegion(){
  int j, k;
  double x, z, y, diameter, zCord;
  int atomStep = 0, moleculeId=1;

  // open the data file for writing
  dataFile = fopen(dataFName.c_str(), "w");

  gatherInfo();

  fprintf(stderr, "Creating Lattice, please wait...\n");

  // Print atom bond and angle info
  fprintf(dataFile,"\n%d atoms\n", atomNum);
  fprintf(dataFile,"\n%d atom types\n", atomType);

  fprintf(dataFile,"%d bonds\n", objNum*(objAtoms-1));
  fprintf(dataFile,"%d bond types\n", bondType);
  
  // Print box bounds
  fprintf(dataFile, "\n%g %g xlo xhi\n", dim[0], dim[1]);
  fprintf(dataFile, "%g %g ylo yhi\n", dim[2], dim[3]);
  fprintf(dataFile, "%g %g zlo zhi\n", dim[4], dim[5]);

  // Print Atoms
  fprintf(dataFile, "\nAtoms\n\n");

  for(x=dim[0]+radius; x<(dim[1]-radius) || areSame(x,(dim[1]-radius)); x+=(2.0*radius+sphereSep)){
    for(y=dim[2]+radius; y<(dim[3]-radius) || areSame(y, (dim[3]-radius)); y+=(2.0*radius+sphereSep)){
      if(FLAG_RAND){
	/*diameter = 2.0*radius;
	zCord = dim[4]+radius;
	while((diameter > (2.0*cutOff)) && (zCord+diameter) < dim[5]){
	  atomStep++;
	  // Print atomNumber atomType xCord yCord zCord diameter density moleculeID
	  fprintf(dataFile, "%d %d %g %g %g %g %g %d\n", atomStep, 1,
		  x, y, zCord, diameter, atomDensity, moleculeId);
	  diameter *= scaleFactor;
	  zCord += sphereSep+diameter;*/
      } else {
	diameter = 2.0*radius;
	zCord = dim[4]+radius;
	while(((diameter > (2.0*cutOff))) && ((zCord+(1.0/2.0)*diameter) < dim[5])){
	  atomStep++;
	   // Print atomNumber atomType xCord yCord zCord diameter density moleculeID
	  if(areSame(zCord, dim[4]+radius))
	    fprintf(dataFile, "%d %d %g %g %g %g %g %d\n", atomStep, 2,
		  x, y, zCord, diameter, atomDensity, moleculeId);
	  else
	    fprintf(dataFile, "%d %d %g %g %g %g %g %d\n", atomStep, 1,
		  x, y, zCord, diameter, atomDensity, moleculeId);

	  zCord += (1.0/2.0)*diameter;
	  diameter *= scaleFactor;
	  zCord += (1.0/2.0)*diameter;
	  //  fprintf(stderr, "d: %g, cutoff: %g, zcord: %g\n", diameter, cutOff, zCord);
	  
	}
      }
    }
  }
  
  // Print Bonds
  fprintf(dataFile, "\nBonds\n\n");
  for(j=0; j < objNum; j++){
    for(k=1; k < objAtoms; k++){
	fprintf(dataFile, "%d %d %d %d\n", j*objAtoms+(k-j), bondType, j*objAtoms+k, j*objAtoms+k+1); 
    }
  }
  
  
  fclose(dataFile);

  printInfo();
  printStderr();
}

int GrassPro::gatherInfo(){
  int k;

  // error if radius for the given dimensions is to small
  if(valRadius(dim, radius)){    
      fprintf(stderr, "ERROR: Radius value of %g is to large. Exiting program now...\n", radius);
      return 1; 
  }
  // Number of objects in simulation domain
  objNum = maxStalks();
  atomNum = objNum*objAtoms;

  return 0;
}

// Print to info file
void GrassPro::printInfo(){
  string temp = dataFName + "_info";
  int k;
  infoFile = fopen(temp.c_str(), "w");

  fprintf(infoFile, "Data file: %s\n", dataFName.c_str());
  fprintf(infoFile, "Atoms produced: %d\n", atomNum);
  fprintf(infoFile, "Stalks produced: %d\n", objNum);
  fprintf(infoFile, "Atom Density: %g\n", atomDensity);
  if(FLAG_SPH){
    fprintf(infoFile, "Objects produced: %d\n", objNum);
    fprintf(infoFile, "Bonds produced: %d\n", (objAtoms-1)*objNum);
  }
  fprintf(infoFile, "Radius: %g\n", radius);
  fprintf(infoFile, "Scaling Factor: %g\n", scaleFactor);
  fprintf(infoFile, "Cutoff: %g\n", cutOff);
  fprintf(infoFile, "Sphere Seperation: %g\n", sphereSep);
  fprintf(infoFile, "Box Boundarys:\n");
  for(k = 0; k < 6; k+=2)
    fprintf(infoFile, "%g %g\n", dim[k], dim[k+1]);
  if(FLAG_TRANSFORM)
    fprintf(infoFile, "Grass was randomly transformed within its bounding box\n");
  if(FLAG_RAND)
    fprintf(infoFile, "Grass was randomly re-sized by as much as %g\n", scaleFactor);
  // Close info File
  fclose(infoFile);
}

// Print to screen
void GrassPro::printStderr(){
  int k; 
  fprintf(stderr, "#######*********************#######\n");
  fprintf(stderr, "#######   LATTICE COMPLETE  #######\n");
  fprintf(stderr, "#######*********************#######\n");
  fprintf(stderr, "Data file: %s\n", dataFName.c_str());
  fprintf(stderr, "Atoms produced: %d\n", atomNum);
  fprintf(stderr, "Stalks produced: %d\n", objNum);
  fprintf(infoFile, "Atom Density: %g\n", atomDensity);
  if(FLAG_SPH){
    fprintf(stderr, "Objects produced: %d\n", objNum);
    fprintf(stderr, "Bonds produced: %d\n", (objAtoms-1)*objNum);
  }
  fprintf(stderr, "Radius: %g\n", radius);
  fprintf(stderr, "Scaling Factor: %g\n", scaleFactor);
  fprintf(stderr, "Cutoff: %g\n", cutOff);
  fprintf(stderr, "Sphere Seperation: %g\n", sphereSep);
  fprintf(stderr, "Box Boundarys:\n");
  for(k = 0; k < 6; k+=2)
    fprintf(stderr, "%g %g\n", dim[k], dim[k+1]);
  if(FLAG_TRANSFORM)
    fprintf(stderr, "Grass was randomly transformed within its bounding box\n");
  if(FLAG_RAND)
    fprintf(stderr, "Grass was randomly re-sized by as much as %g\n", scaleFactor);
}

// find the maximum grass stalks that can be produced in the given region 
int GrassPro::maxStalks(){
  
  double maxXAtoms, maxYAtoms;
  double currentHeight=0;
  
  maxXAtoms = static_cast<int>((width)/(2.0*radius+sphereSep));
  maxYAtoms = static_cast<int>((length)/(2.0*radius+sphereSep));
  
  double x;
  for(x = 2.0*radius; x > 2.0*cutOff && currentHeight < height; x*=scaleFactor){
      currentHeight += x;
      if(currentHeight < height)
	objAtoms++;
  }        

  return static_cast<int>(maxXAtoms*maxYAtoms);
}
    
