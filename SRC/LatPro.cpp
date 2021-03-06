/* Author: Ryan Houlihan
 *
 * Lammps Lattice dump file Producer
 *
 * Usage: LatPro [-FLAG] [ARGUMENT]
 * FLAG: 
 *       [-s fileName]                                 :: read in an sph file to fit to the given lattice 
 *       [-a xMin xMax yMin yMax zMin zMax]            :: read in dimensions 
 *       [-o1 fileName]                                :: output Gran Region 
 *       [-o2 fileName]                                :: output Grass Region
 *       [-r radiusVal]                                :: radius value 
 *       [-v radiHigh radiMid1 radiMid2 radiLow]       :: variable sized radius
 *       [-w] 				               :: wiggle
 *       [-t]                                          :: transform vectors 
 *       [-h]                                          :: hybrid style
 *       [-d density]                                  :: atom density
 *       [-g scalefactor cutoff sphereSep]             :: scalefactor (fraction to reduce diameter by)
 *                                                        cutoff (smallest radius allowed)
 *                                                        sphereSep (sphere seperation distance in meters)
 *       [-R]                                          :: Randomly varry grass size by scaleFactor
 *
 *
 * Grass Example: latPro -a -.05 .05 -.05 .05 0 .05 -r .005 -g .75 .0025 .00001 -d 100 -o2 GrassTest
 * Gran Example:  latPro -s rock.txt -a -.08 .08 -.08 .08 0 .40 -r .0025 -d 2060 -w -t -o1 RigidGranTest
 */

#include "INCLUDE/LatPro.h" 

int main(int argv, char* argc[]){
  double dim[6], radius[RADIUS_CNT], atomDensity = 0, scaleFactor = 0, cutOff = 0, sphereSep = 0;
  string fname, sphName;
  char tmpstr[30];
  int x, y, badRadius = 0;
  int FLAG_SPH = 0, FLAG_DIM = 0, FLAG_OUT = 0, FLAG_RAD = 0, 
    FLAG_GRAN = 0, FLAG_GRASS = 0, FLAG_TRANSFORM = 0, FLAG_VAR = 0, 
    FLAG_WIGGLE = 0, FLAG_RAND = 0;

  // take in user inputs
  for(x = 0; x < argv; x++){
    if(argc[x][0] == (char)'-'){
      // read in sph file
      if(argc[x][1] == (char)'s'){
	FLAG_SPH = 1;
	x++;
	sphName = argc[x];
	// read in dimensions
      } else if(argc[x][1] == (char)'a'){
	FLAG_DIM = 1;
	int z = 0;
	for(z = 0; z < 6; z++){
	  x++;
	  dim[z] = atof(argc[x]);
	}
      } else if(argc[x][1] == (char)'o'){
	  FLAG_OUT = 1;
	  if(argc[x][2] == (char)'1'){
	    x++;
	    FLAG_GRAN = 1;
	    string tmp = argc[x];
	    fname = tmp + ".data";
	  }else if(argc[x][2] == (char)'2'){
	    x++;
	    FLAG_GRASS = 1;
	    string tmp = argc[x];
	    fname = tmp + ".data";
	  }
      } else if(argc[x][1] == (char)'r'){
	FLAG_RAD = 1;
	x++;
	radius[0] = atof(argc[x]);
	radius[1] = '\0';
      } else if(argc[x][1] == (char)'t'){
	FLAG_TRANSFORM = 1;
      } else if(argc[x][1] == (char)'v'){
	FLAG_VAR = 1;
	int j;
	for(j = 0; j < RADIUS_CNT; j++){
	  x++;	
	  radius[j] = atof(argc[x]);
	}	
      }
      else if(argc[x][1] == (char)'w'){
	FLAG_WIGGLE = 1;
      }
      else if(argc[x][1] == (char)'d'){
	x++;
	atomDensity = atof(argc[x]);
      }
      else if(argc[x][1] == (char)'g'){
	x++;
	scaleFactor = atof(argc[x]);
	x++;
	cutOff = atof(argc[x]);
	x++;
	sphereSep = atof(argc[x]);
      }
      else if(argc[x][1] == (char)'R'){
	FLAG_RAND = 1;
      }
    }
  }
  
  fprintf(stderr, "#######************************************************#######\n");
  fprintf(stderr, "#######      LAMMPS Lattice .dump & data. Producer     #######\n");
  fprintf(stderr, "#######************************************************#######\n");

  // Create grass region or granular region depending on set flags
    if(FLAG_GRASS){
      GrassPro* grassRegion;
      grassRegion = new GrassPro(dim, radius[0], fname, sphName, atomDensity, 
				 scaleFactor, cutOff, sphereSep, 
				 FLAG_TRANSFORM, FLAG_SPH, FLAG_RAND);
      grassRegion->createRegion();
    } else if(FLAG_GRAN){
      GranPro* granRegion;
      granRegion = new GranPro(dim, radius, fname, sphName, atomDensity, FLAG_TRANSFORM, 
			       FLAG_SPH, FLAG_VAR, FLAG_WIGGLE);
      granRegion->createRegion();
     } 

  return 0;
}

// Write to dump file based on users inputs and given sph file
int  writeDump(double *dim, double *radius, string dumpFName, 
	       string sphFName, int FLAG_SPH, int FLAG_VAR){
  FILE *dumpFile; 
  int j, k;
  int atomNum, objectNum, objAtoms, atomGroup = 1, atomType = 1, 
    atomStep1 = 1, atomStep2 = 1, radiusCnt[RADIUS_CNT];
  for(j=0;j<RADIUS_CNT;j++)
    radiusCnt[j] = 0;

  double x, z, y, *boundSphere, **sphObj;
   
  fprintf(stderr, "Creating Lattice, please wait...\n");

  if(FLAG_SPH){
    // Find the number of atoms per object    
    objAtoms = atomPerObj(sphFName);
    if(objAtoms < 0){
      fprintf(stderr, "ERROR: Cannot open SPH file, exiting program.\n");
      return 0;
    }

    // Read in the objects individual atoms coordinates
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
  }

  // error if radius for the given dimensions is to small
  if(valRadius(dim, radius[0])){
    fprintf(stderr, "ERROR: Radius value of %g is to large. Exiting program now...\n", boundSphere[3]);
    return 1; 
  }

  if(FLAG_SPH){
    // Number of objects in simulation domain
    objectNum = maxObj(abs(dim[0])+abs(dim[1]), abs(dim[2])+abs(dim[3]), abs(dim[4])+abs(dim[5]), boundSphere[3]);
    atomNum = objectNum*objAtoms;

  }else {
    // Number of objects in simulation domain
    atomNum = maxObj(abs(dim[0])+abs(dim[1]), abs(dim[2])+abs(dim[3]), abs(dim[4])+abs(dim[5]), radius[0]);
  }

  // open the dump file for writing
  dumpFile = fopen(dumpFName.c_str(), "w");

  // Print timestep info
  fprintf(dumpFile,"ITEM: TIMESTEP\n");
  fprintf(dumpFile,"0\n");

  // Print Num atoms
  fprintf(dumpFile,"ITEM: NUMBER OF ATOMS\n");
  fprintf(dumpFile,"%d\n", atomNum);
  
  // Print box bounds
  fprintf(dumpFile, "ITEM: BOX BOUNDS\n");
  for(j = 0; j < 6; j+=2)
    fprintf(dumpFile, "%g %g\n", dim[j], dim[j+1]);

  // Print Atoms
  fprintf(dumpFile, "ITEM: ATOMS id type type x y z radius\n");

  double rad = radius[0], random;
  srand(time(0));
  for(z=dim[4]+radius[0]; z<(dim[5]-radius[0]) || AreSame(z, (dim[5]-radius[0])); z+=(2.0*radius[0])){
    for(y=dim[2]+radius[0]; y<(dim[3]-radius[0]) || AreSame(y, (dim[3]-radius[0])); y+=(2.0*radius[0])){
      for(x=dim[0]+radius[0]; x<(dim[1]-radius[0] ) || AreSame(x,(dim[1]-radius[0])); x+=(2.0*radius[0])){  
	if(FLAG_SPH){
	  for(k=0; k < objAtoms; k++){
	    fprintf(dumpFile, "%d %d %d %g %g %g %g\n", atomStep1, atomGroup, atomType, 
		    x+(sphObj[k][0]-boundSphere[0]), y+(sphObj[k][1]-boundSphere[1]), 
		    z+(sphObj[k][2]-boundSphere[2]), sphObj[k][3]);
    	    atomStep1++;
	  }
	} else{
	  if(FLAG_VAR){
	    random = (double) rand()/RAND_MAX;
	      if(random < (1.0/5.0)){
		rad = radius[2];
		radiusCnt[0]++;
	      }
	      else if((random > (1.0/5.0)) && (random <(1.0/3.0))){
		rad = radius[1];
		radiusCnt[1]++;
	      }
	      else{
		rad = radius[0];
		radiusCnt[2]++;
	      }
	  }
	      fprintf(dumpFile, "%d %d %d %g %g %g %g\n", atomStep2, atomGroup, 
		      atomType, x, y, z, rad);
	}
	atomStep2++;
      }
    }
  }

  // Close dump file
  fclose(dumpFile);

  // Replay info to user
  fprintf(stderr, "#######*********************#######\n");
  fprintf(stderr, "#######   LATTICE COMPLETE  #######\n");
  fprintf(stderr, "#######*********************#######\n");
  fprintf(stderr, "Dump file: %s\n", dumpFName.c_str());
  if(FLAG_SPH)
    fprintf(stderr, "Objects produced: %d\n", objectNum);
  fprintf(stderr, "Atoms produced: %d\n", atomNum);
  if(FLAG_VAR)
    fprintf(stderr, "Radii: %d of size %g\n       %d of size %g\n       %d of size %g\n", 
	    radiusCnt[0], radius[0], radiusCnt[1], radius[1], radiusCnt[2], radius[2]);
  else
    fprintf(stderr, "Radius: %g\n", radius[0]);
  fprintf(stderr, "Box Boundarys:\n");
  for(k = 0; k < 6; k+=2)
    fprintf(stderr, "%g %g\n", dim[k], dim[k+1]);

  return 0;
}


// Write to data file based on users inputs and given sph file
int  writeData(double *dim, double *radius, string dataFName, 
	       string sphFName, int FLAG_TRANSFORM, int FLAG_SPH, int FLAG_VAR, int FLAG_WIGGLE){

  FILE *dataFile, *infoFile; 
  int j, k;
  int atomNum, objNum, objAtoms, atomGroup = 1, atomType = 2,
    bondType = 1, angleType = 0, atomStep = 0, moleculeId=1, radiusCnt[RADIUS_CNT];
  for(j=0;j<RADIUS_CNT;j++){
    radiusCnt[j] = 0;
  }
  double x, z, y, *boundSphere, **sphObj, atomDensity;
 
  fprintf(stderr, "Please enter atom density: ");
  scanf("%lf", &atomDensity);
  
  fprintf(stderr, "Creating Lattice, please wait...\n");

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
  }else{
    // Number of objects in simulation domain
    atomNum = maxObj(abs(dim[0])+abs(dim[1]), abs(dim[2])+abs(dim[3]), abs(dim[4])+abs(dim[5]), radius[0]);
  }
  // open the dump file for writing
  dataFile = fopen(dataFName.c_str(), "w");
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
  for(z=dim[4]+radius[0]; z<(dim[5]-radius[0]) || AreSame(z, (dim[5]-radius[0])); z+=(2.0*radius[0])){
    for(y=dim[2]+radius[0]; y<(dim[3]-radius[0]) || AreSame(y, (dim[3]-radius[0])); y+=(2.0*radius[0])){
      for(x=dim[0]+radius[0]; x<(dim[1]-radius[0]) || AreSame(x,(dim[1]-radius[0])); x+=(2.0*radius[0])){
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
		if(FLAG_WIGGLE){
		  off_x = (radius[0] - radius[1]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		  while(sqrt(pow(off_x, 2)+pow(off_y,2)) > (radius[0]-radius[1]))
		    off_y = (radius[0] - radius[1]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		  while(sqrt(pow(off_z, 2)+pow(off_y,2)) > (radius[0]-radius[1]) || sqrt(pow(off_x, 2)+pow(off_z,2))> (radius[0]-radius[1]))
		  off_z = (radius[0] - radius[1]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		}
	      } else if((random > 0.25) && (random < 0.5)){
		rad = radius[2];
		radiusCnt[2]++;
		if(FLAG_WIGGLE){
		  off_x = (radius[0] - radius[2]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		  while(sqrt(pow(off_x, 2)+pow(off_y,2)) > (radius[0]-radius[2]))
		    off_y = (radius[0] - radius[2]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		  while(sqrt(pow(off_z, 2)+pow(off_y,2))> (radius[0]-radius[2]) || sqrt(pow(off_x, 2)+pow(off_z,2))> (radius[0]-radius[2]))
		  off_z = (radius[0] - radius[2]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		}
	      } else{
		rad = radius[3];
		radiusCnt[3]++;
		if(FLAG_WIGGLE){
		  off_x = (radius[0] - radius[3]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		  while(sqrt(pow(off_x, 2)+pow(off_y,2)) > (radius[0]-radius[3]))
		    off_y = (radius[0] - radius[3]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		  while(sqrt(pow(off_z, 2)+pow(off_y,2))> (radius[0]-radius[3]) || sqrt(pow(off_x, 2)+pow(off_z,2))> (radius[0]-radius[3]))
		  off_z = (radius[0] - radius[3]) * (double) (rand() - RAND_MAX/2)/(RAND_MAX);
		}
	      }
	  }
	  atomStep++;
	  if(AreSame(y+off_y, 0.0) && AreSame(x+off_x, 0.0) && AreSame(z+off_z, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, 0.0, 0.0);
	  else if(AreSame(y+off_y, 0.0) && AreSame(x+off_x,0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, 0.0, z+off_z);
	  else if(AreSame(y+off_y, 0.0) && AreSame(z+off_z,0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, x+off_x, 0.0, 0.0);
	  else if(AreSame(x+off_x, 0.0) && AreSame(z+off_z, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, y+off_y, 0.0);
	  else if(AreSame(x+off_x, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, 0.0, y+off_y, z+off_z);
	  else if(AreSame(y+off_y, 0.0))
	    fprintf(dataFile, "%d %d %g %g %g %g %g\n", atomStep, 
		    1, 2.0*rad , atomDensity, x+off_x, 0.0, z+off_z);
	  else if(AreSame(z+off_z, 0.0))
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

  // Close data file
  fclose(dataFile);

  string temp = dataFName + "_info";
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

  // Print to screen
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
    fprintf(stderr, "Granules translated inside bouding sphere\n");
  if(FLAG_TRANSFORM)
    fprintf(stderr, "Granules were randomly transformed\n");
  fprintf(stderr, "Estimated z height at zero porosity: %g m\n", estimateZHeight(radiusCnt, radius, dim));

  return 0;
}

