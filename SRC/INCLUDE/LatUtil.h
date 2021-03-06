#ifndef LATUTIL_H
#define LATUTIL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define RADIUS_CNT 4

using namespace std;

// Functions
int writeDump(double*, double*, string, string, int, int);
int writeData(double*, double*, string, string, int, int, int, int);
int writeGrass();
int sphObjRead(string, double**);
int maxObj(double, double, double, double);
int atomPerObj(string);
int valRadius(double*, double);
double* boundingSphere(double**, int);
void scaleRad(double*, double**, double, int);
bool areSame(double, double);
double estimateZHeight(int*, double*, double*);

#endif
