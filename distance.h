//#include "histogram.h"

#include "histogram.h"
//#include "readwritedisplay.h"
//#include "minimization.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <ctime>
#include <fstream>
#include <sstream>
using namespace std;

void reset_values2d(long n, double dV[][2], double x[][2], double Box[]);

void getdV2d(double V, double dV[][2], double& minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[]);

long getneighboringcell2d(long xnbins, long ynbins, 
			  long cell, long cellboundarytype[],
			  int Boundary[], long neighboringcells[]);

long getpaircount2d(long cell, long neighboringcells[], long bins[], long bin_dims[]);

void getpairslist2d(long pairs[][3], long paircount, long cellboundarytype[],
		    long cell, long neighboringcells[], long bins[], long bin_dims[]);

void getdist2d(double dist[][2], double x[][2], long pairs [][3], long paircount, double Box []);

void getdV2d_sub(double V, double dV[][2], double& minr, double dist[][2], long pairs[][3], long paircount, double rad[], long rad_type[]);



void reset_values3d(long n, double dV[][3], double x[][3], double Box[]);

void getdV3d(double V, double dV[][3], double& minr, double x[][3], long n, double rad[], long rad_type[], double Box[], int Boundary[]);
long getneighboringcell3d(long xnbins, long ynbins, long znbins, 
			  long cell, long cellboundarytype[],
			  int Boundary[], long neighboringcells[]);
void getdist3d(double dist[][3], double x[][3], long pairs [][3], long paircount, double Box []);
void getdV3d_sub(double V, double dV[][3], double& minr, double dist[][3], long pairs[][3], long paircount, double rad[], long rad_type[]);

long modbin(long a, long b);
void updateneighboringcell3d(long& xbint, long& ybint, long& zbint, 
			   long& xbin_on_boundary, long& ybin_on_boundary, long& zbin_on_boundary,
			   long xnbins, long ynbins, long znbins, int Boundary[],
			   long neighboringcells[], long cellboundarytype[], int& count);
