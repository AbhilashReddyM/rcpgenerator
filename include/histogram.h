// histogram header file

// #include "histogram.h"
// #include "distance.h"
// #include "readwritedisplay.h"
// #include "minimization.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <ctime>
#include <fstream>
#include <sstream>
using namespace std;

long sub2ind(long column, long row, long dim[]);
void khist(double arr[], long length, double xmin, long nbins, double binsize,
					 long bins[], long dim[]);
void khist2d(double arr[][2], long length, double xmin, long xnbins, double xbinsize,
						 double ymin, long ynbins, double ybinsize,
						 long bins[], long dim[], int Boundary[]);
void khist3d(double arr[][3], long length, double xmin, long xnbins, double xbinsize,
						 double ymin, long ynbins, double ybinsize,
						 double zmin, long znbins, double zbinsize,
						 long bins[], long dim[]);
