//minimization

//#include "histogram.h"
#include "distance.h"
//#include "readwritedisplay.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <ctime>
#include <fstream>
#include <sstream>
using namespace std;

void minimizeenergyNLCG2d(double V, double dV[][2], double& minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[], long& k, double minr_tol, double dV_tol);

void minimizeenergySD2d(double V, double dV[][2], double& minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[], long& k);

double normdV2d(double dV[][2], long n);

void vectorize(double x[][2], double x_vec[], double dV[][2], double dV_vec[], long n);

void unvectorize(double x[][2], double x_vec[], double dV[][2], double dV_vec[], long n);

void reset_pos2d(long n, double x[][2], double Box[]);

double dot2d(long n, double v1[], double v2[]);

double get_alpha2d(double dV2[][2], double q[], double V, double dV[][2], double& minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[]);



void minimizeenergyNLCG3d(double V, double dV[][3], double& minr, double x[][3], long n, double rad[], long rad_type[], double Box[], int Boundary[], long& k, double minr_tol, double dV_tol);

double normdV3d(double dV[][3], long n);

void vectorize3d(double x[][3], double x_vec[], double dV[][3], double dV_vec[], long n);

void unvectorize3d(double x[][3], double x_vec[], double dV[][3], double dV_vec[], long n);

void reset_pos3d(long n, double x[][3], double Box[]);

double dot3d(long n, double v1[], double v2[]);

double get_alpha3d(double dV2[][3], double q[], double V, double dV[][3], double& minr, double x[][3], long n, double rad[], long rad_type[], double Box[], int Boundary[]);


