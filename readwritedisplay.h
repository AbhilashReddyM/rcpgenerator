//reading writing and displaying functions

#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;


void SaveConfig2d(string savedir, long n, long dim, double sigma, double num_ratio, double rad[], double Box[], int Boundary[], double phi, double x[][2], long rad_type[], double seed);

void SaveConfig3d(string savedir, long n, long dim, double sigma, double num_ratio, double rad[], double Box[], int Boundary[], double phi, double x[][3], long rad_type[], double seed);

string getdatadir(char Directory[], int length);

string getsavedir(string datadir);


