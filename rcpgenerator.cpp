#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include "twister.h"
#include <ctime>
#include "histogram.h"
#include "distance.h"
#include "readwritedisplay.h"
#include "minimization.h"
#include <fstream>
#include <sstream>
//#include <time>


//will do multidimensional array like idl [column][row]


using namespace std;

//define prototypes
//*********************************************************
long sub2ind(long column, long row, long dim []);
//string getdatadir(char Directory[], int length);
//string getsavedir(string datadir);
//*********************************************************


//Declare variables
 double time1, timedif;
int testvar=1;
const double pi = 3.14159;
string setup_file="GenConfig.txt";
string dirtmp = "Data     ";
string datadir;
int verbose = 0;
string savedir;
int single = 0;
//int Boundary[] = {1, 1, 1};

//Main
//*********************************************************
int main(int argc, char *argv[])
{
  
  //Read in info from command line
  //******************************************************************
  //******************************************************************  
  int c;
  static char optstring[] = "f:o:vi";
  opterr=0;
    
  while ((c=getopt(argc, argv, optstring)) != -1)
    {
      switch(c)
	{
	case 'f':
                           
	  setup_file=string(optarg);
	  break; 
	
	case 'o':
	  savedir = string(optarg); 
	  break; 

	case 'v':
	  verbose = 1;
	  break; 

	case 'i':
	  single = 1;
	  break;

	case '?':
	  cout << "Found an option that was not in optstring";
	  cout << endl;
	}
    }
  //******************************************************************
  //******************************************************************  


 
  //Read in config info from file
  //******************************************************************
  //******************************************************************
  long dim = 2; 

  ifstream setup; // indata is like cin
  int num; // variable for input value
  setup.open(setup_file.c_str()); // opens the file
  if(!setup) { // file couldn't be opened
    cerr << "Error: could not open setup file" << endl;
    exit(1);
  }

  char buffer[250];
  int a; 
  setup.getline (buffer,150); a=sscanf(buffer,"%ld", &dim); //dim

  //******************************************************************
  //******************************************************************


  if (dim == 2)
    {


      long n = 100; double phi = 0.25;
      double delta_phi = 0.005; double sigma = 1.4;
      double num_ratio = 0.5; double Box[] = {0, 1, 0, 8};
      int Boundary[] = {1, 0, 1}; int fixy = 0;
      double minr_tol = 0.9999; double dV_tol = 1E-7;
      double seed = 0;
      double params[] = {1.1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


      setup.getline (buffer,150); a=sscanf(buffer,"%ld", &n); //n
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &phi); //phi
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &delta_phi); //delta_phi
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &sigma); //sigma
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &num_ratio); //num_ratio
      setup.getline (buffer,150); 
      a=sscanf(buffer,"%lf %lf %lf %lf", &Box[0], &Box[1], &Box[2], &Box[3]); //Box
      setup.getline (buffer,150); 
      a=sscanf(buffer,"%d %d", &Boundary[0], &Boundary[1]); //Boundary
      setup.getline (buffer,150); a=sscanf(buffer,"%d", &fixy); //fixy
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &minr_tol); //minr_tol
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &dV_tol); //dV_tol
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &seed); //seed
      setup.getline (buffer,150); 
      //a=sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &params[0], &params[1], &params[2], &params[3], &params[4], &params[5], &params[6], &params[7], &params[8], &params[9]); //params


      double h = 0; double old_h = Box[3];
      double BoxArea = (Box[1] - Box[0])*(Box[3] - Box[2]);
      double small = sqrt(2*phi*BoxArea/(n*pi*(1 + (sigma*sigma))));
      double rad[] = {small*sigma, small};
      double x [n][2];
      long rad_type[n]; double minr = 2;
      double V = 0; double dV [n][2];
      init_genrand(seed);


      //Initialize initial positions of particles
      //******************************************************************
      //******************************************************************  
      double n_small = 0;
      double n_big = 0;
      if (Boundary[0] > -2)
	{
	  for(long i=0; i < n; i++)
	    {
	      x[i][0]=double(genrand_real2())*(Box[1] - Box[0]) + Box[0];
	      x[i][1]=double(genrand_real2())*(Box[3] - Box[2]) + Box[2];
	      if (i < n*num_ratio)
		{
		  rad_type[i] = 1;
		  n_small = n_small + 1;
		}
	      else
		{
		  rad_type[i] = 0;		  
		  n_big = n_big + 1;
		}
	    }
	}
      small = sqrt(phi*BoxArea/(pi*(n_small + (sigma*sigma)*n_big)));
      rad[0] = small*sigma;
      rad[1] = small;
      	  
      if (Boundary[0] <= -2)
	{
	  //usrintializepacking(n, x, rad, rad_type, Box, Boundary, params);
	}

      num_ratio = n_small/(n_small + n_big);

      //******************************************************************
      //****************************************************************** 

      //Set radius and Box based on info from config file
      //******************************************************************
      //******************************************************************  
      
      if (fixy == 0 && Boundary[0] > -1)
	{
	  h = Box[3];
	  small = 2*phi*(Box[1] - Box[0])*h/(pi*(n_small + n_big*sigma*sigma));
	  rad[0] = small*sigma;
	  rad[1] = small;  
	  Box[3] = 2*h*rad[1];
	  BoxArea = (Box[1] - Box[0])*(Box[3] - Box[2]); 
	  phi = pi*(n_big*rad[0]*rad[0] + n_small*rad[1]*rad[1])/BoxArea;
	  old_h = Box[3];

	}
      else if (Boundary[0] == -1)
	{

	  rad[1] = sqrt(phi*(Box[1]/2)*(Box[1]/2)/((n_small + (sigma*sigma)*n_big)));
	  rad[0] = sigma*rad[1];

	  int i=0;
	  while (i < n)
	    {
	      x[i][0]=double(genrand_real2())*(Box[1] - Box[0]) + Box[0];
	      x[i][1]=double(genrand_real2())*(Box[3] - Box[2]) + Box[2];
	      
	      double r_hat[] = {x[i][0] - Box[1]/2, x[i][1] - Box[1]/2};
	      double r = sqrt(r_hat[0]*r_hat[0] + r_hat[1]*r_hat[1]);
	      
	      if (r < (Box[1]/2 - rad[rad_type[i]]))
		{
		  //rad_type[i] = i % 2;
		  i = i + 1;
		}
	    }

	}
      else if (Boundary[0] <= -2)
	{
	  //usrpackingfraction(n, phi, rad, sigma, Box, Boundary, params);
	}
      reset_values2d(n, dV, x, Box);
      //******************************************************************
      //******************************************************************  
      
      
 
      
      
      //Stuff to compute time to run code
      time1 = (double) clock();            /* get initial time */
      time1 = time1 / CLOCKS_PER_SEC;      /*    in seconds    */
      
      
      
      //Expansion contraction loop
      //******************************************************************
      //******************************************************************  
      int phi_type = 1; // 1 for expansion and 0 for contraction
      int count = 0; //counts number of consecutive expansions
      long iter = 0; //total number of iterations
      long tmp_count=0; //number of iteration for an individual energy minimization
      long loop_count = 0; //number of times minimization has been performed
      while (delta_phi > 1E-6)
	{
	  
	  tmp_count=0;
	  minimizeenergyNLCG2d(V, dV, minr, x, n, rad, rad_type,
			       Box, Boundary, tmp_count, minr_tol, dV_tol);
	  iter = iter + tmp_count;
	  
	  //check if expansion or contraction needs to be performed
	  if (minr < minr_tol)
	    {
	      if (phi_type == 1){delta_phi = delta_phi/2;}
	      phi = phi - delta_phi;
	      phi_type = -1;
	      count = 0;
	      
	    }
	  else 
	    {
	      if (phi_type == -1){delta_phi = delta_phi/2;}
	      phi = phi + delta_phi;
	      phi_type = 1;
	      
	      if (phi_type = 1 & phi > .7 & delta_phi < 0.001){count = count + 1;}
	      if (count == 5){delta_phi = 2*delta_phi; count = 0;}
	      
	    }

	  //update radius
	  if (fixy == 0 && Boundary[0] > -1)	    
	    {
	      small = 2*phi*(Box[1] - Box[0])*h/(pi*(n_small + n_big*sigma*sigma));
	      rad[0] = small*sigma;
	      rad[1] = small;  
	      Box[3] = 2*h*rad[1];
	      BoxArea = (Box[1] - Box[0])*(Box[3] - Box[2]); 
	      phi = pi*(n_big*rad[0]*rad[0] + n_small*rad[1]*rad[1])/BoxArea;
	      old_h = Box[3]/old_h;
	      for(long i=0; i < n; i++){x[i][1] = x[i][1]*old_h;};
	      old_h = Box[3];
	    }
	  else
	    {
	      if (Boundary[0] > -1)
		{
		  small = sqrt(phi*BoxArea/(pi*(n_small + (sigma*sigma)*n_big)));
		  rad[0] = small*sigma;
		  rad[1] = small;
		}
	      else if (Boundary[0] == -1)
		{
		  rad[1] = sqrt(phi*(Box[1]/2)*(Box[1]/2)/((n_small + (sigma*sigma)*n_big)));
		  rad[0] = sigma*rad[1];

		}
	      else if (Boundary[0] <= -2)
		{
		  //usrpackingfraction(n, phi, rad, sigma, Box, Boundary, params);
		}
	    }
	  
	  //print info every 20 loop iterations if verbose mode is set
	  if (loop_count % 10 == 0 & verbose == 1)
	    {
	      cout << "phi = " <<  phi << endl;
	      cout << "Delta phi = " <<  delta_phi << endl;
	      cout << "total iteration count = " << iter << endl; 
	    }
	  
	  loop_count = loop_count + 1;


	  //output data at each iteration
	  if (single == 1)
	    {
	      int c = loop_count;
	      string out;
	      
	      for (int i = 0; i < 9 - floor(log10(c))-1; i++)
		{
		  out = out + "0";
		}
	      
	      stringstream ss;
	      ss << loop_count;
	      
	      out = out + ss.str();	  
	      
	      string tmps;
	      tmps = savedir + "/"  + out;
	      const char *tmpc3;
	      tmpc3 = tmps.c_str();
	      
	      ofstream outfile;
	      outfile.open (tmpc3);//,ios::app);
	      
	      outfile << rad[0] << " " << rad[1] << endl;
	      outfile << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << endl;


	      for(long i=0; i < n; i++)
		{
		  outfile << x[i][0] << " " << x[i][1] << " " << rad_type[i]+1 << endl;
		}
	      
	      outfile.close();
	    }
	  
	  //**********
	  
	}
      //******************************************************************
      //******************************************************************  
      
      //if (Boundary[0] <= -2){Box[2] = params[0];}

      //Save data to savedir
      SaveConfig2d(savedir, n, dim, sigma, num_ratio, rad, Box, Boundary, phi, x, rad_type, seed);

      if (verbose == 1)
	{
	  cout << "Saving Data to " + savedir << endl; 
	  cout << "number of iterations needed = " << iter << endl;
	  timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - time1;
	  cout << "time to run in minutes = " << timedif/60 << endl;
	}

    }
  else
    {
   
      long n = 100; double phi = 0.25;
      double delta_phi = 0.005; double sigma = 1.4;
      double num_ratio = 0.5; double Box[] = {0, 1, 0, 1, 0, 1};
      int Boundary[] = {1, 1, 1}; int fixy = 0;
      double minr_tol = 0.9999; double dV_tol = 1E-7;
      double seed = 0;


      setup.getline (buffer,150); a=sscanf(buffer,"%ld", &n); //n
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &phi); //phi
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &delta_phi); //delta_phi
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &sigma); //sigma
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &num_ratio); //num_ratio
      setup.getline (buffer,150); 
      a=sscanf(buffer,"%lf %lf %lf %lf %lf %lf", &Box[0], &Box[1], &Box[2], &Box[3], &Box[4], &Box[5]); //Box
      setup.getline (buffer,150); 
      a=sscanf(buffer,"%d %d %d", &Boundary[0], &Boundary[1], &Boundary[2]); //Boundary
      setup.getline (buffer,150); a=sscanf(buffer,"%d", &fixy); //fixy
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &minr_tol); //minr_tol
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &dV_tol); //dV_tol
      setup.getline (buffer,150); a=sscanf(buffer,"%lf", &seed); //seed

      //Set radius and Box based on info from config file
      //******************************************************************
      //******************************************************************  
      double h = 0; double old_h = Box[5];
      double BoxVol = (Box[1] - Box[0])*(Box[3] - Box[2])*(Box[3] - Box[2]);
      double small = pow(2*phi*(3/(4*pi))*BoxVol/(n*(1 + (sigma*sigma*sigma))),1/3.);
      double rad[] = {small*sigma, small};
      double x [n][3];
      long rad_type[n]; double minr = 2;
      double V = 0; double dV [n][3];
      init_genrand(seed);

      double n_small = 0;
      double n_big = 0;
      for(long i=0; i < n; i++)
	{
	  x[i][0]=double(genrand_real2())*(Box[1] - Box[0]) + Box[0];
	  x[i][1]=double(genrand_real2())*(Box[3] - Box[2]) + Box[2];
	  x[i][2]=double(genrand_real2())*(Box[5] - Box[4]) + Box[4];
	  if (i < n*num_ratio)
	    {
	      rad_type[i] = 1;
	      n_small = n_small + 1;
	    }
	  else
	    {
		  rad_type[i] = 0;		  
		  n_big = n_big + 1;
	    }
	}

      small = pow(phi*3./4.*BoxVol/(pi*(n_small + (sigma*sigma)*n_big)),1./3.);
      rad[0] = small*sigma;
      rad[1] = small;

      //Need to update for 3D case
      if (fixy == 0)
	{
	  h = Box[5];
	  small = sqrt(3*phi*(Box[1] - Box[0])*(Box[3] - Box[2])*h/(2*pi*(n_small + n_big*(sigma*sigma*sigma))));
	  rad[0] = small*sigma;
	  rad[1] = small;  
	  Box[5] = 2*h*rad[1];
	  BoxVol = (Box[1] - Box[0])*(Box[3] - Box[2])*(Box[5] - Box[4]);
	  phi = 4*pi/3*(n_big*rad[0]*rad[0]*rad[0] + n_small*rad[1]*rad[1]*rad[1])/BoxVol;
	  old_h = Box[5];
	}
      //******************************************************************
      //******************************************************************  
      

      //Initialize initial positions of particles
      //******************************************************************
      //******************************************************************  


      //need to create reset_values3d
      //reset_values2d(n, dV, x, Box);
      
      //******************************************************************
      //******************************************************************  
      
      
      //Stuff to compute time to run code
      time1 = (double) clock();            /* get initial time */
      time1 = time1 / CLOCKS_PER_SEC;      /*    in seconds    */


      //Expansion contraction loop
      //******************************************************************
      //******************************************************************  
      int phi_type = 1; // 1 for expansion and 0 for contraction
      int count = 0; //counts number of consecutive expansions
      long iter = 0; //total number of iterations
      long tmp_count=0; //number of iteration for an individual energy minimization
      long loop_count = 0; //number of times minimization has been performed
      while (delta_phi > 1E-6)
	{

	  tmp_count=0;
	  minimizeenergyNLCG3d(V, dV, minr, x, n, rad, rad_type,
			       Box, Boundary, tmp_count, minr_tol, dV_tol);
	  iter = iter + tmp_count;

	  //check if expansion or contraction needs to be performed
	  if (minr < minr_tol)
	    {
	      if (phi_type == 1){delta_phi = delta_phi/2;}
	      phi = phi - delta_phi;
	      phi_type = -1;
	      count = 0;
	      
	    }
	  else 
	    {
	      if (phi_type == -1){delta_phi = delta_phi/2;}
	      phi = phi + delta_phi;
	      phi_type = 1;
	      
	      if (phi_type = 1 & phi > .7 & delta_phi < 0.001){count = count + 1;}
	      if (count == 5){delta_phi = 2*delta_phi; count = 0;}
	      
	    }

	  //update radius
	  if (fixy == 0)
	    {
	      small = sqrt(3*phi*(Box[1] - Box[0])*(Box[3] - Box[2])*h/(2*pi*(n_small + n_big*(sigma*sigma*sigma))));
	      rad[0] = small*sigma;
	      rad[1] = small;  
	      Box[5] = 2*h*rad[1];
	      BoxVol = (Box[1] - Box[0])*(Box[3] - Box[2])*(Box[5] - Box[4]);
	      phi = 4*pi/3*(n_big*rad[0]*rad[0]*rad[0] + n_small*rad[1]*rad[1]*rad[1])/BoxVol;
	      old_h = Box[5];
	    }
	  else
	    {
	      small = pow(2*phi*(3/(4*pi))*BoxVol/(n*(1 + (sigma*sigma*sigma))),1/3.);
	      rad[0] = small*sigma;
	      rad[1] = small;
	    }
	  
	  //print info every 20 loop iterations if verbose mode is set
	  if (loop_count % 5 == 0 & verbose == 1)
	    {
	      cout << "phi = " <<  phi << endl;
	      cout << "Delta phi = " <<  delta_phi << endl;
	      cout << "total iteration count = " << iter << endl; 
	    }
	  
	  loop_count = loop_count + 1;

	  //if (phi > 0.61){delta_phi = 1E-7;}

	}
      //******************************************************************
      //******************************************************************  

	  //output data at each iteration
	  if (single == 1)
	    {
	      int c = loop_count;
	      string out;
	      
	      for (int i = 0; i < 9 - floor(log10(c))-1; i++)
		{
		  out = out + "0";
		}
	      
	      stringstream ss;
	      ss << loop_count;
	      
	      out = out + ss.str();	  
	      
	      string tmps;
	      tmps = savedir + "/"  + out;
	      const char *tmpc3;
	      tmpc3 = tmps.c_str();
	      
	      ofstream outfile;
	      outfile.open (tmpc3);//,ios::app);
	      
	      outfile << rad[0] << " " << rad[1] << endl;
	      outfile << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << " " << Box[4] << " " << Box[5] << endl;


	      for(long i=0; i < n; i++)
		{
		  outfile << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << rad_type[i]+1 << endl;
		}
	      
	      outfile.close();
	    }
	  
	  //**********


      //Save data to savedir
      SaveConfig3d(savedir, n, dim, sigma, num_ratio, rad, Box, Boundary, phi, x, rad_type, seed);

      if (verbose == 1)
	{
	  cout << "Saving Data to " + savedir << endl; 
	  cout << "number of iterations needed = " << iter << endl;
	  timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - time1;
	  cout << "time to run in minutes = " << timedif/60 << endl;
	}

    }
      
  return 0;
}









//Get Cell
//*********************************************************
long getcell2d(long xnbins, long ynbins, long cell, long& xbin, long& ybin)
{
  xbin = cell % xnbins;
  ybin = (cell - xbin)/xnbins;
}

long kmod(long a, long b)
{

  while ( a < 0 ){a = a +b;}

  return a % b;

}


