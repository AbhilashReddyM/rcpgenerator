#include "readwritedisplay.h"

void SaveConfig2d(string savedir, long n, long dim, double sigma, double num_ratio, double rad[], double Box[], int Boundary[], double phi, double x[][2], long rad_type[], double seed)
{

  // date
  // time
  // n
  // dim
  // sigma
  // number ratio
  // rad
  // Box
  // Boundary
  // phi

  string tmps = "mkdir " + savedir;
  const char *tmpc;
  tmpc = tmps.c_str();
  system(tmpc);
  
  time_t date; // Make a time_t object that'll hold the date
  time(&date); //  Set the date variable to the current date
  //cout << "The date and Time: " << ctime(&date) << endl; // Display the date in human-readable format

  ofstream outfile;
  tmps = savedir + "/system";
  const char *tmpc2;
  tmpc2 = tmps.c_str();
  //cout << tmpc2 << endl;
  outfile.open (tmpc2);//,ios::app);

  //date and time
  stringstream ss;
  ss << ctime(&date);

  char buffer[256];
  ss.getline (buffer,100);
  outfile << buffer[8] << buffer[9] << "-" 
       << buffer[4] << buffer[5] << buffer[6] << "-"
       << buffer[20] << buffer[21] << buffer[22] << buffer[23] << endl;

  outfile << buffer[11] << buffer[12] << buffer[13] << buffer[14] 
       << buffer[15] << buffer[16] << buffer[17] << buffer[18] << endl;

  //sys info
  outfile << n << endl;
  outfile << dim << endl;
  outfile << sigma << endl;
  outfile << num_ratio << endl;
  outfile << rad[0] << endl;
  outfile << rad[1] << endl;
  outfile << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3];

  if (dim == 3){outfile << Box[4] << " " << Box[5] << " " << endl;}
  else{outfile << endl;}

  if (Boundary[0] == 0){outfile << "C ";}
  else{if (Boundary[0] == 1) {outfile << "P ";}
    else{outfile << "R ";} }

  if (Boundary[1] == 0){outfile << "C ";}
  else{outfile << "P ";}

  if (dim == 3){
      if (Boundary[1] == 0){outfile << "C ";}
      else{outfile << "P ";}
    }

  outfile << endl;
  outfile << phi << endl;
  outfile << seed << endl;

  outfile.close();


  //save coordinates
  //ofstream outfile;
  tmps = savedir + "/"  + "FinalConfig";
  const char *tmpc3;
  tmpc3 = tmps.c_str();

  outfile.open (tmpc3);//,ios::app);

  for(long i=0; i < n; i++)
    {
      outfile << x[i][0] << " " << x[i][1] << " " << rad_type[i]+1 << endl;
    }

  outfile.close();

}




void SaveConfig3d(string savedir, long n, long dim, double sigma, double num_ratio, double rad[], double Box[], int Boundary[], double phi, double x[][3], long rad_type[], double seed)
{

  // date
  // time
  // n
  // dim
  // sigma
  // number ratio
  // rad
  // Box
  // Boundary
  // phi

  string tmps = "mkdir " + savedir;
  const char *tmpc;
  tmpc = tmps.c_str();
  system(tmpc);
  
  time_t date; // Make a time_t object that'll hold the date
  time(&date); //  Set the date variable to the current date
  //cout << "The date and Time: " << ctime(&date) << endl; // Display the date in human-readable format

  ofstream outfile;
  tmps = savedir + "/system";
  const char *tmpc2;
  tmpc2 = tmps.c_str();
  //cout << tmpc2 << endl;
  outfile.open (tmpc2);//,ios::app);

  //date and time
  stringstream ss;
  ss << ctime(&date);

  char buffer[256];
  ss.getline (buffer,100);
  outfile << buffer[8] << buffer[9] << "-" 
       << buffer[4] << buffer[5] << buffer[6] << "-"
       << buffer[20] << buffer[21] << buffer[22] << buffer[23] << endl;

  outfile << buffer[11] << buffer[12] << buffer[13] << buffer[14] 
       << buffer[15] << buffer[16] << buffer[17] << buffer[18] << endl;

  //sys info
  outfile << n << endl;
  outfile << dim << endl;
  outfile << sigma << endl;
  outfile << num_ratio << endl;
  outfile << rad[0] << endl;
  outfile << rad[1] << endl;
  outfile << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3];

  if (dim == 3){outfile << " " << Box[4] << " " << Box[5] << " " << endl;}
  else{outfile << endl;}

  if (Boundary[0] == 0){outfile << "C ";}
  else{outfile << "P ";}

  if (Boundary[1] == 0){outfile << "C ";}
  else{outfile << "P ";}

  if (dim == 3){
      if (Boundary[1] == 0){outfile << "C ";}
      else{outfile << "P ";}
    }

  outfile << endl;
  outfile << phi << endl;
  outfile << seed << endl;

  outfile.close();


  //save coordinates
  //ofstream outfile;
  tmps = savedir + "/"  + "FinalConfig";
  const char *tmpc3;
  tmpc3 = tmps.c_str();

  outfile.open (tmpc3);//,ios::app);

  for(long i=0; i < n; i++)
    {
      outfile << x[i][0] << " " << x[i][1] << " " << x[i][2] << " " << rad_type[i]+1 << endl;
    }

  outfile.close();

}


string getdatadir(char Directory[], int length)
{

  int str_count = 0;
  for (int i=0; i < length; i++)
    {
      char tmp1[] = {Directory[i],' '};
      char tmp2[] = {Directory[length-1],' '};

      if (Directory[i] != ' ')
	{
	  str_count = str_count+1;
	}  
    }

  string Dir;
  for (int i=0; i < str_count; i++)
    {
      if (Directory[i] != ' ')
	{
	  string s(1, Directory[i]);
	  //strcat(Dir, s);
	  Dir = Dir + s;
	}  
    } 
  
  return Dir;

}


string getsavedir(string datadir)
{


  string getdir = "ls ";
  getdir = getdir + datadir + "/ > tmpfile";
  const char *tmp;
  tmp = getdir.c_str();
  
  //char test1[] = "ls ";
  //char test2[] = "Data       ";
  //char test3[] = " ";
  //strcat( test1, " " );
  //strcat( test1, Directory);
  //strcat(test1, " > tmpfile");
  //cout << test1 << endl;
  system(tmp);


  char buffer[256];
  ifstream indata; // indata is like cin
  int num; // variable for input value
  indata.open("tmpfile"); // opens the file


  //  if(!indata) { // file couldn't be opened
  //  cerr << "Error: file could not be opened" << endl;
  //  exit(1);
  //}

  int line_count = 0;
  int a;
  int b = 0;
  while ( !indata.eof() ) { // keep reading until end-of-file
    indata.getline (buffer,100);
    if ( !indata.eof())
      {
	int j=0;
	char tmparr[256];
	int tmpval;
	int flag = 0;
	for (int i=0; i < 9; i++ )
	  {
	    char testvar[] = {buffer[i]};
	    a=sscanf(testvar,"%d", &tmpval);
	    if (tmpval != 0 | flag != 0)
	      {
		//cout << "here" << endl;
		tmparr[j]=buffer[i];
		j = j + 1;
		flag = 1;
	      }
	    }
	a=sscanf(tmparr,"%d", &b);
      }
  }
  //cout << b << endl;
  //cout << floor(log10(b))+1 << endl; 
  indata.close();

  //cout << b << endl;

  int c = b;
  if (b == 0){c = 1;}

  string out;

  //cout << 9 - floor(log10(c+1))-1 << endl;;

  for (int i = 0; i < 9 - floor(log10(c+1))-1; i++)
    {
      out = out + "0";
    }


  stringstream ss;
  ss << b+1;

  out = out + ss.str();
  system("rm tmpfile");

  return out;

}
