#include "histogram.h"

// sub2ind
//*********************************************************
long sub2ind(long column, long row, long dim[])
{
  long out = dim[1] * column + row;
  return out;
}
//*********************************************************

// 1D histogram
//*********************************************************
void khist(double arr[], long length, double xmin, long nbins, double binsize,
           long bins[], long dim[])
{
  long tmp;
  for (long i = 0; i < length; i++)
  {
    tmp = floor((arr[i] - xmin) / binsize);
    bins[sub2ind(0, tmp, dim)] = bins[sub2ind(0, tmp, dim)] + 1;
    bins[sub2ind(bins[sub2ind(0, tmp, dim)], tmp, dim)] = i;
  }
}
//*********************************************************

// 2D histogram
//*********************************************************
void khist2d(double arr[][2], long length, double xmin, long xnbins, double xbinsize,
             double ymin, long ynbins, double ybinsize,
             long bins[], long dim[], int Boundary[])
{
  long xbin;
  long ybin;
  long actualbin;

  for (long i = 0; i < length; i++)
  {

    xbin = floor((arr[i][0] - xmin) / xbinsize);
    ybin = floor((arr[i][1] - ymin) / ybinsize);
    actualbin = ybin * xnbins + xbin;
    bins[sub2ind(0, actualbin, dim)] = bins[sub2ind(0, actualbin, dim)] + 1;
    bins[sub2ind(bins[sub2ind(0, actualbin, dim)], actualbin, dim)] = i;
  }
}
//*********************************************************

// 3D histogram
//*********************************************************
void khist3d(double arr[][3], long length, double xmin, long xnbins, double xbinsize,
             double ymin, long ynbins, double ybinsize,
             double zmin, long znbins, double zbinsize,
             long bins[], long dim[])
{
  long xbin;
  long ybin;
  long zbin;
  long actualbin;

  for (long i = 0; i < length; i++)
  {
    xbin = floor((arr[i][0] - xmin) / xbinsize);
    ybin = floor((arr[i][1] - ymin) / ybinsize);
    zbin = floor((arr[i][2] - zmin) / zbinsize);
    actualbin = zbin * xnbins * ynbins + ybin * xnbins + xbin;
    bins[sub2ind(0, actualbin, dim)] = bins[sub2ind(0, actualbin, dim)] + 1;
    bins[sub2ind(bins[sub2ind(0, actualbin, dim)], actualbin, dim)] = i;
  }
}
//*********************************************************
