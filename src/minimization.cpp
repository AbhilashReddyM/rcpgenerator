#include "minimization.h"

// Energy Minimization via Nonlinear Conjugate Gradient
void minimizeenergyNLCG2d(double V, double dV[][2], double &minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[], long &k, double minr_tol, double dV_tol)
{

  double x_vec[2 * n];
  double dV_vec[2 * n];
  double dV2[n][2];
  double dV2_vec[2 * n];
  double alpha;
  double beta;
  double q[2 * n];
  double s[2 * n];
  long reset_count = 0;

  // Evaluate Function
  reset_values2d(n, dV, x, Box);
  getdV2d(V, dV, minr, x, n, rad, rad_type, Box, Boundary);
  vectorize(x, x_vec, dV, dV_vec, n);

  // Initialize search direction q = -dV
  for (long i = 1; i < 2 * n; i++)
  {
    q[i] = -dV_vec[i];
  }

  k = 0;
  minr = 0;

  if (normdV2d(dV, n) / n / 2 < dV_tol)
  {
    minr = 1;
  }

  double dV2_vecminusdV_vec[2 * n];

  while (minr<minr_tol & normdV2d(dV, n) / n / 2> dV_tol)
  {

    // Find alpha and return dV2
    alpha = get_alpha2d(dV2, q, V, dV, minr, x, n, rad, rad_type, Box, Boundary);
    vectorize(x, x_vec, dV2, dV2_vec, n);

    // Update position
    for (long i = 0; i < 2 * n; i++)
    {
      x_vec[i] = x_vec[i] + alpha * q[i];
    }
    unvectorize(x, x_vec, dV, dV_vec, n);

    for (long i = 0; i < 2 * n; i++)
    {
      dV2_vecminusdV_vec[i] = dV2_vec[i] - dV_vec[i];
    }

    beta = dot2d(2 * n, dV2_vec, dV2_vecminusdV_vec) / dot2d(2 * n, dV_vec, dV_vec);
    if (beta < 0)
    {
      beta = 0;
    }

    for (long i = 0; i < 2 * n; i++)
    {
      q[i] = -dV2_vec[i] + beta * q[i];
    }

    if (minr<0.53 | reset_count> 25 |
        dot2d(2 * n, dV2_vec, dV_vec) / dot2d(2 * n, dV_vec, dV_vec) > 2.5 | beta == 0)
    {
      for (long i = 0; i < 2 * n; i++)
      {
        q[i] = -dV2_vec[i];
      }
      reset_count = 0;
    }

    // Update info
    reset_values2d(n, dV, x, Box);
    getdV2d(V, dV, minr, x, n, rad, rad_type, Box, Boundary);
    vectorize(x, x_vec, dV, dV_vec, n);
    reset_count = reset_count + 1;
    k = k + 1;
  }

  reset_values2d(n, dV, x, Box);
}

// Energy Minimization via steepest descent
void minimizeenergySD2d(double V, double dV[][2], double &minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[], long &k)
{

  k = 0;
  minr = 0;
  dV[0][0] = 100;

  while (minr<.9999 & normdV2d(dV, n) / n / 2> 0.00001 & k < 25000)
  {
    reset_values2d(n, dV, x, Box);
    getdV2d(V, dV, minr, x, n, rad, rad_type, Box, Boundary);

    for (int i = 0; i < n; i++)
    {
      x[i][0] = x[i][0] - dV[i][0] / 5000;
      x[i][1] = x[i][1] - dV[i][1] / 5000;
    }

    k = k + 1;
  }
  reset_values2d(n, dV, x, Box);
}

double normdV2d(double dV[][2], long n)
{
  double out = 0;
  for (int i = 0; i < n; i++)
  {
    out = out + sqrt(dV[i][0] * dV[i][0] + dV[i][1] * dV[i][1]);
  }

  return out;
}

// Vectorize variables
void vectorize(double x[][2], double x_vec[], double dV[][2], double dV_vec[], long n)
{
  double out = 0;
  long count = 0;
  for (int i = 0; i < n; i++)
  {
    x_vec[count] = x[i][0];
    x_vec[count + 1] = x[i][1];
    dV_vec[count] = dV[i][0];
    dV_vec[count + 1] = dV[i][1];
    count = count + 2;
  }
}

void unvectorize(double x[][2], double x_vec[], double dV[][2], double dV_vec[], long n)
{
  double out = 0;
  long count = 0;
  for (int i = 0; i < n; i++)
  {
    x[i][0] = x_vec[count];
    x[i][1] = x_vec[count + 1];
    dV[i][0] = dV_vec[count];
    dV[i][1] = dV_vec[count + 1];
    count = count + 2;
  }
}

void reset_pos2d(long n, double x[][2], double Box[])
{
  for (long i = 0; i < n; i++)
  {
    while (x[i][0] < Box[0])
      x[i][0] += (Box[1] - Box[0]);
    while (x[i][1] < Box[2])
      x[i][1] += (Box[3] - Box[2]);
    while (x[i][0] > Box[1])
      x[i][0] -= (Box[1] - Box[0]);
    while (x[i][1] > Box[3])
      x[i][1] -= (Box[3] - Box[2]);
  }
}

double dot2d(long n, double v1[], double v2[])
{
  double out = 0;

  for (long i = 0; i < n; i++)
  {
    out = out + v1[i] * v2[i];
  }

  return out;
}

// Get alpha
double get_alpha2d(double dV2[][2], double q[], double V, double dV[][2], double &minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[])
{
  double alpha = 1E-8;
  double tmp[n][2];
  double tmp_vec[2 * n];
  double x_vec[2 * n];
  double dV_vec[2 * n];
  double dV2_vec[2 * n];
  double Hs[2 * n];
  double tmpval = sqrt(dot2d(2 * n, q, q));

  vectorize(x, x_vec, dV, dV_vec, n);
  for (long i = 0; i < 2 * n; i++)
  {
    tmp_vec[i] = x_vec[i] + alpha * q[i] / tmpval;
  }
  unvectorize(tmp, tmp_vec, dV, dV_vec, n);
  reset_values2d(n, dV2, tmp, Box);

  getdV2d(V, dV2, minr, tmp, n, rad, rad_type, Box, Boundary);
  vectorize(tmp, tmp_vec, dV2, dV2_vec, n);

  for (long i = 0; i < 2 * n; i++)
  {
    Hs[i] = (dV2_vec[i] - dV_vec[i]) / alpha * tmpval;
  }
  alpha = -4 * (dot2d(2 * n, dV_vec, q) / dot2d(2 * n, q, Hs));

  double norm_dfq = abs(dot2d(2 * n, dV_vec, q));

  if (alpha < 1E-5)
  {
    alpha = 1E-2;
  }

  int flag = 1;
  while (flag == 1)
  {

    for (long i = 0; i < 2 * n; i++)
    {
      tmp_vec[i] = x_vec[i] + alpha * q[i];
    }
    unvectorize(tmp, tmp_vec, dV, dV_vec, n);
    reset_values2d(n, dV2, tmp, Box);
    getdV2d(V, dV2, minr, tmp, n, rad, rad_type, Box, Boundary);
    vectorize(tmp, tmp_vec, dV2, dV2_vec, n);

    if (abs(dot2d(2 * n, dV2_vec, q)) <= norm_dfq / 2.3)
    {
      flag = 0;
    }
    alpha = alpha / 1.3;

    if (alpha < 5E-7)
    {
      flag = 0;
    }
  }

  return 1.3 * alpha;
}

//*********************************************************
//*********************************************************
void minimizeenergyNLCG3d(double V, double dV[][3], double &minr, double x[][3], long n, double rad[], long rad_type[], double Box[], int Boundary[], long &k, double minr_tol, double dV_tol)
{

  double x_vec[3 * n];
  double dV_vec[3 * n];
  double dV2[n][3];
  double dV2_vec[3 * n];
  double alpha;
  double beta;
  double q[3 * n];
  double s[3 * n];
  long reset_count = 0;

  // Evaluate Function
  reset_values3d(n, dV, x, Box);
  getdV3d(V, dV, minr, x, n, rad, rad_type, Box, Boundary);
  vectorize3d(x, x_vec, dV, dV_vec, n);

  // Initialize search direction q = -dV
  for (long i = 1; i < 3 * n; i++)
  {
    q[i] = -dV_vec[i];
  }

  k = 0;
  minr = 0;

  if (normdV3d(dV, n) / n / 2 < dV_tol)
  {
    minr = 1;
  }

  double dV2_vecminusdV_vec[3 * n];

  while (minr<minr_tol & normdV3d(dV, n) / n / 2> dV_tol) // & k < 25000)
  {

    alpha = get_alpha3d(dV2, q, V, dV, minr, x, n, rad, rad_type, Box, Boundary);
    vectorize3d(x, x_vec, dV2, dV2_vec, n);

    // Update position
    for (long i = 0; i < 3 * n; i++)
    {
      x_vec[i] = x_vec[i] + alpha * q[i];
    }
    unvectorize3d(x, x_vec, dV, dV_vec, n);

    for (long i = 0; i < 3 * n; i++)
    {
      dV2_vecminusdV_vec[i] = dV2_vec[i] - dV_vec[i];
    }
    beta = dot3d(3 * n, dV2_vec, dV2_vecminusdV_vec) / dot3d(3 * n, dV_vec, dV_vec);
    if (beta < 0)
    {
      beta = 0;
    }
    for (long i = 0; i < 3 * n; i++)
    {
      q[i] = -dV2_vec[i] + beta * q[i];
    }

    if (minr<0.53 | reset_count> 25 |
        dot3d(3 * n, dV2_vec, dV_vec) / dot3d(3 * n, dV_vec, dV_vec) > 2.5 | beta == 0)
    {
      for (long i = 0; i < 3 * n; i++)
      {
        q[i] = -dV2_vec[i];
      }
      reset_count = 0;
    }

    // Update info
    reset_values3d(n, dV, x, Box);
    getdV3d(V, dV, minr, x, n, rad, rad_type, Box, Boundary);
    vectorize3d(x, x_vec, dV, dV_vec, n);
    reset_count = reset_count + 1;
    k = k + 1;
  }
}

double normdV3d(double dV[][3], long n)
{
  double out = 0;
  for (int i = 0; i < n; i++)
  {
    out = out + sqrt(dV[i][0] * dV[i][0] + dV[i][1] * dV[i][1] + dV[i][2] * dV[i][2]);
  }

  return out;
}

// Vectorize variables
void vectorize3d(double x[][3], double x_vec[], double dV[][3], double dV_vec[], long n)
{
  double out = 0;
  long count = 0;
  for (int i = 0; i < n; i++)
  {
    x_vec[count] = x[i][0];
    x_vec[count + 1] = x[i][1];
    x_vec[count + 2] = x[i][2];
    dV_vec[count] = dV[i][0];
    dV_vec[count + 1] = dV[i][1];
    dV_vec[count + 2] = dV[i][2];
    count = count + 3;
  }
}

void unvectorize3d(double x[][3], double x_vec[], double dV[][3], double dV_vec[], long n)
{
  double out = 0;
  long count = 0;
  for (int i = 0; i < n; i++)
  {
    x[i][0] = x_vec[count];
    x[i][1] = x_vec[count + 1];
    x[i][2] = x_vec[count + 2];
    dV[i][0] = dV_vec[count];
    dV[i][1] = dV_vec[count + 1];
    dV[i][2] = dV_vec[count + 2];
    count = count + 3;
  }
}

void reset_pos3d(long n, double x[][3], double Box[])
{
  for (long i = 0; i < n; i++)
  {
    while (x[i][0] < Box[0])
      x[i][0] += (Box[1] - Box[0]);
    while (x[i][1] < Box[2])
      x[i][1] += (Box[3] - Box[2]);
    while (x[i][2] < Box[4])
      x[i][2] += (Box[5] - Box[4]);
    while (x[i][0] > Box[1])
      x[i][0] -= (Box[1] - Box[0]);
    while (x[i][1] > Box[3])
      x[i][1] -= (Box[3] - Box[2]);
    while (x[i][2] > Box[5])
      x[i][2] -= (Box[5] - Box[4]);
  }
}

double dot3d(long n, double v1[], double v2[])
{
  double out = 0;

  for (long i = 0; i < n; i++)
  {
    out = out + v1[i] * v2[i];
  }

  return out;
}

// Get alpha
double get_alpha3d(double dV2[][3], double q[], double V, double dV[][3], double &minr, double x[][3], long n, double rad[], long rad_type[], double Box[], int Boundary[])
{
  double alpha = 1E-8;
  double tmp[n][3];
  double tmp_vec[3 * n];
  double x_vec[3 * n];
  double dV_vec[3 * n];
  double dV2_vec[3 * n];
  double Hs[3 * n];
  double tmpval = sqrt(dot3d(3 * n, q, q));

  vectorize3d(x, x_vec, dV, dV_vec, n);

  for (long i = 0; i < 3 * n; i++)
  {
    tmp_vec[i] = x_vec[i] + alpha * q[i] / tmpval;
  }
  unvectorize3d(tmp, tmp_vec, dV, dV_vec, n);
  reset_values3d(n, dV2, tmp, Box);

  getdV3d(V, dV2, minr, tmp, n, rad, rad_type, Box, Boundary);
  vectorize3d(tmp, tmp_vec, dV2, dV2_vec, n);

  for (long i = 0; i < 3 * n; i++)
  {
    Hs[i] = (dV2_vec[i] - dV_vec[i]) / alpha * tmpval;
  }
  alpha = -4 * (dot3d(3 * n, dV_vec, q) / dot3d(3 * n, q, Hs));

  double norm_dfq = abs(dot3d(3 * n, dV_vec, q));

  if (alpha < 1E-5)
  {
    alpha = 1E-2;
  }

  int flag = 1;
  while (flag == 1)
  {

    for (long i = 0; i < 3 * n; i++)
    {
      tmp_vec[i] = x_vec[i] + alpha * q[i];
    }
    unvectorize3d(tmp, tmp_vec, dV, dV_vec, n);
    reset_values3d(n, dV2, tmp, Box);
    getdV3d(V, dV2, minr, tmp, n, rad, rad_type, Box, Boundary);
    vectorize3d(tmp, tmp_vec, dV2, dV2_vec, n);

    if (abs(dot3d(3 * n, dV2_vec, q)) <= norm_dfq / 2.3)
    {
      flag = 0;
    }
    alpha = alpha / 1.3;

    if (alpha < 5E-7)
    {
      flag = 0;
    }
  }

  return 1.3 * alpha;
}
