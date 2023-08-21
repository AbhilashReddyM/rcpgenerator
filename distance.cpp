#include "distance.h"

// Reset values
//*********************************************************
void reset_values2d(long n, double dV[][2], double x[][2], double Box[])
{
	for (long i = 0; i < n; i++)
	{
		dV[i][0] = 0;
		dV[i][1] = 0;

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

// getdV2d
//*********************************************************
void getdV2d(double V, double dV[][2], double &minr, double x[][2], long n, double rad[], long rad_type[], double Box[], int Boundary[])
{

	minr = 1;

	// Set x params for hist
	double xmin = Box[0];
	double xmax = Box[1];
	double xbinsize = rad[0] * 2;
	long xnbins = floor((xmax - xmin) / xbinsize);
	xbinsize = (xmax - xmin) / double(xnbins);

	// Set x params for hist
	double ymin = Box[2];
	double ymax = Box[3];
	double ybinsize = rad[0] * 2;
	long ynbins = floor((ymax - ymin) / ybinsize);
	ybinsize = (ymax - ymin) / double(ynbins);

	long maxh = 15;

	long bins[xnbins * ynbins * maxh];
	long bin_dims[] = {maxh, xnbins * ynbins};

	for (long i = 0; i < maxh * xnbins * ynbins; i++)
	{
		bins[i] = 0;
	}

	khist2d(x, n, xmin, xnbins, xbinsize,
					ymin, ynbins, ybinsize, bins, bin_dims, Boundary);

	for (long i = 0; i < xnbins * ynbins; i++)
	{

		long cell = i;
		long cellboundarytype[5];
		long neighboringcells[] = {0, 0, 0, 0, 0};

		if (bins[sub2ind(0, cell, bin_dims)] > 0)
		{

			getneighboringcell2d(xnbins, ynbins, cell, cellboundarytype,
													 Boundary, neighboringcells);

			// Get pair count
			long paircount = 0;
			paircount = getpaircount2d(cell, neighboringcells,
																 bins, bin_dims);

			// Gen Pairs List 2d
			if (paircount > 0)
			{
				long pairs[paircount][3];
				getpairslist2d(pairs, paircount, cellboundarytype,
											 cell, neighboringcells,
											 bins, bin_dims);

				double dist[paircount][2];
				getdist2d(dist, x, pairs, paircount, Box);
				double minr_tmp = 0;
				getdV2d_sub(V, dV, minr_tmp, dist, pairs, paircount, rad, rad_type);
				if (minr_tmp < minr)
				{
					minr = minr_tmp;
				}
			}

			// Test if cell is on boundary
			long xbin = cell % xnbins;
			long ybin = (cell - xbin) / xnbins;

			if (xbin == 0 & Boundary[0] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][0] - radtmp < Box[0])
					{
						double r = 2 * (x[particle][0] - Box[0]);
						double d = 2 * radtmp;
						dV[particle][0] = dV[particle][0] - (1 - r / d) / d;
					}
				}
			}

			if (ybin == 0 & Boundary[1] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][1] - radtmp < Box[2])
					{
						double r = 2 * (x[particle][1] - Box[2]);
						double d = 2 * radtmp;
						dV[particle][1] = dV[particle][1] - (1 - r / d) / d;
					}
				}
			}

			if (xbin == xnbins - 1 & Boundary[0] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][0] + radtmp > Box[1])
					{
						double r = 2 * (x[particle][0] - Box[1]);
						double d = 2 * radtmp;
						dV[particle][0] = dV[particle][0] + (1 + r / d) / d;
					}
				}
			}

			if (ybin == ynbins - 1 & Boundary[1] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][1] + radtmp > Box[3])
					{
						double r = 2 * (x[particle][1] - Box[3]);
						double d = 2 * radtmp;
						dV[particle][1] = dV[particle][1] + (1 + r / d) / d;
					}
				}
			}

			if (Boundary[0] == -1)
			{
				for (int particle = 0; particle < n; particle++)
				{

					double r_hat[] = {x[particle][0] - Box[1] / 2, x[particle][1] - Box[1] / 2};
					double r = sqrt(r_hat[0] * r_hat[0] + r_hat[1] * r_hat[1]);
					r_hat[0] = r_hat[0] / r;
					r_hat[1] = r_hat[1] / r;
					double d = 2 * rad[rad_type[particle]];

					if (r + d / 2 > Box[1] / 2)
					{
						r = 2 * (Box[1] / 2 - r);
						dV[particle][0] = dV[particle][0] + (1 - r / d) / d * r_hat[0];
						dV[particle][1] = dV[particle][1] + (1 - r / d) / d * r_hat[1];
					}
				}
			}
		}
	}
}
//*********************************************************

// Get Neighboring Cell
//*********************************************************
void getneighboringcell2d(long xnbins, long ynbins,
													long cell, long cellboundarytype[],
													int Boundary[], long neighboringcells[])
{

	long xbin = cell % xnbins;
	long ybin = (cell - xbin) / xnbins;

	long xbin_on_boundary = 0;
	long ybin_on_boundary = 0;

	if (xbin == xnbins - 1)
	{
		xbin_on_boundary = 1;
	}
	if (ybin == ynbins - 1)
	{
		ybin_on_boundary = 1;
	}
	if (ybin == 0)
	{
		ybin_on_boundary = -1;
	}

	if (xbin_on_boundary == 0)
	{
		if (ybin_on_boundary == 0)
		{
			neighboringcells[0] = 4;
			neighboringcells[1] = cell - xnbins;
			cellboundarytype[1] = 0;
			neighboringcells[2] = cell - xnbins + 1;
			cellboundarytype[2] = 0;
			neighboringcells[3] = cell + 1;
			cellboundarytype[3] = 0;
			neighboringcells[4] = cell + xnbins + 1;
			cellboundarytype[4] = 0;
		}
		else if (ybin_on_boundary == -1)
		{
			// On bottom y boundary, but not on x boundary
			if (Boundary[1] == 1)
			{
				neighboringcells[0] = 4;
				neighboringcells[1] = (ynbins - 1) * xnbins + xbin;
				cellboundarytype[1] = 2;
				neighboringcells[2] = (ynbins - 1) * xnbins + xbin + 1;
				cellboundarytype[2] = 2;
				neighboringcells[3] = cell + 1;
				cellboundarytype[3] = 0;
				neighboringcells[4] = cell + xnbins + 1;
				cellboundarytype[4] = 0;
			}
			else
			{
				neighboringcells[0] = 2;
				neighboringcells[1] = cell + 1;
				cellboundarytype[1] = 0;
				neighboringcells[2] = cell + xnbins + 1;
				cellboundarytype[2] = 0;
			}
		}
		else if (ybin_on_boundary == 1)
		{

			if (Boundary[1] == 1)
			{
				neighboringcells[0] = 4;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
				neighboringcells[2] = cell - xnbins + 1;
				cellboundarytype[2] = 0;
				neighboringcells[3] = cell + 1;
				cellboundarytype[3] = 0;
				neighboringcells[4] = xbin + 1;
				cellboundarytype[4] = 2;
			}
			else
			{
				neighboringcells[0] = 3;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
				neighboringcells[2] = cell - xnbins + 1;
				cellboundarytype[2] = 0;
				neighboringcells[3] = cell + 1;
				cellboundarytype[3] = 0;
			}
		}
	}
	else
	{

		if (ybin_on_boundary == 0)
		{
			if (Boundary[0] == 1)
			{
				neighboringcells[0] = 4;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
				neighboringcells[2] = cell - 2 * xnbins + 1;
				cellboundarytype[2] = 1;
				neighboringcells[3] = cell - xnbins + 1;
				cellboundarytype[3] = 1;
				neighboringcells[4] = cell + 1;
				cellboundarytype[4] = 1;
			}
			else
			{
				neighboringcells[0] = 1;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
			}
		}
		else if (ybin_on_boundary == -1)
		{
			if (Boundary[0] == 1 & Boundary[1] == 1)
			{
				neighboringcells[0] = 4;
				neighboringcells[1] = ynbins * xnbins - 1;
				cellboundarytype[1] = 2;
				neighboringcells[2] = ynbins * xnbins - xnbins;
				cellboundarytype[2] = 4;
				neighboringcells[3] = 0;
				cellboundarytype[3] = 1;
				neighboringcells[4] = cell + 1;
				cellboundarytype[4] = 1;
			}
			else if (Boundary[0] == 0 & Boundary[1] == 0)
			{
				neighboringcells[0] = 0;
			}
			else if (Boundary[0] == 1 & Boundary[1] == 0)
			{
				neighboringcells[0] = 2;
				neighboringcells[1] = 0;
				cellboundarytype[1] = 1;
				neighboringcells[2] = cell + 1;
				cellboundarytype[2] = 1;
			}
			else if (Boundary[0] == 0 & Boundary[1] == 1)
			{
				neighboringcells[0] = 1;
				neighboringcells[1] = ynbins * xnbins - 1;
				cellboundarytype[1] = 2;
			}
		}
		else if (ybin_on_boundary == 1)
		{
			if (Boundary[0] == 1 & Boundary[1] == 1)
			{
				neighboringcells[0] = 4;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
				neighboringcells[2] = cell - 2 * xnbins + 1;
				cellboundarytype[2] = 1;
				neighboringcells[3] = cell - xnbins + 1;
				cellboundarytype[3] = 1;
				neighboringcells[4] = 0;
				cellboundarytype[4] = 4;
			}
			else if (Boundary[0] == 0 & Boundary[1] == 0)
			{
				neighboringcells[0] = 1;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
			}
			else if (Boundary[0] == 1 & Boundary[1] == 0)
			{
				neighboringcells[0] = 3;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
				neighboringcells[2] = cell - 2 * xnbins + 1;
				cellboundarytype[2] = 1;
				neighboringcells[3] = cell - xnbins + 1;
				cellboundarytype[3] = 1;
			}
			else if (Boundary[0] == 0 & Boundary[1] == 1)
			{
				neighboringcells[0] = 1;
				neighboringcells[1] = cell - xnbins;
				cellboundarytype[1] = 0;
			}
		}
	}
}
//*********************************************************

// Get pair count
//*********************************************************
long getpaircount2d(long cell, long neighboringcells[], long bins[], long bin_dims[])
{
	long within_cell = bins[sub2ind(0, cell, bin_dims)];
	long paircount = within_cell * within_cell / 2 - within_cell / 2;
	for (long j = 0; j < neighboringcells[0]; j++)
	{
		paircount = paircount + within_cell * bins[sub2ind(0, neighboringcells[j + 1], bin_dims)];
	}

	return paircount;
}
//*********************************************************

// Gen Pairs List 2d
//*********************************************************
void getpairslist2d(long pairs[][3], long paircount, long cellboundarytype[],
										long cell, long neighboringcells[], long bins[], long bin_dims[])
{

	long pc = 0; // pair count
	// Get pairs within cell
	if (bins[sub2ind(0, cell, bin_dims)] > 1)
	{
		for (long p1 = 1; p1 < bins[sub2ind(0, cell, bin_dims)]; p1++)
		{
			for (long p2 = p1 + 1; p2 < bins[sub2ind(0, cell, bin_dims)] + 1; p2++)
			{
				pairs[pc][0] = bins[sub2ind(p1, cell, bin_dims)];
				pairs[pc][1] = bins[sub2ind(p2, cell, bin_dims)];
				pairs[pc][2] = 0;
				pc = pc + 1;
			}
		}
	}

	// Get pairs within neighboring cells
	// loop over neighboring cells
	for (long k = 1; k < neighboringcells[0] + 1; k++)
	{
		for (long p1 = 1; p1 < bins[sub2ind(0, cell, bin_dims)] + 1; p1++)
		{
			for (long p2 = 1;
					 p2 < bins[sub2ind(0, neighboringcells[k], bin_dims)] + 1; p2++)
			{
				pairs[pc][0] = bins[sub2ind(p1, cell, bin_dims)];
				pairs[pc][1] = bins[sub2ind(p2, neighboringcells[k], bin_dims)];
				pairs[pc][2] = cellboundarytype[k];
				pc = pc + 1;
			}
		}
	}
}
//*********************************************************

// Get Cell
//*********************************************************
void getdist2d(double dist[][2], double x[][2], long pairs[][3], long paircount, double Box[])
{

	for (long k = 0; k < paircount; k++)
	{

		// cout << "printing pairs info" << endl;
		// cout << pairs[k][2] << endl;

		if (pairs[k][2] == 0)
		{
			dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0];
			dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1];
		}
		else if (pairs[k][2] == 1)
		{
			dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0] + Box[1];
			dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1];
			if (dist[k][0] > Box[1])
			{
				dist[k][0] = dist[k][0] - 2 * Box[1];
			}
		}
		else if (pairs[k][2] == 2)
		{
			dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0];
			dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1] + Box[3];
			if (dist[k][1] > Box[3])
			{
				dist[k][1] = dist[k][1] - 2 * Box[3];
			}
		}
		else if (pairs[k][2] == 3 | pairs[k][2] == 4)
		{
			dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0] + Box[1];
			if (dist[k][0] > Box[1])
			{
				dist[k][0] = dist[k][0] - 2 * Box[1];
			}
			dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1] + Box[3];
			if (dist[k][1] > Box[3])
			{
				dist[k][1] = dist[k][1] - 2 * Box[3];
			}
		}
	}
}
//*********************************************************

// Get dV 2d
//*********************************************************
void getdV2d_sub(double V, double dV[][2], double &minr, double dist[][2], long pairs[][3], long paircount, double rad[], long rad_type[])
{

	double r2 = 0;
	double r = 0;
	double d = 0;
	double rat = 0;
	double tmp[] = {0.0, 0.0};
	minr = 2;
	V = 0;

	// Get force between particles
	for (long k = 0; k < paircount; k++)
	{
		r2 = dist[k][0] * dist[k][0] + dist[k][1] * dist[k][1];
		d = rad[rad_type[pairs[k][0]]] + rad[rad_type[pairs[k][1]]];
		rat = r2 / (d * d);
		if (rat < 1)
		{
			r = sqrt(r2);
			if (sqrt(rat) < minr)
			{
				minr = sqrt(rat);
			}
			tmp[0] = (1 / d - 1 / r) / d * dist[k][0];
			tmp[1] = (1 / d - 1 / r) / d * dist[k][1];

			dV[pairs[k][0]][0] = dV[pairs[k][0]][0] + tmp[0];
			dV[pairs[k][0]][1] = dV[pairs[k][0]][1] + tmp[1];

			dV[pairs[k][1]][0] = dV[pairs[k][1]][0] - tmp[0];
			dV[pairs[k][1]][1] = dV[pairs[k][1]][1] - tmp[1];

			V = V + 0.5 * (1 - rat) * (1 - rat);
		}
	}

	if (minr == 2)
	{
		minr = 1;
	}
}
//*********************************************************

// getdV3d
//*********************************************************
void getdV3d(double V, double dV[][3], double &minr, double x[][3], long n, double rad[], long rad_type[], double Box[], int Boundary[])
{

	minr = 1;

	// Set x params for hist
	double xmin = Box[0];
	double xmax = Box[1];
	double xbinsize = rad[0] * 2;
	long xnbins = floor((xmax - xmin) / xbinsize);
	xbinsize = (xmax - xmin) / double(xnbins);

	// Set y params for hist
	double ymin = Box[2];
	double ymax = Box[3];
	double ybinsize = rad[0] * 2;
	long ynbins = floor((ymax - ymin) / ybinsize);
	ybinsize = (ymax - ymin) / double(ynbins);

	// Set z params for hist
	double zmin = Box[4];
	double zmax = Box[5];
	double zbinsize = rad[0] * 2;
	long znbins = floor((zmax - zmin) / zbinsize);
	zbinsize = (zmax - zmin) / double(znbins);

	long maxh = 30;

	long *bins = NULL;
	bins = new long[xnbins * ynbins * znbins * maxh];
	long bin_dims[] = {maxh, xnbins * ynbins * znbins};

	for (long i = 0; i < maxh * xnbins * ynbins; i++)
	{
		bins[i] = 0;
	}

	khist3d(x, n, xmin, xnbins, xbinsize,
					ymin, ynbins, ybinsize,
					zmin, znbins, zbinsize, bins, bin_dims);

	long tmp2 = 0;
	for (long i = 0; i < xnbins * ynbins * znbins; i++)
	{

		long cell = i;
		long cellboundarytype[14];
		long neighboringcells[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

		if (bins[sub2ind(0, cell, bin_dims)] > 0)
		{
			getneighboringcell3d(xnbins, ynbins, znbins, cell, cellboundarytype,
													 Boundary, neighboringcells);

			// Get pair count
			long paircount = 0;
			paircount = getpaircount2d(cell, neighboringcells,
																 bins, bin_dims);

			// Gen Pairs List 2d
			if (paircount > 0)
			{
				long pairs[paircount][3];
				getpairslist2d(pairs, paircount, cellboundarytype,
											 cell, neighboringcells,
											 bins, bin_dims);

				long xbin = cell % xnbins;
				long zbin = floor(cell / (xnbins * ynbins));
				long ybin = (cell - zbin - xbin) / xnbins;

				double dist[paircount][3];
				getdist3d(dist, x, pairs, paircount, Box);

				double minr_tmp = 0;
				getdV3d_sub(V, dV, minr_tmp, dist, pairs, paircount, rad, rad_type);
				if (minr_tmp < minr)
				{
					minr = minr_tmp;
				}
			}

			// Test if cell is on boundary
			long xbin = cell % xnbins;
			long zbin = floor(cell / (xnbins * ynbins));
			long ybin = (cell - zbin * xnbins * ynbins - xbin) / xnbins;

			if (xbin == 0 & Boundary[0] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][0] - radtmp < Box[0])
					{
						double r = 2 * (x[particle][0] - Box[0]);
						double d = 2 * radtmp;
						dV[particle][0] = dV[particle][0] - (1 - r / d) / d;
					}
				}
			}

			if (ybin == 0 & Boundary[1] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][1] - radtmp < Box[2])
					{
						double r = 2 * (x[particle][1] - Box[2]);
						double d = 2 * radtmp;
						dV[particle][1] = dV[particle][1] - (1 - r / d) / d;
					}
				}
			}

			if (zbin == 0 & Boundary[2] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][2] - radtmp < Box[4])
					{
						double r = 2 * (x[particle][2] - Box[4]);
						double d = 2 * radtmp;
						dV[particle][2] = dV[particle][2] - (1 - r / d) / d;
					}
				}
			}

			if (xbin == xnbins - 1 & Boundary[0] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][0] + radtmp > Box[1])
					{
						double r = 2 * (x[particle][0] - Box[1]);
						double d = 2 * radtmp;
						dV[particle][0] = dV[particle][0] + (1 + r / d) / d;
					}
				}
			}

			if (ybin == ynbins - 1 & Boundary[1] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][1] + radtmp > Box[3])
					{
						double r = 2 * (x[particle][1] - Box[3]);
						double d = 2 * radtmp;
						dV[particle][1] = dV[particle][1] + (1 + r / d) / d;
					}
				}
			}

			if (zbin == znbins - 1 & Boundary[2] == 0)
			{
				for (int j = 1; j < bins[sub2ind(0, cell, bin_dims)] + 1; j++)
				{
					int particle = bins[sub2ind(j, cell, bin_dims)];
					double radtmp = rad[rad_type[particle]];
					if (x[particle][2] + radtmp > Box[5])
					{
						double r = 2 * (x[particle][2] - Box[5]);
						double d = 2 * radtmp;
						dV[particle][2] = dV[particle][2] + (1 + r / d) / d;
					}
				}
			}
		}
	}

	delete[] bins;
	bins = NULL;
}

// Reset values
//*********************************************************
void reset_values3d(long n, double dV[][3], double x[][3], double Box[])
{
	for (long i = 0; i < n; i++)
	{
		dV[i][0] = 0;
		dV[i][1] = 0;
		dV[i][2] = 0;

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

void getneighboringcell3d(long xnbins, long ynbins, long znbins,
													long cell, long cellboundarytype[],
													int Boundary[], long neighboringcells[])
{

	long xbin = cell % xnbins;
	long zbin = floor(cell / (xnbins * ynbins));
	long ybin = (cell - zbin * xnbins * ynbins - xbin) / xnbins;

	long xbin_on_boundary = 0;
	long ybin_on_boundary = 0;
	long zbin_on_boundary = 0;

	long xbint = 0;
	long zbint = 0;
	long ybint = 0;

	int count = 1;

	// bottom first row on left
	xbint = xbin - 1;
	ybint = ybin - 1;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom first row in middle
	xbint = xbin;
	ybint = ybin - 1;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom first row on right
	xbint = xbin + 1;
	ybint = ybin - 1;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom second row on left
	xbint = xbin - 1;
	ybint = ybin;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom second row in middle
	xbint = xbin;
	ybint = ybin;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom second row on right
	xbint = xbin + 1;
	ybint = ybin;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom third row on left
	xbint = xbin - 1;
	ybint = ybin + 1;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom third row in middle
	xbint = xbin;
	ybint = ybin + 1;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// bottom third row on right
	xbint = xbin + 1;
	ybint = ybin + 1;
	zbint = zbin - 1;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// middle down
	xbint = xbin;
	ybint = ybin - 1;
	zbint = zbin;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// middle down and right
	xbint = xbin + 1;
	ybint = ybin - 1;
	zbint = zbin;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// middle right
	xbint = xbin + 1;
	ybint = ybin;
	zbint = zbin;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	// middle uo and right
	xbint = xbin + 1;
	ybint = ybin + 1;
	zbint = zbin;
	updateneighboringcell3d(xbint, ybint, zbint, xbin_on_boundary, ybin_on_boundary,
													zbin_on_boundary, xnbins, ynbins, znbins, Boundary,
													neighboringcells, cellboundarytype, count);

	neighboringcells[0] = count - 1;
}

void getdist3d(double dist[][3], double x[][3], long pairs[][3], long paircount, double Box[])
{

	for (long k = 0; k < paircount; k++)
	{

		if (pairs[k][2] == 0)
		{
			dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0];
			dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1];
			dist[k][2] = x[pairs[k][0]][2] - x[pairs[k][1]][2];
		}
		else if (1)
		{
			if (pairs[k][2] == 1 || pairs[k][2] == 3 || pairs[k][2] == 5 || pairs[k][2] == 7)
			{
				dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0] + Box[1];
				if (dist[k][0] > Box[1])
				{
					dist[k][0] = dist[k][0] - 2 * Box[1];
				}
			}
			else
			{
				dist[k][0] = x[pairs[k][0]][0] - x[pairs[k][1]][0];
			}
			if (pairs[k][2] == 2 || pairs[k][2] == 3 || pairs[k][2] == 6 || pairs[k][2] == 7)
			{
				dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1] + Box[3];
				if (dist[k][1] > Box[3])
				{
					dist[k][1] = dist[k][1] - 2 * Box[3];
				}
			}
			else
			{
				dist[k][1] = x[pairs[k][0]][1] - x[pairs[k][1]][1];
			}
			if (pairs[k][2] == 4 || pairs[k][2] == 5 || pairs[k][2] == 6 || pairs[k][2] == 7)
			{
				dist[k][2] = x[pairs[k][0]][2] - x[pairs[k][1]][2] + Box[5];
				if (dist[k][2] > Box[5])
				{
					dist[k][2] = dist[k][2] - 2 * Box[5];
				}
			}
			else
			{
				dist[k][2] = x[pairs[k][0]][2] - x[pairs[k][1]][2];
			}
		}
	}
}

// Get dV 3d
//*********************************************************
void getdV3d_sub(double V, double dV[][3], double &minr, double dist[][3], long pairs[][3], long paircount, double rad[], long rad_type[])
{

	double r2 = 0;
	double r = 0;
	double d = 0;
	double rat = 0;
	double tmp[] = {0.0, 0.0, 0.0};
	minr = 2;
	V = 0;

	// Get force between particles
	for (long k = 0; k < paircount; k++)
	{
		r2 = dist[k][0] * dist[k][0] + dist[k][1] * dist[k][1] + dist[k][2] * dist[k][2];
		d = rad[rad_type[pairs[k][0]]] + rad[rad_type[pairs[k][1]]];
		rat = r2 / (d * d);
		if (rat < 1)
		{
			r = sqrt(r2);
			if (sqrt(rat) < minr)
			{
				minr = sqrt(rat);
			}
			tmp[0] = (1 / d - 1 / r) / d * dist[k][0];
			tmp[1] = (1 / d - 1 / r) / d * dist[k][1];
			tmp[2] = (1 / d - 1 / r) / d * dist[k][2];

			dV[pairs[k][0]][0] = dV[pairs[k][0]][0] + tmp[0];
			dV[pairs[k][0]][1] = dV[pairs[k][0]][1] + tmp[1];
			dV[pairs[k][0]][2] = dV[pairs[k][0]][2] + tmp[2];

			dV[pairs[k][1]][0] = dV[pairs[k][1]][0] - tmp[0];
			dV[pairs[k][1]][1] = dV[pairs[k][1]][1] - tmp[1];
			dV[pairs[k][1]][2] = dV[pairs[k][1]][2] - tmp[2];

			V = V + 0.5 * (1 - rat) * (1 - rat);
		}
	}

	if (minr == 2)
	{
		minr = 1;
	}
}
//*********************************************************

long modbin(long a, long b)
{

	long out = a;

	if (a == b)
	{
		out = 0;
	}
	else if (a < 0)
	{
		out = b + a;
	}

	return out;
}

void updateneighboringcell3d(long &xbint, long &ybint, long &zbint,
														 long &xbin_on_boundary, long &ybin_on_boundary, long &zbin_on_boundary,
														 long xnbins, long ynbins, long znbins, int Boundary[],
														 long neighboringcells[], long cellboundarytype[], int &count)
{

	xbin_on_boundary = 0;
	ybin_on_boundary = 0;
	zbin_on_boundary = 0;

	if (modbin(xbint, xnbins) != xbint)
	{
		xbin_on_boundary = 1;
	}
	xbint = modbin(xbint, xnbins);
	if (modbin(ybint, ynbins) != ybint)
	{
		ybin_on_boundary = 1;
	}
	ybint = modbin(ybint, ynbins);
	if (modbin(zbint, znbins) != zbint)
	{
		zbin_on_boundary = 1;
	}
	zbint = modbin(zbint, znbins);

	if (zbin_on_boundary == 0 | (zbin_on_boundary == 1 & Boundary[2] == 1))
	{
		if (ybin_on_boundary == 0 | (ybin_on_boundary == 1 & Boundary[1] == 1))
		{
			if (xbin_on_boundary == 0 | (xbin_on_boundary == 1 & Boundary[0] == 1))
			{
				neighboringcells[count] = zbint * xnbins * ynbins + ybint * xnbins + xbint;
				cellboundarytype[count] = xbin_on_boundary + 2 * ybin_on_boundary + 4 * zbin_on_boundary;
				count = count + 1;
			}
		}
	}
}
