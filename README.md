# RCPGENERATOR
This is a code developed by Kenneth Desmond. I found it [here](https://physics.emory.edu/faculty/weeks/ken/RCPAlgorithm.html). I put it here on Github for my convenience.


>I have developed an algorithm based on one proposed by Xu[1] to generate random close packing of disks in 2D and spheres in 3D. The source code can be downloaded from here. The procedure involves slowly swelling and contracting disks or spheres with time until an rcp configuration is generated. We start by initializing a random arrangement of particles. Then we swell the disks until some of the disks overlap, at which point we apply an energy minimization technique to move the disks so they no longer overlap. If energy minimization produces a non-overlapping state then we begin swelling the particles again, otherwise we contract the particles until energy minimization will give a non-overlapping state. By repeating this swelling/contraction process while slowly decreasing the swelling/contracting rate an rcp configuration can be generated. Below is a movie demonstrating our algorithm along with a flow chart outlining the algorithm.
>
>[1] N. Xu, J. Blawzdziewicz, and C. S. O'Hern, Random close packing revisited: Ways to pack frictionless disks, Phys Rev E 71, 061306 (2005)



# Building and Usage

Running `make` should build the executable. I tested it on Windows using mingw-gcc and linux. There is [documentation](./Doucmentation.pdf) available.


# Credit
If you use this code please cite this paper:

**"Random Close Packing of Disks and Spheres in Confined Geometries"
Kenneth W. Desmond & Eric R. Weeks, Phys. Rev. E 80, 051305 (2009)**