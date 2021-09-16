Simple test for the GPU Matrix-free FEM algorithm of Ljunkvist being applied to a 2D hyperbolic system. In essence, it solves the shock tube problem of Sod in a 2D or 3D environment. To test it, do the following:

   - Create a "mesh" folder containing the typical alya mesh files (*.dims.dat, *.geo.dat, *.fix.bou);
   - Generate the initial fields (VELOC.alya, DENSI.alya, PRESS.alya), PRESS can be changed for TEMPE.alya;
   - Copy the executable to your case folder, one behind "mesh";
   - During runtime, provide the name of the mesh files to be read;
