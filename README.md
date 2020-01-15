#fdaPDE_dev

This repository contains the development version of fdaPDE package (future 2.0 version).

New features wrt CRAN: smooth regression for manifold and volumetric domains, also with areal data. Smooth fPCA over 2D, 2.5D and 3D domains, also with areal data.

Bugs in fPCA, boundary conditions and space-varying regression have been fixed. Issues of point location in 2.5D have been fixed. Areal data are still undergoing tests.

Compiled in both Win RStudio and Ubuntu 18.04 using g++ compiler. If using a Linux machine, it is advisable to install rgl, plot3D and plot3Drgl before fdaPDE.

Subfolder structure:
/src contains all C++ code and a special file named Makevars necessary to build and install the R package, 
/R contains the R functions that wrap the C++ calls

Remarks:

1) the meshes 2.5D and 3D contained in /data are outdated (nodes are stored in a vector instead of a matrix). To create a consistent mesh    from those, do:

   data(meshName)

   vertici <- matrix(data=meshName$nodes,ncol=3,byrow=T)
   triangoli <- matrix(data=meshName$triangles, ncol=3, byrow=T) # tetrahedrons with 4 columns in 3D

   mesh <- fdaPDE::create.MESH.2.5D(nodes = vertici,triangles = triangoli)
   
2) the test scripts in /tests are outdated, and will be updated as soon as possible. Scripts in /RScripts have not been verified yet.

3) the shift of indexes from R to C++ is done within the R functions smooth.FEM.basis and smooth.FEM.FPCA. Do not use C++ scripts directly on the R mesh objects, unless you take care of shifing indexes by yourself.
 
