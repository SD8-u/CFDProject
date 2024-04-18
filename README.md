---Compiling---
To compile the simulator natively, create a build directory and execute compile.sh as follows:
*mkdir build*
*bash compile.sh*
---Apptainer---
Using apptainer will containerise the simulator and dependencies. To generate the .sif file:
*apptainer build cfd_apptainer.def  cfd_apptainer.sif*
Executing the .sif file both compiles and runs the simulator. 
