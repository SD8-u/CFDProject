<h1>Parallel Fluid Simulator using Finite Element Method (FEM) with PETSc </h1>
Solves the incompressible Navier Stokes equation using a semi-implicit approach with BDF2 time integration and Chroin-Temam Projection. Simulates laminar 2D fluid flow. 
<h2>Compiling</h2>
<p>To compile the simulator natively, create a build directory and execute compile.sh as follows: </p>

*mkdir build* \
*bash compile.sh* 

<h2>Apptainer</h2>
<p>Using apptainer will containerise the simulator and dependencies. Executing the .sif file both compiles and runs the simulator. To generate the .sif file:</p>

*apptainer build cfd_apptainer.def  cfd_apptainer.sif* 

<h2>Examples</h2>

![1000re05s](https://github.com/uol-feps-soc-comp3931-2324-classroom/final-year-project-SD8-u/assets/69975595/91c98c87-ed3c-4400-af4f-db61debae93d)

![backstep1](https://github.com/uol-feps-soc-comp3931-2324-classroom/final-year-project-SD8-u/assets/69975595/ac883fa0-3fa3-42c4-bf19-0d525ff9764d)
