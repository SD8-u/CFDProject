import sys
import cfd
import plot
from mpi4py import MPI

#comm = MPI.Comm.Get_parent()
comm = MPI.COMM_WORLD

#Start CFD module and compute flow for given parameters with MPI
cfd.startUp()
x, y, u, v, p = cfd.computeFlow(int(sys.argv[1]), int(sys.argv[2]), 
                                      float(sys.argv[3]), float(sys.argv[4]), 
                                      float(sys.argv[5]), sys.argv[6])
cfd.cleanUp()

#Use 1st processor to generate the resulting visual plot
if(comm.Get_rank() == 0):
    plot.generate_plot(x, y, u, v, p)

comm.barrier()