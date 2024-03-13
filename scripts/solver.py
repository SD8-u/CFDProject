import sys
import bloodflow
import plot
from mpi4py import MPI

comm = MPI.Comm.Get_parent()
#comm = MPI.COMM_WORLD

bloodflow.startUp()
x, y, u, v, p = bloodflow.computeFlow(int(sys.argv[1]), int(sys.argv[2]), 
                                      float(sys.argv[3]), float(sys.argv[4]), 
                                      float(sys.argv[5]), sys.argv[6])
bloodflow.cleanUp()
if(comm.Get_rank() == 0):
    plot.generate_plot(x, y, u, v, p)

comm.barrier()