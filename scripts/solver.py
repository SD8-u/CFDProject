import sys
import bloodflow
import plot
from mpi4py import MPI
import numpy as np

comm = MPI.Comm.Get_parent()

N = np.array(comm.Get_rank(), dtype='i')

bloodflow.startUp()
x, y, u, v, p = bloodflow.computeFlow(int(sys.argv[1]), int(sys.argv[2]), 
                                      float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
bloodflow.cleanUp()
if(comm.Get_rank() == 0):
    plot.generate_plot(x, y, u, v, p)

comm.Reduce([N, MPI.INT], None, op=MPI.SUM, root=0)