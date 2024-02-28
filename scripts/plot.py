import bloodflow
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

def generate_plot(refinement, timesteps, velocity, dt, viscosity):
    x, y, u, v, p = bloodflow.computeFlow(refinement, timesteps, velocity, dt, viscosity)
    
    tri = Triangulation(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    contour_plot = ax.tricontourf(tri, p, cmap='viridis') 
    fig.colorbar(contour_plot, label='Pressure', ax=ax)
    ax.quiver(x, y, u, v, color='red')
    ax.set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Fluid Vector Field')
    fig.savefig('vector_field.png', dpi=500)
    plt.close() 

    u = u.reshape(len(np.unique(x)), len(np.unique(y)))
    v  = v.reshape(len(np.unique(x)), len(np.unique(y)))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    contour_plot = ax.tricontourf(tri, p, cmap='viridis')
    fig.colorbar(contour_plot, label='Pressure')
    ax.streamplot(np.unique(x), np.unique(y), u.transpose(), 
                    v.transpose(), density=1.3, linewidth=0.8, arrowsize=0.8, color='red')
    ax.set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Fluid Stream Plot')
    fig.savefig('stream_plot.png', dpi=500)
    plt.close()
