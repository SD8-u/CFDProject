import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

def generate_plot(x, y, u, v, p): 
    p = p.reshape(len(np.unique(x)), len(np.unique(y)))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    contour_plot = ax.contourf(np.unique(x), np.unique(y), p.transpose(), 70, cmap='turbo') 
    fig.colorbar(contour_plot, label='Pressure', ax=ax)
    ax.quiver(x, y, u, v, color='white')
    ax.set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Fluid Vector Field')
    fig.savefig('vector_field.png', dpi=500)
    plt.close() 

    u = u.reshape(len(np.unique(x)), len(np.unique(y)))
    v = v.reshape(len(np.unique(x)), len(np.unique(y)))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    contour_plot = ax.contourf(np.unique(x), np.unique(y), p.transpose(), 70, cmap='turbo')
    fig.colorbar(contour_plot, label='Pressure')
    ax.streamplot(np.unique(x), np.unique(y), u.transpose(), 
                    v.transpose(), density=1, linewidth=1, arrowsize=0.5, color='white')
    ax.set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Fluid Stream Plot')
    fig.savefig('stream_plot.png', dpi=500)
    plt.close()
