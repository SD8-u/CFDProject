import matplotlib.pyplot as plt
import numpy as np

#Generate plots for visualsing fluid
def generate_plot(x, y, u, v, p): 
    #Reshape velocity/pressure fields for visualisation
    p = p.reshape(len(np.unique(x)), len(np.unique(y)))
    u = u.reshape(len(np.unique(x)), len(np.unique(y)))
    v = v.reshape(len(np.unique(x)), len(np.unique(y)))

    #Derive magnitude for velocity distribution
    m = np.sqrt(np.square(u) + np.square(v))

    #Plot contour for pressure distribution and fluid stream
    fig = plt.figure()
    ax = fig.add_subplot(111)
    contour_plot = ax.contourf(np.unique(x), np.unique(y), p.transpose(), 70, cmap='turbo')
    fig.colorbar(contour_plot, label='Pressure')
    ax.streamplot(np.unique(x), np.unique(y), u.transpose(), 
                    v.transpose(), density=2, linewidth=0.9, arrowsize=0.5, color='white')
    ax.set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Fluid Stream Plot')
    fig.savefig('stream_plot.png', dpi=500)
    plt.close()

    #Plot contour for velocity magnitude
    fig = plt.figure()
    ax = fig.add_subplot(111)
    contour_plot = ax.contourf(np.unique(x), np.unique(y), m.transpose(), 70, cmap='turbo')
    fig.colorbar(contour_plot, label='Velocity')
    ax.set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Velocity Magnitude Plot')
    fig.savefig('velocity_plot.png', dpi=500)
    plt.close()

