import bloodflow
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

x, y, u, v, p = bloodflow.computeFlow(4, 50, 3, 0.0001, 0.001)
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

contour_plot = plt.tricontourf(tri, p, cmap='viridis') 
plt.colorbar(contour_plot, label='Pressure')
plt.streamplot(np.unique(x), np.unique(y), u.transpose(), 
                v.transpose(), density=1.3, linewidth=1, color='red')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Stream Plot')
plt.savefig('stream_plot.png', dpi=500)
plt.close() 