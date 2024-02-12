import bloodflow
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

x, y, u, v, p = bloodflow.computeFlow(4, 500, 1, 0.0001, 100)
tri = Triangulation(x, y)

contour_plot = plt.tricontourf(tri, p, cmap='viridis') 
plt.colorbar(contour_plot, label='Pressure')
plt.quiver(x, y, u, v, color='red')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Vector Field')
plt.grid(True)
plt.savefig('vector_field.png', dpi=500)
plt.close() 

u = u.reshape(34, 34)
v  = v.reshape(34, 34)

contour_plot = plt.tricontourf(tri, p, cmap='viridis') 
plt.colorbar(contour_plot, label='Pressure')
plt.streamplot(np.unique(x), np.unique(y), u.transpose(), 
                v.transpose(), density=1.3, linewidth=1, color='red')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Stream Plot')
plt.grid(True)
plt.savefig('stream_plot.png', dpi=500)
plt.close() 