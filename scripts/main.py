import bloodflow
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

x, y, u, v, p = bloodflow.computeFlow(4, 10, 20, 0.0001, 1000)

x1 = []
y1 = []
p1 = []
for i in range(len(p)):
    if p[i] != -1 and p[i] < 1000000:
        x1.append(x[i])
        y1.append(y[i])
        p1.append(p[i])

tri = Triangulation(x1, y1)
contour_plot = plt.tricontourf(tri, p1, cmap='viridis') 
plt.colorbar(contour_plot, label='Pressure')

plt.quiver(x, y, u, v, color='red')

#for i in range(len(p1)):
    #print(p1[i])

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Vector Field')
plt.grid(True)
plt.savefig('fluid_plot2.png', dpi=500)
plt.close() 