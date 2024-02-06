import bloodflow
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

x, y, u, v, p = bloodflow.computeFlow(4, 30, 1, 0.0001, 100)

x1 = []
y1 = []
p1 = []
for i in range(len(p)):
    if p[i] != -1 and p[i] < 1000:
        x1.append(x[i])
        y1.append(y[i])
        p1.append(p[i])

tri = Triangulation(x1, y1)
contour_plot = plt.tricontourf(tri, p1, cmap='viridis') 
plt.colorbar(contour_plot, label='Pressure')

plt.quiver(x, y, u, v)

for i in range(len(p1)):
    print(p1[i])

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Vector Field')
plt.grid(True)
plt.savefig('fluid_plot4.png', dpi=500)
plt.close() 