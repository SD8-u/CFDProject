import bloodflow
import matplotlib.pyplot as plt
import numpy as np

x, y, u, v = bloodflow.computeFlow(3, 100, 10, 0.1, 10)

for i in range(0, len(x)):
    if x[i] < 0.15 or y[i] < 0.15 or x[i] > 2.85 or y[i] > 2.85:
        u[i] = 0
        v[i] = 0

# Plot vector field
plt.figure()
plt.quiver(x, y, u, v, scale=1)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Vector Field')
plt.grid(True)
plt.savefig('fluid_plot.png', dpi=500)
plt.close()