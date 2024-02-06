import bloodflow
import matplotlib.pyplot as plt
import numpy as np

x, y, u, v = bloodflow.computeFlow(4, 10, 1, 0.0001, 100)

plt.figure()
plt.quiver(x, y, u, v)

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fluid Vector Field')
plt.grid(True)
plt.savefig('fluid_plot1.png', dpi=500)
plt.close() 