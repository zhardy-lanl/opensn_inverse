import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

path = os.path.abspath(os.path.dirname(__file__))

data = np.loadtxt(f"{path}/data.txt")
X, Y, F = data[:, 0], data[:, 1], data[:, 2]
X, Y = np.meshgrid(np.unique(X), np.unique(Y))
F = F.reshape(X.shape)

plt.figure()
plt.xlabel("Background Density")
plt.ylabel("Target Density")
plt.title("Objective Function")
norm = LogNorm(F.min(), F.max())
plt.pcolormesh(X, Y, F, cmap='jet', norm=norm)
plt.colorbar()
plt.tight_layout()

plt.figure()
plt.xlabel("Background Density")
plt.ylabel("Target Density")
plt.title("Objective Function")
plt.pcolormesh(X, Y, F, cmap='jet')
plt.colorbar()
plt.tight_layout()

plt.show()
