import os
import numpy as np
from matplotlib import pyplot as plt

obj = []
alpha = []
densities = []
with open("res.txt", "r") as f:
    for line in f:
        entries = line.strip().split()
        for i, entry in enumerate(entries):
            if entry == "Func:":
                obj.append(float(entries[i + 1]))
            if entry == "Alpha:":
                alpha.append(float(entries[i + 1]))
            if entry == "Densities:":
                vals = []
                for j in range(i + 1, len(entries)):
                    if entries[j] != "CONVERGED":
                        val_str = entries[j].strip("[").strip("]")
                        vals.append(float(val_str))
                densities.append(vals)

plt.figure()
plt.xlabel("Iteration", fontsize=12)
plt.ylabel("Convergence", fontsize=12)
plt.grid(True)

iters = list(range(1, len(obj) + 1))
plt.semilogy(iters, obj, "-k*", label="Objective Function")

densities = np.array(densities).T
styles = ["-*b", "-*r"]
for i in range(len(densities)):
    ref = 5.0 if i == 0 else 2.0
    y = np.abs(abs(densities[i] - ref))
    color = "Red" if i == 0 else "Gray"
    plt.semilogy(iters, y, styles[i], label=f"{color} Material L1-Error")

plt.legend()
plt.tight_layout()
plt.show()
