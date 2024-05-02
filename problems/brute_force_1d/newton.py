import os
import numpy as np
import matplotlib.pyplot as plt

path = os.path.abspath(os.path.dirname(__file__))

data = np.loadtxt(f"{path}/data.txt")
X, Y = data[:, 0], data[:, 1]


def f(x):
    if isinstance(x, float):
        idx = np.searchsorted(X, x)
        x_low, y_low = X[idx], Y[idx]
        return y_low + df(x) * (x - x_low)
    elif isinstance(x, (list, np.ndarry)):
        return np.array([f(xi) for xi in x])
    else:
        TypeError()


def df(x):
    if isinstance(x, float):
        idx = np.searchsorted(X, x)
        x_low, x_high = X[idx], X[idx + 1]
        y_low, y_high = Y[idx], Y[idx + 1]
        dx, dy = x_high - x_low, y_high - y_low
        return dy / dx
    elif isinstance(x, (list, np.ndarray)):
        return np.array([df(xi) for xi in x])
    else:
        TypeError()


def plot_iteration(x0, label=""):
    y0, m = f(x0), df(x0)
    x_int = x0 - y0 / m

    x = [x0 - 0.5 if m < 0.0 else x_int - 0.5,
         x_int + 0.5 if m < 0.0 else x0 + 0.5]
    y = [y0 + m * (x[0] - x0), y0 + m * (x[1] - x0)]

    line, = plt.plot(x0, y0, marker="o", label=label)
    plt.plot(x, y, "-.", color=line.get_color())


plt.figure()
plt.xlabel("Density")
plt.ylabel("Objective Function")

plt.plot(X, Y, "-k", label="$F(x)$")
plt.axhline(y=0, color="r")

r = np.random.random()
x_n = min(X) + (max(X) - min(X)) * (1 - r)
plot_iteration(x_n, "$x_0$")
print(f"Initial Guess: x_n {x_n:8.3g}")

for nit in range(100):
    delta = f(x_n) / df(x_n)
    x_n -= delta

    plot_iteration(x_n, fr"$solution_{nit + 1}, \Delta x {delta:.3e}$")
    print(f"Iteration {nit + 1:3d}: x_n {x_n:8.3g}  delta {delta:12.3e}")

    if np.abs(delta) < 1.0e-10:
        break

plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig(f"{path}/newtons_method.pdf")
plt.show()
