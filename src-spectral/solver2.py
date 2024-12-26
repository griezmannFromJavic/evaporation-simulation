import numpy as np
from scipy.linalg import solve

def chebyshev_points(N):
    """Generate Chebyshev collocation points and weights."""
    x = np.cos(np.pi * np.arange(N) / (N - 1))
    return x

def chebyshev_diff_matrix(N):
    """Generate the Chebyshev differentiation matrix."""
    x = chebyshev_points(N)
    c = np.ones(N)
    c[0], c[-1] = 2, 2
    c = c * (-1)**np.arange(N)

    X = np.tile(x, (N, 1))
    dX = X - X.T
    D = (c[:, None] / c) / (dX + np.eye(N))
    D -= np.diag(np.sum(D, axis=1))
    return D, x

def solve_poisson_2d(N, f_func):
    """Solve the 2D Poisson equation using the Chebyshev method."""
    D, x = chebyshev_diff_matrix(N)
    D2 = D @ D  # Second derivative matrix

    # Grid in (x, y) space
    X, Y = np.meshgrid(x, x)
    f = f_func(X, Y)

    # Construct the Laplacian operator
    I = np.eye(N)
    L = np.kron(I, D2) + np.kron(D2, I)

    # Apply boundary conditions: u = 0 on boundaries
    u = np.zeros((N, N))
    f = f.flatten()
    for i in [0, N - 1]:
        L[i * N:(i + 1) * N, :] = 0  # Top/Bottom rows
        L[:, i * N:(i + 1) * N] = 0
        L[i * N + np.arange(N), i * N + np.arange(N)] = 1
        f[i * N:(i + 1) * N] = 0

    # Solve
    u_flat = solve(L, f)
    u = u_flat.reshape(N, N)
    return X, Y, u

# Define the right-hand side of the Poisson equation
def f_func(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)

# Solve and visualize
N = 32
X, Y, U = solve_poisson_2d(N, f_func)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, U, cmap='viridis')
plt.title('Solution to Poisson Equation')
plt.show()
