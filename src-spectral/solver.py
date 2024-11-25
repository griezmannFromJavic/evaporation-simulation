import numpy as np
from velicine import l_c

N = 64 # no. of distribution points

def gaussLobato(N):
    xj = np.array([np.cos(np.pi * j / N) for j in range(N + 1)])
    return xj

nodes = gaussLobato(N)

def chebyshev_polynomial(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        T0, T1 = 1, x
        for _ in range(2, n + 1):
            Tn = 2 * x * T1 - T0
            T0, T1 = T1, Tn
        return Tn

def stable_polynomial_evaluation(a, x):
    """
    Evaluates the polynomial u_N(x) = a_0 T_0(x) + a_1 T_1(x) + ... + a_N T_N(x)
    using a stable recursive algorithm, where T_k are Chebyshev polynomials.

    Parameters:
    - a : list or array of coefficients [a_0, a_1, ..., a_N]
    - x : the value at which to evaluate the polynomial

    Returns:
    - u_N(x) : the evaluated polynomial at x
    """
    N = len(a) - 1
    B_N_plus_1 = 0
    B_N = a[N]

    for k in range(N - 1, -1, -1):
        B_k = a[k] + 2 * x * B_N - B_N_plus_1
        B_N_plus_1 = B_N
        B_N = B_k

    # Final calculation of u_N(x)
    u_N = a[0] - B_N_plus_1 + B_N * x
    return u_N

def chebyshev_derivative_matrix(N):
    '''
    Note that the computation of D is very sensitive to round-off errors. Therefore, it is rec-
    ommended for practical computations to use specificially adapted versions of the above
    formula, e.g. chebdif.m from “A Matlab Differentiation Matrix Suite” by J.A.C. Wei-
    deman and S.C. Reddy.
    '''
    # Compute Chebyshev nodes x_j
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Initialize the matrix D with zeros
    D = np.zeros((N + 1, N + 1))

    # Chebyshev coefficients
    c = np.ones(N + 1)
    c[0] = 2
    c[N] = 2

    # Fill the matrix entries
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                # Off-diagonal entries
                D[i, j] = ((-1) ** (i + j)) / (c[i] * c[j] * (x[i] - x[j]))
            else:
                # Diagonal entries
                if i == 0:
                    D[i, i] = (2 * N ** 2 + 1) / 6
                elif i == N:
                    D[i, i] = -(2 * N ** 2 + 1) / 6
                else:
                    D[i, i] = -x[i] / (2 * (1 - x[i] ** 2))

    return D

def chebyshev_derivative_matrix_stable(N):
    '''
    You're right that calculating the Chebyshev derivative matrix directly is sensitive to round-off errors,
    particularly due to division by small differences in the Chebyshev nodes. Here are some strategies and
    modifications to reduce these round-off errors.
    '''
    # Chebyshev nodes x_j
    x = np.cos(np.pi * np.arange(N + 1) / N)

    # Initialize the matrix D with zeros
    D = np.zeros((N + 1, N + 1))

    # Chebyshev coefficients
    c = np.ones(N + 1)
    c[0], c[N] = 2, 2  # Set c_0 and c_N to 2 as required

    # Calculate the off-diagonal entries to use symmetry for stability
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                # Off-diagonal entries
                D[i, j] = (c[i] * (-1) ** (i + j)) / (c[j] * (x[i] - x[j]))

    # Adjust diagonal entries based on the row sum to maintain precision
    for i in range(N + 1):
        D[i, i] = -np.sum(D[i, :])

    return D

def xCoordinate(l):
    'x ranges from -1 to 1 while l ranges from 0 to L.'
    x = 2 * l / l_c - 1
    return x
