#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 16  // Number of Chebyshev nodes for space
#define M 16  // Number of Chebyshev nodes for time

// Function to compute Chebyshev nodes
double* chebyshev_nodes(int n) {
    double* x = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        x[i] = cos(M_PI * i / (n - 1));
    }
    return x;
}

// Function to compute Chebyshev differentiation matrix
void chebyshev_diff_matrix(int n, double* x, double** D) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                if (i == 0 || i == n - 1) {
                    D[i][j] = (2.0 * (n - 1) * (n - 1) + 1.0) / 6.0;
                } else {
                    D[i][j] = -x[i] / (2.0 * (1.0 - x[i] * x[i]));
                }
            } else {
                double c_i = (i == 0 || i == n - 1) ? 2.0 : 1.0;
                double c_j = (j == 0 || j == n - 1) ? 2.0 : 1.0;
                D[i][j] = c_i / c_j * pow(-1.0, i + j) / (x[i] - x[j]);
            }
        }
    }
}

// Function to initialize the solution
void initialize(double** u, int n, int m, double* x, double* t) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            u[i][j] = exp(-10.0 * x[i] * x[i]) * sin(M_PI * t[j]); // Example initial condition
        }
    }
}

// Main function
int main() {
    // Allocate memory for nodes, differentiation matrices, and solution
    double* x = chebyshev_nodes(N);
    double* t = chebyshev_nodes(M);

    double** Dx = (double**)malloc(N * sizeof(double*));
    double** Dt = (double**)malloc(M * sizeof(double*));
    for (int i = 0; i < N; i++) {
        Dx[i] = (double*)malloc(N * sizeof(double));
    }
    for (int i = 0; i < M; i++) {
        Dt[i] = (double*)malloc(M * sizeof(double));
    }

    chebyshev_diff_matrix(N, x, Dx);
    chebyshev_diff_matrix(M, t, Dt);

    double** u = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        u[i] = (double*)malloc(M * sizeof(double));
    }

    initialize(u, N, M, x, t);

    // Compute the 2D Laplacian in Chebyshev space
    double** Lu = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        Lu[i] = (double*)calloc(M, sizeof(double));
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < N; k++) {
                Lu[i][j] += Dx[i][k] * u[k][j];
            }
            for (int k = 0; k < M; k++) {
                Lu[i][j] += Dt[j][k] * u[i][k];
            }
        }
    }

    // Open a file for writing
    FILE *file = fopen("solution.csv", "w");
    if (file == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

    // Output solution (Laplacian result as example)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            fprintf(file, "%.2f", Lu[i][j]); // Print the number with 2 decimals
            if (j < M - 1) {
                fprintf(file, ","); // Add a comma between values
            }
        }
        fprintf(file, "\n"); // Newline after each row
    }

    // Free memory
    free(x);
    free(t);
    for (int i = 0; i < N; i++) {
        free(Dx[i]);
        free(u[i]);
        free(Lu[i]);
    }
    for (int i = 0; i < M; i++) {
        free(Dt[i]);
    }
    free(Dx);
    free(Dt);
    free(u);
    free(Lu);

    return 0;
}
