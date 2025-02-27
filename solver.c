#include <stdio.h>
#include <stdlib.h>
#include "Lab3IO.h"
#include "solver_util.h"
#include "timer.h"
#include <omp.h>
#include <math.h>

typedef struct
{
    double value; // Stores the maximum absolute value
    int index;    // Stores the row index where max value was found
} Pivot;

#pragma omp declare reduction(max_pivot:Pivot : omp_out = (omp_in.value > omp_out.value) ? omp_in : omp_out) \
    initializer(omp_priv = omp_orig)

double **GausianElimination(double **G, double **U, int n)
{
    Pivot kp;
#pragma omp parallel
    {
        for (int k = 0; k < n - 1; k++)
        {
            kp.index = k;
            kp.value = fabs(U[k][k]);
#pragma omp for reduction(max_pivot : kp) schedule(static)
            for (int i = k + 1; i < n; i++)
            {
                if (fabs(U[i][k]) > kp.value)
                {
                    kp.index = i;
                    kp.value = fabs(U[i][k]);
                }
            }

            if (k != kp.index)
            {
#pragma omp for
                for (int j = 0; j < n + 1; j++)
                {
                    double tmp = U[k][j];
                    U[k][j] = U[kp.index][j];
                    U[kp.index][j] = tmp;
                }
            }

            // Elimination
            double elimination;
#pragma omp for
            // for (int j = k; j < n + 1; j++) { // Loop over columns first
            //     double pivot_col_val = U[k][j]; // Extract the pivot row value for column j
            //     for (int i = k + 1; i < n; i++) { // Loop over rows below the current row k
            //         U[i][j] -= (U[i][k] / U[k][k]) * pivot_col_val;
            //     }
            // }
            for (int i = k + 1; i < n; i++)
            { // Loop over rows below currnet row k
                double factor = U[i][k] / U[k][k];
                for (int j = k; j < n + 1; j++)
                { // loop over columns    diagonal
                    U[i][j] -= factor * U[k][j];
                }
            }
        }
    }
    if (DEBUG)
    {
        printf("\nGausian Elmination\n");
        PrintMat(U, n, n + 1);
    }

    return U;
}

double **JordanElimination(double **U, double *x, int n)
{
#pragma omp parallel
    {
        for (int k = n - 1; k > 0; k--)
        {
            double kn = U[k][n];
#pragma omp for
            for (int i = 0; i < k; i++)
            {
                U[i][n] -= (U[i][k] / U[k][k]) * U[k][n];
                U[i][k] = 0;
            }
        }
        for (int i = 0; i < n; i++)
        {
            x[i] = U[i][n] / U[i][i];
        }
    }
    if (DEBUG)
    {
        printf("\nJordan Elimination\n");
        PrintMat(U, n, n + 1);
    }
    return U;
}

double *solveForX(double **U, double *x, int n)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        x[i] = U[i][n] / U[i][i];
    }
    if (DEBUG)
    {
        printf("\nOutput Vector\n");
        PrintVec(x, n);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: ./main numThreads");
        return 1;
    }
    int numThreads = atoi(argv[1]);
    omp_set_num_threads(numThreads);
    double **A, **U;
    double *x;
    int size;

    Lab3LoadInput(&A, &size);
    U = CreateMat(size, size + 1); // copy for modifications
    copyMatrix(A, U, size);
    x = CreateVec(size);

    if (DEBUG)
    {
        printf("Input Matrix\n");
        PrintMat(A, size, size + 1);
    }
    double start, end, gStart, gEnd, jStart, jEnd, xStart, xEnd;

    GET_TIME(start);

    GausianElimination(A, U, size);
    // GET_TIME(gStart);
    JordanElimination(U, x, size);

    GET_TIME(end);
    // printf("Gausian: %10.6f, Jordan: %10.6f\n", start - gStart, gStart - end);

    Lab3SaveOutput(x, size, end - start);

    return 0;
}