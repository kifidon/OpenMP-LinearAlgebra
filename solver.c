#include <stdio.h>
#include <stdlib.h>
#include "Lab3IO.h"
#include "solver_util.h"
#include <omp.h>
#include <math.h>

#define THREADS 4

double ** GausianElimination(double ** G, int n){
    // Copy Matrix for saftey
    double  ** U;
    // #pragma omp parallel num_threads(THREADS) shared(U)
    // {
        // #pragma omp single
        // {
            U = CreateMat(n, n+1);
            copyMatrix(&G, &U, n);
        // }
        for (int k = 0; k < n-1; k++){
            int kp = k; //ins Pivot row 
            for (int i = k+1; i<n; i++) {
                if (fabs(U[i][k]) > fabs(U[kp][k])){
                    kp = i;
                }
            }
            if (kp != k) {
                for (int j = 0; j < n+1; j++) {
                    double tmp = U[k][j];
                    U[k][j] = U[kp][j];
                    U[kp][j] = tmp;
                }
            }
            // Elimination
            for (int i = k+1; i < n; i++) {
                double factor = U[i][k] / U[k][k];
                for (int j = k; j < n+1; j++) {
                    U[i][j] -= factor * U[k][j];
                }
            }
        } 
        PrintMat(U, n, n+1);
        return U;
    // }
}



int main(int argc, char* argv[]) {
    double ** A;
    int size; 
    Lab3LoadInput(&A, &size);
    return 0;
}