#ifndef SOLVER_UTIL
#define SOLVER_UTIL 

#include <stdio.h>


void copyMatrix(double*** G, double *** U, int n) {

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n+1; j++) {  // N+1 because it's an augmented matrix
            U[i][j] = G[i][j];  // Copy each element
        }
    }
}

#endif