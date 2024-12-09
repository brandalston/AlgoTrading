//
// Created by Brandon Alston on 11/2/24.
//

#ifndef MATRIXUTILITY_H
#define MATRIXUTILITY_H

#include "constants.h"
#include "statUtility.h"
#include "/Users/brandonalston/newmat10/newmat.h"
#include "/Users/brandonalston/newmat10/newmatap.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
// #include <ql>
// using QuantLib;
using namespace std;

class MatrixUtil {
public:
    double* genCorrelatedDeviates(const SymmetricMatrix& R, double dt, double z[])
    {
        int i, j;
        double sum[4] = {0.0};
        double deviate = 0.0; // standard normal deviate
        int m = R.Nrows(); // number of rows in correlation matrix
        vector<double> dz; // vector of correlated deviates
        vector<double> eigenValue; // vector of eigenvalues
        vector<double> eigenVector[4]; // array of vector of eigenvectors
        vector<double>::iterator eigenVecIter; // vector iterator
        double lambda[4] = {0.0}; // stores eigenvalues of correlation matrix R
        double dw[4] = {0.0}; // stores correlated deviates
        DiagonalMatrix D(m); // diagonal matrix
        Matrix V(m,m); // m x n matrix
        D = EigenValues(R); // get eigenvalues
        V = genEigenVectors(R); // get eigenvectors
        // store eigenvalues
        for (i = 0; i < m; i++) {
            eigenValue.push_back(D.element(i,i));
            lambda[i] = D.element(i,i);
        }
        // stores rows of eigenvectors so that we can compute
        // dz[i] = v[i][1]*sqrt(eigenvalue[1])*dw1 + v[i][2]*sqrt(eigenvalue[2])*dw2
        // + . . .
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                eigenVector[i].push_back(V.element(i,j));
            }
        }
        srand(0); long seed = (long) rand() % 100; long *idum = &seed;
        // generate uncorrelated deviates
        for (i = 0; i < m; i++) {
            deviate = util.NormalDeviate(idum);
            dw[i] = deviate*sqrt(dt);
        }

        // generate correlated deviates
        for (i = 0; i < m; i++) {
            eigenVecIter = eigenVector[i].begin();
            for (j = 0; j < m; j++) {
                sum[i] += (*eigenVecIter)*sqrt(lambda[j])*dw[j];
                eigenVecIter++;
                }
            z[i] = sum[i];
        }
    return z;
    }

    double* genCorrelatedDeviatesCholesky(const SymmetricMatrix& R, double dt, double z[]) {
        int m = R.Nrows(), n = R.Ncols();
        Matrix lb(m,n);
        StatUtility statUtil;
        double deviate = 0.0, dw[4] = {0.0}, sum = 0.0;
        long seed = 0; long* idum = 0; int i, j;
        // number of rows
        // number of columns
        // lower-banded (lb) matrix
        // Statistical utility class
        // standard normal deviate
        // stores deviate*sqrt(dt)
        // seed for RNG
        // stores address of seed
        lb = Cholesky(R); // calls Cholesky routine in NEWMAT library
        srand(time(0)); // initialize RNG
        seed = long(rand()) % 100; // generate seed
        idum = &seed; // store address of seed
        // generate uncorrelated deviates
        for (i = 0; i < m; i++) {
            deviate = statUtil.gasdev(idum); // generate normal (gaussian) deviate
            dw[i] = deviate*sqrt(dt);
        }
        // generate correlated deviates
        for (i = 0; i < m; i++) {
            sum = 0;
            for (j = 0; j < m; j++) {
                sum += lb.element(i,j)*dw[j];
            }
            z[i] = sum;
        }
        return z;
    }

};


#endif //MATRIXUTILITY_H
