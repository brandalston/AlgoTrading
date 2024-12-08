//
// Created by Brandon Alston on 11/2/24.
//

#ifndef STATUTILITY_H
#define STATUTILITY_H

#include <ctime>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <queue>
#include "/Users/brandonalston/newmat10/newmat.h"
#include "/Users/brandonalston/newmat10/newmatap.h"

#include <quant>
using namespace QuantLib;
using namespace std;

// compute a Sobol sequence
#define GRAY(n) (n ^ ( n >> 1 )) // for Sobol sequence
#define MAXDIM 5
#define VMAX 30
struct sobolp {
    double sequence[MAXDIM];
    int x[MAXDIM];
    int v[MAXDIM][VMAX];
    double RECIPD;
    int _dim; int _skip;
    unsigned long _nextn;
    unsigned long cur_seed;
    // dimension of the sample space
};

class StatUtility {
public:
    /************************************************************************
    sobolp_generateSamples : generates a Sobol sequence
    [in]: struct sobolp* config : pointer to Sobol structure
    double* samples : pointer to sample values
    [out]: void
    ************************************************************************/
    void sobolp_generateSamples(struct sobolp* config, double* samples) {
        int i;
        nextSobol(config, config->cur_seed);
        config->cur_seed++;
        for(i = 0; i < config->_dim; i++ )
            samples[i] = config->sequence[i];
    }

    /******************************************************************************
    nextSobolNoSeed : generates the next Sobol seed number to
    generate the next Sobol value
    : pointer to Sobol structure
    [in]: struct sobolp* config [out]: void
    ******************************************************************************/
    static void nextSobolNoSeed(struct sobolp* config) {
        int c = 1;
        int i;
        int save = config->_nextn;
        while((save %2) == 1) {
            c += 1;
            save = save /2;
        }
        for(i=0;i<config->_dim;i++) {
            config->x[i] = config->x[i]^(config->v[i][c-1]<< (VMAX-c));
            config->sequence[i] = config->x[i]*config->RECIPD;
        }
        config->_nextn += 1;
    }

    /******************************************************************************
    sobolp_init : initializes the Sobol algorithm
    [in]: sobolp* config : pointer to Sobol
    int dim : dimension of the sample spaces
    unsigned long seed : seed for Sobol number generator
    [out] : void
    ******************************************************************************/
    void sobolp_init(struct sobolp* config, int dim, unsigned long seed) {
        int d[MAXDIM], POLY[MAXDIM];
        int save;
        int m,i,j,k;
        config->_dim = dim;
        config->_nextn = 0;
        config->RECIPD = 1.0 / pow( 2.0, VMAX );
        config->cur_seed = seed;
        POLY[0] = 3; d[0] = 1;  // x + 1
        POLY[1] = 7; d[1] = 2;  // x^2 + x + 1
        POLY[2] = 11; d[2] = 3; // x^3 + x + 1
        POLY[3] = 19; d[3] = 4; // x^4 + x + 1
        POLY[4] = 37; d[4] = 5; // x^5 + x^2 + 1

        for(i = 0; i < config->_dim; i++ ) {
            for(j = 0; j < d[i]; j++ ) {
                config->v[i][j] = 1;
            }
        }
        for( i = 0; i < config->_dim; i++ ) {
            for( j = d[i]; j < VMAX; j++ ) {
                config->v[i][j] = config->v[i][j-d[i]];
                save = POLY[i];
                m = pow( 2, d[i] );
                for( k = d[i]; k > 0; k-- ) {
                    config->v[i][j] = config->v[i][j] ^ m*(save%2)*config->v[i][j-k];
                    save = save/2;
                    m = m/2;
                }
            }
        }
        for( i = 0; i < config->_dim; i++ )
            config->x[i]=0;
        config->_skip = pow( 2, 6 );
        for( i = 1; i <= config->_skip; i++ )
            nextSobolNoSeed(config);
    }

    /*****************************************************************************
    generateFaure M : generates a Faure sequence of length M
    [in] long N : number of time steps
    long M : number of simulations
    [out]: vector<double> X : the Faure sequence
    ******************************************************************************/
    vector<double> generateFaure(long N, long M) {
        int p = generatePrime(N);
        int l, q, k;
        long v1, v2, v3;
        long value = 0;
        long a[250][250] = {0};
        int m = (int) (log(M)/log(p));
        if (m == 0)
            m = 1;
        long x[] = {0};
        unsigned long fact = 0;
        for (k = 1; k <= N; k++) {
            for (l = 0; l <= m; l++)
            {
                value = pow(p,l+1);
                a[0][l] = (int)((M % value)/p);
                for (q = l; q <= m; q++)
                {
                    v1 = factorial(q);
                    v2 = factorial(q-l);
                    v3 = factorial(l);
                    fact = v1/(v2*v3);
                    value = fact*a[k-1][q] % p;
                    a[k][l] = a[k][l] + value;
                }
                x[k] = x[k] + a[k][l]/pow(p,l+1);
            }
            x.push_back((double)x[k]);
        }
        return x;
        // number of time steps
    }

    /******************************************************************************
    factorial : computes the factorial of a number
    [in]: N : number to factorialize
    [out]: N!
    ******************************************************************************/
    long factorial(long N) {
        if ((N == 1) || (N == 0))
            return 1;
        else
            return N*factorial(N-1);
    }

    /*****************************************************************************
    generatePrime: This function computes the smallest prime greater than or equal to N
    [in]: long N : find prime >= N
    [out]: prime >= N
    /*****************************************************************************/
    long generatePrime(long N) {
        long i = N;
        bool flag = false;
        do
        {
            // check if number is prime
            if ((i % 2 != 0) && (i % 3 != 0) && (i % 4 != 0) && (i % 5 != 0)
            && (i % 7 != 0) && (i % 8 != 0) && (i % 9 != 0))
                flag = true;
            else
                i++;
        }
        while (flag != true);
        return i;
    }

    /*****************************************************************************
    polarRejection
    This function computes two standard deviates using polar rejection
    (transformation) method Returns the first deviate and stores the second
    deviates in a vector Y so that is can be used for another call rather than
    throwing it away.
    [in]: double y : seed value
    int i : ith standard deviate
    [out]: Y[i] : ith standard normal deviate in Y
    ******************************************************************************/
    double polarRejection(double y, int i)
    {
        double w = 0.0;
        double x1, x2, z1, z2, c;
        double temp = 0.0;
        double *idum = &y;
        do {
            x1 = gasdev((long*)idum);
            x2 = gasdev((long*)idum);
            w = x1*x1 + x2*x2;
        }
        while (w >= 1);
        c = sqrt(-2*(log(w)/w));
        z1 = c*x1;
        Y.push_back(z1);
        z2 = c*x2;
        Y.push_back(z2);
        return Y[i];
    }

    // Hull’s approximation of the cumulative normal distribution
    double normalCalc(double d)
    {
        const double a1 = 0.319381530;
        const double a2 = -0.356563782;
        const double a3 = 1.781477937;
        const double a4 = -1.821255978;
        const double a5 = 1.330274429;
        const double gamma = 0.2316419;
        const double k1 = 1/(1 + gamma*d);
        const double k2 = 1/(1 – gamma*d);
        const double normalprime = (1/(sqrt(2*PI)))*exp(-d*d/2);
        double value = 0.0;
        double h = 0.0;
        if (d >= 0)
            value = 1- normalprime*(a1*k1 + a2*pow(k1,2) + a3*pow(k1,3) + a4*pow(k1,4) +
            a5*pow(k1,5));
        else
            value = normalprime*(a1*k2 + a2*pow(k2,2) + a3*pow(k2,3) + a4*pow(k2,4) +
            a5*pow(k2,5));
        return value;
    }
};



#endif //STATUTILITY_H
