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
#include <cassert>
#include <queue>
#include <random>
// #include <quant>
#include "/Users/brandonalston/newmat10/newmat.h"
#include "/Users/brandonalston/newmat10/newmatap.h"
#include "nrutil.h"
using namespace std;

#define GRAY(n) (n ^ ( n >> 1 )) // for Sobol sequence
#define MAXBIT 30
#define MAXDIM 5
#define MAXDIM_NR 6
class StatUtility {
public:
    /************************************************************************
    Compute a Sobol sequence from Numerical Recipes in C by Press et al. (1992)
    When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM
    different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[1..n]
    the next values from n of these sequences. (n must not be changed between initializations.)
    ************************************************************************/

    void sobseq(int *n, float x[]) {
        int j,k,l;
        unsigned long i,im,ipp;
        static float fac;
        static unsigned long in,ix[MAXDIM_NR+1],*iu[MAXBIT+1];
        static unsigned long mdeg[MAXDIM_NR+1]={0,1,2,3,3,4,4};
        static unsigned long ip[MAXDIM_NR+1]={0,0,1,1,2,1,4};
        static unsigned long iv[MAXDIM_NR*MAXBIT+1]={0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
        if (*n < 0) { // Initialize, don’t return a vector.
            for (k=1;k<=MAXDIM_NR;k++) {
                ix[k]=0;
            }
            in=0;
            if (iv[1] != 1) return;
            fac=1.0/(1L << MAXBIT);
            for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM_NR) {
                iu[j] = &iv[k];
            }
            // To allow both 1D and 2D addressing.
            for (k=1;k<=MAXDIM_NR;k++) {
                for (j=1;j<=mdeg[k];j++) {
                    iu[j][k] <<= (MAXBIT-j);
                }
                // Stored values only require normalization.
                for (j=mdeg[k]+1;j<=MAXBIT;j++) { // Use the recurrence to get other valipp=ip[k]; ues.
                    i=iu[j-mdeg[k]][k];
                    i ^= (i >> mdeg[k]);
                    for (l=mdeg[k]-1;l>=1;l--) {
                        if (ipp & 1) i ^= iu[j-l][k];
                        ipp >>= 1;
                    }
                    iu[j][k]=i;
                }
            }
        } else { // Calculate the next vector in the seim=in++; quence.
            for (j=1;j<=MAXBIT;j++) { // Find the rightmost zero bit.
                if (!(im & 1)) break;
                im >>= 1;
            }
            if (j > MAXBIT) nrerror("MAXBIT too small in sobseq");
            im=(j-1)*MAXDIM_NR;
            for (k=1;k<=IMIN(*n,MAXDIM_NR);k++) { // XOR the appropriate direction number into each component of the vector and convert to a floating number.
                ix[k] ^= iv[im+k];
                x[k]=ix[k]*fac;
            }
        }
    }

    struct sobolp {
        double sequence[MAXDIM];
        int x[MAXDIM];
        int v[MAXDIM][MAXBIT];
        double RECIPD;
        int _dim;
        int _skip;
        unsigned long _nextn;
        unsigned long cur_seed;
        // dimension of the sample space
    };

    void nextSobol(struct sobolp* sobolp, unsigned long cur_seed) {
        /*
         *int c = 1;
        int i;
        int save = sobolp->_nextn;
        while((save %2) == 1) {
            c += 1;
            save = save /2;
        }
        for(i=0;i<sobolp->_dim;i++) {
            sobolp->x[i] = sobolp->x[i]^(sobolp->v[i][c-1]<< (MAXBIT-c));
            sobolp->sequence[i] = sobolp->x[i]*sobolp->RECIPD;
        }
        sobolp->_nextn += 1;
        */
    };

    void sobolp_generateSamples(struct sobolp* config, double* samples) {
        int i;
        nextSobol(config, config->cur_seed);
        config->cur_seed++;
        for(i = 0; i < config->_dim; i++ )
            samples[i] = config->sequence[i];
    }

    /******************************************************************************
    nextSobolNoSeed : generates the next Sobol seed number to
    generate the next Sobol value pointer to Sobol structure
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
            config->x[i] = config->x[i]^(config->v[i][c-1]<< (MAXBIT-c));
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
        config->RECIPD = 1.0 / pow( 2.0, MAXBIT );
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
            for( j = d[i]; j < MAXBIT; j++ ) {
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
        long p = generatePrime(N); p = int(p);
        int l, q, k;
        long v1, v2, v3;
        long value = 0;
        long a[250][250] = {0};
        int m = (int) (log(M)/log(p));
        if (m == 0)
            m = 1;
        vector<double> x = {0};
        long fact = 0;
        for (k = 1; k <= N; k++) {
            for (l = 0; l <= m; l++) {
                value = pow(p,l+1);
                a[0][l] = (int)((M % value)/p);
                for (q = l; q <= m; q++) {
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
    double polarRejection(double y, int i) {
        double w = 0.0;
        double x1, x2, z1, z2, c;
        // double temp = 0.0;
        double *idum = &y;
        vector<double> Y;
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
    double normalCalc(double d) {
        const double a1 = 0.319381530;
        const double a2 = -0.356563782;
        const double a3 = 1.781477937;
        const double a4 = -1.821255978;
        const double a5 = 1.330274429;
        const double gamma = 0.2316419;
        const double k1 = 1/(1 + gamma*d);
        const double k2 = 1/(1 - gamma*d);
        const double normalprime = (1/(sqrt(2*M_PI)))*exp(-d*d/2);
        double value = 0.0;
        // double h = 0.0;
        if (d >= 0)
            value = 1- normalprime*(a1*k1 + a2*pow(k1,2) + a3*pow(k1,3) + a4*pow(k1,4) + a5*pow(k1,5));
        else
            value = normalprime*(a1*k2 + a2*pow(k2,2) + a3*pow(k2,3) + a4*pow(k2,4) + a5*pow(k2,5));
        return value;
    }
    /************************************************
    Box-Muller algorithm to generate a Gaussian (normal) deviate from U(O,1)
    ************************************************/
    double gasdev(long *seed) {
        mt19937_64 rng;
        rng.seed(seed);
        uniform_real_distribution<double> unif(0, 1);
        // initialize a uniform distribution between 0 and 1
        double W=2.0, u1=0.0, u2=0.0, v1=0.0, v2=0.0;
        while (W>1) {
            u1 = unif(rng);
            u2 = unif(rng);
            v1 = 2*u1 -1, v2 = 2*u2 - 1;
            W = pow(v1,2) + pow(v2,2);
        }
        double N1 = v1*sqrt(-2*log(W)/W), N2 = v2*sqrt(-2*log(W)/W);
        return sqrt(pow(N1,2)+pow(N2,2));
    }

    int poisson(double lambda) {
        /*
         *assert (lambda > 0. );
        double a = exp( -lambda );
        double b = 1;
        // initialize random number generator
        srand(0);
        long seed = (long) rand() % 100;
        long* idum = &seed;
        for (int i = 0; b >= a; i++ )
            b *= gasdev(idum);
        return i - 1;
        */
        double L = exp(-lambda);
        double p = 1;
        double k = 0;
        srand(0);
        long seed = (long) rand() % 100;
        long* idum = &seed;
        do {
            k++;
            p *= gasdev(idum);
        } while( p > L);
        return (k-1);
    }
};


#endif //STATUTILITY_H