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

    double MonteCarloMethod::MonteCarloSobol(double price, double strike, double vol, double rate, double div, double T, char type, long N, long M) {
        int i, j;
        double sum1 = 0.0; // sum of payoffs
        double sum2 = 0.0; // sum of squared payoffs
        double value = 0.0; // stores value of option for each simulation
        double S1 = price; // stock price for +deviate
        double S2 = price; // stock price for -deviate
        double lnS1 = log(price); // log of the initial stock price for +deviate
        double lnS2 = log(price); // log of the initial stock price for -deviate
        double SD; // standard deviation
        double SE; // standard error
        long dim = N; // dimension of Sobol sequence
        double dt = T/N; // time step
        double mu = rate - div - 0.5*vol*vol;// drift
        double rands[5]; // stores random variables
        cout.precision(4); // output decimal format precision
        int cnt = 0; // counter
        struct sobolp sp; // Sobol sequence structure
        srand(time(0)); // initialize RNG
        long seed = (long) rand() % 100; // generate seed
        // initialize Sobol sequnce
        sobolp_init(&sp,dim,seed);
        for (i = 0; i < M; i++) {
            // initalize stock price for the next simulation
            lnS1 = log(price);
            lnS2 = log(price);
            for (j = 0; j < N; j++) {
                // generate Sobol samples
                sobolp_generateSamples(&sp,rands);
                // generate path and antithetic path
                lnS1 = lnS1 + mu*dt + vol*sqrt(dt)*rands[cnt];
                lnS2 = lnS2 = mu*dt + vol*sqrt(dt)*(-rands[cnt]);
                // keep track of Sobol number to use
                if ((cnt + 1) % N == 0)
                    cnt = 0;
                else
                    cnt++;
            }
            // convert back to lognormal random variables
            S1 = exp(lnS1);
            S2 = exp(lnS2);
            if (type == 'C') {
                value = 0.5*(max(0, S1 - strike) + max(0, S2 - strike));
            } else {
                value = 0.5*(max(0, strike - S1) + max(0,strike - S2));
            }
            sum1 = sum1 + value;
            sum2 = sum2 + value*value;
        }
        // compute standard deviation
        SD = sqrt((exp(-2*rate*T)/(M-1))*(sum2 - (sum1*sum1)/M));
        cout << "stddev " << " " < SD < endl;
        // compute standard error
        SE = SD/sqrt(M);
        cout << "stdderr " << " " << SE << endl;
        return exp(-rate*T)*sum1/M;
    }
};

#endif //STATUTILITY_H
