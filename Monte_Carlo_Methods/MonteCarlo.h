//
// Created by Brandon Alston on 11/13/24.
//

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include "diffusionProcesses.h"
#include "statUtility.h"
#include "constants.h"
#include "matrixUtility.h"
using namespace std;

class MonteCarloMethod {
    /*********************************************************************************
    MonteCarloSobol : values a European call option using Faure sequence for variance reduction
    [in]: double price : asset price
    double strike : strike price
    double vol : volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : option maturity
    char type : type of option
    long N : number of time steps
    long M : number of simulations
    [out] : double : call price
    **********************************************************************************/
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

    /*********************************************************************************
    MonteCarloFaureQuasiRandom : values a European call option using Faure sequence for variance reduction
    [in]: double S : asset price
    double X : strike price
    double vol : volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : option maturity
    long N : number of time steps
    long M: : number of simulations
    [out] : double : callValue
    **********************************************************************************/
    double MonteCarloMethod::MonteCarloFaureQuasiRandom (double S, double X, double vol, double rate, double div, double T, long N, long M) {
        int i, j, k;
        double dt = T/N; double mudt = (rate - div - 0.5*vol*vol)*dt; double voldt = vol*sqrt(dt); double sum = 0.0;
        double sum1 = 0.0;
        double lnSt, lnSt1, St, St1;
        double lnS = log(S);
        double deviate = 0.0;
        double callValue = 0.0;I
        double SD = 0.0; double SE = 0.0; vector<double> x; cout.setf(ios::showpoint);
        cout.precision(3);
        // step step
        // drift
        // diffusion term
        // standard deviation
        // standard error
        // stores Faure sequence
        k = 0;
        for (i = 1; i <= M; i++)
        {
            // generate Faure sequence
            x = generateFaure(N,M);
            // initialize log asset prices for next simulation path
            lnSt = lnS;
            lnSt1 = lnS;
            for (j = 0; j < N; j++)
            {
                // get standard normal deviate using polar rejection method
                deviate = util.polarRejection(x[j],k);
                nSt = lnSt + mudt + voldt*deviate;
                // compute antithetic
                lnSt1 = lnSt1 + mudt + voldt*(-deviate);
                // increment index to retrieve deviate stored in vector Y in polar rejection
                method
                k++;
            }
            St = exp(lnSt);
            St1 = exp(lnSt1);
            callValue = 0.5*(max(0, St - X) + max(0,St1-X));
            sum = sum + callValue;
            sum1 = sum1 + callValue*callValue;
        }
        callValue = exp(-rate*T)*(sum/M)
        SD = sqrt(exp(-2*rate*T)*(sum1/M) - callValue*callValue);
        cout << "stdev = " << SD << endl;
        SE = SD/sqrt(M-1);
        cout << "stderr = " << SE << endl;
        return callValue;
    }


    /*********************************************************************************
     *MonteCarloAntithetic : values a European call option using antithetic variates
    [in]: double price : asset price
    double strike : strike price
    double vol : volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : option maturity
    long N : number of time steps
    long M : number of simulations
    [out] double value call value
    **********************************************************************************/
    double MonteCarloMethod::MonteCarloAntithetic (double price, double strike, double vol, double rate, double div, double T, long M, long N, char type) {
        int i, j;
        double deviate; double sum1 = 0.0; double sum2 = 0.0; double value = 0.0; double S1 = price; double S2 = price; double lnS1 = log(price); // standard normal deviate
        // sum of payoffs
        // sum of squared payoffs
        // value of option
        // stock price for +deviate
        // stock price for -deviate
        // log of the initial stock price for
        // +deviate
        double lnS2 = log(price); // log of the initial stock price for
        // -deviate
        double SD; // standard deviation
        double SE; // standard error
        double deltat = (double) T/N; // time step
        double mu = rate - div - 0.5*vol*vol; // drift
        srand(time(0)); // initialize RNG
        long seed = (long) rand() % 100; // generate seed
        long* idum = &seed; // store seed address
        cout.setf(ios::showpoint);
        cout.precision(4);
        for (i = 0; i < M; i++)
        {
            // initalize stock price for the next simulation
            lnS1 = log(price);
            lnS2 = log(price);
            for (j = 0; j < N; j++)
            {
                deviate = util.gasdev(idum);
                // simulate paths
                72 MONTE CARLO SIMULATION
                lnS1 = lnS1 + mu*deltat + vol*sqrt(deltat)*deviate;
                lnS2 = lnS2 + mu*deltat + vol*sqrt(deltat)*(-deviate);
            }
            // convert back to lognormal random variables
            S1 = exp(lnS1);
            S2 = exp(lnS2);
            if (type == ‘C’)
                value = 0.5*(max(0, S1 - strike) + max(0,S2 - strike));
            else // if put
                value = 0.5*(max(0, strike - S1) + max(0, strike - S2));
            sum1 = sum1 + value;
            sum2 = sum2 + value*value;
        }
        value = exp(-rate*T)*sum1/M
        cout << "value = " << value << endl;
        // compute standard deviation
        SD = sqrt((exp(-2*rate*T)/(M-1))*(sum2 - (sum1*sum1)/M));
        cout << " stdev = " << SD << endl;
        // compute standard error
        SE = SD/sqrt(M);
        cout << " stderr = " << SE << endl;
        return value;
    }
    /**********************************************************************************
    dynamicReplication : synthetically replicates option using stock and money
    market account
    [in]: double price : stock price
    double strike : strike price
    double vol : volatility
    double rate : interest rate
    double div : dividend yield
    double T : option maturity
    char type : ‘C’all or ‘P’ut
    long M : number of simulations
    long N : number of time steps
    [out] double : synthetic option price
    **********************************************************************************/
    double MonteCarloMethod::dynamicReplication(double price, double strike, double vol, double rate, double div, double T, char type, long M, long N) {
        // initialize variables
        int i, j;
        double S = 0.0; // stock price
        double lnS; // log of S
        double delta; // delta of option
        double totalStockShares = 0; // total shares of stock
        double totalMMAShares = 0.0; // total number of MMA shares
        double numShares = 1.0; // number of shares bought or sold at time t
        double numMMA = 0.0; // number of MMA shares bought at time t
        double MMAValue = 0.0; // value of money market account at time t
        double totalMMAValue; // = MMAValue*totalMMAShares
        double d1 = 0.0; // used to calculate delta
        double portValue = 0.0; // portfolio value
        double deviate = 0.0; // normal deviates used for Monte Carlo
        double temp = 0.0; // temp variable to hold delta value
        double totalStockValue = 0.0; // total stock value
        long seed = -1; // initial seed value
        long* idum = 0; // used for gasdev function
        double dt = T/M; // step size
        double mu = 0.0; // drift
        StatUtility util;
        // initial states
        d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*(T))/(vol*sqrt(T));
        delta = util.normalCalc(d1);
        numShares = delta; // number of shares
        totalStockValue = numShares*price;
        MMAValue = numShares*price; // initialize value of money market account
        numMMA = numShares;
        totalMMAValue = numMMA*price;
        totalMMAShares = numMMA;
        totalStockShares = numShares;
        temp = delta;
        portValue = totalStockValue - totalMMAValue;
        srand(unsigned(0));
        seed = (long) rand() % 100;
        idum = &seed;
        for (i = 0; i < M; i++) {
            // initialize starting price
            lnS = log(price);
            // do simulations on each path
            for (j = 0; j < N; j++) {
                deviate = util.gasdev(idum);
                lnS = lnS + (rate - div - 0.5*vol*vol)*dt + vol*sqrt(dt)*deviate;
            }
            S = exp(lnS);
            MMAValue = MMAValue*exp(rate*dt);
            // compute current delta
            if (i != M-1) {
                d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*(T-i*dt))/(vol*sqrt(T-
                i*dt));
                delta = util.normalCalc(d1);
            }
            else
                delta = 1.0;
            // adjust total delta
            delta = delta - temp;
            if (S >= price) {
                // buy shares
                temp = delta;
                numShares = delta;
                totalStockShares = totalStockShares + numShares;
                totalStockValue = totalStockShares*price;
                // finance purchase of stock by selling shares (borrowing) from MMA
                numMMA = numShares;
                totalMMAShares = totalMMAShares + numMMA;
                MMAValue = MMAValue + numMMA*price;
                totalMMAValue = MMAValue*totalMMAShares;
                portValue = totalStockValue - totalMMAValue;
            }
            else {
                // sell shares
                temp = delta;
                numShares = delta;
                totalStockShares = totalStockShares - numShares;
                totalStockValue = totalStockShares*price;
                // buy back the money market shares shorted
                numMMA = -numShares;
                totalMMAShares = totalMMAShares + numMMA;
                MMAValue = MMAValue + numMMA*price;
                totalMMAValue = MMAValue*totalMMAShares;
                portValue = totalStockValue - totalMMAValue;
            }
        }
        std::cout << "final cost: " << totalMMAValue - totalStockValue << endl;
        return totalMMAValue – totalStockValue;
    }

    /**********************************************************************************
MonteCarloADG : values a European Call option using Monte Carlo with antithetic,
delta, and gamma control variates. Adapted from Clewlow and Strickland (1998a)
[in]: double S : asset price
double X : strike price
double vol : volatility
double rate : risk-free rate
double div : dividend yield
double T : option maturity
long N : number of time steps
long M : number of simulations
[out] : double callValue
**********************************************************************************/
    double MonteCarloMethod::MonteCarloADG(double S, double X, double vol, double rate,
    double div, double T, long N, long M)
    {
        int i, j;
        double dt = T/N;
        double mudt = (rate - div - 0.5*vol*vol)*dt;
        double voldt = vol*sqrt(dt);
        double erddt = exp((rate - div)*dt); // helps compute E[Si] efficiently
        double egamma = exp((2*(rate - div) // helps compute gamma control variate
        + vol*vol)*dt)-2*erddt + 1;
        double beta1 = -1; // fixed beta coefficent on delta control
        // variate
        double beta2 = -0.5; // fixed gamma coefficient of gamma control
        // variate
        double sum = 0.0; // summation of call values
        double sum1 = 0.0; // summation of squared call values
        double t; // current time
        double St, St1; // stock prices at current time
        double Stn, Stn1; // stock prices at next time step
        double CT; // call value at maturity at end of
        // simulation path
        double cv1; // delta control variate
        double cv2; // gamma control variate
        double delta, gamma; // delta and gamma of positive antithetic
        double delta1, gamma1; // delta and gamma of negative antithetic
        double deviate; // standard deviate
        double SD; // standard deviation
        double SE; // standard error
        srand(time(0)); // initialize RNG
        long seed = (long) rand() % 100; // seed for random number generator
        long *idum = &seed; // used for generating standard normal
        // deviate
        double callValue; // call value
        cout.setf(ios::showpoint); // output format
        cout.precision(4); // set output decimal precision
        for (i = 1; i <= M; i++)
        {
            // initialize variables for simulation
            St = S;
            St1 = S;
            cv1 = 0;
            cv2 = 0;
            for (j = 1; j <= N; j++)
            {
                // compute hedge sensitivities
                t = (j-1)*dt;
                delta = og.calcDelta(St,X,rate,div,vol,T,t);
                delta1 = og.calcDelta(St1,X,rate,div,vol,T,t);
                gamma = og.calcGamma(St,X,rate,div,vol,T,t);
                gamma1 = og.calcGamma(St1,X,rate,div,vol,T,t);
                // generate gaussian deviate
                deviate = util.gasdev(idum);
                // evolve asset price
                Stn = St*exp(mudt + voldt*deviate);
                Stn1 = St1*exp(mudt + voldt*(-deviate));
                // accumulate control deviates
                cv1 = cv1 + delta*(Stn - St*erddt) + delta1*(Stn1 - St1*erddt);
                cv2 = cv2 + gamma*((Stn - St)*(Stn - St) - pow(St,2*egamma))
                + gamma1*((Stn1 - St1)*(Stn1 - St1) - pow(St1,2*egamma));
                St = Stn;
                St1 = Stn1;
            }
            // compute value with antithetic and control variates
            CT = 0.5*(max(0,St - X) + max(0, St1 - X) + beta1*cv1 + beta2*cv2);
            sum = sum + CT;
            sum1 = sum1 + CT*CT;
        }
        callValue = exp(-rate*T)*(sum/M);
        cout << “value = ” << callValue << endl;
        SD = sqrt((sum1 - sum1*sum1/M)*exp(-2*rate*T)/(M-1));
        cout << “stddev = ” << SD << endl;
        2.7 Hedge Control Variates 81
        SE = SD/sqrt(M);
        cout << “stderr = ” << SE << endl;
        return callValue;
    }


    /*********************************************************************************
    JumpModel : values a European call option using a jump-diffusion
    process
    [in]: double price: : asset price
    double strike : strike price
    double vol : volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : option maturity
    int N : number of time steps
    int M : number of simulations
    double lambda : rate (intensity) of jumps
    double kappa : average jump sized measured as a proportional increase in
    the stock price
    [out] : double callValue
    **********************************************************************************/
    double MonteCarlo::JumpDiffusion(double price, double strike, double vol, double
    rate, double div, double T, int M, int N, double lambda, double kappa)
    {
        int i, j;
        double dt; double deviate = 0.0; double deviate1 = 0.0; double payoff = 0.0; double sum = 0.0; double S = price; double mu = rate - div - lambda*kappa; long seed; long* idum = 0; StatUtility util; // time increment, (yrs)
        // standard normal deviate
        // Poisson deviate
        // option payoff
        // sum payoffs
        // store stock price
        // expected return
        // seed for random number generator
        // identifies address of seed
        // statistic utility class
        srand(time(0)); seed = (long) rand() % 100; idum = &seed;
        dt = T/N; // initialize random number generator
        // generate seed
        // time step
        for (i = 0; i < M; i++)
        {
            // initialize stock price for each simulation
            S = price;
            for(j = 0; j < N; j++)
            {
                deviate = util.gasdev(idum); deviate1 = util.poisson(lambda); // generate gaussian deviate
                // generate Poisson deviate
                S = S*exp(mu*dt+ vol*sqrt(dt)*deviate + sqrt(dt)*deviate1);
            }
            payoff = max(S – strike, 0);
            sum += payoff;
        }
        return exp(-rate*T)*(sum/M);
    }











};
#endif //MONTECARLO_H
