//
// Created by Brandon Alston on 12/9/24.
//

#ifndef TREEMODEL_H
#define TREEMODEL_H

#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include "diffusionProcess.h"
#include "statUtility.h"
#include "constants.h"
#include "matrixUtility.h"
#include "genericUtility.h"
using namespace std;


class TreeModel : public Instrument {
public:
    enum Exercise { European = 'E', American = 'A' };
    enum Type { Call = 'C', Put = 'P' };
    TreeModel();
    virtual TreeModel() {}
    BinomialTreeCRRAmerican(double price, double strike, double rate, double div, double vol, double T, int N, char type);
    TwoVarBinomialTree(double S1, double S2, double strike, double rate, double div1, double div2, double rho, double vol1, double  vol2, double T, int N, char exercise, char type);
    friend class OptionGreeks;
    void setPricingEngine(const Handle<PricingEngine>& engine);
    virtual void performCalculations() const;
    virtual void setupEngine() const = 0; // set up pricing engine
    virtual double calculate() const = 0; // compute price
    /**********************************************************************************
    buildBinomialTreeCRRAmerican : computes the value of an American option using backwards induction in a Cox-Ross-Rubenstein binomial tree model
    [in]: double price : asset price
    double strike : strike price
    double rate : risk-free interest rate
    double div : dividend yield
    double vol : volatility
    double T : time to maturity
    int N : number of time steps
    char exercise : ‘E’uropean or ‘A’merican
    char type : ‘C’all or ‘P’ut
    [out]: double : value of American option
    **********************************************************************************/
    double TreeModel::BinomialCRRAmerican(double price, double strike, double rate, double div, double vol, double T, int N, char type) {
        int i,j;
        double prob; // probability of up movement
        double S[200][200] = {0.0}; // stock price at node i,j
        double c[200][200] = {0.0}; // call price at node i, j
        double a;
        double num = 0.0;
        double up = 0.0;
        double down = 0.0;
        double dt = 0.0;
        dt = T/N; // time step size
        up = exp(vol*sqrt(dt)); // up movement
        down = 1/up; // down movement
        a = exp((rate-div)*dt); // growth rate in prob
        prob = (a - down)/(up - down); // compute stock price at each node
        // initialize call prices
        for (i = 0; i <= N; i++) {
            for (j = 0; j <= i; j++) {
                S[i][j] = price*(pow(up,j))*(pow(down,i-j));
                c[i][j] = 0;
            }
        }
        // compute terminal payoffs
        for (j = N; j >= 0; j--) {
            if (type == ‘C’)
                c[N][j] = max(S[N][j]-strike,0);
            else
                c[N][j] = max(strike-S[N][j],0);
        }
        // work backwards
        for (i = N-1; i >= 0; i--) {
            for (j = i; j >= 0; j--) {
                c[i][j] = exp(-rate*dt)*(prob*(c[i+1][j+1]) + (1- prob)*(c[i+1][j]));
                if (type == ‘C’)
                    c[i][j] = max(S[i][j] – strike,c[i][j]);
                else
                    c[i][j] = max(strike – S[i][j],c[i][j]);
            }
        }
        return c[0][0];
    }

    /**********************************************************************************
    buildTwoVarBinomialTree : computes the value of an American spread option using a 2
    variable binomial tree
    [in]: double S1 : asset price 1
    double S2 : asset price 2
    double strike : strike price of spread option
    double rate : risk-free interest rate
    double div1 : dividend yield of asset 1
    double div2 : dividend yield of asset 2
    double rho : correlation of asset 1 and asset 2
    double vol1 : volatility of asset 1
    double vol2 : volatility of asset 2
    double T : time to maturity
    int N : number of time steps
    char exercise : ‘E’uropean or ‘A’merican
    char type : ‘C’all or ‘P’ut
    [out]: double : value of spread option
    **********************************************************************************/
    double TreeModel::TwoVarBinomialTree (double S1, double S2, double strike, double rate, double div1, double div2, double rho, double vol1, double  vol2, double T, int N, char exercise, char type) {
        double dt = T/N; // time step
        double mu1 = rate – div1 – 0.5*vol1*vol1; // drift for stock 1
        double mu2 = rate – div2 – 0.5*vol2*vol2; // drift for stock 2
        double dx1 = vol1*sqrt(dt); // state step for stock 1
        double dx2 = vol2*sqrt(dt); // state step for stock 2
        double puu, pud, pdu, pdd, dx; // probabilities
        double S1t[100] = { 0.0 }; // array of stock price 1
        double S2t[100] = { 0.0 }; // array of stock price 2
        double C[100][100] = { 0.0 }; int i,j,k; // call price at time step i and node j
        // compute probabilities
        puu = ((dx1*dx2 + (dx2*mu1 + dx1*mu2 + rho*vol1*vol2)*dt)/(4*dx1*dx2));
        pud = ((dx1*dx2 + (dx2*mu1 – dx1*mu2 – rho*vol1*vol2)*dt)/(4*dx1*dx2));
        pdu = ((dx1*dx2 + (-dx2*mu1 + dx1*mu2 – rho*vol1*vol2)*dt)/(4*dx1*dx2));
        pdd = ((dx1*dx2 + (-dx2*mu1 – dx1*mu2 + rho*vol1*vol2)*dt)/(4*dx1*dx2));
        // initialize asset prices at maturity
        S1t[-N] = S1*exp(-N*dx1);
        S2t[-N] = S2*exp(-N*dx2);
        // compute stock prices at each node
        for (j = -N+1; j <= N; j++) {
            S1t[j] = S1t[j-1]*exp(dx1);
            S2t[j] = S2t[j-1]*exp(dx2);
        }
        // compute early exercise payoff at each node
        for (j = -N; j <= N; j += 2) {
            for (k = -N; k <= N; k += 2) {
                if (type == ‘C’)
                    C[j][k] = max(0.0, S1t[j] – S2t[k] – strike);
                else
                    C[j][k] = max(0.0, strike – St1[j] + St2[k]);
            }
        }
        // step back through the tree applying early exercise
        for (i = N-1; i >= 0; i--) {
            for (j = -i; j <= i; j +=2 ) {
                for (k = -i; k <= i; k += 2) {
                    // compute risk-neutral price
                    C[j][k] = exp(-rate*T)*(pdd*C[j-1][k-1] + pud*C[j+1][k-1] + pdu*C[j-1][k+1] + puu*C[j+1][k+1]);
                    if (exercise == ‘A’) {
                        if (type == ‘C’)
                            C[j][k] = max(C[j][k], S1t[j] – S2t[k] – strike);
                        else
                            C[j][k] = max(C[j][k], strike – St1[j] + St2[k]);
                        }
                }
            }
        }
        return C[0][0];
    }

    double TreeModel::TrinomialCRRAmerican(double price, double strike, double vol, double rate, double div, double T, long N, char type) {
        int i, j;
        double pd; // down probability
        double pm; // middle probability
        double pu; // up probability
        double S[250][250]; // stock price at node i, j
        double c[250][250]; // call price at node i,j
        double up = 0.0; // up movement
        double down =0.0; // down movement
        double dt = T/N; // time step
        double drift = rate - div - 0.5*vol*vol; // drift
        pu = 0.33333 + (drift/vol)*sqrt(dt/6);
        pd = 0.33333 - (drift/vol)*sqrt(dt/6);
        pm = 0.33333;
        up = exp(vol*sqrt(3*dt/2));
        down = 1/up;
        // compute stock prices at each node
        for (i = N; i >= 0; i--) {
            for (j = -i; j <= i; j++)
            {
                S[i][j] = price*pow(up,j);
            }
        }
        // compute payoffs at the final time step
        for (j = N; j >= -N; j--) {
            if (type == ‘C’)
                c[N][j] = max(S[N][j] - strike,0);
            else
                c[N][j] = max(strike - S[N][j],0);
        }
        // backwards induction
        for (i=N-1; i >= 0; i--) {
            for (j=i; j >= -i; j--) {
                if (type == ‘C’)
                    c[i][j] = max(exp(-rate*dt)*(pu*c[i+1][j+1] + pm*c[i+1][j] + pd*c[i+1][j-
                    1]), S[i][j] – strike);
                else
                    c[i][j] = max(exp(-rate*dt)*(pu*c[i+1][j+1] + pm*c[i+1][j] + pd*c[i+1][j-
                    1]), strike - S[i][j]);
            }
        }
        return c[0][0];
    }

};


class ConvertibleBond {
public:
    ConvertibleBond();
    virtual ConvertibleBond();
    double calcConvertibleBond(double price, double vol, double rate, double dividend, double T, double principal, double coupon, double frequency, int N, double conversionRatio, double conversionPrice, double creditSpread);
private:
    double S[20][20]; // value of stock price at node i,j
    double V[20][20]; // value of convertible bond at node i,j
    double cp[20][20]; // conversion probability at node i,j
    double creditAdjustedRate; // credit spread at each node i,j
    double call[20][20]; // callable value

    /*********************************************************************************
    calcConvertibleBond
    computes the value of convertible bond with callable provisions
    [in]: double price : stock price
    double vol : stock volatility
    vector<double> rates : contains zero-curve rates
    double dividend : stock dividend yield
    double T : time to maturity of convertible bond
    double principal : par value of bond
    double couponRate : coupon rate of bond
    double frequency : frequency of coupon payments
    int N : number of time steps
    double conversionRatio : conversion ratio
    double conversionPrice : conversion price
    double creditSpread : credit spread of issuer
    map<int,double> callSchedule : call schedule map of times to call prices
    [out] double
    *********************************************************************************/
    double ConvertibleBond::calcConvertibleBond(double price, double vol, vector<double> rates, double dividend, double T, double principal, double couponRate, double frequency, int N, double conversionRatio, double conversionPrice, double creditSpread, map<int,double> callSchedule) {
        int i,j;
        double up = 0.0; // up movement
        double down = 0.0; // down movement
        double interest = 0.0; // interest
        double H = 0.0; // holding value
        double rate = rates[rates.size()-1]; // initial short rate
        double dt = T/N;  // compute time step
        up = exp(vol*sqrt(dt)); // up movement
        down = 1/up; // down movement

        // build CRR stock tree
        for (i = 0; i <= N; i++) {
            for (j = 0; j <= i; j++) {
                S[i][j] = price*(pow(up,j))*(pow(down,i-j));
            }
        }
        interest = principal*coupon*dt; // interest payment
        for (j = N; j >= 0; j--) {
            double payment = principal + principal*coupon*dt;
            if (S[N][j] >= conversionPrice)
                V[N][j] = max(conversionRatio*S[N][j],payment);
            else
                V[N][j] = payment;
            if (V[N][j] == conversionRatio*S[N][j])
                cp[N][j] = 1.0;
            else
                cp[N][j] = 0.0;
        }
        // work backwards
        for (i = N-1; i >= 0; i--) {
            for (j = i; j >= 0; j--) {
                // compute call schedule price at current node
                // in practice, we would want to check that the call date coincides exactly
                // with the time step on the tree. However, we are not using enough time
                // steps for them to coincide (N would have to be substantially increased)
                // and just assume that the bond is callable at each time step
                call[i][j] = callSchedule[i];
                // compute conversion probability
                cp[i][j] = 0.5*(cp[i+1][j+1] + cp[i+1][j]);
                // compute credit adjusted discount rate
                creditAdjustedRate = cp[i][j]*rates[i] + (1- cp[i][j])*creditSpread;
                // compute holding value
                H = 0.5*((V[i+1][j+1] + interest)/(1 + creditAdjustedRate*dt) + (V[i+1][j]
                + interest)/(1 + creditAdjustedRate*dt));
                // check that stock price exceeds conversion price
                if (S[i][j] >= conversionPrice)
                    V[i][j] = max(conversionRatio*S[i][j],min(H,call[i][j] + interest));
                else
                    V[i][j] = min(H,call[i][j] + interest);
            }
        }
        return V[0][0];
    }

};
#endif //TREEMODEL_H

