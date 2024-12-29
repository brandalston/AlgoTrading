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


class Tree : public Instrument
{
public:
    enum Exercise { European = 'E', American = 'A' };
    enum Type { Call = 'C', Put = 'P' };
    Tree();
    virtual Tree() {}
};
class BinomialTree : public Tree {
public:
    BinomialTree() { }
    BinomialTree(const Handle<DiffusionProcess>& process,
    const TimeGrid& timeGrid,
    bool isPositive = false);
    double dx(Size i) const { return dx_[i]; }
    double underlying(Size i, Size index) const;
    const TimeGrid& timeGrid() const { return timeGrid_; }
    inline int descendant(int i, int index, int branch) const {
        return branchings_[i]->descendant(index, branch);
    }
    inline double probability(int i, int j, int b) const {
        return branchings_[i]->probability(j, b);
    }
    inline int size(int i) const {
        if (i==0)
            return 1;
        const std::vector<int>& k = branchings_[i-1]->k_;
        int jMin = *std::min_element(k.begin(), k.end()) - 1;
        int jMax = *std::max_element(k.begin(), k.end()) + 1;
        return jMax - jMin + 1;
    }
    double underlying(int i, int index) const {
        if (i==0) return x0_;
        const std::vector<int>& k = branchings_[i-1]->k_;
        int jMin = *std::min_element(k.begin(), k.end()) - 1;
        return x0_ + (jMin*1.0 + index*1.0)*dx(i);
    }
protected:
    std::vector<Handle<TrinomialBranching> > branchings_;
    double x0_;
    std::vector<double> dx_; // vector of step sizes
    TimeGrid timeGrid_;

    BinomialTreeCRRAmerican(double price, double strike, double rate, double div, double vol, double T, int N, char type);
    TwoVarBinomialTree(double S1, double S2, double strike, double rate, double div1, double div2, double rho, double vol1, double  vol2, double T, int N, char exercise, char type);
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
    double BinomialTree::BinomialCRRAmerican(double price, double strike, double rate, double div, double vol, double T, int N, char type) {
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
            if (type == 'C')
                c[N][j] = max(S[N][j]-strike,0);
            else
                c[N][j] = max(strike-S[N][j],0);
        }
        // work backwards
        for (i = N-1; i >= 0; i--) {
            for (j = i; j >= 0; j--) {
                c[i][j] = exp(-rate*dt)*(prob*(c[i+1][j+1]) + (1- prob)*(c[i+1][j]));
                if (type == 'C')
                    c[i][j] = max(S[i][j] - strike, c[i][j]);
                else
                    c[i][j] = max(strike - S[i][j], c[i][j]);
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
    double BinomialTree::TwoVarBinomialTree (double S1, double S2, double strike, double rate, double div1, double div2, double rho, double vol1, double  vol2, double T, int N, char exercise, char type) {
        double dt = T/N; // time step
        double mu1 = rate - div1 - 0.5*vol1*vol1; // drift for stock 1
        double mu2 = rate - div2 - 0.5*vol2*vol2; // drift for stock 2
        double dx1 = vol1*sqrt(dt); // state step for stock 1
        double dx2 = vol2*sqrt(dt); // state step for stock 2
        double puu, pud, pdu, pdd, dx; // probabilities
        double S1t[100] = { 0.0 }; // array of stock price 1
        double S2t[100] = { 0.0 }; // array of stock price 2
        double C[100][100] = { 0.0 }; int i,j,k; // call price at time step i and node j
        // compute probabilities
        puu = ((dx1*dx2 + (dx2*mu1 + dx1*mu2 + rho*vol1*vol2)*dt)/(4*dx1*dx2));
        pud = ((dx1*dx2 + (dx2*mu1 - dx1*mu2 - rho*vol1*vol2)*dt)/(4*dx1*dx2));
        pdu = ((dx1*dx2 + (-dx2*mu1 + dx1*mu2 - rho*vol1*vol2)*dt)/(4*dx1*dx2));
        pdd = ((dx1*dx2 + (-dx2*mu1 - dx1*mu2 + rho*vol1*vol2)*dt)/(4*dx1*dx2));
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
                if (type == 'C')
                    C[j][k] = max(0.0, S1t[j] - S2t[k] - strike);
                else
                    C[j][k] = max(0.0, strike - S1t[j] + S2t[k]);
            }
        }
        // step back through the tree applying early exercise
        for (i = N-1; i >= 0; i--) {
            for (j = -i; j <= i; j +=2 ) {
                for (k = -i; k <= i; k += 2) {
                    // compute risk-neutral price
                    C[j][k] = exp(-rate*T)*(pdd*C[j-1][k-1] + pud*C[j+1][k-1] + pdu*C[j-1][k+1] + puu*C[j+1][k+1]);
                    if (exercise == 'A') {
                        if (type == 'C')
                            C[j][k] = max(C[j][k], S1t[j] - S2t[k] - strike);
                        else
                            C[j][k] = max(C[j][k], strike - S1t[j] + S2t[k]);
                        }
                }
            }
        }
        return C[0][0];
    }

};

class TrinomialTree : public Tree
{
public:
    TrinomialTree() { }
    TrinomialTree(const Handle<DiffusionProcess>& process,
    const TimeGrid& timeGrid,
    bool isPositive = false);
    double dx(Size i) const { return dx_[i]; }
    double underlying(Size i, Size index) const;
    const TimeGrid& timeGrid() const { return timeGrid_; }
    inline int descendant(int i, int index, int branch) const {
        return branchings_[i]->descendant(index, branch);
    }
    inline double probability(int i, int j, int b) const {
        return branchings_[i]->probability(j, b);
    }
    inline int size(int i) const {
        if (i==0)
            return 1;
        const std::vector<int>& k = branchings_[i-1]->k_;
        int jMin = *std::min_element(k.begin(), k.end()) - 1;
        int jMax = *std::max_element(k.begin(), k.end()) + 1;
        return jMax - jMin + 1;
    }
    double underlying(int i, int index) const {
        if (i==0) return x0_;
        const std::vector<int>& k = branchings_[i-1]->k_;
        int jMin = *std::min_element(k.begin(), k.end()) - 1;
        return x0_ + (jMin*1.0 + index*1.0)*dx(i);
    }
protected:
    std::vector<Handle<TrinomialBranching> > branchings_;
    double x0_;
    std::vector<double> dx_; // vector of step sizes
    TimeGrid timeGrid_;

    TrinomialTree::TrinomialTree(const Handle<DiffusionProcess>& process, const TimeGrid& timeGrid, bool isPositive) : Tree(timeGrid.size()), dx_(1, 0.0), timeGrid_(timeGrid) {
        x0_ = process->x0();
        int nTimeSteps = timeGrid.size() - 1;
        int jMin = 0;
        int jMax = 0;
        for (int i = 0; i < nTimeSteps; i++)
        {
            Time t = timeGrid[i];
            Time dt = timeGrid.dt(i);
            // variance must be independent of x
            double v2 = process->variance(t, 0.0, dt);
            double v = sqrt(v2);
            dx_.push_back(v*sqrt (3.0));
            Handle<TrinomialBranching> branching(new TrinomialBranching());
            for (int j = jMin; j <= jMax; j++)
            {
                double x = x0_ + j*dx_[i];
                double m = process->expectation(t, x, dt);
                int temp = (int)floor ((m-x0_)/dx_[i+1] + 0.5);
                if (isPositive)
                {
                    while (x0_+(temp-1)*dx_[i+1] <= 0)
                        temp++;
                }
                branching->k_.push_back(temp);
                double e = m - (x0_ + temp*dx_[i+1]);
                double e2 = e*e;
                double e3 = e*sqrt (3.0);
                branching->probs_[0].push_back((1.0 + e2/v2 - e3/v)/6.0);
                branching->probs_[1].push_back((2.0 - e2/v2)/3.0);
                branching->probs_[2].push_back((1.0 + e2/v2 + e3/v)/6.0);
            }
            branchings_.push_back(branching);
            const std::vector<int>& k = branching->k_;
            jMin = *std::min_element(k.begin(), k.end()) - 1;
            jMax = *std::max_element(k.begin(), k.end()) + 1;
        }
    }
};

/*****************************************************************************
class TrinomialBranching : Recombining trinomial tree class
This class defines a recombining trinomial tree approximating a diffusion. The
diffusion term of the SDE must be independent of the underlying process.
*****************************************************************************/
class TrinomialBranching {
public:
    TrinomialBranching() : probs_(3) {}
    virtual TrinomialBranching() {}
    inline Size descendant(Size index, Size branch) const {
        return (k_[index] - jMin()) - 1 + branch;
    }
    inline double probability(Size index, Size branch) const {
        return probs_[branch][index];
    }
    inline int jMin() const {
        return *std::min_element(k_.begin(), k_.end()) - 1;
    }
private:
    friend class TrinomialTree;
    std::vector<int> k_; // branch k
    std::vector<std::vector<double> > probs_; // branching probabilities
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
                H = 0.5*((V[i+1][j+1] + interest)/(1 + creditAdjustedRate*dt) + (V[i+1][j]+ interest)/(1 + creditAdjustedRate*dt));
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

