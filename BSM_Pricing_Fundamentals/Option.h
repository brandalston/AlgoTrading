//
// Created by Brandon Alston on 11/2/24.
//

#ifndef OPTION_H
#define OPTION_H

#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include "diffusionProcesses.h"
#include "statUtility.h"
#include "constants.h"
#include "matrixUtility.h"
using namespace std;

// Hull’s approximation of the cumulative normal distribution
double StatUtility::normalCalc(double d) {
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

};

class Option : public Instrument {
public:
    enum Exercise { European = 'E', American = 'A' };
    enum Type { Call = 'C', Put = 'P' };
    Option();
    Option(double price, double strike, double vol, double rate, double div, double T, char type, char exercise);
    Option(const Handle<PricingEngine>& engine);
    virtual∼Option() {}
    friend class OptionGreeks;
    void setPricingEngine(const Handle<PricingEngine>& engine);
    virtual void performCalculations() const;
    virtual void setupEngine() const = 0; // set up pricing engine
    virtual double calculate() const = 0; // compute price

    // option greeks
    class OptionGreeks {
    public:
        StatUtility util; // statistical utility class
        OptionGreeks() {}
        double calcVega(double price, double strike, double rate, double div, double vol, double T, double t);
        double calcDelta(double price, double strike, double rate, double div, double vol, double T, char type);
        double calcGamma(double price, double strike, double rate, double div, double vol, double T);
        double calcRho(double price, double strike, double rate, double div, double vol, double T, char type);
        double calcTheta(double price, double strike, double rate, double div, double vol, double T, char type);
    private:
        // Greek sensitivities
        double delta; // delta
        double gamma; // gamma
        double theta; // theta
        double vega; // vega
        double rho; // rho
    };
protected:
    double strike_; // strike price
    double rate_; // interest rate
    double T_; // maturity
    double price_; // underlying asset
    double vol_; // volatility
    double dividend_; // dividend yield
    char type_; // option type ‘C’all or ‘P’ut
    char exercise_; // exercise type ‘E’uropean and ‘A’merican
    Handle<PricingEngine> engine_; // pricing engine
    OptionGreeks og; // option greeks
    StatUtility util; // statistical utility class
    struct MatrixUtil* mu; // matrix utility class
};

/*
// default constructor
Option::Option() : price_(50.0), strike_(50.0), rate_(0.06), dividend_(0.0), T_(1), type_('C'), exercise_('E'){}
// overloaded constructor
Option::Option(double price, double strike, double vol, double rate, double div, double T, char type, char exercise)
: price_(price), strike_(strike), vol_(vol), rate_(rate), dividend_(div), T_(T), type_(type), exercise_(exercise) {}
*/

double Option::OptionGreeks::calcDelta(double price, double strike, double rate, double div, double vol, double T, char type) {
    double d1, delta;
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    if (type == 'C') {
        delta = exp(-div*T)*util.normalCalcPrime(d1);
    } else {
        delta = exp(-div*T)*(util.normalCalc(d1) - 1);
    }
    return delta;
}

double Option::OptionGreeks::calcVega(double price, double strike, double rate, double div, double vol, double T, double t) {
    double d1, vega, normalPrime;
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    normalPrime = util.normalCalcPrime(d1);
    vega = (normalPrime*exp(-div*T))*price*sqrt(T);
    return vega;
}

double Option::OptionGreeks::calcGamma(double price, double strike, double rate, double div, double vol, double T) {
    double d1, gamma, normalPrime;
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    normalPrime = util.normalCalcPrime(d1);
    gamma = (normalPrime*exp(-div*T))/(price*vol*sqrt(T));
    return gamma;
}

double Option::OptionGreeks::calcRho(double price, double strike, double rate, double div, double vol, double T, char type) {
    double d1, d2;
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    d2 = d1 - vol*sqrt(T);
    double rho = 0.0;
    if (type == 'C') {
        rho = strike*T*exp(-rate*T)*util.normalCalc(d2);
    } else {
        rho = -strike*T*exp(-rate*T)*util.normalCalc(-d2);
    }
    return rho;
}

double Option::OptionGreeks::calcTheta(double price, double strike, double rate, double div, double vol, double T, char type) {
    double d1, d2
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    d2 = d1 - vol*sqrt(T);
    double theta = 0.0;
    if (type == 'C') {
        theta = (-price*util.normalCalc(d1)*vol*exp(-div*T))/(2*sqrt(T)) +
        div*price*util.normalCalc(d1)*exp(-div*T) -
        rate*strike*exp(-rate*T)*util.normalCalc(d2);
    } else {
        theta = (-price*util.normalCalc(d1)*vol*exp(-div*T))/(2*sqrt(T)) -
        div*price*util.normalCalc(-d1)*exp(-div*T) +
        rate*strike*exp(-rate*T)*util.normalCalc(-d2);
    }
    return theta;
}

// overloaded constructor
Option::Option(const Handle<PricingEngine>& engine) : engine_(engine) {
        QL_REQUIRE(!engine_.isNull(), "Option::Option : null pricing engine not allowed");
    }

void Option::setPricingEngine(const Handle<PricingEngine>& engine) {
    QL_REQUIRE(!engine.isNull(), "Option::setPricingEngine : null pricing engine not allowed");
    engine_ = engine;
    // this will trigger recalculation and notify observers
    update();
    setupEngine();
}

void Option::performCalculations() const {
    setupEngine();
    engine_->calculate();
    const OptionValue* results = dynamic_cast<const OptionValue*>(engine_->results());
    QL_ENSURE(results != 0, "Option::performCalculations : no results returned from option pricer");
    NPV_ = results->value;
}

// Vanilla option (no discrete dividends, no barriers) on a single asset
class VanillaOption : public Option {
public:
    VanillaOption() { }
    VanillaOption(double price, double strike, double rate, double div, double vol, double T, Type type, Exercise exercise, const Handle<PricingEngine>&engine);
    double impliedVolatility(double targetValue, double accuracy = 1.0e-4, Size maxEvaluations = 100, double minVol = 1.0e-4, double maxVol = 4.0)const;
    double delta() const; //get delta
    double gamma() const; //get gamma
    double theta() const; //get theta
    double vega() const;  //get vega
    double rho() const;   //get rho
protected:
    void setupEngine() const;
    void performCalculations() const;
    virtual double calculate() const { return NPV_; }
    Date exerciseDate_; // exercise Date
    RelinkableHandle<TermStructure> riskFreeRate; // spot rate term structure
    // results
    mutable double delta_, gamma_, theta_, vega_, rho_, dividendRho_;
    // arguments
    Type type_;
    Exercise exercise_;
    double underlying_; // underlying price
    double strike_; // strike price
    double dividendYield_; // dividend yield
    double riskFreeRate_; // spot risk-free rate
    double maturity_; // time to maturity (years)
    double volatility_; // volatility
private:
    // helper class for implied volatility calculation
    class ImpliedVolHelper : public ObjectiveFunction {
    public:
        StatUtility util;
        ImpliedVolHelper(const Handle<PricingEngine>& engine, double targetValue);
        map<int,double> calcImpliedVols(double price, vector<double> opPrices, vector<int>strikes, double rate, double dividend, double T, Type type);
        map<pair<double,int>,double> calcImpliedSurface(double price, vector<double> opPrices, vector<int>strikes, vector<double> T, map<double,double> rates, double dividend, Type type);
        double operator()(double x) const;
    private:
        Handle<PricingEngine> engine_;
        double targetValue_;
        const OptionValue* results_;
    };
};

// Black-Scholes-Merton Option
class BSMOption : public VanillaOption {
public:
    BSMOption() { }
    BSMOption(Option::Type type, double underlying, double strike, double dividendYield, double riskFreeRate, double residualTime, double volatility);
    virtual∼ BSMOption() {}
    // modifiers
    virtual void setVolatility(double newVolatility);
    virtual void setRiskFreeRate(double newRate);
    virtual void setDividendYield(double newDividendYield);
    double calcBSCallPrice(double price, double strike, double vol, double rate, double div, double T);
    double calcBSPutPrice(double vol, double rate, double div, double strike, double price, double T);
protected:
    Option::Type type_;
    Option::Exercise exercise_;
    double underlying_;
    double strike_;
    double dividendYield_;
    double riskFreeRate_;
    double residualTime_;
    double volatility_;
    double value_;
};

double BSMOption::calcBSCallPrice(double vol, double rate, double div, double strike, double price, double T) {
    double prob1;
    double prob2;
    double d1, d2;
    double callprice;
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    d2 = d1 – vol*sqrt(T);
    prob1 = util.normalCalc(d1);
    prob2 = util.normalCalc(d2);
    callprice = price*exp(-div*T)*prob1 - strike*exp(-rate*T)*prob2;
    return callprice;
}

double BSMOption::calcBSPutPrice(double vol, double rate, double div, double strike, double price, double T) {
    double prob1;
    double prob2;
    double putprice;
    double d1, d2;
    d1 = (log(price/strike) + (rate - div + (vol)*(vol)/2)*T)/(vol*sqrt(T));
    d2 = d1 - vol*sqrt(T);
    prob1 = util.normalCalc(-d1);
    prob2 = util.normalCalc(-d2);
    putprice = strike*exp(-rate*T)*prob2 - price*exp(-div*T)*prob1;
    return putprice;
}

#endif //OPTION_H
