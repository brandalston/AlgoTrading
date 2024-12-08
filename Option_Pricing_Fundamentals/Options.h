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
    char type_; // option type 'C’all or 'P’ut
    char exercise_; // exercise type 'E’uropean and 'A’merican
    Handle<PricingEngine> engine_; // pricing engine
    OptionGreeks og; // option greeks
    StatUtility util; // statistical utility class
    struct MatrixUtil* mu; // matrix utility class


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
};


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
class BSMOption : public Option {
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
};

// Spread Option
class SpreadOption : public Option {
public:
    SpreadOption() { }
    SpreadOption(Option::Type type, double price1, double price2, double strike, double rate, double vol1, double vol2, double div1, double div2, double  rho, double T, char type, int M, int N);
    virtual~ SpreadOption() {}
    double calcSpreadCallPrice(double price, double strike, double vol, double rate, double div, double T);
    double calcSpreadPutPrice(double vol, double rate, double div, double strike, double price, double T);
    virtual void setupEngine() const { }
    virtual double calculate() const { return value_ ; }
    /*********************************************************************************
        calcMCEuroSpreadOption : computes the value of a European spread option
        [in] double price1 : price of asset 1
        double price2 : price of asset 2
        double strike : strike price
        double rate : risk-free rate
        double vol1 : volatility of asset 1
        double vol2 : volatility of asset 2
        double div1 : dividend yield of asset 1
        double div2 : dividend yield of asset 2
        double rho : correlation of dz1 and dz2
        double T : maturity of option
        char type : option type (C)all or (P)ut
        long M : number of simulations
        long N : number of time steps
        [out] : double : price of spread option
        **********************************************************************************/
    double SpreadOption::calcMCEuroSpreadOption(double price1, double price2, double strike, double rate, double vol1, double vol2, double div1, double div2, double  rho, double T, char type, int M, int N) {
        int i, j;
        double dt = T/N; // size of time step
        double mu1 = (rate - div1 - 0.5*vol1*vol1); // drift for stock price 1
        double mu2 = (rate - div2 - 0.5*vol2*vol2); // drift for stock price 2
        double srho = sqrt(1 - rho*rho); // square root of 1 – rho*rho
        double sum1 = 0.0; // sum of all the call values on stock 1 at time T
        double sum2 = 0.0; // sum of all the call values on stock 2 at time T
        double S1 = 0.0; // stock price 1 
        double S2 = 0.0; // stock price 2
        double deviate1 = 0.0;  // deviate for stock price 1
        double deviate2 = 0.0  // deviate for stock price 2
        double z1 = 0.0; // correlated deviate for stock price 1
        double z2 = 0.0; // correlated deviate for stock price 2
        double CT = 0.0; // option price at maturity
        double SD = 0.0; // standard deviate of price
        double SE = 0.0; // standard error of price
        double value = 0.0; // spread option price
        
        srand(time(0)); // initialize RNG
        long seed = (long) rand() % 100; // generate seed
        long* idum = &seed;
        N = 1; // no path dependency
        for (i = 0; i < M; i++) {
            // initialize prices for each simulation
            S1 = price1;
            S2 = price2;
            for (j = 0; j < N; j++) {
                // generate deviates
                deviate1 = util.gasdev(idum);
                deviate2 = util.gasdev(idum);
                // calculate correlated deviates
                z1 = deviate1;
                z2 = rho*deviate1 + srho*deviate2;
                S1 = S1*exp(mu1*dt + vol1*z1*sqrt(dt));
                S2 = S2*exp(mu2*dt + vol2*z2*sqrt(dt));
            }
            if (type == 'C’)
                CT = max(S1 - S2 - strike, 0);
            else
                CT = max(strike - S1 + S2,0);
            sum1 = sum1 + CT;
            sum2 = sum2 + CT*CT;
        }
        value = exp(-rate*T)*(sum1/M);
        SD = sqrt((sum2 - sum1*sum1/M)*exp(-2*rate*T)/(M-1));
        SE = SD/sqrt(M);
        return value;
    }
};

class AsianOption : public Option {
public:
    AsianOption(double price, double strike, double vol, double rate, double div, double T);
    AsianOption() : value_(0.0) {}
    virtual~ AsianOption() {}
    // modified Black Scholes pricing formula
    double calcBSAsianPrice(double price, double strike, double vol, double rate, double div, double T, char type);
    // calculate arithemic ave. Asian option using Monte Carlo (MCA)
    double calcMCAAsianPrice(double price, double strike, double vol, double rate, double div, double T, char type, int M, int N);
    // calculate geometric ave. Asian option using Monte Carlo (MCG)
    double calcMCGAsianPrice(double price, double strike, double vol, double rate, double div, double T, char type, int M, int N);
    virtual void setupEngine() const { }
    virtual double calculate() const { return value_ ; }
private:
    double volA; double qA; // Arithmetic ave. volatility for modified Black-Scholes formula
    // Arithmetic ave. dividend yield for modified Black-Scholes
    // formula
    double value_; // Asian option price

    /**********************************************************************************
    calcMCGAsianPrice: computes the price of a geometric Asian option using Monte Carlo
    simulation
    [in]: double price : initial stock price
    double strike : strike price
    double vol : volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : time to maturity
    char type : (C)all or (P)ut
    int M : number of simulations
    int N : number of time steps
    [out] double : price of geometric Asian option
    **********************************************************************************/
    double AsianOption::calcMCGAsianPrice(double price, double strike, double vol, double rate, double div, double T, char type, long M, long N) {
        // initialize variables
        int i, j;
        double G = 0.0; // price of geometric average Asian option
        double mu = 0.0; // drift
        double deviate; // normal deviate
        double S = 0.0; // stock price
        double sum = 0.0; // sum of payoffs
        double sum2 = 0.0; // sum of squared payoffs
        double product = 0.0; // product of stock prices
        double payoff = 0.0; // option payoff
        double deltat = 0.0; // step size
        double stddev = 0.0; // standard deviation
        double stderror = 0.0; // standard error
        double deltat = T/N; // compute change in step size
        double mu = rate - div - 0.5*vol*vol; 
        // set output decimal format
        cout.precision(4); 
        srand(time(0)); // initialize RNG
        long seed = (long) rand() % 100; // generate random number generator
        long *idum = &seed; // store address of seed
        // for each simulation
        for (i = 0; i <= M; i++) {
            S = price;
            product = 1;
            for (j = 0; j < N; j++)
            {
                deviate = util.gasdev(idum);
                S = S*exp(mu*deltat + vol*sqrt(deltat)*deviate);
                product *= S;
            }
            // compute geometric average
            G = pow(product,(double)1/N);
            if (type == 'C’)
                payoff = max(G – strike,0);
            else
                payoff = max(strike – G,0);
            sum += payoff;
            sum2 += payoff*payoff;
        }
        value_ = exp(-rate*T)*(sum/M);
        stddev = sqrt((sum2 - sum*sum/M)*exp(-2*rate*T)/(M-1));
        stderror = stddev/sqrt(M);
        return value_;
    }

    /*********************************************************************************
    calcMCAsianPrice : computes the price of an arithmetic Asian option using Monte Carlo simulation
    [in]: double price : initial stock price
    double strike : strike price
    double vol : stock volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : time to maturity
    char type : (C)all or (P)ut
    int M : number of simulations
    int N : number of time steps
    [out]: double : price of arithmetic Asian option
    **********************************************************************************/
    double AsianOption::calcMCAAsianPrice(double price, double strike, double vol, double rate, double div, double T, char type, int M, int N) {
        // initialize variables
        double A = 0.0; double mu = 0.0; // arithmetic average
        // drift
        int i, j;
        double deviate; // normal deviate
        double stddev = 0.0; // standard deviation
        double stderror = 0.0; // standard error
        double S = 0.0; // stock price
        double sum = 0.0; // sum of payoffs
        double sum1 = 0.0; // sum of stock prices
        double sum2 = 0.0; // sum of squared payoffs
        double payoff = 0.0; // payoff of option
        deltat = T/N; // step size
        mu = rate - div - 0.5*vol*vol; // compute drift
        cout.precision(4); // set output decimal format
        srand(time(0)); // initializer RNG
        long seed = (long) rand() % 100; // generate seed
        long *idum = &seed;
        // for each simulation
        for (i = 0; i <= M; i++)
        {
            // reinitialize for each simulation
            S = price;
            sum1 = 0;
            for (j = 0; j <N; j++)
            {
            }
            deviate = util.gasdev(idum);
            S = S*exp(mu*deltat + vol*sqrt(deltat)*deviate);
            sum1 += S;
            A = sum1/N;
            if (type == 'C’)
                payoff = max(A - strike, 0);
            else
                payoff = max(strike – A,0);
            sum += payoff;
            sum2 += payoff*payoff;
        }
        value_= exp(-rate*T)*(sum/M);
        cout << value = "  << value_ <<endl;
        stddev = sqrt((sum2 - sum*sum/M)*exp(-2*rate*T)/(M-1));
        stderror = stddev/sqrt(M);
        cout << "  stddev = " << stddev << "  " << " stderror " << stderror << endl;
        return value_;
    }

    /*********************************************************************************
    calcMCGAsianPrice : computes the price of an geometric Asian option with a
    control variate using Monte Carlo simulation
    [in]: double price : initial stock price
    double strike : strike price
    double vol : stock volatility
    double rate : risk-free rate
    double div : dividend yield
    double T : time to maturity
    char type : (C)all or (P)ut
    int M : number of simulations
    int N : number of time steps
    [out]: double : price of geometic Asian option with a control variate
    **********************************************************************************/
    double AsianOption::calcMCAAsianGCV(double price, double strike, double vol, double
    rate, double div, double T, char type, int M, int N)
    {
        // initialize variables
        int i, j;
        double geo = 0.0; // geometric average
        double ave = 0.0;  // arithmetic average
        double mu = 0.0; // drift
        double stddev = 0.0; // standard deviation
        double stderror = 0.0; // standard error
        double deviate; // standard deviate
        double S = 0.0; // stock price
        double sum = 0.0; // sum of payoffs
        double sum1 = 0.0; // sum of squared payoffs
        double product = 0.0; // product of stock prices
        double payoff = 0.0; // option payoff
        double dt = T/N; // step size
        cout.precision(4); // set output decimal format
        srand(time(0)); // initialize RNG
        long seed = (long) rand() % 100; // generate seed
        long* idum = &seed; // store address of seed
        mu = rate - div - 0.5*vol*vol; // drift
        // simulation
        for (i = 0; i <= M; i++) {
            // initialize for each simulation
            S = price;
            product = 1;
            sum = 0;
            sum1 = 0;
            for (j = 0; j < N; j++) {
                deviate = util.gasdev(idum);
                S = S*exp(mu*deltat + vol*sqrt(dt)*deviate);
                sum = sum + S;
                product *= S;
            }
            ave = sum/N; geo = pow(product,(double)1/N);
            if (type == 'C') {
                payoff = max(0, (ave - strike) - (geo - strike));
            } else {
                payoff = max(0, (strike - ave) - (strike - geo));
            }sum += payoff;
            sum1 += payoff*payoff;
            // calculate arithmetic average
            // calculate geometric average
        }
        value_ = exp(-rate*T)*(sum/M) + calcMCGAsianPrice(price,strike,vol,rate,div,T,'C');
        cout << "value = " << value_ <<endl;
        stddev = sqrt((sum1 - sum*sum/M)*exp(-2*rate*T)/(M-1));
        stderror = stddev/sqrt(M);
        cout << "  stddev = " << stddev << "  " << " stderror = " << stderror << endl;
        return value_;
    }
};


#endif //OPTION_H
