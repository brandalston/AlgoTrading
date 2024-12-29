//
// Created by Brandon Alston on 11/2/24.
//

#ifndef DIFFUSIONPROCESSES_H
#define DIFFUSIONPROCESSES_H

#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <cmath>
#include "diffusionProcess.h"
#include "statUtility.h"
#include "constants.h"
#include "matrixUtility.h"
using namespace std;

typedef double Time;
typedef double Rate;
/**********************************************************************************
Abstract Instrument class
This class is purely abstract and defines the interface of concrete instruments
which will be derived from this one. It implements the Observable interface
**********************************************************************************/
class Instrument : public Patterns::Observer, public Patterns::Observable {
public:
    Instrument(const string& isinCode, const string& description) : NPV_(0.0), isExpired_(false), isinCode_(isinCode), description_(description), calculated(false) {}
    virtual Instrument() {}
    // inline definitions
    // returns the ISIN code of the instrument, when given.
    inline string isinCode() const {
        return isinCode_;
    }
    // returns a brief textual description of the instrument.
    inline string description() const {
        return description_;
    }
    // returns the net present value of the instrument.
    inline double NPV() const {
        calculate();
        return (isExpired_ ? 0.0 : NPV_);
    }
    // returns whether the instrument is still tradable.
    inline bool isExpired() const {
        calculate();
        return isExpired_;
    }
    // updates dependent instrument classes
    inline void update() {
        calculated = false;
        notifyObservers();
    }
    /*
    This method forces the recalculation of the instrument value and other results
    which would otherwise be cached. It is not declared as const since it needs to
    call the non-const notifyObservers method. Explicit invocation of this method
    is not necessary if the instrument registered itself as observer with the
    structures on which such results depend.
    */
    inline void recalculate() {
        performCalculations();
        calculated = true;
        notifyObservers();
    }
    /*
    method.
    This method performs all needed calculations by calling the performCalculations
    Instruments cache the results of the previous calculation. Such results will be
    returned upon later invocations of calculate. The results depend on arguments
    such as term structures which could change between invocations; the instrument
    must register itself as observer of such objects for the calculations to be
    performed again when they change.
    This method should not be redefined in derived classes. The method does not
    modify the structure of the instrument and is therefore declared as constant.
    Temporary variables are declared as mutable.
    */
    inline double calculate() const {
        if (!calculated)
            performCalculations();
        calculated = true;
        return 0.0;
    }
protected:
    // This method must implement any calculations which must be
    // (re)done in order to calculate the NPV of the instrument.
    virtual void performCalculations() const = 0;
    // The value of these attributes must be set in the body of the
    // performCalculations method.
    mutable double NPV_;
    mutable bool isExpired_;
private:
    string isinCode_, description_; // description of instrument
    mutable bool calculated; // tracks if instrument was calculated
};


/**********************************************************************************
General diffusion process classes
This class describes a stochastic process governed by dx(t) = mu(t, x(t))dt +
sigma(t, x(t))dz(t).
**********************************************************************************/
class DiffusionProcess {
public:
    DiffusionProcess(double x0) : x0_(x0) {}
    virtual DiffusionProcess() {}
    double x0() const { return x0_; }
    // returns the drift part of the equation, i.e. mu(t, x_t)
    virtual double drift(Time t, double x) const = 0;
    // returns the diffusion part of the equation, i.e. sigma(t,x_t)
    virtual double diffusion(Time t, double x) const = 0;
    // returns the expectation of the process after a time interval
    // returns E(x_{t_0 + delta t} | x_{t_0} = x_0) since it is Markov.
    // By default, it returns the Euler approximation defined by
    // x_0 + mu(t_0, x_0) delta t.
    virtual double expectation(Time t0, double x0, Time dt) const {
        return x0 + drift(t0, x0)*dt;
    }
    // returns the variance of the process after a time interval
    // returns Var(x_{t_0 + Delta t} | x_{t_0} = x_0).
    // By default, it returns the Euler approximation defined by
    // sigma(t_0, x_0)^2 \Delta t .
    virtual double variance(Time t0, double x0, Time dt) const {
        double sigma = diffusion(t0, x0);
        return sigma*sigma*dt;
    }
private:
    double x0_;
};


/**********************************************************************************
Black-Scholes diffusion process class
1.4 Black-Scholes and Diffusion Process Implementation 19
This class describes the stochastic process governed by dS = (r – 0.5{sigma^2}) dt
+ sigmadz(t).
**********************************************************************************/
class BlackScholesProcess : public DiffusionProcess {
public:
    BlackScholesProcess(Rate rate, double volatility, double s0 = 0.0)
    : DiffusionProcess(s0), r_(rate), sigma_(volatility) {}
    double drift(Time t, double x) const {
        return r_ - 0.5*sigma_*sigma_;
    }
    double diffusion(Time t, double x) const {
        return sigma_;
    }
private:
    double r_, sigma_;
};

/**********************************************************************************
Ornstein-Uhlenbeck process class
This class describes the Ornstein-Uhlenbeck process governed by dx = -a x(t) dt +
sigma dz(t).
**********************************************************************************/
class OrnsteinUhlenbeckProcess : public DiffusionProcess {
public:
    OrnsteinUhlenbeckProcess(double speed, double vol, double x0 = 0.0)
    : DiffusionProcess(x0), speed_(speed), volatility_(vol) {}
    double drift(Time t, double x) const {
        return -speed_*x;
    }
    double diffusion(Time t, double x) const {
        return volatility_;
    }
    double expectation(Time t0, double x0, Time dt) const {
        return x0*exp(-speed_*dt);
    }
    double variance(Time t0, double x0, Time dt) const {
        return 0.5*volatility_*volatility_/speed_*(1.0 - exp(-2.0*speed_*dt));
    }
private:
    double speed_, volatility_;
};

/**********************************************************************************
Square-root process class
This class describes a square-root process governed by dx = a (b – x_t) dt + \sigma
sqrt{x_t} dW_t.
**********************************************************************************/
class SquareRootProcess : public DiffusionProcess {
public:
    SquareRootProcess(double b, double a, double sigma, double x0 = 0)
    : DiffusionProcess(x0), mean_(b), speed_(a), volatility_(sigma) {}
    double drift(Time t, double x) const {
        return speed_*(mean_ - x);
    }
    double diffusion(Time t, double x) const {
        return volatility_*sqrt(x);
    }
private:
    double mean_, speed_, volatility_;
};



#endif //DIFFUSIONPROCESSES_H
