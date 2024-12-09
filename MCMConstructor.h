//
// Created by Brandon Alston on 11/2/24.
//

#ifndef PACKAGES_H
#define PACKAGES_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "handle.h"
#include "ql/MonteCarlo/path.h"
#include "ql/RandomNumbers/randomarraygenerator.h"
#include "ql/diffusionprocess.h"
using namespace std;

size_t Size;
double Time;

namespace QuantLib {
namespace MonteCarlo {
    // General purpose Monte Carlo model for path samples
    /* Any Monte Carlo which uses path samples has three main components,
    namely,
    - S, a sample accumulator,
    - PG, a path generator,
    - PP, a path pricer.
    MonteCarloModel<S, PG, PP> puts together these three elements.
    The constructor accepts two safe references, i.e. two smart
    pointers, one to a path generator and the other to a path pricer.
    In case of control variate technique the user should provide the
    additional control option, namely the option path pricer and the
    option value.
    The minimal interfaces for the classes S, PG, and PP are:
    class S {
        void add(VALUE_TYPE sample, double weight) const;
    };
    class PG {
        Sample<PATH_TYPE> next() const;
    };
    class PP :: unary_function<PATH_TYPE, VALUE_TYPE> {
        VALUE_TYPE operator()(PATH_TYPE &) const;
    };
    */
    template<class S, class PG, class PP>
    class MonteCarloModel
    {
    public:
        typedef typename PG::sample_type sample_type;
        typedef typename PP::result_type result_type;
        MonteCarloModel(const Handle<PG>& pathGenerator,
        const Handle<PP>& pathPricer,
        const S& sampleAccumulator,
        const Handle<PP>& cvPathPricer = Handle<PP>(),
        result_type cvOptionValue = result_type());
        void addSamples(Size samples);
        const S& sampleAccumulator(void) const;
    private:
        Handle<PG> pathGenerator_; // path generator
        Handle<PP> pathPricer_; // path pricer
        S sampleAccumulator_; // sample accumulator
        Handle<PP> cvPathPricer_; // control variate path price
        result_type cvOptionValue_; // control variate option value
        bool isControlVariate_;

        // inline definitions
        template<class S, class PG, class PP>
        inline MonteCarloModel<S, PG, PP>::MonteCarloModel( const Handle<PG>& pathGenerator,
            const Handle<PP>& pathPricer, const S& sampleAccumulator, const Handle<PP>& cvPathPricer,
            MonteCarloModel<S, PG, PP>::result_type cvOptionValue: pathGenerator_(pathGenerator), pathPricer_(pathPricer),
            sampleAccumulator_(sampleAccumulator), cvPathPricer_(cvPathPricer),cvOptionValue_(cvOptionValue)) {
            if (cvPathPricer_.isNull())
                isControlVariate_= false; // no control variates
            else
                isControlVariate_= true; // use control variates
        }

        template<class S, class PG, class PP>
        inline void MonteCarloModel<S, PG, PP>::addSamples(Size samples) {
            for(Size j = 1; j <= samples; j++) {
                sample_type path = pathGenerator_->next();
                result_type price = (*pathPricer_)(path.value);
                if (isControlVariate_)
                    price += cvOptionValue_-(*cvPathPricer_)(path.value);
                sampleAccumulator_.add(price, path.weight);
            }
        }

        template<class S, class PG, class PP>
        inline const S& MonteCarloModel<S, PG, PP>::sampleAccumulator() const {
            return sampleAccumulator_;
        }
    };



}
}
#endif //PACKAGES_H
