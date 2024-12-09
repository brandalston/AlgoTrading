//
// Created by Brandon Alston on 12/7/24.
//

#ifndef PATHDEPENDENTPROCESSES_H
#define PATHDEPENDENTPROCESSES_H
#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include "statUtility.h"
#include "matrixUtility.h"
#include "genericUtility.h"
#include "constants.h"

class BrownianBridge {
public:
    BrownianBridge(unsigned long numberOfSteps);
    virtual BrownianBridge() {};
    void buildPath(vector<double>& theWienerProcessPath, const vector<double>&
    gaussianVariates);
    void generateDeviates(unsigned long numberOfSteps);
private:
    unsigned long numberOfSteps;
    vector<unsigned long> leftIndex;
    vector<unsigned long> rightIndex;
    vector<unsigned long> bridgeIndex;
    vector<double> leftWeight;
    vector<double> rightWeight;
    vector<double> stddev;
    vector<double> normalVariates;
    StatUtility util;

    /**********************************************************************************
    BrownianBridge : Constructor , initializes Brownian bridge
    [in] numberOfSteps: number of steps on path
    [out] none
    **********************************************************************************/
    BrownianBridge::BrownianBridge(unsigned long numberOfSteps):numberOfSteps(numberOfSteps), leftIndex(numberOfSteps), rightIndex(numberOfSteps), bridgeIndex(numberOfSteps), leftWeight(numberOfSteps), rightWeight(numberOfSteps), stddev(numberOfSteps) {
        vector<unsigned long> map(numberOfSteps);
        unsigned long i, j, k, l;
        // map is used to indicated which points are already constructed. If map[i] is zero, path point i is yet
        // unconstructed. map[i] - 1 is the index of the variate that constructs the path point i.
        map[numberOfSteps-1] = 1; // the first point in the construction is the global step
        bridgeIndex[0] = numberOfSteps - 1; // bridge index
        stddev[0] = sqrt(numberOfSteps); // the standard deviation of the global step
        leftWeight[0] = rightWeight[0] = 0; // global step to the last point in time
        for (j = 0, i = 0; i < numberOfSteps; ++i) {
            while (map[j])
                ++j; // find the next unpopulated entry in the map
            k = j;
            while ((!map[k])) // find the next unpopulated entry in the map from there
                ++k;
            l = j + ((k - 1 - j) >> 1); // l is now the index of the point to be constructed next
            map[l] = i;
            bridgeIndex[i] = l; // the ith gaussian variate to be used to set point l
            leftIndex[i] = j;  // point j-1 is the left strut of the bridge for point l
            rightIndex[i] = k; // point k is the right strut of the bridge
            leftWeight[i] = (k - l)/(k + 1- j);
            rightWeight[i] = (l + 1 - j)/(k+ 1- j);
            stddev[i] = sqrt(((l + 1 - j)*(k-1))/(k+1-j));
            j = k + 1;
            if (j >= numberOfSteps)
                j = 0; // wrap around
        }
    }

    /**********************************************************************************
    buildPath: builds a path using a Brownian bridge
    [in] path : simulated Brownian path
    [in] normalVariates : vector of normal deviates
    [out] none
    **********************************************************************************/
    void BrownianBridge::buildPath(vector<double>& path, const vector<double>&normalVariates) {
        assert(normalVariates.size() == numberOfSteps && path.size() == numberOfSteps);
        unsigned long i, j, k, l;
        path[numberOfSteps - 1] = stddev[0]*normalVariates[0];
        for (i = 1; i < numberOfSteps; i++) {
            j = leftIndex[i];
            k = rightIndex[i];
            l = bridgeIndex[i];
            if (j) {
                path[l] = leftWeight[i]*path[j-1] + rightWeight[i]*path[k] +
                stddev[i]*normalVariates[i];
            } else {
                path[l] = rightWeight[i]*path[k] + stddev[i]*normalVariates[i];
            }
        }
    }
    /**********************************************************************************
    generateDeviates: generates a sequences of normal random deviates
    [in] numberOfSteps: number of steps per path (= number of deviates needed per path)
    [out] none
    **********************************************************************************/
    void BrownianBridge::generateDeviates(unsigned long numberOfSteps) {
        double deviate = 0.0;
        srand(time(0)); // initialize RNG
        long seed = (long) rand(); // generate random seed
        long* idum = &seed;
        for (int i=0; i < numberOfSteps; i++) {
            deviate = util.gasdev(idum);
            normalVariates.push_back(deviate);
        }
    }
};



#endif //PATHDEPENDENTPROCESSES_H

