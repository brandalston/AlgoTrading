//
// Created by Brandon Alston on 12/7/24.
//

#ifndef GENERICUTILITY_H
#define GENERICUTILITY_H

#include <vector>
#include <set>
#include <iostream>
#include <cmath>

class genericUtility {
    inline int poisson(double lambda) {
        assert (lambda > 0. );
        double a = exp( -lambda );
        double b = 1;
        // initialize random number generator
        srand(0);
        long seed = (long) rand() % 100;
        long* idum = &seed;
        for (int i = 0; b >= a; i++ )
            b *= gasdev(idum);
        return i - 1;
    }

};

#endif //GENERICUTILITY_H
