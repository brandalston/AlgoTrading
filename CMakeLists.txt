cmake_minimum_required(VERSION 3.29)
project(AlgoTrading)

set(CMAKE_CXX_STANDARD 20)

add_executable(AlgoTrading main.cpp
        main.cpp
        MCMConstructor.h
        Utility/statUtility.h
        Utility/matrixUtility.h
        Utility/genericUtility.h
        Utility/constants.h
        Utility/handle.h
        Option_Pricing_Fundamental/Option.h
        Option_Pricing_Fundamental/diffusionProcess.h
        Monte_Carlo_Method/MonteCarlo.h
        Monte_Carlo_Method/pathDependentProcess.h
        Tree_Models/Tree.h
        Utility/nrutil.cpp
        Utility/nrutil.h
)
target_include_directories (AlgoTrading PUBLIC Option_Pricing_Fundamental Monte_Carlo_Method Tree_Models Utility)
