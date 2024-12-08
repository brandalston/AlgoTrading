cmake_minimum_required(VERSION 3.29)
project(AlgoTrading)

set(CMAKE_CXX_STANDARD 20)

add_executable(AlgoTrading main.cpp
        main.cpp
        packages.h
        Utility/statUtility.h
        Utility/matrixUtility.h
        constants.h
        Monte_Carlo_Methods/MonteCarlo.h
        BSM_Pricing_Fundamentals/pathDependentProcesses.h
        Utility/genericUtility.h
)
target_include_directories (AlgoTrading PUBLIC BSM_Pricing_Fundamentals Monte_Carlo_Methods Utility)
