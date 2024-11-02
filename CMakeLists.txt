cmake_minimum_required(VERSION 3.29)
project(AlgoTrading)

set(CMAKE_CXX_STANDARD 20)

add_executable(AlgoTrading main.cpp
        main.cpp
        packages.h
        statUtility.h
        matrixUtility.h
        constants.h
)
target_include_directories (AlgoTrading PUBLIC BSM_Pricing_Fundamentals)
