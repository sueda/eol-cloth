#pragma once

#ifdef WIN32
    #define EIGEN_GUROBI_DLLIMPORT __declspec(dllimport)
    #define EIGEN_GUROBI_DLLEXPORT __declspec(dllexport)
#else
    #define EIGEN_GUROBI_DLLIMPORT
    #define EIGEN_GUROBI_DLLEXPORT
#endif

#ifdef EIGEN_GUROBI_EXPORT
    #define EIGEN_GUROBI_API EIGEN_GUROBI_DLLEXPORT
#else
    #define EIGEN_GUROBI_API EIGEN_GUROBI_DLLIMPORT
#endif
