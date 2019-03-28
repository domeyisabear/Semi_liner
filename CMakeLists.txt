cmake_minimum_required(VERSION 3.12)
project(semi_liner)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

include_directories(/Library/gurobi702/mac64/include)

link_libraries(/Library/gurobi702/mac64/lib/libgurobi70.so)
link_libraries(/Library/gurobi702/mac64/lib/libgurobi_g++4.2.a)

add_executable(semi_liner
        basic.cpp
        basic.h
        cJSON.h
        cJSON.c
        CJsonObject.hpp
        CJsonObject.cpp
        InputParam.cpp
        InputParam.h
        main.cpp
        MasterProblem.h
        MP_benders.cpp
        MP_benders.h
        MP_cutplane.cpp
        MP_cutplane.h
        MP_cutplane_multi.cpp
        MP_cutplane_multi.h
        MP_extensive.cpp
        MP_extensive.h
        #Problem.cpp
        ReadCsv.cpp
        ReadCsv.h
        SLBPCallback.cpp
        SLBPCallback.h
        SMP_extensive.cpp
        SMP_extensive.h
        SubLBProblem.cpp
        SubLBProblem.h
        SubLBProblem_benders.h
        SubLBProblem_benders.cpp
        SubLPProblem.cpp
        SubLPProblem.h
        SubMIPProblem.h
        threadpool.h
        tools.cpp
        tools.h)
