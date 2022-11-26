cmake_minimum_required(VERSION 3.15)
project(NF)

set(CMAKE_CXX_STANDARD 14)

add_executable(NF adv_ann.cpp
        adv_ann.h
        ann_ext.cpp
        ann_ext.h
        arbfn.cpp
        arbfn.h
        eigen.cpp
        eigen.h
        filesystem.cpp
        filesystem.h
        gconfig.cpp
        gconfig.h
        gng.cpp
        gng.h
        indiv.h
        main.cpp
        matrix.cpp
        matrix.h
        metamodel.cpp
        metamodel.h
        myblas.cpp
        myblas.h
        mymath.h
        param.cpp
        param.h
        rng.cpp
        rng.h
        svd.cpp
        svd.h
        util.cpp
        util.h)