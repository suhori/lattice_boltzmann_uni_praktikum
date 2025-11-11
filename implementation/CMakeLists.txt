cmake_minimum_required(VERSION 4.0)
project(lattice_boltzmann_uni_praktikum)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(.)

add_executable(lattice_boltzmann_uni_praktikum
        LBM.cpp
        LBM.h
        main.cpp
        seconds.cpp
        seconds.h)
