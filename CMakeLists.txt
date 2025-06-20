cmake_minimum_required(VERSION 3.30)
project(Bio_dpt)

set(CMAKE_CXX_STANDARD 20)

add_executable(Bio_dpt main.cpp
        src/pcr_lib.cpp
        src/pcr_lib.h)
