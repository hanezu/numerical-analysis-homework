cmake_minimum_required(VERSION 3.6)
project(programs_elliptic_pde)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    000-Poisson-Rect/poisson-rect.c
        001-Poisson-Param-Region/poisson-param-region.c
    002-Poisson-Sublevelset/poisson-levelset.c 004-Poisson-Neumann-Region/poisson-neumann-region.h test/004.c test/0042.c 004-Poisson-Neumann-Region/test.c 004-Poisson-Neumann-Region/test.c)

add_executable(programs_elliptic_pde ${SOURCE_FILES})