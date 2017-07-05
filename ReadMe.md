# RTFEM
RTFEM stands for Real-Time Finite Element Method.
The library allows to compute the entire process of Finite Element Method for Solid Mechanics in real time, namely:
 1. Pre-Processing:
    * 3D Finite Element Meshing of given Mesh.
 1. FEM Solver:
    * [CPU] Linear elastic, linear deformation, static.
    * TODO
 1. Post-Processing:
    * Output all Solid Mechanics FEM artifacts.

# Language
Library is written in c++11.

# Build
Project is build using cmake.
Run "cmake /." command in order to build the library.

## Unittests 
You can disable/enable building unittests in CMakeLists.txt.
Unittests are written using googletest

# Example application
TODO