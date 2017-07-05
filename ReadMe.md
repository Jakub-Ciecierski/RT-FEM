# RTFEM
RTFEM stands for Real-Time Finite Element Method. <br/>
The library computes the entire process of Finite Element Method for solid mechanics in real time, namely:
 1. Pre-Processing:
    * 3D Finite Element Meshing of given geometry.
 1. FEM Solver:
    * [CPU] Linear elastic, linear deformation, static.
    * TODO
 1. Post-Processing:
    * Output all Solid Mechanics FEM artifacts.

# Language
Library is written in c++11.
TODO

# Build
RTFEM can be built on any platform supporting cmake. <br/>
To build the library, run the following command in root directory:
```
cmake .
```

### Unittests
Building Unittests can be disabled/enabled in CMakeLists.txt. <br/>
Unittests are written using googletest

### Benchmarks
TODO

# Example application
TODO