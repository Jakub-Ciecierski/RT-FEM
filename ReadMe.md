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

# Supported Platforms
Library is written in c++11.
TODO
* Linux 32/64bit (Tested)
* Windows 32/64bit (Not Tested TODO)

# Build
RTFEM can be built on any platform supporting cmake. <br/>
To build the library, run the following command in root directory:
```
cmake .
```

### Requirements
The following software is required to be installed: 
* Cmake 3.3

### External Dependencies
External Dependencies can be found in ./external and tests/external directories. <br/>
They are installed automatically when using cmake installation.
* TetGen (external/tetgen1.5.0)
    * 3D Tetrahedron meshing algorithm
* Eigen (external/eigen)
* googletest (tests/external/googletest)
    * Unit testing

### Unittests
Building Unittests can be disabled/enabled in CMakeLists.txt. <br/>
Unittests are written using googletest

### Benchmarks
TODO

# Demo application
TODO