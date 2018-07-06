# RTFEM
RTFEM stands for Real-Time Finite Element Method. <br/>
The library computes the entire process of Finite Element Method for solid mechanics in real time, namely:
 1. Pre-Processing:
    * 3D Finite Element Meshing of given geometry.
 1. FEM Solver:
    * CPU and GPU Linear elastic, linear deformation, dynamic.
Sample video:
https://drive.google.com/file/d/1A7MO1FB9vLsI6-7j8mlbYS4NfCRSP8pY/view?usp=sharing

# Supported Platforms
Library is written in c++11.

* Linux 32/64bit
* Windows 32/64bit

# Build
RTFEM can be built on any platform supporting cmake. <br/>
To build the library, run the following command in root directory:
```
cmake .
```

### Requirements
The following software is required to be installed: 
* Cmake 3.3
* CUDA 8
    * Set environment variable CUDA_BIN_PATH to root directory (e.g. /opt/cuda/)
### External Dependencies
External Dependencies can be found in ./external and tests/external directories. <br/>
These are installed automatically when using cmake installation.
* TetGen (external/tetgen1.5.0)
    * 3D Tetrahedron meshing algorithm
* Eigen (external/eigen)
    * Matrix operations and System of Linear Equation Solver
* googletest (tests/external/googletest)
    * Unit testing

### Unittests
Building Unittests can be disabled/enabled in CMakeLists.txt. <br/>
Unittests are written using googletest

# Other Info
Code formatted using CLion built-in Google settings with modifications:
* All indent is 4
* Indent of class/struct visiblity is 0
