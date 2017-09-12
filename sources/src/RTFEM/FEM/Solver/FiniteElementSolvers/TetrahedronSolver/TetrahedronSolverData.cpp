#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver/TetrahedronSolverData.h>

namespace rtfem {

template
struct TetrahedronShapeFunctionCoefficients<double>;
template
struct TetrahedronShapeFunctionCoefficients<float>;

template
struct Edges<double>;
template
struct Edges<float>;

template
struct FacesArea<double>;
template
struct FacesArea<float>;

}
