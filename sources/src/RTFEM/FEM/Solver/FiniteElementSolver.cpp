#include <RTFEM/FEM/Solver/FiniteElementSolver.h>

namespace rtfem {

template
class FiniteElementSolver<double>;
template
class FiniteElementSolver<float>;

template
struct TractionForces<double>;
template
struct TractionForces<float>;

template
struct FiniteElementSolverData<double>;
template
struct FiniteElementSolverData<float>;

}