#include <RTFEM/FEM/Solver/FiniteElementSolver.h>

namespace rtfem {

template<class T>
FiniteElementSolver<T>::FiniteElementSolver() :
        status_(FiniteElementSolverStatus::OK){}

template<class T>
const FiniteElementSolverStatus& FiniteElementSolver<T>::GetStatus(){
    return status_;
}

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