#ifndef PROJECT_GPUSPARSELINEARSOLVER_H
#define PROJECT_GPUSPARSELINEARSOLVER_H

namespace rtfem {

template<class T>
class GPUSparseLinearSolver {
public:
    GPUSparseLinearSolver() = default;
    ~GPUSparseLinearSolver() = default;

    void Solve();
private:

    void genTridiag(int *I, int *J, double *val, int N, int nz);
};

}

#endif //PROJECT_GPUSPARSELINEARSOLVER_H
