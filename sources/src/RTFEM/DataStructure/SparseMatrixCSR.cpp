#include "RTFEM/DataStructure/SparseMatrixCSR.h"

namespace rtfem {

template<class T>
SparseMatrixCSR<T>::SparseMatrixCSR(std::vector<T> values,
                                    std::vector<int> row_extents,
                                    std::vector<int> columns_indices,
                                    unsigned int m,
                                    unsigned int n,
                                    MatrixType type) :
        values_(values),
        row_extents_(row_extents),
        columns_indices_(columns_indices),
        m_(m), n_(n),
        type_(type){}

template
class SparseMatrixCSR<double>;
template
class SparseMatrixCSR<float>;

}