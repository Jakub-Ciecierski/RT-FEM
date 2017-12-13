#include "RTFEM/DataStructure/Dense2SparseMatrix.h"

namespace rtfem {

template<class T>
SparseMatrixCSR<T> Dense2SparseMatrix<T>::Transform(
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
        dense_matrix){
    unsigned int m = dense_matrix.rows();
    unsigned int n = dense_matrix.cols();

    std::vector<T> values;
    std::vector<int> row_extents(m+1);
    std::vector<int> columns_indices;

    row_extents[0] = 0;
    unsigned int last_row_non_zero_elemen_count = 0;
    for(unsigned int i = 0; i < m; i++){
        unsigned int row_non_zero_elemen_count = 0;
        for(unsigned int j = 0; j < n; j++){
            auto value = dense_matrix(i, j);
            if(value != 0){
                values.push_back(value);
                columns_indices.push_back(j);
                row_non_zero_elemen_count++;
            }
        }
        if(i > 0){
            row_extents[i] = row_extents[i - 1] +
                             last_row_non_zero_elemen_count;
        }
        last_row_non_zero_elemen_count = row_non_zero_elemen_count;
    }
    row_extents[m] = row_extents[m - 1] +
                     last_row_non_zero_elemen_count;

    return SparseMatrixCSR<T>(values, row_extents, columns_indices, m, n);
}

template
class Dense2SparseMatrix<double>;
template
class Dense2SparseMatrix<float>;

}