#ifndef PROJECT_SPARSEMATRIXCSR_H
#define PROJECT_SPARSEMATRIXCSR_H

#include<vector>

namespace rtfem {

/**
 * m x n
 * @tparam T
 */
template<class T>
class SparseMatrixCSR {
public:
    SparseMatrixCSR(std::vector<T> values,
                    std::vector<unsigned int> row_extents,
                    std::vector<unsigned int> columns_indices,
                    unsigned int m,
                    unsigned int n);
    ~SparseMatrixCSR() = default;

    const std::vector<T>& values(){return values_;}
    const std::vector<unsigned int>& row_extents(){return row_extents_;}
    const std::vector<unsigned int>& columns_indices(){return columns_indices_;}
    unsigned int m(){return m_;}
    unsigned int n(){return n_;}

private:
    std::vector<T> values_;
    std::vector<unsigned int> row_extents_;
    std::vector<unsigned int> columns_indices_;

    unsigned int m_;
    unsigned int n_;
};
}


#endif //PROJECT_SPARSEMATRIXCSR_H
