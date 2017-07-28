#include "RTFEM/DataStructure/Matrix.h"

#include <stdexcept>

namespace rtfem {

Matrix::Matrix(unsigned int row_count, unsigned int column_count) :
        dimensions_{row_count, column_count}{
    InitData(0.0);
}

Matrix::Matrix(double value, unsigned int row_count, unsigned int column_count) :
        dimensions_{row_count, column_count}{
    InitData(value);
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> lists) :
        dimensions_{0,0} {
    unsigned int max_size = 0;
    for(const auto& list : lists){
        if(max_size < list.size())
            max_size = list.size();
    }
    dimensions_ = MatrixDimension{(unsigned int)lists.size(), max_size};
    InitData(0.0);
    unsigned int i = 0;

    for(const auto& list : lists){
        unsigned int j = 0;
        for(const auto& value : list){
            data_[i][j] = value;
            j++;
        }
        i++;
    }
}

Matrix::Matrix(const MatrixDimension &&matrix_dimension) :
        dimensions_(matrix_dimension) {
    InitData(0.0);
}

Matrix::~Matrix() {}

std::vector<double>& Matrix::operator[](unsigned int i) {
    return GetRow(i);
}

const std::vector<double>& Matrix::operator[] (unsigned int i) const {
    return GetRow(i);
}

std::vector<double>& Matrix::GetRow(unsigned int i){
    if(i >= dimensions_.row_count)
        throw std::out_of_range("Matrix Index out of Range");

    return data_[i];
}

const std::vector<double>& Matrix::GetRow(unsigned int i) const {
    if(i >= dimensions_.row_count)
        throw std::out_of_range("Matrix Index out of Range");

    return data_[i];
}

void Matrix::InitData(double value){
    data_.resize(dimensions_.row_count);

    for(unsigned int i = 0; i < dimensions_.row_count; i++){
        data_[i].resize(dimensions_.column_count);
        for(unsigned int j = 0; j < dimensions_.column_count; j++){
            data_[i][j] = 0;
        }
    }
}

bool Matrix::operator==(const Matrix &rhs) const {
    if(dimensions().column_count != rhs.dimensions().column_count ||
            dimensions().row_count != rhs.dimensions().row_count)
        return false;

    unsigned int n = dimensions().row_count;
    unsigned int m = dimensions().column_count;
    for(unsigned int i = 0; i < n; i++){
        for(unsigned int j = 0; j < m; j++){
            if((*this)[i][j] != rhs[i][j])
                return false;
        }
    }
    return true;
}

bool Matrix::operator!=(const Matrix &rhs) const {
    return !(rhs == *this);
}

Matrix& Matrix::operator+=(const Matrix& rhs){
    if(this->dimensions().row_count != rhs.dimensions().row_count ||
            this->dimensions().column_count != rhs.dimensions().column_count)
        throw std::invalid_argument("Matrix Addition: Dimensions are not valid");

    for(unsigned int i = 0; i < this->dimensions().row_count; i++){
        for(unsigned int j = 0; j < this->dimensions().column_count; j++){
            (*this)[i][j] = (*this)[i][j] + rhs[i][j];
        }
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m){
    const std::string seperator = "   ";
    for(unsigned int i = 0; i < m.dimensions().row_count; i++){
        for(unsigned int j = 0; j < m.dimensions().column_count; j++){
            os << std::to_string(m[i][j]) << seperator;
        }
        os << std::endl;
    }
    return os;
}

Matrix operator*(double x, const Matrix &m) {
    Matrix new_matrix(std::move(m.dimensions()));

    for(unsigned int i = 0; i < m.dimensions().row_count; i++){
        for(unsigned int j = 0; j < m.dimensions().column_count; j++){
            new_matrix[i][j] = m[i][j] * x;
        }
    }
    return new_matrix;
}

Matrix operator*(const Matrix &m, double x) {
    Matrix new_matrix(std::move(m.dimensions()));

    for(unsigned int i = 0; i < m.dimensions().row_count; i++){
        for(unsigned int j = 0; j < m.dimensions().column_count; j++){
            new_matrix[i][j] = m[i][j] * x;
        }
    }
    return new_matrix;
}

Matrix operator/(const Matrix &m, double x) {
    Matrix new_matrix(std::move(m.dimensions()));

    for(unsigned int i = 0; i < m.dimensions().row_count; i++){
        for(unsigned int j = 0; j < m.dimensions().column_count; j++){
            new_matrix[i][j] = m[i][j] / x;
        }
    }
    return new_matrix;
}

Matrix operator*(const Matrix& m1, const Matrix& m2){
    if(m1.dimensions().column_count != m2.dimensions().row_count)
        throw std::invalid_argument("Matrix multiplication: Dimensions are not valid");
    unsigned int n = m1.dimensions().column_count;

    Matrix m(m1.dimensions().row_count, m2.dimensions().column_count);

    for(unsigned int i = 0; i < m1.dimensions().row_count; i++){ // 3
        for(unsigned int j = 0; j < m2.dimensions().column_count; j++){ // 4
            for(unsigned int k = 0; k < n; k++){
                m[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return m;
}

Matrix operator+(const Matrix& m1, const Matrix& m2){
    if(m1.dimensions().row_count != m2.dimensions().row_count ||
            m1.dimensions().column_count != m2.dimensions().column_count)
        throw std::invalid_argument("Matrix Addition: Dimensions are not valid");
    Matrix m(m1.dimensions().row_count, m1.dimensions().column_count);

    for(unsigned int i = 0; i < m1.dimensions().row_count; i++){
        for(unsigned int j = 0; j < m1.dimensions().column_count; j++){
            m[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return m;
}

}