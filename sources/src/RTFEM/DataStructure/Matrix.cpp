#include "RTFEM/DataStructure/Matrix.h"

#include <stdexcept>

namespace rtfem {

Matrix::Matrix(UInt row_count, UInt column_count) :
        dimensions_{row_count, column_count}{
    InitData();
}

Matrix::Matrix(const MatrixDimension &&matrix_dimension) :
        dimensions_(matrix_dimension) {
    InitData();
}

Matrix::~Matrix() {}

std::vector<Float>& Matrix::operator[](UInt i) {
    return GetRow(i);
}

const std::vector<Float>& Matrix::operator[] (UInt i) const {
    return GetRow(i);
}

std::vector<Float>& Matrix::GetRow(UInt i){
    if(i >= dimensions_.row_count)
        throw std::out_of_range("Matrix Index out of Range");

    return data_[i];
}

const std::vector<Float>& Matrix::GetRow(UInt i) const {
    if(i >= dimensions_.row_count)
        throw std::out_of_range("Matrix Index out of Range");

    return data_[i];
}

void Matrix::InitData(){
    data_.resize(dimensions_.row_count);

    for(UInt i = 0; i < dimensions_.row_count; i++){
        data_[i].resize(dimensions_.column_count);
        for(UInt j = 0; j < dimensions_.column_count; j++){
            data_[i][j] = 0;
        }
    }
}

}