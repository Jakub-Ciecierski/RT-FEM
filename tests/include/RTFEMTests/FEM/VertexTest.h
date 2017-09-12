#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

#include <Eigen/Core>

namespace rtfem {
template<class T>
class Vertex;
}

struct VertexPack {
    unsigned int id;
    Eigen::Matrix<double, 3, 1> coordinates;
    std::unique_ptr<rtfem::Vertex<double>> vertex;
};

class VertexTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    VertexPack default_vertex_;
};