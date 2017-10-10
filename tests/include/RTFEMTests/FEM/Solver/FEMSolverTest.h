#include "gtest/gtest.h"

#include <memory>

namespace rtfem {

template<class T>
class FEMModel;

}
class FEMSolverTest : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

    std::shared_ptr<rtfem::FEMModel<double>> fem_model_;
};