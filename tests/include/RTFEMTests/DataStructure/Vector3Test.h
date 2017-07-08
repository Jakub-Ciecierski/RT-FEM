#ifndef PROJECT_VECTOR3TEST_H
#define PROJECT_VECTOR3TEST_H

#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>

#include <memory>

namespace rtfem {
struct Vector3;
}

class Vector3Test : public ::testing::Test {
protected:
    virtual void SetUp() override;

    virtual void TearDown() override;

protected:
    std::unique_ptr<rtfem::Vector3> vector_;

};


#endif //PROJECT_VECTOR3TEST_H
