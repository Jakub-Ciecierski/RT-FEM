#include "RTFEMTests/GPU/DeviceProperties/DevicePropertiesTest.h"

#include "RTFEM/GPU/DeviceProperties/DevicesProperties.cuh"

void DevicePropertiesTest::SetUp() {
}

void DevicePropertiesTest::TearDown() {
}

TEST_F(DevicePropertiesTest, OStreamOperator){
    rtfem::DevicesProperties devices_properties;
    devices_properties.Update();

    for(const auto& device_properties : devices_properties){
        std::cout << device_properties << std::endl;
    }
}