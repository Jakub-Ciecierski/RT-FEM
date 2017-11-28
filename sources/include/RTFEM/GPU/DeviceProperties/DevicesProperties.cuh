#ifndef PROJECT_DEVICESPROPERTIES_H
#define PROJECT_DEVICESPROPERTIES_H

#include "RTFEM/GPU/DeviceProperties/DeviceProperties.cuh"

#include <vector>

namespace rtfem {

class DevicesProperties {
public:
    DevicesProperties() = default;
    ~DevicesProperties() = default;

    const std::vector<DeviceProperties>& device_properties(){
        return device_properties_;
    }

    void Update();

    DeviceProperties* begin();
    DeviceProperties* end();

    const DeviceProperties* begin() const ;
    const DeviceProperties* end() const ;

private:
    std::vector<DeviceProperties> device_properties_;
};

}

#endif //PROJECT_DEVICESPROPERTIES_H
