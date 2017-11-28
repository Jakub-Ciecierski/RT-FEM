#include "RTFEM/GPU/DeviceProperties/DevicesProperties.cuh"

#include "cuda_runtime.h"

namespace rtfem {

void DevicesProperties::Update(){
    int device_count;
    cudaGetDeviceCount(&device_count);

    device_properties_.clear();
    for(int i = 0; i < device_count; i++){
        cudaDeviceProp properties;
        cudaGetDeviceProperties(&properties, i);

        device_properties_.push_back(DeviceProperties{properties});
    }
}

DeviceProperties* DevicesProperties::begin(){
    return &(device_properties_[0]);
}

DeviceProperties* DevicesProperties::end(){
    return &(device_properties_[device_properties_.size()]);
}

const DeviceProperties* DevicesProperties::begin() const{
    return &(device_properties_[0]);
}

const DeviceProperties* DevicesProperties::end() const{
    return &(device_properties_[device_properties_.size()]);
}

}