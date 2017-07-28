#include "RTFEMTests/Builder/FEMModelSampleBuilder.h"

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/FiniteElement.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/Memory/UniquePointer.h>

FEMModelSampleBuilder::FEMModelSampleBuilder() {}

FEMModelSampleBuilder::~FEMModelSampleBuilder() {}

std::shared_ptr<rtfem::FEMModel<double>>
FEMModelSampleBuilder::CreateRandomFEMModel() {
    auto finite_elements = std::vector<std::shared_ptr<rtfem::FiniteElement<double>>>(finite_element_count_);
    auto vertices = std::vector<std::shared_ptr<rtfem::Vertex<double>>>(vertex_count_);

    vertices[0] = std::make_shared<rtfem::Vertex<double>>(0, Eigen::Vector3<double>(0,0,0));
    vertices[1] = std::make_shared<rtfem::Vertex<double>>(1, Eigen::Vector3<double>(1,0,0));
    vertices[2] = std::make_shared<rtfem::Vertex<double>>(2, Eigen::Vector3<double>(1,1,0));
    vertices[3] = std::make_shared<rtfem::Vertex<double>>(3, Eigen::Vector3<double>(1,1,1));
    vertices[4] = std::make_shared<rtfem::Vertex<double>>(4, Eigen::Vector3<double>(0,1,1));
    vertices[5] = std::make_shared<rtfem::Vertex<double>>(5, Eigen::Vector3<double>(2,2,2));
    vertices[6] = std::make_shared<rtfem::Vertex<double>>(6, Eigen::Vector3<double>(1,0,1));
    vertices[7] = std::make_shared<rtfem::Vertex<double>>(7, Eigen::Vector3<double>(0,1,0));
    vertices[8] = std::make_shared<rtfem::Vertex<double>>(8, Eigen::Vector3<double>(0,0,2));

    finite_elements[0] = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertices[0],
                                                                                         vertices[1],
                                                                                         vertices[7],
                                                                                         vertices[8]);
    finite_elements[1] = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertices[1],
                                                                                         vertices[2],
                                                                                         vertices[3],
                                                                                         vertices[8]);
    finite_elements[2] = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertices[3],
                                                                                         vertices[4],
                                                                                         vertices[5],
                                                                                         vertices[8]);
    finite_elements[3] = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertices[5],
                                                                                         vertices[6],
                                                                                         vertices[0],
                                                                                         vertices[8]);

    return std::make_shared<rtfem::FEMModel<double>>(finite_elements, vertices, rtfem::Material<double>{80000, 0.3});
}
