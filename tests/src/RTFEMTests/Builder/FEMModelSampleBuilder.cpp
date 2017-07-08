#include "RTFEMTests/Builder/FEMModelSampleBuilder.h"

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/FiniteElement.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/DataStructure/Vector3.h>
#include <RTFEM/Memory/UniquePointer.h>

FEMModelSampleBuilder::FEMModelSampleBuilder() {}

FEMModelSampleBuilder::~FEMModelSampleBuilder() {}

std::shared_ptr<rtfem::FEMModel>
FEMModelSampleBuilder::CreateRandomFEMModel() {
    auto finite_elements = std::vector<std::shared_ptr<rtfem::FiniteElement>>(finite_element_count_);
    auto vertices = std::vector<std::shared_ptr<rtfem::Vertex>>(vertex_count_);

    vertices[0] = std::make_shared<rtfem::Vertex>(0, rtfem::Vector3(0,0,0));
    vertices[1] = std::make_shared<rtfem::Vertex>(1, rtfem::Vector3(1,0,0));
    vertices[2] = std::make_shared<rtfem::Vertex>(2, rtfem::Vector3(1,1,0));
    vertices[3] = std::make_shared<rtfem::Vertex>(3, rtfem::Vector3(1,1,1));
    vertices[4] = std::make_shared<rtfem::Vertex>(4, rtfem::Vector3(0,1,1));
    vertices[5] = std::make_shared<rtfem::Vertex>(5, rtfem::Vector3(2,2,2));
    vertices[6] = std::make_shared<rtfem::Vertex>(6, rtfem::Vector3(1,0,1));
    vertices[7] = std::make_shared<rtfem::Vertex>(7, rtfem::Vector3(0,1,0));
    vertices[8] = std::make_shared<rtfem::Vertex>(8, rtfem::Vector3(0,0,2));

    finite_elements[0] = std::make_shared<rtfem::TetrahedronFiniteElement>(vertices[0],
                                                                           vertices[1],
                                                                           vertices[7],
                                                                           vertices[8]);
    finite_elements[1] = std::make_shared<rtfem::TetrahedronFiniteElement>(vertices[1],
                                                                           vertices[2],
                                                                           vertices[3],
                                                                           vertices[8]);
    finite_elements[2] = std::make_shared<rtfem::TetrahedronFiniteElement>(vertices[3],
                                                                           vertices[4],
                                                                           vertices[5],
                                                                           vertices[8]);
    finite_elements[3] = std::make_shared<rtfem::TetrahedronFiniteElement>(vertices[5],
                                                                           vertices[6],
                                                                           vertices[0],
                                                                           vertices[8]);

    return std::make_shared<rtfem::FEMModel>(finite_elements, vertices, rtfem::Material{80000, 0.3});
}
