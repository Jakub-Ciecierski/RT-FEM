#include "RTFEMTests/Builder/FEMModelBoxBuilder.h"

#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/BoundaryConditionContainer.h>
#include <RTFEM/FEM/BoundaryCondition.h>
#include <RTFEM/Memory/UniquePointer.h>

std::shared_ptr<rtfem::FEMModel<double>> FEMModelBoxBuilder::Create(){
    rtfem::FEMGeometry<double> fem_geometry;

    fem_geometry.finite_elements =
        std::vector<std::shared_ptr<rtfem::FiniteElement<double>>>(
            6);

    fem_geometry.vertices =
        std::vector<std::shared_ptr<rtfem::Vertex<double>>>(8);

    fem_geometry.vertices[0] = std::make_shared<rtfem::Vertex<double>>(0,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            1,
                                                                            -1,
                                                                            -1));
    fem_geometry.vertices[1] = std::make_shared<rtfem::Vertex<double>>(1,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            -1,
                                                                            1,
                                                                            -1));
    fem_geometry.vertices[2] = std::make_shared<rtfem::Vertex<double>>(2,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            -1,
                                                                            -1,
                                                                            1));
    fem_geometry.vertices[3] = std::make_shared<rtfem::Vertex<double>>(3,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            1,
                                                                            -1,
                                                                            1));
    fem_geometry.vertices[4] = std::make_shared<rtfem::Vertex<double>>(4,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            -1,
                                                                            -1,
                                                                            -1));
    fem_geometry.vertices[5] = std::make_shared<rtfem::Vertex<double>>(5,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            1,
                                                                            1,
                                                                            1));
    fem_geometry.vertices[6] = std::make_shared<rtfem::Vertex<double>>(6,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            1,
                                                                            1,
                                                                            -1));
    fem_geometry.vertices[7] = std::make_shared<rtfem::Vertex<double>>(7,
                                                                        Eigen::Vector3<
                                                                            double>(
                                                                            -1,
                                                                            1,
                                                                            1));

    fem_geometry.finite_elements[0]
        = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(5, 1, 6, 0);
    fem_geometry.finite_elements[1] =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(4, 1, 2, 0);
    fem_geometry.finite_elements[2] =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(3, 2, 5, 0);
    fem_geometry.finite_elements[3] =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(7, 1, 5, 0);
    fem_geometry.finite_elements[4] =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(1, 7, 2, 0);
    fem_geometry.finite_elements[5] =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(2, 7, 5, 0);

    return std::make_shared<rtfem::FEMModel<double>>(
        fem_geometry,
        rtfem::Material<double>{80000, 0.3});
}