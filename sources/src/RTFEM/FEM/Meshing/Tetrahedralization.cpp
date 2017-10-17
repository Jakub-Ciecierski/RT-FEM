#include "RTFEM/FEM/Meshing/Tetrahedralization.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/FEM/Vertex.h>

#include <tetgen.h>

#include <iostream>
#include <RTFEM/Memory/UniquePointer.h>

namespace rtfem {

const unsigned int DIMENSION_COUNT = 3;
const unsigned int TETRAHEDRA_DOF_COUNT = 4;

template<class T>
void Tetrahedralization<T>::SetOptions(
    const TetrahedralizationOptions &options) {
    options_ = options;
}

template<class T>
std::unique_ptr<FEMGeometry<T>> Tetrahedralization<T>::Compute(
    const TriangleMeshIndexed<T> &triangle_mesh) {
    tetgenio tetgen_input, tetgen_output;
    tetgenbehavior tetgen_options;

    SetupInput(triangle_mesh, tetgen_input, tetgen_options);
    tetrahedralize(&tetgen_options, &tetgen_input, &tetgen_output);
    return FetchOutput(tetgen_output);
}

template<class T>
void Tetrahedralization<T>::SetupInput(
    const TriangleMeshIndexed<T> &triangle_mesh,
    tetgenio &tetgen_input,
    tetgenbehavior &tetgen_options) {
    SetupInputPoints(triangle_mesh, tetgen_input);
    SetupInputFacets(triangle_mesh, tetgen_input);
    SetupInputOptions(tetgen_options);
}

template<class T>
void Tetrahedralization<T>::SetupInputPoints(
    const TriangleMeshIndexed<T> &triangle_mesh,
    tetgenio &tetgen_input) {
    tetgen_input.numberofpoints = triangle_mesh.points.size();
    tetgen_input.pointlist =
        new REAL[tetgen_input.numberofpoints * DIMENSION_COUNT];
    for (unsigned int i = 0; i < (unsigned int) tetgen_input.numberofpoints;
         i++) {
        tetgen_input.pointlist[(i * DIMENSION_COUNT) + 0] =
            triangle_mesh.points[i][0];
        tetgen_input.pointlist[(i * DIMENSION_COUNT) + 1] =
            triangle_mesh.points[i][1];
        tetgen_input.pointlist[(i * DIMENSION_COUNT) + 2] =
            triangle_mesh.points[i][2];
    }
}

template<class T>
void Tetrahedralization<T>::SetupInputFacets(
    const TriangleMeshIndexed<T> &triangle_mesh,
    tetgenio &tetgen_input) {
    tetgen_input.numberoffacets = triangle_mesh.triangles.size();
    tetgen_input.facetlist = new tetgenio::facet[tetgen_input.numberoffacets];
    tetgen_input.facetmarkerlist = new int[tetgen_input.numberoffacets];
    for (unsigned int i = 0; i < (unsigned int) tetgen_input.numberoffacets;
         i++) {
        auto f = &tetgen_input.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];

        auto p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = triangle_mesh.triangles[i].v1;
        p->vertexlist[1] = triangle_mesh.triangles[i].v2;
        p->vertexlist[2] = triangle_mesh.triangles[i].v3;

        f->numberofholes = 0;
        f->holelist = NULL;

        tetgen_input.facetmarkerlist[i] = 1;
    }
}

template<class T>
void Tetrahedralization<T>::SetupInputOptions(
    tetgenbehavior &tetgen_options) {
    std::vector<char> commandline;
    commandline.push_back('p');
    commandline.push_back('q');

    if (options_.maximum_volume > 0) {
        commandline.push_back('a');
        auto maximum_volume_str = std::to_string(options_.maximum_volume);
        for (unsigned int i = 0; i < maximum_volume_str.size(); i++) {
            commandline.push_back(maximum_volume_str[i]);
        }
    }
    commandline.push_back('\0');

    if (!tetgen_options.parse_commandline(commandline.data())) {
        throw std::invalid_argument("parse_commandline incorrect commandline");
    }
}

template<class T>
std::unique_ptr<FEMGeometry<T>> Tetrahedralization<T>::FetchOutput(
    tetgenio &tetgen_output) {
    auto fem_geometry = rtfem::make_unique<FEMGeometry<T>>();

    FetchPoints(*fem_geometry, tetgen_output);
    FetchTetrahedra(*fem_geometry, tetgen_output);
    FetchFaces(*fem_geometry, tetgen_output);

    return fem_geometry;
}

template<class T>
void Tetrahedralization<T>::FetchPoints(FEMGeometry<T> &fem_geometry,
                                        tetgenio &tetgen_output) {
    for (int i = 0; i < tetgen_output.numberofpoints; i++) {
        Eigen::Vector3<T> coordinates{
            static_cast<T>(tetgen_output.pointlist[(i * DIMENSION_COUNT) + 0]),
            static_cast<T>(tetgen_output.pointlist[(i * DIMENSION_COUNT) + 1]),
            static_cast<T>(tetgen_output.pointlist[(i * DIMENSION_COUNT) + 2])
        };
        auto vertex = std::make_shared<Vertex<T>>(i, coordinates);
        fem_geometry.vertices.push_back(vertex);
    }
}

template<class T>
void Tetrahedralization<T>::FetchTetrahedra(FEMGeometry<T> &fem_geometry,
                                            tetgenio &tetgen_output) {
    if (tetgen_output.numberofcorners != TETRAHEDRA_DOF_COUNT) {
        throw std::invalid_argument("tetgen_output.numberofcorners must be 4");
    }

    for (int i = 0; i < tetgen_output.numberoftetrahedra; i++) {
        auto start_index = i * tetgen_output.numberofcorners;

        auto finite_element = std::make_shared<TetrahedronFiniteElement<T>>(
            tetgen_output.tetrahedronlist[start_index + 0],
            tetgen_output.tetrahedronlist[start_index + 1],
            tetgen_output.tetrahedronlist[start_index + 2],
            tetgen_output.tetrahedronlist[start_index + 3]);

        fem_geometry.finite_elements.push_back(finite_element);
    }
}

template<class T>
void Tetrahedralization<T>::FetchFaces(FEMGeometry<T> &fem_geometry,
                                       tetgenio &tetgen_output) {
    // TODO Faces
    std::cout << "Number of trifaces: "
              << tetgen_output.numberoftrifaces
              << std::endl;
    for(int i = 0 ; i < tetgen_output.numberoftrifaces; i++){
        std::cout << tetgen_output.trifacemarkerlist[i] << std::endl;
    }

}

template
class Tetrahedralization<float>;
template
class Tetrahedralization<double>;

}