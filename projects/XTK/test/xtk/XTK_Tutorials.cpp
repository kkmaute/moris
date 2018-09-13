/*
 * XTK_Tutorials.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: doble
 */
#include "catch.hpp"

/*
 * XTK at a minimum needs:
 * 1.)An XTK model
 *
 */
#include "xtk/cl_XTK_Model.hpp"

/*
 * 4.)The Matrix class
 */
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

/*
 * 3.) The mesh
 */
#include "mesh/cl_Mesh_Data.hpp" // Mesh Data API
#include "mesh/cl_Mesh_Builder_Stk.hpp" // Method to build a mesh
#include "mesh/cl_Mesh_Enums.hpp" // convenient enums for mesh (i.e. entity ranks)

/*
 * 5.) A geometry and geometry engine
 */
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

using namespace mesh;
using namespace moris;
namespace xtk
{

TEST_CASE("XTK Tutorial","[XTK_Tutorial]"){
    // Setting up the background mesh

    std::cout<<"Setting up background"<<std::endl;

    /*
     *  Specify the type of mesh building implementation to use, in this case,
     *  a STK mesh builder is used meaning a STK mesh will be the implementation
     *  for the background mesh
     */
    mesh::Mesh_Builder_Stk<real, size_t, DDRMat, DDSTMat> tMeshBuilder;

    /*
     * Next a file name (absolute path) needs to be specified)
     * XTKROOT is an environment variable specifying a directory where
     * some commonly used meshes are stored
     */
    std::string tPrefix = std::getenv("XTKROOT");
    std::string tMeshFileName = tPrefix + "/TestExoFiles/sandwich.e";

    /*
     * Specify the nodal fields to declare on mesh.
     * note: these do note need to be on the exodus file.
     */
    Cell<std::string> tNodalFieldNames = {"NODEFIELD1"} ;

    /*
     * Specify whether or not to create faces and edges
     * Note: XTK needs faces and edges, so in general
     * this should be set to true
     */
    bool tCreateFaces = true;

    /*
     * Load the mesh into the Mesh_Data
     */
    std::shared_ptr<Mesh_Data<real, size_t, DDRMat, DDSTMat>> tBackgroundMesh
    = tMeshBuilder.build_mesh_from_string( tMeshFileName,
                                           tNodalFieldNames,
                                           tCreateFaces);


    /*
     * Setup the geometries and geometry engine
     *
     * In this case, a problem with two geometries is set-up
     */
    size_t tNumGeometries = 2;
    /*
     * The first sphere parameters
     */
    real tRadius1  = 5.1;
    real tXCenter1 = 0.0;
    real tYCenter1 = 0.0;
    real tZCenter1 = 0.0;
    Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere1
    = Sphere<real, size_t, DDRMat, DDSTMat>(tRadius1,tXCenter1,tYCenter1,tZCenter1);

    /*
     * The second sphere parameters
     */
    real tRadius2 = 3.1;
    real tXCenter2 = 0.0;
    real tYCenter2 = 0.0;
    real tZCenter2 = 0.0;
    Sphere<real, size_t, DDRMat, DDSTMat> tLevelSetSphere2
    = Sphere<real, size_t, DDRMat, DDSTMat>(tRadius2,tXCenter2,tYCenter2,tZCenter2);

    /*
     * Place pointers to geometries in a vector
     */
    Cell<Geometry<real,size_t,DDRMat,DDSTMat>*> tGeometryVector = {&tLevelSetSphere1, &tLevelSetSphere2};

    /*
     * The phase table is used to interpret an locations,
     * inside vs outside with respect to a given geometry
     * and interpret that as a material phase. This uses a
     * 2^n material phase table
     */

    Phase_Table<size_t, DDSTMat> tPhaseTable (tNumGeometries,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, DDRMat, DDSTMat> tGeometryEngine
    = Geometry_Engine<real, size_t, DDRMat, DDSTMat> (tGeometryVector,tPhaseTable);

    /*
     * Specify the isocontour of the level set field.
     */
    tGeometryEngine.mThresholdValue = 0.0;

    /*
     * Specify to compute sensitivities
     */
    tGeometryEngine.mComputeDxDp = true;


    /*
     * Specify the model dimension(Note: XTK currently only works in 3D)
     */
    size_t tModelDimension = 3;

    /*
     * Tell the XTK model about:
     *  the model dimension,
     *  the background mesh
     *  the geometry engine
     */
    Model<real, size_t, DDRMat, DDSTMat> tXTKModel(tModelDimension, tBackgroundMesh, tGeometryEngine);

    /*
     * Specify the method to use for decomposing the mesh
     * Note: only the following two methods are implemented
     *
     * A Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 method
     * means the hex8 background cells intersected by the geometries
     * will be regularly subdivided into 24 TET4s
     *
     * A Subdivision_Method::C_HIERARCHY_TET4 method
     * takes the regularly subdivided mesh from method 1
     * and generates a conformal interface mesh using a node
     * ids as a metric for subdivision
     */
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                           Subdivision_Method::C_HIERARCHY_TET4};

    /*
     * Tell the model to decompose the mesh with the specified methods
     */
    tXTKModel.decompose(tDecompositionMethods);

    /*
     * Perform Enrichment
     */

    /*
     * Working with XTK Cut Mesh
     */


}
}
