/*
 * cl_XTK_Tet_Decomposition.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: ktdoble
 */

#include <memory>
#include <mpi.h>

#include "catch.hpp"
#include "paths.hpp"

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "cl_Mesh_Enums.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes
#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"


#include "geometry/cl_Mesh_Field_Geometry.hpp"
#include "geometry/cl_Multi_Cylinder.hpp"
#include "cl_Sphere.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"




namespace xtk
{
TEST_CASE("Tet background mesh analytic sphere","[TET_START_SPHERE]")
{
    /*
     * Set up:
     * MatrixFactory,
     * Geometry,
     * Geometry Engine,
     * Mesh
     */

    real tRadius  = 3;
    real tXCenter = 0.0;
    real tYCenter = 0.0;
    real tZCenter = 0.0;
    Sphere tLevelSetSphere(tRadius,tXCenter,tYCenter,tZCenter);

    tRadius  = 3;
    tXCenter = 0.0;
    tYCenter = 0.0;
    tZCenter = 1.0;
    Sphere tLevelSetSphere2(tRadius,tXCenter,tYCenter,tZCenter);
    Geometry_Engine tGeometryEngine({&tLevelSetSphere,&tLevelSetSphere2});
    tGeometryEngine.mComputeDxDp = false;
    // Load the mesh
    std::string tPrefix = moris::get_base_moris_dir();
    std::string tMeshFileName = tPrefix + "/TestExoFiles/tet_cube_mesh.e";
    moris::Cell<std::string> tFieldNames;
    mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName,tFieldNames,true);

    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);


    /*
     * Convert tet4 to tet10
     */

//    tXTKModel.convert_mesh_tet4_to_tet10();

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    // Specify there are 2 possible phases
    size_t tNumPhases = 4;

    // Output phase 0 and 1
    Cell<size_t> tPhasesToOutput = {0,1};

    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

    tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/tet_cube_mesh_output.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);
}



TEST_CASE("Tet background mesh analytic cylinder","[TET_START_CYLINDER]")
{

    /*
     * Set up:
     * MatrixFactory,
     * Geometry,
     * Geometry Engine,
     * Mesh
     */
    real tCordLength = 10;
    Cell<Cell<real>> tCenter = {{0,0,0}};
    Cell<real> tRadius = {tCordLength/4};
    Cell<real> tLength = {1.1*tCordLength};
    Cell<Cell<real>> tAxis   = {{1,0,0}};

    Multi_Cylinder<real, size_t, moris::DDRMat, moris::DDSTMat> tMultiCylinder(tCenter,tRadius,tLength, tAxis);
    Geometry_Engine tGeometryEngine(tMultiCylinder);

    // Load the mesh
    std::string tPrefix = moris::get_base_moris_dir();
    std::string tMeshFileName = tPrefix + "/TestExoFiles/tet_cube_mesh.e";
    moris::Cell<std::string> tFieldNames;
    mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName,tFieldNames,true);

    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);


    /*
     * Get the output mesh and write to exodus file
     */


    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    // Specify there are 2 possible phases
    size_t tNumPhases = 2;

    // Say I only want to output phase 1
    Cell<size_t> tPhasesToOutput = {0,1};

    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

    tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/tet_cube_mesh_cylinder_discrete_output.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);
}

TEST_CASE("Tet background mesh discrete cylinder","[TET_START_CYLINDER_DISCRETE]")
{

    /*
     * Set up:
     * MatrixFactory,
     * Geometry,
     * Geometry Engine,
     * Mesh
     */

    // Specify the mesh information and how to build it
    std::string tPrefix = moris::get_base_moris_dir();
    std::string tMeshFileName = tPrefix + "/TestExoFiles/tet_cube_mesh.e";
    moris::Cell<std::string> tFieldNames;
    mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;

    // Setup Analytic LSF, then Discretize it

    // Cylinder Parameters
    real tCordLength = 10;
    Cell<Cell<real>> tCenter = {{0,0,0}};
    Cell<real> tRadius = {tCordLength/4};
    Cell<real> tLength = {1.1*tCordLength};
    Cell<Cell<real>> tAxis   = {{1,0,0}};
    Multi_Cylinder<real, size_t, moris::DDRMat, moris::DDSTMat> tMultiCylinder(tCenter,tRadius,tLength, tAxis);

    // Discretize the cylinder on the mesh field "lsf"
    Cell<std::string> tScalarFieldNames = {"lsf"};
    Cell<xtk::Geometry<real, size_t, moris::DDRMat, moris::DDSTMat>*> tLevelSetFunctions = {&tMultiCylinder};
    Mesh_Field_Geometry tLevelSetMesh(tLevelSetFunctions,tMeshFileName,tScalarFieldNames,tMeshBuilder);

    // Tell the geometry engine about the discrete field mesh and how to interpret phases
    Geometry_Engine tGeometryEngine(tLevelSetMesh);

    // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGeometryEngine);

    // Do the cutting
    tXTKModel.decompose(tDecompositionMethods);


    // Specfiy the output options
    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    // Specify there are 2 possible phases
    size_t tNumPhases = 2;

    // Say I only want to output phase 0 (inside the cylinder)
    Cell<size_t> tPhasesToOutput = {0};

    // Give this information to the output options
    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    // Tell XTK to construct the output meshd ata
    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

    // Write output mesh
    tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/tet_cube_mesh_cylinder_output.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);
}



}
