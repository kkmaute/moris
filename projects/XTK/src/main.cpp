/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */



// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"

// XTKL: Geometry  Include
#include "ios/cl_Logger.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Discrete_Level_Set.hpp"

#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "xtk/fn_write_element_ownership_as_field.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "geometry/cl_Composite_Fiber.hpp"
#include "geometry/cl_Gyroid.hpp"
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

#include "linalg_typedefs.hpp"

#include "cl_Param_List.hpp" // CON/src

// MPI Header
#include <mpi.h>


#include <iostream>
#include <ctime>
// ---------------------------------------------------------------------


using namespace xtk;
int
main( int    argc,
      char * argv[] )
{

    MPI_Init(&argc,&argv);

    std::clock_t startf;

    startf = std::clock();

    mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;


    // Set up parameters
    moris::Param_List< boost::variant<real, std::string> > tParams;

    tParams.insert( "rad"  , 21.1 );
    tParams.insert( "xcenter", 20.0 );
    tParams.insert( "ycenter", 20.0 );
    tParams.insert( "zcenter", 20.0 );
    tParams.insert( "mesh_in", "generated:50x50x50" );


    // Geometry Engine Setup -----------------------
    // Using a Level Set Sphere as the Geometry
    real tRadius = tParams.get< real >( "rad" );
    real tXCenter = tParams.get< real >( "xcenter" );
    real tYCenter = tParams.get< real >( "ycenter" );
    real tZCenter = tParams.get< real >( "zcenter" );
    Sphere<real, size_t, moris::DDRMat, moris::DDSTMat> tLevelSetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    std::string tLevelSetMeshFileName = tParams.get< std::string >( "mesh_in" );;
    xtk::Cell<std::string> tScalarFieldNames = {"LEVEL_SET_SPHERE"};
    xtk::Cell<xtk::Geometry<xtk::real, xtk::size_t, moris::DDRMat, moris::DDSTMat>*> tLevelSetFunctions = {&tLevelSetSphere};
    xtk::Discrete_Level_Set<xtk::real, xtk::size_t, moris::DDRMat, moris::DDSTMat> tLevelSetMesh(tLevelSetFunctions,tLevelSetMeshFileName,tScalarFieldNames,tMeshBuilder);
    Phase_Table<size_t, moris::DDSTMat> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, moris::DDRMat, moris::DDSTMat> tGeometryEngine(tLevelSetMesh,tPhaseTable);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, moris::DDRMat, moris::DDSTMat> tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGeometryEngine);

    tXTKModel.mSameMesh = true;

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                           Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;
    tOutputOptions.change_phases_to_output(2,Cell<size_t>(1,1));

    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/unit_sandwich_sphere_05_threshold_discrete.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);


    std::cout << "Total Time: " << (std::clock() - startf) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;


    MPI_Finalize();

    return 0;
}
