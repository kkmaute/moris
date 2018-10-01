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
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"

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
// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src

moris::Comm_Manager gMorisComm;

using namespace xtk;
int
main( int    argc,
      char * argv[] )
{

    // Initialize Moris global communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);

    // Geometry Engine Setup -----------------------
    // Using a Levelset Sphere as the Geometry

    real tRadius = 4.25;
    real tXCenter = 5.0;
    real tYCenter = 5.0;
    real tZCenter = 5.0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:10x10x10";
    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

    std::cout<<"Number of nodes = " <<     tCutMeshData->get_num_nodes()<<std::endl;

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/test_new_mtk.e";

    std::cout<<"tMeshOutputFile = "<< tMeshOutputFile<<std::endl;
    tCutMeshData->create_output_mesh(tMeshOutputFile);


    delete tMeshData;

    // finalize moris global communication manager
    gMorisComm.finalize();

    return 0;
}
