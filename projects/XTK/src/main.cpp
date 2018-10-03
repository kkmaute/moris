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
#include "cl_Mesh_Tools.hpp"
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

// HMR
#include "cl_HMR_Parameters.hpp"
#define private public
#define protected public
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#undef private
#undef protected
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

    //------------------------------------------------------------------------------
    // Geometry Engine Setup
    //------------------------------------------------------------------------------


    // Using a Levelset Sphere as the Geometry
    real tRadius =  0.1;
    real tXCenter = 0.125;
    real tYCenter = 0.125;
    real tZCenter = 0.125;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

    //------------------------------------------------------------------------------
    // Create Mesh
    //------------------------------------------------------------------------------

    /*!
     * This setup creates a minimal mesh and assumes an element edge length
     * of 1.
     *
     * HMR also supports parameter lists ( see tutorial )
     * For internal tests, using the parameter object is more convenient
     */
    moris::hmr::Parameters tParameters;

    // create a 2x2x2 meshB
    tParameters.set_number_of_elements_per_dimension( Matrix< DDLUMat >{ {2}, {2}, {2} } );

    // create mesh order 1 ( XTK does only support 1st order so far
    tParameters.set_mesh_order( 1 );

    // we do not want tructation for this test
    tParameters.set_bspline_truncation( false );

    // we do not need a buffer if there is no truncation
    tParameters.set_buffer_size( 0 );

    // do not print debug information during test
    tParameters.set_verbose( false );

    //------------------------------------------------------------------------------

    // create an HMR object and perform an arbitrary refinement
    moris::hmr::HMR tHMR( tParameters );
    tHMR.flag_element( 0 );

    // the optional flag resets the pattern
    // ( if it is set, we HMR always starts with a tensor mesh )

    tHMR.perform_refinement( true );
    tHMR.flag_element( 0 );
    tHMR.perform_refinement( false );

    //------------------------------------------------------------------------------

    // These are private functions that I use until MTK can write STK again
    tHMR.mLagrangeMeshes( 1 )->save_faces_to_vtk("Faces.vtk");

    if( tParameters.get_number_of_dimensions() == 3)
    {
        tHMR.mLagrangeMeshes( 1 )->save_edges_to_vtk("Edges.vtk");
    }
    tHMR.mLagrangeMeshes( 1 )->save_to_vtk("Elements.vtk");

    //------------------------------------------------------------------------------
    // create mesh interface
    mtk::Mesh * tMesh = tHMR.create_output_mesh();

    //------------------------------------------------------------------------------


    /*!
     * Setup xtk model with HMR MTK mesh
     */
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMesh,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                Subdivision_Method::C_HIERARCHY_TET4};

    // Decompose the mesh
    tXTKModel.decompose(tDecompositionMethods);

    // Get the XTK mesh as an MTK mesh
    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();


    /*!
      * Write the output mesh
      * \code{.cpp}
      * std::string tPrefix = std::getenv("XTKOUTPUT");
      * std::string tMeshOutputFile = tPrefix + "/hmr_to_xtk_intersected.e";
      * tCutMeshData->create_output_mesh(tMeshOutputFile);
      * \endcode
      */
    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/hmr_to_xtk_intersected.e";

    tCutMeshData->create_output_mesh(tMeshOutputFile);


//
//    // --------------------------------------------
//    // Repeat with STK
//    // --------------------------------------------
//
//    // geometry engine using same sphere/phase table
//    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngineSTK(tLevelsetSphere,tPhaseTable);
//
//    // Create Mesh ---------------------------------
//    std::string tMeshFileName = "generated:1x1x1";
//    moris::mtk::Mesh* tSTKMesh = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );
//
//    // Setup XTK Model -----------------------------
//    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModelSTK(tModelDimension,tSTKMesh,tGeometryEngineSTK);
//
//    // Decompose the mesh
//    tXTKModelSTK.decompose(tDecompositionMethods);
//
//    // Output decomposed STK mesh
//    moris::mtk::Mesh* tSTKCut = tXTKModelSTK.get_output_mesh();
//    tMeshOutputFile = tPrefix + "/stk_cut_mesh.e";
//    tSTKCut->create_output_mesh(tMeshOutputFile);


    // Clean up
//    delete tSTKCut;
//    delete tMesh;
//    delete tCutMeshData;

    // finalize moris global communication manager
    gMorisComm.finalize();

    return 0;
}
