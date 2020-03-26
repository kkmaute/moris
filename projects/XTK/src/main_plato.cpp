/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */
//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// XTK
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_Sphere.hpp"
#include "cl_Discrete_Level_Set.hpp"
#include "cl_MGE_Geometry_Engine.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src
#include "fn_compute_xtk_model_volumes.hpp"

//------------------------------------------------------------------------------
#include "cl_GEN_Discrete_Level_Set.hpp"

// select namespaces
using namespace moris;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;


int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    // Generate mesh from string
    std::string tMeshInputFileName  = argv[1];
    std::string tMeshOutputFileName = argv[2];
    std::string tLSFName            = argv[3];

    std::cout<<"Input: "<<tMeshInputFileName<<std::endl;
    std::cout<<"Output: "<<tMeshOutputFileName<<std::endl;
    std::cout<<"Field: "<<tLSFName<<std::endl;

    // Create mesh with supplementary data
    moris::mtk::Interpolation_Mesh* tMeshData   = moris::mtk::create_interpolation_mesh(MeshType::STK, tMeshInputFileName );


    moris::ge::Discrete_Level_Set tLevelSetMesh(tMeshData,{tLSFName});

    // Tell the geometry engine about the discrete field mesh and how to interpret phases
//    Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//    Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);
    moris::ge::GEN_Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);


    // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
    tXTKModel.mSameMesh = true;
    tXTKModel.mVerbose  =  false;

    // Do the cutting
    tXTKModel.decompose(tDecompositionMethods);

    // Sensitivity
    tXTKModel.compute_sensitivity();


    // Specfiy the output options
    Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;
    tOutputOptions.mHaveInterface = false;
    tOutputOptions.mSeparateInterfaceBlock = false;
    tOutputOptions.mPackageDxDpSparsely = true;
    // Specify there are 2 possible phases
    size_t tNumPhases = 2;

    //Say I only want to output phase 1
    Cell<size_t> tPhasesToOutput = {1};

    // Give this information to the output options
    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    // Tell XTK to construct the output meshd ata
    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

    // Write output mesh
    tCutMeshData->create_output_mesh(tMeshOutputFileName);

    // Compute volume
    moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();
    moris::real tParentPhase0Vol = compute_non_intersected_parent_element_volume_by_phase(0,tNodeCoords,tXTKModel);
    moris::real tParentPhase1Vol = compute_non_intersected_parent_element_volume_by_phase(1,tNodeCoords,tXTKModel);
    moris::real tChildPhase0Vol  = compute_child_element_volume_by_phase(0,tNodeCoords,tXTKModel);
    moris::real tChildPhase1Vol  = compute_child_element_volume_by_phase(1,tNodeCoords,tXTKModel);

    std::cout<<"XTK Volumes: "<<std::endl;
    std::cout<<"Phase 0: "<<tParentPhase0Vol + tChildPhase0Vol<<std::endl;
    std::cout<<"Phase 1: "<<tParentPhase1Vol + tChildPhase1Vol<<std::endl;

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    delete tMeshData;
    delete tCutMeshData;

    return 0;

}
