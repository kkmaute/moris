/*
 * main.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: doble
 */


// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "HDF5_Tools.hpp"

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
#include "fn_compute_interface_surface_area.hpp"

//------------------------------------------------------------------------------

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

    Cell<real> tRadii = {4,8,12,16,19};

    Matrix<DDRMat> tFullTimes(tRadii.size(),1);
    Matrix<DDRMat> tDecompTimes(tRadii.size(),1);
    Matrix<DDRMat> tSurfaceArea(tRadii.size(),1);

    for(moris::uint i = 0 ; i < tRadii.size(); i++)
    {
        // start timing on this decomposition
        std::clock_t start = std::clock();

        std::string tMeshInputFileName  = "generated:40x40x40";
        moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, tMeshInputFileName );

        real tRadius = tRadii(i);
        real tXCenter = 20.0;
        real tYCenter = 20.0;
        real tZCenter = 20.0;
        Sphere tLevelSetSphere(tRadius,tXCenter,tYCenter,tZCenter);
        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelSetSphere,tPhaseTable);


        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = true;

        // Do the cutting
        std::clock_t startdecomp = std::clock();
        tXTKModel.decompose(tDecompositionMethods);
        tDecompTimes(i) = (std::clock() - startdecomp);
//        tXTKModel.unzip_interface();

        tFullTimes(i) = (std::clock() - start);

        //    // Compute volume
        moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();
        moris::real tParentPhase0Vol = compute_non_intersected_parent_element_volume_by_phase(0,tNodeCoords,tXTKModel);
        moris::real tParentPhase1Vol = compute_non_intersected_parent_element_volume_by_phase(1,tNodeCoords,tXTKModel);
        moris::real tChildPhase0Vol  = compute_child_element_volume_by_phase(0,tNodeCoords,tXTKModel);
        moris::real tChildPhase1Vol  = compute_child_element_volume_by_phase(1,tNodeCoords,tXTKModel);
        moris::real tMySurfaceArea   = compute_interface_surface_area(tNodeCoords,tXTKModel);

        moris::real tMyPhase0Volume = tParentPhase0Vol + tChildPhase0Vol;
        moris::real tMyPhase1Volume = tParentPhase1Vol + tChildPhase1Vol;

        // Collect all volumes
        moris::real tGlbVolume0 = 0.0;
        sum_all(tMyPhase0Volume,tGlbVolume0);

        // Collect all volumes
        moris::real tGlbVolume1 = 0.0;
        sum_all(tMyPhase1Volume,tGlbVolume1);


        // Collect all volumes
        moris::real tGlbSurf = 0.0;
        sum_all(tMySurfaceArea,tGlbSurf);
        tSurfaceArea(i) = tGlbSurf;

        if(par_rank() == 0)
        {
            std::cout<<"Phase 0 Volume: "<<tGlbVolume0<<std::endl;
            std::cout<<"Phase 1 Volume: "<<tGlbVolume1<<std::endl;
            std::cout<<"Surface Area: "<<tGlbSurf<<std::endl;
        }



        delete tMeshData;

    }


    print(tFullTimes,"tTimes");

    std::string tHdfTime = "./surface_area_times_"+std::to_string(par_size()) + ".hdf5";
    hid_t tHDFId = create_hdf5_file(tHdfTime.c_str());

    herr_t tErr;
    save_matrix_to_hdf5_file(tHDFId,"full_time",tFullTimes,tErr);
    save_matrix_to_hdf5_file(tHDFId,"decomp_time",tDecompTimes,tErr);
    save_matrix_to_hdf5_file(tHDFId,"surf_area",tSurfaceArea,tErr);
    // finalize MORIS global communication manager
    gMorisComm.finalize();





    //    // Generate mesh from string
    //    std::string tMeshInputFileName  = argv[1];
    //    std::string tMeshOutputFileName = argv[2];
    //    std::string tLSFName            = argv[3];
    //
    //    std::cout<<"Input: "<<tMeshInputFileName<<std::endl;
    //    std::cout<<"Output: "<<tMeshOutputFileName<<std::endl;
    //    std::cout<<"Field: "<<tLSFName<<std::endl;
    //
    //    // Create mesh with supplementary data
    //    moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, tMeshInputFileName );
    //
    //
    //    xtk::Discrete_Level_Set tLevelSetMesh(tMeshData,{tLSFName});
    //
    //    // Tell the geometry engine about the discrete field mesh and how to interpret phases
    //    Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    //    Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);
    //
    //    tGeometryEngine.mThresholdValue = 45;
    //
    //
    //    // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
    //    size_t tModelDimension = 3;
    //    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    //    Model tXTKModel(tModelDimension,tLevelSetMesh.get_level_set_mesh(),tGeometryEngine);
    //    tXTKModel.mSameMesh = true;
    //    tXTKModel.mVerbose = true;
    //
    //    // Do the cutting
    //    tXTKModel.decompose(tDecompositionMethods);
    //
    //    tXTKModel.unzip_interface();

    // Sensitivity
    //    tXTKModel.compute_sensitivity();

    //    // Tell XTK to construct the output meshd ata
    //    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
    //
    //    // Write output mesh
    //    tCutMeshData->create_output_mesh(tMeshOutputFileName);
    //
    //    delete tCutMeshData;

    return 0;

}

