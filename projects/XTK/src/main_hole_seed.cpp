/*
 * main_hole_seed.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: doble
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
#include "cl_XTK_Hole_Seeder.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src
#include "fn_compute_xtk_model_volumes.hpp"

#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

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
    std::string tMeshInputFileName  = argv[1];
    std::string tMeshOutputFileName = argv[2];
    std::string tLSFName            = argv[3];
    moris::real tRadius             = 2;
    moris::uint tNumX               = 6;
    moris::uint tNumY               = 6;
    moris::uint tNumZ               = 6;

    std::cout<<"Input: "<<tMeshInputFileName<<std::endl;
    std::cout<<"Output: "<<tMeshOutputFileName<<std::endl;
    std::cout<<"Field: "<<tLSFName<<std::endl;
    std::cout<<"Radius: "<<tRadius<<std::endl;
    std::cout<<"Num in x: "<<tNumX<<std::endl;
    std::cout<<"Num in y: "<<tNumY<<std::endl;
    std::cout<<"Num in z: "<<tNumZ<<std::endl;

    // Declare scalar node field
    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
    std::string tFieldName1 = tLSFName;
    tNodeField1.set_field_name(tFieldName1);
    tNodeField1.set_field_entity_rank(EntityRank::NODE);


    // Initialize field information container
    moris::mtk::MtkFieldsInfo tFieldsInfo;

    // Place the node field into the field info container
    moris::mtk::add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

    // Declare some supplementary fields
    moris::mtk::MtkMeshData tSuppMeshData;
    tSuppMeshData.FieldsInfo = &tFieldsInfo;


    // Load up mesh
    moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, tMeshInputFileName,&tSuppMeshData ,false );
    tMeshData->mVerbose = true;

    // construct hole seeder
    xtk::Hole_Seeder tHoleSeeder(tMeshData,tRadius,tRadius,tRadius,2,tNumX,tNumY,tNumZ);

    // perform hole seeding
    tHoleSeeder.seed_field();


    // write to mesh
    tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName,EntityRank::NODE, tHoleSeeder.get_seeded_field());

    // write output mesh to exodus
    tMeshData->create_output_mesh(tMeshOutputFileName);

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    delete tMeshData;

    return 0;


}


