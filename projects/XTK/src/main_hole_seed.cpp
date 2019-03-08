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
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_Sphere.hpp"
#include "cl_Discrete_Level_Set.hpp"
#include "cl_MGE_Geometry_Engine.hpp"
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

    moris::uint tNumNodes = tMeshData->get_num_entities(EntityRank::NODE);
    moris::Matrix< moris::DDRMat > tCoords(tNumNodes,3);
    moris::Matrix< moris::DDRMat > tCoordsX(tNumNodes,1);
    moris::Matrix< moris::DDRMat > tCoordsY(tNumNodes,1);
    moris::Matrix< moris::DDRMat > tCoordsZ(tNumNodes,1);
    for(moris::moris_index i = 0; i <(moris::moris_index)tNumNodes; i++)
    {
        moris::Matrix< moris::DDRMat > tNodeCoord =  tMeshData->get_node_coordinate(i);
        tCoordsX(i) = tNodeCoord(0);
        tCoordsY(i) = tNodeCoord(1);
        tCoordsZ(i) = tNodeCoord(2);
        tCoords.get_row(i) = tNodeCoord.matrix_data();
    }

    // figure out bounding box
    moris::real tMaxX = tCoordsX.max();
    moris::real tMinX = tCoordsX.min();
    moris::real tMaxY = tCoordsY.max();
    moris::real tMinY = tCoordsY.min();
    moris::real tMaxZ = tCoordsZ.max();
    moris::real tMinZ = tCoordsZ.min();

    std::cout<<"Bounding Box:"<<std::endl;
    std::cout<<"    x in ["<<tMinX<<","<<tMaxX<<"]"<<std::endl;
    std::cout<<"    y in ["<<tMinY<<","<<tMaxY<<"]"<<std::endl;
    std::cout<<"    z in ["<<tMinZ<<","<<tMaxZ<<"]"<<std::endl;

    moris::real lx = tMaxX-tMinX;
    moris::real ly = tMaxY-tMinY;
    moris::real lz = tMaxZ-tMinZ;

    moris::real xOffset = lx/(tNumX+1);
    moris::real yOffset = ly/(tNumY-1);
    moris::real zOffset = lz/(tNumZ-1);

    std::cout<<"Offsets:"<<std::endl;
    std::cout<<"    x offset ="<<xOffset<<std::endl;
    std::cout<<"    y offset ="<<yOffset<<std::endl;
    std::cout<<"    z offset ="<<zOffset<<std::endl;

    moris::uint tNumSpheres = tNumX*tNumY*tNumZ;

    moris::Cell<moris::Matrix<moris::DDRMat>> tCenters(tNumSpheres);

    moris::uint tCount = 0;
    for(moris::uint i = 0; i <tNumX; i++)
    {
        for(moris::uint j = 0; j <tNumY; j++)
        {
            for(moris::uint k = 0; k<tNumZ; k++)
            {
                moris::Matrix<moris::DDRMat> tSphereCenter = {{tMinX + xOffset +  i*xOffset, tMinY+j*yOffset, tMinZ + k*zOffset  }};
                tCenters(tCount) = tSphereCenter;
                tCount++;
            }
        }
    }

    // construct spheres
    moris::Cell<xtk::Sphere> tSpheres(tNumSpheres);
    for(moris::uint i = 0; i <tNumSpheres; i++)
    {
        tSpheres(i) = Sphere(tRadius,tCenters(i)(0),tCenters(i)(1),tCenters(i)(2));
    }

    // iterate through node to compute a level set value at each node
    Cell<moris::Matrix<moris::DDRMat>> tSphereLSV(tNumNodes);
    moris::Matrix<moris::DDRMat> tMaxSphereLSV(tNumNodes,1);
//    moris::real tNodeVal = 0;
    for(moris::uint i =0; i <tNumNodes; i++)
    {
        tSphereLSV(i) = moris::Matrix<moris::DDRMat>(1,tNumSpheres);

        // iterate through all spheres
        for(moris::uint iSphere =0; iSphere<tNumSpheres; iSphere++)
        {
            tSphereLSV(i)(iSphere) = tSpheres(iSphere).evaluate_field_value_with_coordinate(i,tCoords);
        }

        tMaxSphereLSV(i) = tSphereLSV(i).min();
    }

    // write to mesh
    tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName,EntityRank::NODE, tMaxSphereLSV);

    // write output mesh to exodus
    tMeshData->create_output_mesh(tMeshOutputFileName);





    // finalize MORIS global communication manager
    gMorisComm.finalize();

    delete tMeshData;

    return 0;

}


