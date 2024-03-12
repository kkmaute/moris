/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_SpaceTime_Test.cpp
 *
 */

#include "catch.hpp"
#include "moris_typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#define protected public
#define private   public
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#undef protected
#undef private

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set.hpp"                      //FEM/INT/src

#include "fn_PRM_MSI_Parameters.hpp"                      //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"

namespace moris
{
namespace MSI
{
TEST_CASE( "MSI_SPace_Time", "[moris],[MSI],[MSI_Space_Time]" )
{
    if(par_size() == 1 )
    {
        // Define the Interpolation Mesh
        std::string tInterpString = "generated:1x1x1";

        // construct the interpolation and the integration meshes
        mtk::Interpolation_Mesh* tInterpMesh = mtk::create_interpolation_mesh( mtk::MeshType::STK, tInterpString, NULL );
        mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( mtk::MeshType::STK, tInterpMesh );

        // construct a mesh manager
        mtk::Mesh_Manager tMesh;
        tMesh.register_mesh_pair( tInterpMesh, tIntegMesh );

        //std::cout<<" Create the fem nodes "<<std::endl;
        // Create the fem nodes ------------------------------------------------------

        // number of mesh nodes
        uint tNumOfNodes = tInterpMesh->get_num_nodes();

        //create a cell of fem nodes
        Vector< fem::Node_Base* > tNodes( tNumOfNodes, nullptr );

        // loop over the mesh nodes
        for( uint iNode = 0; iNode < tNumOfNodes; iNode++ )
        {
            // create a fem node for each mesh node
            tNodes( iNode ) = new fem::Node( &tInterpMesh->get_mtk_vertex( iNode ) );
        }

        //std::cout<<" Create the IWGs "<<std::endl;
        // Create the IWGs -----------------------------------------------------------
        // create a L2 IWG
        fem::IWG_Factory tIWGFactory;
        std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::HJTEST );
        tIWG->set_residual_dof_type( { { MSI::Dof_Type::LS1 } } );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::LS1 }}, mtk::Leader_Follower::LEADER );

        //std::cout<<" Create the elements "<<std::endl;
        // Create the elements -------------------------------------------------------

        // create a cell of FEM sets
        Vector< MSI::Equation_Set * > tFEMSets( 1, nullptr );

        // get the mesh set from the integration mesh
        moris::mtk::Set * tMeshSet = tIntegMesh->get_set_by_index( 0 );

        // define set info
         Vector< fem::Set_User_Info > tSetInfo( 1 );
         tSetInfo( 0 ).set_mesh_index( 0 );
         tSetInfo( 0 ).set_IWGs( { tIWG } );

        // create a set
         tFEMSets( 0 ) = new fem::Set( nullptr,
                                       tMeshSet,
                                       tSetInfo( 0 ),
                                       tNodes );

        //std::cout<<" Create the model solver interface "<<std::endl;
        // Create the model solver interface -----------------------------------------

        // FIXME force the communication table
        Matrix< IdMat > tCommunicationTable( 1, 1, 0 );

        // FIXME force the coeff map
        map< moris_id, moris_index > tCoefficientsMap;

        // FIXME number of coeff
        uint tNumCoeff = 1000000;

        sint tDofOrder = 1;
        moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
        tMSIParameters.set( "VX", tDofOrder );
        tMSIParameters.set( "LS1", tDofOrder );

        moris::MSI::Model_Solver_Interface* tModelSolverInterface
            = new moris::MSI::Model_Solver_Interface( tMSIParameters,
                                                      tFEMSets,
                                                      tCommunicationTable,
                                                      tCoefficientsMap,
                                                      tNumCoeff,
                                                      tInterpMesh );

        tModelSolverInterface->get_dof_manager()->set_time_levels_for_type( MSI::Dof_Type::LS1, 2 );

        tFEMSets( 0 )->finalize( tModelSolverInterface );

        tModelSolverInterface->finalize();

        CHECK( equal_to( tModelSolverInterface->mDofMgn.mPdofHostList.size(), 8 ) );
        CHECK( equal_to( tModelSolverInterface->mDofMgn.mAdofList.size(), 16 ) );

        moris::uint tCounter = 0;
        moris::sint tCounter_2 = 0;
        for( uint k = 0; k < tModelSolverInterface->mDofMgn.mAdofList.size(); k++ )
        {
            CHECK( equal_to( tModelSolverInterface->mDofMgn.mAdofList( k )->mAdofId, k ) );
            CHECK( equal_to( tModelSolverInterface->mDofMgn.mAdofList( k )->mAdofExternalInd, tCounter++ ) );
            CHECK( equal_to( tModelSolverInterface->mDofMgn.mAdofList( k )->mAdofTypeTimeIdentifier, tCounter_2 ) );

            if( tCounter == tModelSolverInterface->mDofMgn.mPdofHostList.size() )
            {
                tCounter = 0;
                tCounter_2 = 1;
            }
        }

        //std::cout<<" Clean up "<<std::endl;
        // Clean up -----------------------------------------------------------------

        // delete mesh pointer
        //delete tInterpMesh;
        delete tIntegMesh;

        // delete node pointers
        for( fem::Node_Base* tNode : tNodes )
        {
            delete tNode;
        }

        // delete model solver interface pointer
        delete tModelSolverInterface;

    }/* if( par_size() */
}/* TEST_CASE */
}
}

