/*
 * cl_MSI_Multigrid.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */

#include "catch.hpp"
#include "typedefs.hpp"
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

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_Mesh_Enums.hpp"
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

namespace moris
{
namespace MSI
{
TEST_CASE( "MSI_SPace_Time", "[moris],[MSI],[MSI_Space_Time]" )
{
    if(par_size() == 1 )
    {
        //std::cout<<" Create a 3D mesh of HEX8 "<<std::endl;
        // Create a 3D mesh of HEX8 using MTK ------------------------------------------

        // Define the Interpolation Mesh
        std::string tInterpString = "generated:1x1x1";

        // construct the interpolation and the integration meshes
        mtk::Interpolation_Mesh* tInterpMesh = mtk::create_interpolation_mesh( MeshType::STK, tInterpString, NULL );
        mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh );

        // construct a mesh manager
        mtk::Mesh_Manager tMesh;
        tMesh.register_mesh_pair( tInterpMesh, tIntegMesh );

        //std::cout<<" Create the fem nodes "<<std::endl;
        // Create the fem nodes ------------------------------------------------------

        // number of mesh nodes
        uint tNumOfNodes = tInterpMesh->get_num_nodes();

        //create a cell of fem nodes
        moris::Cell< fem::Node_Base* > tNodes( tNumOfNodes, nullptr );

        // loop over the mesh nodes
        for( uint iNode = 0; iNode < tNumOfNodes; iNode++ )
        {
            // create a fem node for each mesh node
            tNodes( iNode ) = new fem::Node( &tInterpMesh->get_mtk_vertex( iNode ) );
        }

        //std::cout<<" Create the IWGs "<<std::endl;
        // Create the IWGs -----------------------------------------------------------

        // input a cell of IWG types to be created
        Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::HJTEST };

        // number of IWGs to be created
        uint tNumOfIWGs = tIWGTypeList.size();

        // a factory to create the IWGs
        fem::IWG_Factory tIWGFactory;

        // create a cell of IWGs for the problem considered
        moris::Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );

        // loop over the IWG types
        for( uint iIWG = 0; iIWG < tNumOfIWGs; iIWG++)
        {
            // create an IWG with the factory for the ith IWG type
            tIWGs( iIWG ) = tIWGFactory.create_IWGs( tIWGTypeList( iIWG ) );

            // set residual dof type
            tIWGs( iIWG )->set_residual_dof_type( { MSI::Dof_Type::LS1 } );

            // set active dof type
            tIWGs( iIWG )->set_dof_type_list( {{ MSI::Dof_Type::LS1 }} );
        }

        //std::cout<<" Create the elements "<<std::endl;
        // Create the elements -------------------------------------------------------

        // create a cell of FEM sets
        moris::Cell< MSI::Equation_Set * > tFEMSets( 1, nullptr );

        // get the mesh set from the integration mesh
        moris::mtk::Set * tMeshSet = tIntegMesh->get_block_by_index( 0 );

        // create a property user defined info storage
        moris::Cell< moris::Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 1 );

        // create a constitutive user defined info storage
        moris::Cell< moris::Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 1 );

        // create a set
        tFEMSets( 0 ) = new fem::Set( tMeshSet,
                                     fem::Element_Type::BULK,
                                     tIWGs,
                                     tPropertyUserDefinedInfo,
                                     tConstitutiveUserDefinedInfo,
                                     tNodes );

        //std::cout<<" Create the model solver interface "<<std::endl;
        // Create the model solver interface -----------------------------------------

        // FIXME force the communication table
        Matrix< IdMat > tCommunicationTable( 1, 1, 0 );

        // FIXME force the coeff map
        map< moris_id, moris_index > tCoefficientsMap;

        // FIXME number of coeff
        uint tNumCoeff = 1000000;

        moris::MSI::Model_Solver_Interface* tModelSolverInterface
            = new moris::MSI::Model_Solver_Interface( tFEMSets,
                                                      tCommunicationTable,
                                                      tCoefficientsMap,
                                                      tNumCoeff,
                                                      tInterpMesh );

        uint tDofOrder = 1;
        tModelSolverInterface->set_param( "VX" )   = (sint)tDofOrder;
        tModelSolverInterface->set_param( "LS1" )  = (sint)tDofOrder;

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
        delete tInterpMesh;
        delete tIntegMesh;

        // delete node pointers
        for( fem::Node_Base* tNode : tNodes )
        {
            delete tNode;
        }

        // delete IWG pointers
        for( fem::IWG* tIWG : tIWGs )
        {
            delete tIWG;
        }

        // delete model solver interface pointer
        delete tModelSolverInterface;

    }/* if( par_size() */
}/* TEST_CASE */
}
}


