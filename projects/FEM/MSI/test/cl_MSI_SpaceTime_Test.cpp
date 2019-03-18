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
#include "cl_MTK_Mesh.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#include "cl_MSI_Model_Solver_Interface.hpp"


namespace moris
{
    namespace MSI
    {
    TEST_CASE( "MSI_SPace_Time", "[moris],[MSI],[MSI_Space_Time]" )
            {
                if(par_size() == 1 )
                {

                    // Create a 3D mesh of HEX8 using MTK ------------------------------------------
                    std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
                    //------------------------------------------------------------------------------
                    uint aNumElemTypes = 1; // only 1 element type ( quad )
                    uint aNumDim = 3;       // number of spatial dimensions

                    // element connectivity
                    Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4, 5, 6, 7, 8 }};

                    // local to global element map
                    Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};

                    // node coordinates
                    Matrix< DDRMat > aCoords = {{ 0.0, 0.0, 0.0 },
                                                { 1.0, 0.0, 0.0 },
                                                { 1.0, 1.0, 0.0 },
                                                { 0.0, 2.0, 0.0 },
                                                { 0.0, 0.0, 1.0 },
                                                { 1.0, 0.0, 1.0 },
                                                { 1.0, 1.0, 1.0 },
                                                { 0.0, 2.0, 1.0 }};

                    // specify the local to global map
                    Matrix< IdMat > aNodeLocalToGlobal = {{ 1,  2,  3,  4,  5,  6,  7,  8 }};

                    // create mesh MTK database
                    mtk::MtkMeshData tMeshData( aNumElemTypes );
                    tMeshData.CreateAllEdgesAndFaces  = true;
                    tMeshData.SpatialDim              = & aNumDim;
                    tMeshData.ElemConn( 0 )           = & aElementConnQuad;
                    tMeshData.NodeCoords              = & aCoords;
                    tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
                    tMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

                    mtk::Mesh* tMesh = create_mesh( MeshType::STK, tMeshData );

                    //1) Create the fem nodes ------------------------------------------------------
                    std::cout<<" Create the fem nodes "<<std::endl;
                    //------------------------------------------------------------------------------

                    // number of mesh nodes
                    uint tNumOfNodes = tMesh->get_num_nodes();

                    //create a cell of fem nodes
                    moris::Cell< fem::Node_Base* > tNodes( tNumOfNodes, nullptr );

                    // loop over the mesh nodes
                    for( uint k = 0; k < tNumOfNodes; k++ )
                    {
                        // create a fem node for each mesh node
                        tNodes( k ) = new fem::Node( & tMesh->get_mtk_vertex( k ) );
                    }

                    //2) Create the IWGs -----------------------------------------------------------
                    std::cout<<" Create the IWGs "<<std::endl;
                    //------------------------------------------------------------------------------

                    // input a cell of IWG types to be created
                    Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::HJTEST };

                    // number of IWGs to be created
                    uint tNumOfIWGs = tIWGTypeList.size();

                    // a factory to create the IWGs
                    fem::IWG_Factory tIWGFactory;

                    // create a cell of IWGs for the problem considered
                    moris::Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );

                    // loop over the IWG types
                    for( uint i = 0; i < tNumOfIWGs; i++)
                    {
                        // create an IWG with the factory for the ith IWG type
                        tIWGs( i ) = tIWGFactory.create_IWGs( tIWGTypeList( i ) );
                    }

                    //3) Create the elements -------------------------------------------------------
                    std::cout<<" Create the elements "<<std::endl;
                    //------------------------------------------------------------------------------

                    // a factory to create the elements
                    fem::Element_Factory tElementFactory;

                    // ask mesh about number of elements
                    uint tNumOfElements = tMesh->get_num_elems();

                    // create equation objects
                    Cell< MSI::Equation_Object* > tElements( tNumOfElements, nullptr );

                    // loop over the mesh elements
                    for( uint k = 0; k < tNumOfElements; k++ )
                    {
                        // create the element
                        tElements( k ) = tElementFactory.create_element( fem::Element_Type::BULK,
                                                                         & tMesh->get_mtk_cell( k ),
                                                                         tIWGs,
                                                                         tNodes );

                        tElements( k )->mTimeSteps( 0, 0 )= 2;
                    }

                    //4) Create the model solver interface -----------------------------------------
                    std::cout<<" Create the model solver interface "<<std::endl;
                    //------------------------------------------------------------------------------

                    //FIXME force the communication table
                    Matrix< IdMat > tCommunicationTable( 1, 1, 0 );

                    // FIXME: get map from mesh
                    uint tDofOrder = 1;
                    map< moris_id, moris_index > tCoefficientsMap;
                    //tMesh->get_adof_map( tDofOrder, tCoefficientsMap );

                    uint tNumCoeff = 1000000;
                    //= tMesh->get_num_coeffs( 1 )

                    moris::MSI::Model_Solver_Interface* tModelSolverInterface
                        = new moris::MSI::Model_Solver_Interface( tElements,
                                                                  tCommunicationTable,
                                                                  tCoefficientsMap,
                                                                  tNumCoeff,
                                                                  tMesh );

                    tModelSolverInterface->set_param( "VX" )  = (sint)tDofOrder;
                    tModelSolverInterface->set_param( "LS1" )  = (sint)tDofOrder;

                    tModelSolverInterface->finalize();


                    // 8) Clean up -----------------------------------------------------------------
                    std::cout<<" Clean up "<<std::endl;
                    //------------------------------------------------------------------------------
                    delete tMesh;

                    for( uint i = 0; i < tNumOfNodes; i++ )
                    {
                        delete tNodes( i );
                    }

                    for( uint i = 0; i < tNumOfIWGs; i++)
                    {
                        delete tIWGs( i );
                    }

                    for( uint i = 0; i < tNumOfElements; i++ )
                    {
                        delete tElements( i );
                    }

                    delete tModelSolverInterface;
                }/* if( par_size() */
            }/* TEST_CASE */
    }
}


