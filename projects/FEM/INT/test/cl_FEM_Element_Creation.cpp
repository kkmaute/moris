/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Creation.cpp
 *
 */

#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Element_Creation", "[moris],[fem],[ElementCreate]" )
        {
            //create a 2D MORIS mesh of quad4's using MTK database ( from Noah )
            //------------------------------------------------------------------------------
            uint aNumElemTypes = 1; // quad
            uint aNumDim = 2;       // specify number of spatial dimensions

            // specify element connectivity of quad for mesh
            Matrix< IdMat > aElementConnQuad = {{ 1, 2, 5, 8 },
                                                { 2, 3, 4, 5 },
                                                { 8, 5, 6, 7 },
                                                { 5, 4, 9, 6 }};

            // specify the local to global element map for quads
            Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1, 2, 3, 4 }};

            // Node coordinate matrix
            Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
                                        { 1.0, 0.0 },
                                        { 2.0, 0.0 },
                                        { 2.0, 1.0 },
                                        { 1.0, 1.0 },
                                        { 1.0, 2.0 },
                                        { 0.0, 2.0 },
                                        { 0.0, 1.0 },
                                        { 2.0, 2.0 }};

           // specify the local to global map
           Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};

           //------------------------------------------------------------------------------
           // create MORIS mesh using MTK database
           mtk::MtkMeshData aMeshData( aNumElemTypes );
           aMeshData.CreateAllEdgesAndFaces  = true;
           aMeshData.SpatialDim              = & aNumDim;
           aMeshData.ElemConn( 0 )           = & aElementConnQuad;
           aMeshData.NodeCoords              = & aCoords;
           aMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
           aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

           mtk::Mesh* tMesh2D_Quad4 = create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

            //1) Create the fem nodes -----------------------------------------------------
           std::cout<<" Create the fem nodes "<<std::endl;
            //------------------------------------------------------------------------------
            // number of mesh nodes
            uint tNumOfNodes = tMesh2D_Quad4->get_num_nodes();

            //create a cell of fem nodes
            moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

            // loop over the mesh nodes
            for( uint k = 0; k < tNumOfNodes; k++ )
            {
                // create a fem node for each mesh node
                tNodes( k ) = new fem::Node( & tMesh2D_Quad4->get_mtk_vertex( k ) );
            }

            //2) Create the IWGs -----------------------------------------------------------
            std::cout<<" Create the IWGs "<<std::endl;
            //------------------------------------------------------------------------------
            // input a cell of IWG types to be created
            Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::L2 };

            // number of IWGs to be created
            uint tNumOfIWGs = tIWGTypeList.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs for the problem considered
            Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );

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
            Element_Factory tElementFactory;

            // ask mesh about number of elements
            uint tNumOfElements = tMesh2D_Quad4->get_num_elems();

            // create equation objects
            Cell< MSI::Equation_Object* > tElements( tNumOfElements, nullptr );

            // loop over the mesh elements
            for( uint k = 0; k < tNumOfElements; k++ )
            {
                // create the element
                tElements( k ) = tElementFactory.create_cluster( Element_Type::BULK,
                                                                 & tMesh2D_Quad4->get_mtk_cell( k ),
                                                                 tIWGs,
                                                                 tNodes );
            }

            //4) Create the model solver interface -----------------------------------------
            std::cout<<" Create the model solver interface "<<std::endl;
            //------------------------------------------------------------------------------
            //FIXME force the communication table
            Matrix< IdMat > tCommunicationTable( 1, 1, 0 );

            // FIXME: get map from mesh
            uint tDofOrder = 1;
            map< moris_id, moris_index > tCoefficientsMap;
            //tMesh2D_Quad4->get_adof_map( tDofOrder, tCoefficientsMap );

            uint tNumCoeff = 1000000;
            //= tMesh2D_Quad4->get_num_coeffs( 1 )

            moris::MSI::Model_Solver_Interface* tModelSolverInterface
                = new moris::MSI::Model_Solver_Interface( tElements,
                                                          tCommunicationTable,
                                                          tCoefficientsMap,
                                                          tNumCoeff,
                                                          tMesh2D_Quad4 );

            tModelSolverInterface->set_param("L2") = (sint)tDofOrder;

            tModelSolverInterface->finalize();

            // calculate AdofMap
            Matrix< DDUMat > tAdofMap = tModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            //4) Create solver interface ---------------------------------------------------
            std::cout<<" Create solver interface "<<std::endl;
            //------------------------------------------------------------------------------

            MSI::MSI_Solver_Interface * tSolverInterface
                = new moris::MSI::MSI_Solver_Interface( tModelSolverInterface );

            // 5) Create Nonlinear Problem -------------------------------------------------
            std::cout<<" Create Nonlinear Problem "<<std::endl;
            //------------------------------------------------------------------------------

            NLA::Nonlinear_Problem* tNonlinearProblem
                = new NLA::Nonlinear_Problem( tSolverInterface );

            // 6) Create Solvers and solver manager ----------------------------------------
            std::cout<<" Create Solvers and solver manager "<<std::endl;
            //------------------------------------------------------------------------------

            // create factory for nonlinear solver
            NLA::Nonlinear_Solver_Factory tNonlinFactory;

            // create nonlinear solver
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
                = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            // create factory for linear solver
            dla::Solver_Factory  tSolFactory;

            // create linear solver
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
                = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            // set default parameters for linear solver
            tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

            // create solver manager
            dla::Linear_Solver*    tLinSolver       = new dla::Linear_Solver();
            NLA::Nonlinear_Solver* tNonlinearSolver = new NLA::Nonlinear_Solver();

            // set manager and settings
            tNonlinearSolverAlgorithm->set_linear_solver( tLinSolver );

            // set first solver
            tLinSolver->set_linear_algorithm( 0, tLinearSolverAlgorithm );

            tNonlinearSolver->set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // 7) Solve --------------------------------------------------------------------
            std::cout<<" Solve "<<std::endl;
            //------------------------------------------------------------------------------
            Matrix<DDRMat> tSolution1;

            // call solver
            tNonlinearSolver->solve( tNonlinearProblem );

            // temporary array for solver
            Matrix< DDRMat > tSolution;
            tNonlinearSolverAlgorithm->get_full_solution( tSolution );

            // get length of array
            uint tLength = tSolution.length();

            // rearrange data into output
            tSolution1.set_size( tLength, 1 );

            for( uint k = 0; k < tLength; k++ )
            {
                tSolution1( k ) = tSolution( tAdofMap( k ) );
            }

            //Matrix< DDRMat > tEqnObjRHS;
            //Dist_Vector * tSolutionVector;

            //tElements( 0 )->get_equation_obj_residual( tEqnObjRHS, tSolutionVector );

            //clean up----------------------------------------------------------------------
            std::cout<<" Clean up "<<std::endl;
            //------------------------------------------------------------------------------
           delete tMesh2D_Quad4;

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

           delete tSolverInterface;

           delete tNonlinearProblem;

           delete tLinSolver;

           delete tNonlinearSolver;

//            //------------------------------------------------------------------------------
//            // create nodes
//            mtk::Vertex* tVertex1_1 = new NodeProxy( 0.0, -1.0, 0 );
//            mtk::Vertex* tVertex1_2 = new NodeProxy( 1.0, -1.0, 1 );
//            mtk::Vertex* tVertex1_3 = new NodeProxy( 1.0,  0.0, 2 );
//            mtk::Vertex* tVertex1_4 = new NodeProxy( 0.0,  0.0, 3 );
//
//            mtk::Vertex* tVertex2_3 = new NodeProxy( 1.0,  1.0, 4 );
//            mtk::Vertex* tVertex2_4 = new NodeProxy( 0.0,  1.0, 5 );
//
//            mtk::Vertex* tVertex3_3 = new NodeProxy( 1.0,  2.0, 6 );
//            mtk::Vertex* tVertex3_4 = new NodeProxy( 0.0,  2.0, 7 );
//
//            mtk::Vertex* tVertex4_3 = new NodeProxy( 1.0,  3.0, 8 );
//            mtk::Vertex* tVertex4_4 = new NodeProxy( 0.0,  3.0, 9 );
//
//            moris::Cell< mtk::Vertex* > allNodes( 10 );
//            allNodes( 0 ) = tVertex1_1;
//            allNodes( 1 ) = tVertex1_2;
//            allNodes( 2 ) = tVertex1_3;
//            allNodes( 3 ) = tVertex1_4;
//            allNodes( 4 ) = tVertex2_3;
//            allNodes( 5 ) = tVertex2_4;
//            allNodes( 6 ) = tVertex3_3;
//            allNodes( 7 ) = tVertex3_4;
//            allNodes( 8 ) = tVertex4_3;
//            allNodes( 9 ) = tVertex4_4;
//
//            //------------------------------------------------------------------------------
//            // create cells of nodes for element
//
//            // cell containing nodes for element 1
//            moris::Cell< mtk::Vertex* > Name1(4);
//            Name1(0) = tVertex1_1;
//            Name1(1) = tVertex1_2;
//            Name1(2) = tVertex1_3;
//            Name1(3) = tVertex1_4;
//
//            // cell containing nodes for element 2
//            moris::Cell< mtk::Vertex* > Name2(4);
//            Name2(0) = tVertex1_4;
//            Name2(1) = tVertex1_3;
//            Name2(2) = tVertex2_3;
//            Name2(3) = tVertex2_4;
//
//            // cell containing nodes for element 1
//            moris::Cell< mtk::Vertex* > Name3(4);
//            Name3(0) = tVertex2_4;
//            Name3(1) = tVertex2_3;
//            Name3(2) = tVertex3_3;
//            Name3(3) = tVertex3_4;
//
//            // cell containing nodes for element 1
//            moris::Cell< mtk::Vertex* > Name4(4);
//            Name4(0) = tVertex3_4;
//            Name4(1) = tVertex3_4;
//            Name4(2) = tVertex4_3;
//            Name4(3) = tVertex4_4;
//
//            //------------------------------------------------------------------------------
//            // create element as mtk::Cell
//            mtk::Cell* tElement1 = new ElementProxy( Name1, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
//            mtk::Cell* tElement2 = new ElementProxy( Name2, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
//            mtk::Cell* tElement3 = new ElementProxy( Name3, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
//            mtk::Cell* tElement4 = new ElementProxy( Name4, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR);
//
//            //------------------------------------------------------------------------------
//            // create a mesh
//            moris::Cell< mtk::Cell* > Elems(4);
//            Elems(0) = tElement1;
//            Elems(1) = tElement2;
//            Elems(2) = tElement3;
//            Elems(3) = tElement4;
//
//            // 1) Create the fem nodes -----------------------------------------------------
//            std::cout<<" Create the fem nodes "<<std::endl;
//            //------------------------------------------------------------------------------
//            // number of mesh nodes
//            uint tNumOfNodes = allNodes.size();
//
//            //create a celle of fem nodes
//            moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );
//
//            // loop obver the mesh nodes
//            for( uint k = 0; k < tNumOfNodes; k++ )
//            {
//                // create a fem node for each mesh node
//                tNodes( k ) = new fem::Node( allNodes( k ) );
//            }
//
//            // 2) create the IWGs -----------------------------------------------------------
//            std::cout<<" Create the IWGs "<<std::endl;
//            //-------------------------------------------------------------------------------
//            // input a cell of IWG types to be created
////            Cell< IWG_Type > tIWGTypeList = { IWG_Type::HELMHOLTZ,
////                                              IWG_Type::HJ,
////                                              IWG_Type::LSNORMAL,
////                                              IWG_Type::OLSSON };
//            Cell< IWG_Type > tIWGTypeList = { IWG_Type::SPATIALDIFF };
//
//
//            // number of IWGs to be created
//            uint tNumOfIWGs = tIWGTypeList.size();
//
//            // a factory to create the IWGs
//            IWG_Factory tIWGFactory;
//
//            // create a cell of IWGs for the problem considered
//            Cell< IWG* > tIWGs( tNumOfIWGs , nullptr );
//
//            // loop over the IWG types
//            for( uint i = 0; i < tNumOfIWGs; i++)
//            {
//                // create an IWG with the factory for the ith IWG type
//                tIWGs( i ) = tIWGFactory.create_IWGs( tIWGTypeList( i ) );
//            }
//
//////            std::cout<<static_cast< uint >(tIWGs2( 0 )->get_residual_dof_type())<<std::endl;
//////            std::cout<<static_cast< uint >(tIWGs2( 1 )->get_residual_dof_type())<<std::endl;
//////            std::cout<<"---------"<<std::endl;
//////
//////            Cell< MSI::Dof_Type > tIWGActivedofTypes02 = tIWGs2( 0 )->get_active_dof_types();
//////            std::cout<<static_cast< uint >(tIWGActivedofTypes02(0))<<std::endl;
//////            std::cout<<static_cast< uint >(tIWGActivedofTypes02(1))<<std::endl;
//////            Cell< MSI::Dof_Type > tIWGActivedofTypes12 = tIWGs2( 1 )->get_active_dof_types();
//////            std::cout<<static_cast< uint >(tIWGActivedofTypes12(0))<<std::endl;
//////            std::cout<<"---------"<<std::endl;
////
//            // 3) create the elements------------------------------------------------------------
//            std::cout<<" Create the elements "<<std::endl;
//            //-------------------------------------------------------------------------------
//            // a factory to create the elements
//            Element_Factory tElementFactory;
//
//            // get the number of mesh elements
//            uint tNumOfElements = Elems.size();
//
//            // create a cell of fem elements
//            Cell< MSI::Equation_Object * > tListOfElements( tNumOfElements, nullptr );
//
//            // loop over the mesh elements
//            for ( uint i = 0; i < tNumOfElements; i++ )
//            {
//                // create a fem element for the ith mesh element
//                std::cout<<" Create an element "<<std::endl;
//                tListOfElements( i ) = tElementFactory.create_cluster( Element_Type::UNDEFINED,
//                                                                       Elems( i ),
//                                                                       tIWGs,
//                                                                       tNodes );
//
//                // evaluate the residual and jacobian for the ith fem element
//                std::cout<<" Compute residual "<<std::endl;
//                tListOfElements( i )->compute_residual();
//                tListOfElements( i )->compute_jacobian();
////                tListOfElements( i )->compute_jacobian_and_residual();
//            }

////
//            //-------------------------------------------------------------------------------
//
////            // create cell of Dof_Type
////            Cell< Cell < MSI::Dof_Type > > tCellDofType( 2 );
////            tCellDofType( 0 ) = {MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ};
////            tCellDofType( 1 ) = {MSI::Dof_Type::TEMP};
////
////            //std::cout<<static_cast< int >( tCellDofType( 0 )( 0 ) )<<std::endl;
////            //std::cout<<static_cast< int >( tCellDofType( 0 )( 1 ) )<<std::endl;
////            //std::cout<<static_cast< int >( tCellDofType( 0 )( 2 ) )<<std::endl;
////            //std::cout<<static_cast< int >( tCellDofType( 1 )( 0 ) )<<std::endl;
////
////            Cell< MSI::Dof_Type > tCellDofType0 = tCellDofType( 0 );
////            //std::cout<<tCellDofType0.size()<<std::endl;
////
////            Cell< MSI::Dof_Type > tCellDofType1 = tCellDofType( 1 );
////            //std::cout<<tCellDofType1.size()<<std::endl;
//
//            //clean up ----------------------------------------------------------------------
//            //-------------------------------------------------------------------------------
//            for( uint i = 0; i < tNumOfNodes; i++ )
//            {
//                delete allNodes( i );
//                delete tNodes( i );
//            }
//
//            for( uint i = 0; i < tNumOfIWGs; i++)
//            {
//                delete tIWGs( i );
//            }
//
//            for( uint i = 0; i < tNumOfElements; i++ )
//            {
//                delete Elems( i );
//                delete tListOfElements( i );
//            }

        }/* TEST_CASE */
    }/* namespace fem */
}/* namespace moris */

