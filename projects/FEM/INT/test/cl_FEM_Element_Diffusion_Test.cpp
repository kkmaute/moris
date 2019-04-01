
#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

namespace moris
{
    namespace fem
    {

        TEST_CASE( "Element_Diffusion", "[moris],[fem],[ElemDiff]" )
        {
            if(par_size() == 1 )
            {

                // Create a 3D mesh of HEX8 using MTK ------------------------------------------
                std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
                //------------------------------------------------------------------------------
//                uint aNumElemTypes = 1; // only 1 element type ( hex )
//                uint aNumDim = 3;       // number of spatial dimensions
//
//                // element connectivity
//                Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4, 5, 6, 7, 8 }};
//
//                // local to global element map
//                Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};
//
//                // node coordinates
//                Matrix< DDRMat > aCoords = {{ 0.0, 0.0, 0.0 },
//                                            { 1.0, 0.0, 0.0 },
//                                            { 1.0, 1.0, 0.0 },
//                                            { 0.0, 1.0, 0.0 },
//                                            { 0.0, 0.0, 1.0 },
//                                            { 1.0, 0.0, 1.0 },
//                                            { 1.0, 1.0, 1.0 },
//                                            { 0.0, 1.0, 1.0 }};
//
//                // specify the local to global map
//                Matrix< IdMat > aNodeLocalToGlobal = {{ 1,  2,  3,  4,  5,  6,  7,  8 }};
//
//                // create mesh MTK database
//                mtk::MtkMeshData tMeshData( aNumElemTypes );
//                tMeshData.CreateAllEdgesAndFaces  = true;
//                tMeshData.SpatialDim              = & aNumDim;
//                tMeshData.ElemConn( 0 )           = & aElementConnQuad;
//                tMeshData.NodeCoords              = & aCoords;
//                tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
//                tMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;
//
//                mtk::Mesh* tMesh = create_mesh( MeshType::STK, tMeshData );

//                uint aNumElemTypes = 1; // quad
//                uint aNumDim = 2;       // specify number of spatial dimensions
//
//                // specify element connectivity of quad for mesh
//                Matrix< IdMat > aElementConnQuad = { { 1, 2, 3, 4 } };
//
//                // specify the local to global element map for quads
//                Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};
//
//                // Node coordinate matrix
//                Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
//                                            { 1.0, 0.0 },
//                                            { 1.0, 1.0 },
//                                            { -1.0, 1.0 }};
//
//                // specify the local to global map
//                Matrix< IdMat > aNodeLocalToGlobal = { { 1, 2, 3, 4 } };
//
//                //------------------------------------------------------------------------------
//                // create MORIS mesh using MTK database
//                mtk::MtkMeshData aMeshData( aNumElemTypes );
//                aMeshData.CreateAllEdgesAndFaces  = true;
//                aMeshData.SpatialDim              = & aNumDim;
//                aMeshData.ElemConn( 0 )           = & aElementConnQuad;
//                aMeshData.NodeCoords              = & aCoords;
//                aMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
//                aMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;
//
//                mtk::Mesh* tMesh = create_mesh( MeshType::STK, aMeshData );

                uint aNumElemTypes = 1; // only 1 element type ( quad )
                uint aNumDim = 3;       // number of spatial dimensions

                // element connectivity
                Matrix< IdMat > aElementConnQuad = {{ 1, 2, 5, 4, 10, 11, 14, 13 },
                                                    { 2, 3, 6, 5, 11, 12, 15, 14 },
                                                    { 4, 5, 8, 7, 13, 14, 17, 16 },
                                                    { 5, 6, 9, 8, 14, 15, 18, 17 },
                                                    { 10, 11, 14, 13, 19, 20, 23, 22 },
                                                    { 11, 12, 15, 14, 20, 21, 24, 23 },
                                                    { 13, 14, 17, 16, 22, 23, 26, 25 },
                                                    { 14, 15, 18, 17, 23, 24, 27, 26 }};

                // local to global element map for quads
                Matrix< IdMat > aElemLocalToGlobalQuad = { { 1, 2, 3, 4, 5, 6, 7, 8 } };

                // node coordinates
                Matrix< DDRMat > aCoords = {{ 0.0, 0.0, 0.0 },
                                            { 1.0, 0.0, 0.0 },
                                            { 2.0, 0.0, 0.0 },
                                            { 0.0, 1.0, 0.0 },
                                            { 1.0, 1.0, 0.0 },
                                            { 2.0, 1.0, 0.0 },
                                            { 0.0, 2.0, 0.0 },
                                            { 1.0, 2.0, 0.0 },
                                            { 2.0, 2.0, 0.0 },
                                            { 0.0, 0.0, 1.0 },
                                            { 1.0, 0.0, 1.0 },
                                            { 2.0, 0.0, 1.0 },
                                            { 0.0, 1.0, 1.0 },
                                            { 1.0, 1.0, 1.0 },
                                            { 2.0, 1.0, 1.0 },
                                            { 0.0, 2.0, 1.0 },
                                            { 1.0, 2.0, 1.0 },
                                            { 2.0, 2.0, 1.0 },
                                            { 0.0, 0.0, 2.0 },
                                            { 1.0, 0.0, 2.0 },
                                            { 2.0, 0.0, 2.0 },
                                            { 0.0, 1.0, 2.0 },
                                            { 1.0, 1.0, 2.0 },
                                            { 2.0, 1.0, 2.0 },
                                            { 0.0, 2.0, 2.0 },
                                            { 1.0, 2.0, 2.0 },
                                            { 2.0, 2.0, 2.0 } };

                // specify the local to global map
                Matrix< IdMat > aNodeLocalToGlobal = {{ 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                       19, 20, 21, 22, 23, 24, 25, 26, 27 }};

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

//                Cell< fem::IWG_Type >tIWGType( 1, fem::IWG_Type::SPATIALDIFF_BULK );
//
//                // create model
//                mdl::Model * tModel = new mdl::Model( tMesh, 1, tIWGType );
                //------------------------------------------------------------------------------

                // number of mesh nodes
                uint tNumOfNodes = tMesh->get_num_nodes();

                //create a cell of fem nodes
                moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

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
                Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
                                                       fem::IWG_Type::SPATIALDIFF_SIDESET,
                                                       fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                       fem::IWG_Type::SPATIALDIFF_NEUMANN };

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
                // nodal weak bc
                Matrix< DDRMat > tNodalValues( tNumOfNodes, 1, 0.0 );
                real tTempValue = 5.0;
                tNodalValues( 0 ) = tTempValue;
                tNodalValues( 3 ) = tTempValue;
                tNodalValues( 6 ) = tTempValue;
                tNodalValues( 9 ) = tTempValue;
                tNodalValues( 12 ) = tTempValue;
                tNodalValues( 15 ) = tTempValue;
                tNodalValues( 18 ) = tTempValue;
                tNodalValues( 21 ) = tTempValue;
                tNodalValues( 24 ) = tTempValue;

                Matrix< DDRMat > tHeatNodalValues( tNumOfNodes, 1, 0.0 );
                real tHeatValue = 20.0;
                tHeatNodalValues( 2 ) = tHeatValue;
                tHeatNodalValues( 5 ) = tHeatValue;
                tHeatNodalValues( 8 ) = tHeatValue;
                tHeatNodalValues( 11 ) = tHeatValue;
                tHeatNodalValues( 14 ) = tHeatValue;
                tHeatNodalValues( 17 ) = tHeatValue;
                tHeatNodalValues( 20 ) = tHeatValue;
                tHeatNodalValues( 23 ) = tHeatValue;
                tHeatNodalValues( 26 ) = tHeatValue;

                // a factory to create the elements
                Element_Factory tElementFactory;

                // ask mesh about number of elements
//                uint tNumOfElements = tMesh->get_num_elems();

                moris::Cell<std::string> tBlockSetsNames = tMesh->get_set_names( EntityRank::ELEMENT);
                Matrix< IndexMat > tBlockSetElementInd = tMesh->get_set_entity_loc_inds(EntityRank::ELEMENT, tBlockSetsNames(0));
                luint tNumOfElements = tBlockSetElementInd.numel();

                // create equation objects
                Cell< MSI::Equation_Object* > tElements( 8, nullptr );

                // init the number of elements
                uint tNumElementsCount = 0;

                // select the IWG list for the blockset
                Cell< fem::IWG* > tIWGsBlockset = { tIWGs( 0 ) };

                // loop over the mesh elements to create bulk elements
                for( uint k = 0; k < tNumOfElements; k++ )
                {
                    // create a bulk element
                    tElements( k ) = tElementFactory.create_element(   Element_Type::BULK,
                                                                     & tMesh->get_mtk_cell( k ),
                                                                       tIWGsBlockset,
                                                                       tNodes );
                }
                // update the total number of element
                tNumElementsCount = tNumElementsCount + tNumOfElements;

//                // select the IWG list for the sideset------------------------------------------
//                //------------------------------------------------------------------------------
//                Cell< fem::IWG* > tIWGsSideset = { tIWGs( 1 ) };
//
//                //sideset ordinals
//                Cell< Matrix< IndexMat > > tListsOfSideOrdinals = { {{ 0, 3, 4 }},
//                                                                    {{ 0, 1, 4 }},
//                                                                    {{ 2, 3, 4 }},
//                                                                    {{ 1, 2, 4 }},
//                                                                    {{ 0, 3, 5 }},
//                                                                    {{ 0, 1, 5 }},
//                                                                    {{ 2, 3, 5 }},
//                                                                    {{ 1, 2, 5 }} };
//
//
//                // loop over the elements to create sideset elements
//                for( uint k = 0; k < tNumOfElements; k++ )
//                {
//
//                    // create a sideset element
//                    tElements( tNumElementsCount + k ) = tElementFactory.create_element(   Element_Type::SIDESET,
//                                                                                         & tMesh->get_mtk_cell( k ),
//                                                                                           tIWGsSideset,
//                                                                                           tNodes );
//                    // create and set the list of side ordinals
//                    Matrix< IndexMat > tListOfSideOrdinals = tListsOfSideOrdinals( k );
//                    tElements( tNumOfElements + k )->set_list_of_side_ordinals( tListOfSideOrdinals );
//                }
//                // update the total number of element
//                tNumElementsCount = tNumElementsCount + tNumOfElements;

//                // select the IWG list for the Dirichlet sideset
//                Cell< fem::IWG* > tIWGsDirichletSideset = { tIWGs( 2 ) };
//
//                // elements included in the Dirichlet sideset
//                Cell< moris_index > tListOfDirichletElements = { 0, 2, 4, 6 };
//
//                // side ordinal for the Dirichlet sideset
//                Matrix< IndexMat > tListOfSideOrdinals2 = { { 3 } };
//
//                for( uint iDirichlet = 0; iDirichlet < tListOfDirichletElements.size(); iDirichlet++ )
//                {
//                    moris_index tTreatedMeshElement = tListOfDirichletElements( iDirichlet );
//
//                    tElements( tNumElementsCount + iDirichlet )
//                        = tElementFactory.create_element(   Element_Type::SIDESET,
//                                                          & tMesh->get_mtk_cell( tTreatedMeshElement ),
//                                                            tIWGsDirichletSideset,
//                                                            tNodes );
//
//                    tElements( tNumElementsCount + iDirichlet )->set_list_of_side_ordinals( tListOfSideOrdinals2 );
//
//                    // get the nodal weak bcs of the element
//                    Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumElementsCount + iDirichlet )->get_weak_bcs();
//
//                    // get the element number of nodes
//                    uint tNumberOfNodes = tElements( tNumElementsCount + iDirichlet )->get_num_nodes();
//
//                    // set size of the element nodal weak bc
//                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );
//
//                    // loop over the element nodes
//                    Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tTreatedMeshElement ).get_vertex_inds();
//
//                    for( uint l = 0; l < tNumberOfNodes; l++ )
//                    {
//                        // copy weak bc into element
//                        tNodalWeakBCs( l ) = tNodalValues( tNodeIndices( l ) );
//                    }
//                }
//                // update the total number of element
//                tNumElementsCount = tNumElementsCount + tListOfDirichletElements.size();
//
//                // select the IWG list for the Neumann sideset
//                Cell< fem::IWG* > tIWGsNeumannSideset = { tIWGs( 3 ) };
//
//                // elements included in the Neumann sideset
//                Cell< moris_index > tListOfNeumannElements = { 1, 3, 5, 7 };
//
//                // side ordinal for the Neumann sideset
//                Matrix< IndexMat > tListOfSideOrdinals3 = { { 1 } };
//
//                for( uint iNeumann = 0; iNeumann < 4; iNeumann++ )
//                {
//                    moris_index tTreatedMeshElement = tListOfNeumannElements( iNeumann );
//
//                    tElements( tNumElementsCount + iNeumann )
//                        = tElementFactory.create_element(   Element_Type::SIDESET,
//                                                          & tMesh->get_mtk_cell( tTreatedMeshElement ),
//                                                            tIWGsNeumannSideset,
//                                                            tNodes );
//
//                    tElements( tNumElementsCount + iNeumann )->set_list_of_side_ordinals( tListOfSideOrdinals3 );
//
//                    // get the nodal weak bcs of the element
//                    Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumElementsCount + iNeumann )->get_weak_bcs();
//
//                    // get the element number of nodes
//                    uint tNumberOfNodes = tElements( tNumElementsCount + iNeumann )->get_num_nodes();
//
//                    // set size of the element nodal weak bc
//                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );
//
//                    // loop over the element nodes
//                    Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tTreatedMeshElement ).get_vertex_inds();
//
//                    for( uint l = 0; l < tNumberOfNodes; l++ )
//                    {
//                        // copy weak bc into element
//                        tNodalWeakBCs( l ) = tHeatNodalValues( tNodeIndices( l ) );
//                    }
//                }


                //4) Create the model solver interface -----------------------------------------
                std::cout<<" Create the model solver interface "<<std::endl;
                //------------------------------------------------------------------------------

                //FIXME force the communication table
                Matrix< IdMat > tCommunicationTable( 1, 1, 0 );

                // FIXME: get map from mesh
                uint tDofOrder = 1;
                map< moris_id, moris_index > tCoefficientsMap;
                //tMesh->get_adof_map( tDofOrder, tCoefficientsMap );

                uint tNumCoeff = 100;
                //= tMesh->get_num_coeffs( 1 )

                moris::MSI::Model_Solver_Interface* tModelSolverInterface
                    = new moris::MSI::Model_Solver_Interface( tElements,
                                                              tCommunicationTable,
                                                              tCoefficientsMap,
                                                              tNumCoeff,
                                                              tMesh );

                tModelSolverInterface->set_param( "TEMP" )  = (sint)tDofOrder;

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
                    = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

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

                // 8) Postprocessing
                // dof type list for the solution to write on the mesh
                moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

                // create a matrix to be filled  with the solution
                Matrix< DDRMat > tTempSolutionField( tNumOfNodes, 1 );

                // loop over the nodes
                for( uint i = 0; i < tNumOfNodes; i++ )
                {
                    // get a list of elements connected to the ith node
                    Matrix<IndexMat> tConnectedElements =
                        tMesh->get_entity_connected_to_entity_loc_inds( static_cast< moris_index >( i ),
                                                                        EntityRank::NODE,
                                                                        EntityRank::ELEMENT );

                    // number of connected element
                    uint tNumConnectElem = tConnectedElements.numel();

                    // reset the nodal value
                    real tNodeVal = 0.0;

                    // loop over the connected elements
                    for( uint j = 0; j < tNumConnectElem; j++ )
                    {
                        // extract the field value at the ith node for the jth connected element
                        real tElemVal
                            = tElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i,
                                                                                                  tDofTypeList);
                        // add up the contribution of each element to the node value
                        tNodeVal = tNodeVal + tElemVal;
                    }
                    // fill the solution matrix with the node value
                    tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
                }
                print( tTempSolutionField, "tTempSolutionField" );

//                CHECK( equal_to( tTempSolutionField( 0, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 1, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 2, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 3, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 4, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 5, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 6, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 7, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 8, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 9, 0 ), 5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 10, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 11, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 12, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 13, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 14, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 15, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 16, 0 ),5.00, 1.0e+08 ) );
//                CHECK( equal_to( tTempSolutionField( 17, 0 ),5.00, 1.0e+08 ) );

                // initialize Scalar_Field_Info structure with a DDRMat as the template type
                moris::mtk::Scalar_Field_Info<DDRMat> tTempField;

                // set the fields names
                std::string tTempFieldName = "Temp_Field";
                tTempField.set_field_name( tTempFieldName );

                // set the entity rank associated with the fields
                tTempField.set_field_entity_rank( EntityRank::NODE );

                // initialize field information container
                moris::mtk::MtkFieldsInfo tFieldsInfo;

                // place the node field into the field info container
                add_field_for_mesh_input( &tTempField, tFieldsInfo );

                // declare some supplementary fields
                tMeshData.FieldsInfo = &tFieldsInfo;

                // create the mesh
                mtk::Mesh* tMeshForOutput = create_mesh( MeshType::STK, tMeshData );

                // add field to the mesh
                tMeshForOutput->add_mesh_field_real_scalar_data_loc_inds( tTempFieldName,
                                                                          EntityRank::NODE,
                                                                          tTempSolutionField );

                // create output mesh
                std::string tOutputFile = "./int_ElemDiff_test.exo";
                tMeshForOutput->create_output_mesh( tOutputFile );

                // 8) Clean up -----------------------------------------------------------------
                std::cout<<" Clean up "<<std::endl;
                //------------------------------------------------------------------------------
                delete tMesh;
                delete tMeshForOutput;

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

            }/* if( par_size() */
        }/* TEST_CASE */


        TEST_CASE( "Element_Diffusion_1", "[moris],[fem],[ElemDiff_1]" )
        {
        if(par_size() == 1 )
        {
            // Create a 3D mesh of HEX8 using MTK ------------------------------------------
            std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
            //------------------------------------------------------------------------------

            std::string tPrefix = std::getenv("MORISROOT");
            std::string tMeshFileName = tPrefix + "projects/FEM/INT/test/data/Cube_with_side_sets.g";
//            std::string tMeshFileName = "generated:2x2x2|sideset:xXyYzZ";

            std::cout<<"Mesh input name = "<< tMeshFileName<<std::endl;

            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
            std::string tFieldName1 = "Temp_Field";
            tNodeField1.set_field_name(tFieldName1);
            tNodeField1.set_field_entity_rank(EntityRank::NODE);

            // Initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;

            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

            // Declare some supplementary fields
            mtk::MtkMeshData tMeshData;
            tMeshData.FieldsInfo = &tFieldsInfo;

            moris::mtk::Mesh* tMesh = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, &tMeshData );

            //1) Create the fem nodes ------------------------------------------------------
            std::cout<<" Create the fem nodes "<<std::endl;
            //------------------------------------------------------------------------------
//                Cell< fem::IWG_Type >tIWGType( 1, fem::IWG_Type::SPATIALDIFF_BULK );
//
//                // create model
//                mdl::Model * tModel = new mdl::Model( tMesh, 1, tIWGType );
            //------------------------------------------------------------------------------

            // number of mesh nodes
            uint tNumOfNodes = tMesh->get_num_nodes();

            //create a cell of fem nodes
            moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

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
            Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
                                                   fem::IWG_Type::SPATIALDIFF_SIDESET,
                                                   fem::IWG_Type::SPATIALDIFF_DIRICHLET };

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
            // nodal weak bc
            Matrix< DDRMat > tNodalValues( tNumOfNodes, 1, 0.0 );
            tNodalValues( 23 ) = 5.0;
            tNodalValues( 12 ) = 5.0;
            tNodalValues( 25 ) = 5.0;
            tNodalValues( 15 ) = 5.0;
            tNodalValues( 4 ) = 5.0;
            tNodalValues( 17 ) = 5.0;
            tNodalValues( 19 ) = 5.0;
            tNodalValues( 8 ) = 5.0;
            tNodalValues( 21 ) = 5.0;

            // a factory to create the elements
            Element_Factory tElementFactory;

            // ask mesh about number of elements
            // uint tNumOfElements = tMesh->get_num_elems();

            moris::Cell<std::string> tBlockSetsNames = tMesh->get_set_names( EntityRank::ELEMENT);
            Matrix< IndexMat > tBlockSetElementInd = tMesh->get_set_entity_loc_inds(EntityRank::ELEMENT, tBlockSetsNames(0));
            luint tNumOfElements = tBlockSetElementInd.numel();

            // create equation objects
            Cell< MSI::Equation_Object* > tElements( 2 * tNumOfElements + 4, nullptr );

            Cell< fem::IWG* > tIWGs1 = { tIWGs( 0 ) };
            Cell< fem::IWG* > tIWGs2 = { tIWGs( 1 ) };
            Cell< fem::IWG* > tIWGs3 = { tIWGs( 2 ) };

            Cell< Matrix< IndexMat > > tListsOfSideOrdinals = { {{ 0, 3, 4 }},
                                                                {{ 0, 1, 4 }},
                                                                {{ 2, 3, 4 }},
                                                                {{ 1, 2, 4 }},
                                                                {{ 0, 3, 5 }},
                                                                {{ 0, 1, 5 }},
                                                                {{ 2, 3, 5 }},
                                                                {{ 1, 2, 5 }} };

            // loop over the mesh elements
            for( uint k = 0; k < tNumOfElements; k++ )
            {
                // create a bulk element-----------------------------------------
                tElements( k )
                    = tElementFactory.create_element(   Element_Type::BULK,
                                                      & tMesh->get_mtk_cell( k ),
                                                        tIWGs1,
                                                        tNodes );

                // create a sideset element---------------------------------------
                tElements( tNumOfElements + k )
                    = tElementFactory.create_element(   Element_Type::SIDESET,
                                                      & tMesh->get_mtk_cell( k ),
                                                        tIWGs2,
                                                        tNodes );
                // create and set the list of side ordinals
                Matrix< IndexMat > tListOfSideOrdinals = tListsOfSideOrdinals( k );

                tElements( tNumOfElements + k )->set_list_of_side_ordinals( tListOfSideOrdinals );
            }

            Cell< moris_index > tListOfBCElements = { 0, 1, 2, 3 };
            for( uint iDirichlet = 0; iDirichlet < 4; iDirichlet++ )
            {
                moris_index tTreatedMeshElement = tListOfBCElements( iDirichlet );

                tElements( 2 * tNumOfElements + iDirichlet )
                    = tElementFactory.create_element(   Element_Type::SIDESET,
                                                      & tMesh->get_mtk_cell( tTreatedMeshElement ),
                                                        tIWGs3,
                                                        tNodes );
                Matrix< IndexMat > tListOfSideOrdinals2 = { { 4 } };
                tElements( 2 * tNumOfElements + iDirichlet )->set_list_of_side_ordinals( tListOfSideOrdinals2 );

                // get the nodal weak bcs of the element
                Matrix< DDRMat > & tNodalWeakBCs = tElements( 2 * tNumOfElements + iDirichlet )->get_weak_bcs();

                // get the element number of nodes
                uint tNumberOfNodes = tElements( 2 * tNumOfElements + iDirichlet )->get_num_nodes();

                // set size of the element nodal weak bc
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                // loop over the element nodes
                Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tTreatedMeshElement ).get_vertex_inds();

                for( uint l = 0; l < tNumberOfNodes; l++ )
                {
                    // copy weak bc into element
                    tNodalWeakBCs( l ) = tNodalValues( tNodeIndices( l ) );
                }
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

            uint tNumCoeff = 100;
            //= tMesh->get_num_coeffs( 1 )

            moris::MSI::Model_Solver_Interface* tModelSolverInterface
                = new moris::MSI::Model_Solver_Interface( tElements,
                                                          tCommunicationTable,
                                                          tCoefficientsMap,
                                                          tNumCoeff,
                                                          tMesh );

            tModelSolverInterface->set_param( "TEMP" )  = (sint)tDofOrder;

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
                = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

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
//            print( tSolution1, "tSolution1" );

            // 8) Postprocessing
            // dof type list for the solution to write on the mesh
            moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

            // create a matrix to be filled  with the solution
            Matrix< DDRMat > tTempSolutionField( tNumOfNodes, 1 );

            // loop over the nodes
            for( uint i = 0; i < tNumOfNodes; i++ )
            {
                // get a list of elements connected to the ith node
                Matrix<IndexMat> tConnectedElements =
                    tMesh->get_entity_connected_to_entity_loc_inds( static_cast< moris_index >( i ),
                                                                    EntityRank::NODE,
                                                                    EntityRank::ELEMENT );

                // number of connected element
                uint tNumConnectElem = tConnectedElements.numel();

                // reset the nodal value
                real tNodeVal = 0.0;

                // loop over the connected elements
                for( uint j = 0; j < tNumConnectElem; j++ )
                {
                    // extract the field value at the ith node for the jth connected element
                    real tElemVal
                        = tElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i,
                                                                                              tDofTypeList);
                    // add up the contribution of each element to the node value
                    tNodeVal = tNodeVal + tElemVal;
                }
                // fill the solution matrix with the node value
                tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
            }
            print( tTempSolutionField, "tTempSolutionField" );

            CHECK( equal_to( tTempSolutionField( 0, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 1, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 2, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 3, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 4, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 5, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 6, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 7, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 8, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 9, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 10, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 11, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 12, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 13, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 14, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 15, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 16, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 17, 0 ),5.00, 1.0e+08 ) );


            // add field to the mesh
            tMesh->add_mesh_field_real_scalar_data_loc_inds( tFieldName1,
                                                                      EntityRank::NODE,
                                                                      tTempSolutionField );

            // create output mesh
            std::string tOutputFile = "./int_ElemDiff_test_1.exo";
            tMesh->create_output_mesh( tOutputFile );

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
            delete tSolverInterface;
            delete tNonlinearProblem;
            delete tLinSolver;
            delete tNonlinearSolver;
        }/* if( par_size() */
    }

        TEST_CASE( "Element_Diffusion_2", "[moris],[fem],[ElemDiff_2]" )
        {
        if(par_size() == 1 )
        {
            // Create a 3D mesh of HEX8 using MTK ------------------------------------------
            std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
            //------------------------------------------------------------------------------

            std::string tPrefix = std::getenv("MORISROOT");
            std::string tMeshFileName = tPrefix + "projects/FEM/INT/test/data/Cube_with_side_sets.g";
//            std::string tMeshFileName = "generated:2x2x2|sideset:xXyYzZ";

            std::cout<<"Mesh input name = "<< tMeshFileName<<std::endl;

            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
            std::string tFieldName1 = "Temp_Field";
            tNodeField1.set_field_name(tFieldName1);
            tNodeField1.set_field_entity_rank(EntityRank::NODE);

            // Initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;

            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

            // Declare some supplementary fields
            mtk::MtkMeshData tMeshData;
            tMeshData.FieldsInfo = &tFieldsInfo;

            moris::mtk::Mesh* tMesh = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, &tMeshData );

            //1) Create the fem nodes ------------------------------------------------------
            std::cout<<" Create the fem nodes "<<std::endl;
            //------------------------------------------------------------------------------
//                Cell< Cell< fem::IWG_Type > >tIWGType( 3 );
//                tIWGType( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK);
//                tIWGType( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_SIDESET);
//                tIWGType( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET);
//
//                // create model
//                mdl::Model * tModel = new mdl::Model( tMesh, 1, tIWGType );
//
//                // set weakBC
//
//                //solve
            //------------------------------------------------------------------------------

            // number of mesh nodes
            uint tNumOfNodes = tMesh->get_num_nodes();

            //create a cell of fem nodes
            moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

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
            Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
                                                   fem::IWG_Type::SPATIALDIFF_SIDESET,
                                                   fem::IWG_Type::SPATIALDIFF_DIRICHLET };

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
            // nodal weak bc
            Matrix< DDRMat > tNodalValues( tNumOfNodes, 1, 0.0 );
            tNodalValues( 23 ) = 5.0;
            tNodalValues( 12 ) = 5.0;
            tNodalValues( 25 ) = 5.0;
            tNodalValues( 15 ) = 5.0;
            tNodalValues( 4 )  = 5.0;
            tNodalValues( 17 ) = 5.0;
            tNodalValues( 19 ) = 5.0;
            tNodalValues( 8 )  = 5.0;
            tNodalValues( 21 ) = 5.0;

            // a factory to create the elements
            Element_Factory tElementFactory;

            // ask mesh about number of elements
            // uint tNumOfElements = tMesh->get_num_elems();

            moris::Cell<std::string> tBlockSetsNames = tMesh->get_set_names( EntityRank::ELEMENT);
            Matrix< IndexMat > tBlockSetElementInd = tMesh->get_set_entity_loc_inds(EntityRank::ELEMENT, tBlockSetsNames(0));
            luint tNumOfElements = tBlockSetElementInd.numel();

            // create equation objects
            Cell< MSI::Equation_Object* > tElements( 2 * tNumOfElements + 4, nullptr );

            Cell< fem::IWG* > tIWGs1 = { tIWGs( 0 ) };
            Cell< fem::IWG* > tIWGs2 = { tIWGs( 1 ) };
            Cell< fem::IWG* > tIWGs3 = { tIWGs( 2 ) };

            Cell< Matrix< IndexMat > > tListsOfSideOrdinals = { {{ 0, 3, 4 }},
                                                                {{ 0, 1, 4 }},
                                                                {{ 2, 3, 4 }},
                                                                {{ 1, 2, 4 }},
                                                                {{ 0, 3, 5 }},
                                                                {{ 0, 1, 5 }},
                                                                {{ 2, 3, 5 }},
                                                                {{ 1, 2, 5 }} };

            // loop over the mesh elements
            for( uint k = 0; k < tNumOfElements; k++ )
            {
                // create a bulk element-----------------------------------------
                tElements( k )
                    = tElementFactory.create_element(   Element_Type::BULK,
                                                      & tMesh->get_mtk_cell( k ),
                                                        tIWGs1,
                                                        tNodes );

                // create a sideset element---------------------------------------
                tElements( tNumOfElements + k )
                    = tElementFactory.create_element(   Element_Type::SIDESET,
                                                      & tMesh->get_mtk_cell( k ),
                                                        tIWGs2,
                                                        tNodes );
                // create and set the list of side ordinals
                Matrix< IndexMat > tListOfSideOrdinals = tListsOfSideOrdinals( k );

                tElements( tNumOfElements + k )->set_list_of_side_ordinals( tListOfSideOrdinals );
            }

            Cell< moris_index > tListOfBCElements = { 0, 1, 2, 3 };
            for( uint iDirichlet = 0; iDirichlet < 4; iDirichlet++ )
            {
                moris_index tTreatedMeshElement = tListOfBCElements( iDirichlet );

                tElements( 2 * tNumOfElements + iDirichlet )
                    = tElementFactory.create_element(   Element_Type::SIDESET,
                                                      & tMesh->get_mtk_cell( tTreatedMeshElement ),
                                                        tIWGs3,
                                                        tNodes );
                Matrix< IndexMat > tListOfSideOrdinals2 = { { 4 } };
                tElements( 2 * tNumOfElements + iDirichlet )->set_list_of_side_ordinals( tListOfSideOrdinals2 );

                // get the nodal weak bcs of the element
                Matrix< DDRMat > & tNodalWeakBCs = tElements( 2 * tNumOfElements + iDirichlet )->get_weak_bcs();

                // get the element number of nodes
                uint tNumberOfNodes = tElements( 2 * tNumOfElements + iDirichlet )->get_num_nodes();

                // set size of the element nodal weak bc
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                // loop over the element nodes
                Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tTreatedMeshElement ).get_vertex_inds();

                for( uint l = 0; l < tNumberOfNodes; l++ )
                {
                    // copy weak bc into element
                    tNodalWeakBCs( l ) = tNodalValues( tNodeIndices( l ) );
                }
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

            uint tNumCoeff = 100;
            //= tMesh->get_num_coeffs( 1 )

            moris::MSI::Model_Solver_Interface* tModelSolverInterface
                = new moris::MSI::Model_Solver_Interface( tElements,
                                                          tCommunicationTable,
                                                          tCoefficientsMap,
                                                          tNumCoeff,
                                                          tMesh );

            tModelSolverInterface->set_param( "TEMP" )  = (sint)tDofOrder;

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
                = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

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
//            print( tSolution1, "tSolution1" );

            // 8) Postprocessing
            // dof type list for the solution to write on the mesh
            moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

            // create a matrix to be filled  with the solution
            Matrix< DDRMat > tTempSolutionField( tNumOfNodes, 1 );

            // loop over the nodes
            for( uint i = 0; i < tNumOfNodes; i++ )
            {
                // get a list of elements connected to the ith node
                Matrix<IndexMat> tConnectedElements =
                    tMesh->get_entity_connected_to_entity_loc_inds( static_cast< moris_index >( i ),
                                                                    EntityRank::NODE,
                                                                    EntityRank::ELEMENT );

                // number of connected element
                uint tNumConnectElem = tConnectedElements.numel();

                // reset the nodal value
                real tNodeVal = 0.0;

                // loop over the connected elements
                for( uint j = 0; j < tNumConnectElem; j++ )
                {
                    // extract the field value at the ith node for the jth connected element
                    real tElemVal
                        = tElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i,
                                                                                              tDofTypeList);
                    // add up the contribution of each element to the node value
                    tNodeVal = tNodeVal + tElemVal;
                }
                // fill the solution matrix with the node value
                tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
            }
            print( tTempSolutionField, "tTempSolutionField" );

            CHECK( equal_to( tTempSolutionField( 0, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 1, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 2, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 3, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 4, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 5, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 6, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 7, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 8, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 9, 0 ), 5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 10, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 11, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 12, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 13, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 14, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 15, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 16, 0 ),5.00, 1.0e+08 ) );
            CHECK( equal_to( tTempSolutionField( 17, 0 ),5.00, 1.0e+08 ) );


            // add field to the mesh
            tMesh->add_mesh_field_real_scalar_data_loc_inds( tFieldName1,
                                                                      EntityRank::NODE,
                                                                      tTempSolutionField );

            // create output mesh
            std::string tOutputFile = "./int_ElemDiff_test.exo";
            tMesh->create_output_mesh( tOutputFile );

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
            delete tSolverInterface;
            delete tNonlinearProblem;
            delete tLinSolver;
            delete tNonlinearSolver;
        }/* if( par_size() */
    }

    }/* namespace fem */
}/* namespace moris */
