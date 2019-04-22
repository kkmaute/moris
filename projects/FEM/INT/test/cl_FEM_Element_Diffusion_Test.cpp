#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src
#include "fn_norm.hpp"
#include "fn_unique.hpp"

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Element_Block.hpp"            //FEM/INT/src

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

// four tests in this .cpp file:
// 1) Line   51: "ElemDiff" - 2x2x2 Block with Neumann + Dirichlet BCs
// 2) Line  720: "ElemDiff_Bar" - 5x1x1 Block with Ghost applied
// 3) Line 1270: "ElemDiff_Bar_quad" - 2x1x1 Block of quadratic elements with Ghost applied
// 4) Line 1836: "ElemDiff_Bar_cut" - 1x1x4 Block with Ghost applied and a cut last element

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//1//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

        TEST_CASE( "Element_Diffusion", "[moris],[fem],[ElemDiff]" )
        {
            if(par_size() == 1 )
            {
                // Create a 3D mesh of HEX8 using MTK ------------------------------------------
                std::cout<<" Create a 3D mesh of HEX8 using MTK "<<std::endl;
                //------------------------------------------------------------------------------
                uint aNumElemTypes = 1; // only 1 element type ( hex )
                uint aNumDim = 3;       // number of spatial dimensions

                // element connectivity
                Matrix< IdMat > aElementConn = {{ 1, 2, 5, 4, 10, 11, 14, 13 },
                                                    { 2, 3, 6, 5, 11, 12, 15, 14 },
                                                    { 4, 5, 8, 7, 13, 14, 17, 16 },
                                                    { 5, 6, 9, 8, 14, 15, 18, 17 },
                                                    { 10, 11, 14, 13, 19, 20, 23, 22 },
                                                    { 11, 12, 15, 14, 20, 21, 24, 23 },
                                                    { 13, 14, 17, 16, 22, 23, 26, 25 },
                                                    { 14, 15, 18, 17, 23, 24, 27, 26 }};

                // local to global element map for hexs
                Matrix< IdMat > aElemLocalToGlobal = { { 1, 2, 3, 4, 5, 6, 7, 8 } };

                // node coordinates
                Matrix< DDRMat > aCoords = {{ 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 },
                                            { 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 2.0, 1.0, 0.0 },
                                            { 0.0, 2.0, 0.0 }, { 1.0, 2.0, 0.0 }, { 2.0, 2.0, 0.0 },
                                            { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 2.0, 0.0, 1.0 },
                                            { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 2.0, 1.0, 1.0 },
                                            { 0.0, 2.0, 1.0 }, { 1.0, 2.0, 1.0 }, { 2.0, 2.0, 1.0 },
                                            { 0.0, 0.0, 2.0 }, { 1.0, 0.0, 2.0 }, { 2.0, 0.0, 2.0 },
                                            { 0.0, 1.0, 2.0 }, { 1.0, 1.0, 2.0 }, { 2.0, 1.0, 2.0 },
                                            { 0.0, 2.0, 2.0 }, { 1.0, 2.0, 2.0 }, { 2.0, 2.0, 2.0 }};

                // specify the local to global map
                Matrix< IdMat > aNodeLocalToGlobal = {{ 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                       19, 20, 21, 22, 23, 24, 25, 26, 27 }};

                // create mesh MTK database
                mtk::MtkMeshData tMeshData( aNumElemTypes );
                tMeshData.CreateAllEdgesAndFaces  = true;
                tMeshData.SpatialDim              = & aNumDim;
                tMeshData.ElemConn( 0 )           = & aElementConn;
                tMeshData.NodeCoords              = & aCoords;
                tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobal;
                tMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;

                mtk::Mesh* tMesh = create_mesh( MeshType::STK, tMeshData );

                //1) Create the fem nodes ------------------------------------------------------
                //std::cout<<" Create the fem nodes "<<std::endl;
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
                //std::cout<<" Create the IWGs "<<std::endl;
                //------------------------------------------------------------------------------

                // input a cell of IWG types to be created
                Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
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

                //3) Create element blocks -------------------------------------------------------
                //std::cout<<" Create element blocks "<<std::endl;
                //------------------------------------------------------------------------------
                // get the number of elements on proc from the mesh
                luint tNumOfElements = tMesh->get_num_elems();

                // create equation objects
                Cell< MSI::Equation_Object* > tElements;
                tElements.reserve( tNumOfElements + 4 + 4 );

                // get the block names from the mesh
                moris::Cell<std::string> tBlockSetsNames = tMesh->get_set_names( EntityRank::ELEMENT);

                // Cell containing the block mesh cell ( a cell of mesh cells )
                moris::Cell<mtk::Cell const *> tBlockSetElement( tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( 0 ) ).numel(), nullptr );

                // loop on the blocks
                for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
                {
                    Matrix< IndexMat > tBlockSetElementInd = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( Ik ) );

                    // loop on the elements in a block
                    for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                    {
                        // Cell containing the block mesh cell with mesh cells
                        tBlockSetElement( k ) = & tMesh->get_mtk_cell( k );
                    }
                }

                // create a fem element block
                moris::Cell< fem::Element_Block * > tElementBlocks( 3, nullptr );
                moris::Cell< IWG* > tIWGBulk = { tIWGs( 0 ) };
                tElementBlocks( 0 ) = new fem::Element_Block( tBlockSetElement,
                                                              fem::Element_Type::BULK,
                                                              tIWGBulk,
                                                              tNodes );

                // put the equation object of block 0 in the global list of equation objects
                tElements.append( tElementBlocks( 0 )->get_equation_object_list() );

                //4) Create side blocks -------------------------------------------------------
                //std::cout<<" Create Dirichlet side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Dirichlet sideset
                moris::Cell<mtk::Cell const *> tSideSetDirichletElement( 4, nullptr );
                tSideSetDirichletElement( 0 ) = & tMesh->get_mtk_cell( 0 );
                tSideSetDirichletElement( 1 ) = & tMesh->get_mtk_cell( 2 );
                tSideSetDirichletElement( 2 ) = & tMesh->get_mtk_cell( 4 );
                tSideSetDirichletElement( 3 ) = & tMesh->get_mtk_cell( 6 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGDirichlet = { tIWGs( 1 ) };
                tElementBlocks( 1 ) = new fem::Element_Block( tSideSetDirichletElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGDirichlet,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 1 )->get_equation_object_list() );

                //5) Create side blocks -------------------------------------------------------
                //std::cout<<" Create Neumann side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Neumann sideset
                moris::Cell<mtk::Cell const *> tSideSetNeumannElement( 4, nullptr );
                tSideSetNeumannElement( 0 ) = & tMesh->get_mtk_cell( 1 );
                tSideSetNeumannElement( 1 ) = & tMesh->get_mtk_cell( 3 );
                tSideSetNeumannElement( 2 ) = & tMesh->get_mtk_cell( 5 );
                tSideSetNeumannElement( 3 ) = & tMesh->get_mtk_cell( 7 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGNeumann = { tIWGs( 2 ) };
                tElementBlocks( 2 ) = new fem::Element_Block( tSideSetNeumannElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGNeumann,
                                                              tNodes );

                // put the equation object of block 2 in the global list of equation objects
                tElements.append( tElementBlocks( 2 )->get_equation_object_list() );

                //6) Impose boundary conditions -----------------------------------------------
                //std::cout<<" Impose boundary conditions "<<std::endl;
                //-----------------------------------------------------------------------------
                // nodal weak bc for Dirichlet
                real tTempValue = 5.0;
                Cell< moris_index > tDirichletElements = { 0, 2, 4, 6 };
                for( uint iDirichlet = 0; iDirichlet < 4; iDirichlet++ )
                {
                    // set the list of face ordinals
                    tElements( tNumOfElements + iDirichlet )->set_list_of_side_ordinals( {{ 3 }} );

                    //get the nodal weak bcs of the element
                    Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumOfElements + iDirichlet )->get_weak_bcs();

                    // get the element number of nodes
                    uint tNumberOfNodes = tElements( tNumOfElements + iDirichlet )->get_num_nodes();

                    // set size of the element nodal weak bc
                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                    // loop over the element nodes
                    Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tDirichletElements( iDirichlet ) ).get_vertex_inds();

                    for( uint l = 0; l < tNumberOfNodes; l++ )
                    {
                        // copy weak bc into element
                        tNodalWeakBCs( l ) = tTempValue;
                    }
                }

                // nodal weak bc for Neumann
                real tFluxValue = 20.0;
                Cell< moris_index > tNeumannElements = { 0, 2, 4, 6 };
                for( uint iNeumann = 0; iNeumann < 4; iNeumann++ )
                {
                    // set the list of face ordinals
                    tElements( tNumOfElements + tDirichletElements.size() + iNeumann )->set_list_of_side_ordinals( {{ 1 }} );

                    //get the nodal weak bcs of the element
                    Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumOfElements + tDirichletElements.size() + iNeumann )->get_weak_bcs();

                    // get the element number of nodes
                    uint tNumberOfNodes = tElements( tNumOfElements + tDirichletElements.size() + iNeumann )->get_num_nodes();

                    // set size of the element nodal weak bc
                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                    // loop over the element nodes
                    Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tDirichletElements( iNeumann ) ).get_vertex_inds();

                    for( uint l = 0; l < tNumberOfNodes; l++ )
                    {
                        // copy weak bc into element
                        tNodalWeakBCs( l ) = tFluxValue;
                    }
                }

//                //3) Create elements -----------------------------------------------
//                std::cout<<" Create elements "<<std::endl;
//                //-----------------------------------------------------------------------------
//                // a factory to create the elements
//                Element_Factory tElementFactory;
//
//                // create equation objects
//                Cell< MSI::Equation_Object* > tElements( 16, nullptr );
//
//                // init the number of elements
//                uint tNumElementsCount = 0;
//
//                // select the IWG list for the blockset
//                Cell< fem::IWG* > tIWGsBlockset = { tIWGs( 0 ) };
//
//                // loop over the mesh elements to create bulk elements
//                for( uint k = 0; k < tNumOfElements; k++ )
//                {
//                    // create a bulk element
//                    tElements( tNumElementsCount + k ) = tElementFactory.create_element(   Element_Type::BULK,
//                                                                                         & tMesh->get_mtk_cell( k ),
//                                                                                           tIWGsBlockset,
//                                                                                           tNodes );
//                }
//                // update the total number of element
//                tNumElementsCount = tNumElementsCount + tNumOfElements;
//
//                // select the IWG list for the Dirichlet sideset
//                Cell< fem::IWG* > tIWGsDirichletSideset = { tIWGs( 1 ) };
//
//                // elements included in the Dirichlet sideset
//                Cell< moris_index > tListOfDirichletElements = { 0, 2, 4, 6 };
//
//                // loop over the elements included in the Dirichlet sideset
//                for( uint iDirichlet = 0; iDirichlet < tListOfDirichletElements.size(); iDirichlet++ )
//                {
//                    // get the treated element index
//                    moris_index tTreatedMeshElement = tListOfDirichletElements( iDirichlet );
//
//                    // create a sideset element
//                    tElements( tNumElementsCount + iDirichlet )
//                        = tElementFactory.create_element(   Element_Type::SIDESET,
//                                                          & tMesh->get_mtk_cell( tTreatedMeshElement ),
//                                                            tIWGsDirichletSideset,
//                                                            tNodes );
//                    // set the list of face ordinals
//                    tElements( tNumElementsCount + iDirichlet )->set_list_of_side_ordinals( {{ 3 }} );
//
//                    // impose nodal weak bcs
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
//                Cell< fem::IWG* > tIWGsNeumannSideset = { tIWGs( 2 ) };
//
//                // elements included in the Neumann sideset
//                Cell< moris_index > tListOfNeumannElements = { 1, 3, 5, 7 };
//
//                // loop over the elements included in the Neumann sideset
//                for( uint iNeumann = 0; iNeumann < 4; iNeumann++ )
//                {
//                    // get the treated element index
//                    moris_index tTreatedMeshElement = tListOfNeumannElements( iNeumann );
//
//                    // create a sideset element
//                    tElements( tNumElementsCount + iNeumann )
//                        = tElementFactory.create_element(   Element_Type::SIDESET,
//                                                          & tMesh->get_mtk_cell( tTreatedMeshElement ),
//                                                            tIWGsNeumannSideset,
//                                                            tNodes );
//
//                    // set the list of face ordinals
//                    tElements( tNumElementsCount + iNeumann )->set_list_of_side_ordinals( {{ 1 }} );
//
//                    // impose nodal weak bcs
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

                //7) Create the model solver interface -----------------------------------------
                //std::cout<<" Create the model solver interface "<<std::endl;
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

                tElementBlocks( 0 )->finalize( tModelSolverInterface );
                tElementBlocks( 1 )->finalize( tModelSolverInterface );
                tElementBlocks( 2 )->finalize( tModelSolverInterface );

                tModelSolverInterface->finalize();

                // calculate AdofMap
                Matrix< DDUMat > tAdofMap = tModelSolverInterface->get_dof_manager()->get_adof_ind_map();

                //8) Create solver interface ---------------------------------------------------
                //std::cout<<" Create solver interface "<<std::endl;
                //------------------------------------------------------------------------------

                MSI::MSI_Solver_Interface * tSolverInterface
                    = new moris::MSI::MSI_Solver_Interface( tModelSolverInterface );

                // 9) Create Nonlinear Problem -------------------------------------------------
                //std::cout<<" Create Nonlinear Problem "<<std::endl;
                //------------------------------------------------------------------------------

                NLA::Nonlinear_Problem* tNonlinearProblem
                    = new NLA::Nonlinear_Problem( tSolverInterface );

                // 10) Create Solvers and solver manager ----------------------------------------
                //std::cout<<" Create Solvers and solver manager "<<std::endl;
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

                // 11) Solve --------------------------------------------------------------------
                //std::cout<<" Solve "<<std::endl;
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

                // 12) postprocessing------------------------------------------------------------
                //std::cout<<" Postprocessing "<<std::endl;
                //------------------------------------------------------------------------------

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

                // checking the solution--------------------------------------------------------
                //------------------------------------------------------------------------------
                // Expected solution
                Matrix< DDRMat > tExpectedSolution = {{ 5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0,
                                                        5.0, 25.0, 45.0 }};

                // define an epsilon environment
                real tEpsilon = 1E-12;

                // define a bool for solution check
                bool tCheckNodalSolution = true;

                // loop over the node and chyeck solution
                for ( uint i = 0; i < tNumOfNodes; i++ )
                {
                    // check solution
                    tCheckNodalSolution = tCheckNodalSolution
                                       && ( std::abs( tTempSolutionField( i ) - tExpectedSolution( i ) ) < tEpsilon );
                }
                // check bool is true
                REQUIRE( tCheckNodalSolution );

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

                // 13) Clean up -----------------------------------------------------------------
                //std::cout<<" Clean up "<<std::endl;
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


        real
		heat_flux_function(Matrix<DDRMat> aNodeCoordinate)
        {
        	return 1;
        }

        Matrix<IndexMat>
        get_linear_tet_side_vertex_map()
        {
        	return {{0, 1, 3},
        		    {2, 1, 3},
					{0, 2, 3},
					{0, 2, 1}};
        }

        Matrix<IndexMat>
        get_linear_hex_side_vertex_map()
		{
        	return {{0,1,5,4},
        			{1,2,6,5},
					{2,3,7,6},
					{0,4,7,3},
					{0,3,2,1},
					{4,5,6,7}};
		}

        Matrix<IndexMat>
        get_nodes_on_side_set(moris::mtk::Mesh* aMesh,
        					  std::string aSideSetName)
		{
            moris::Cell<mtk::Cell*> tSideSetCells;
            Matrix<IndexMat> tSideOrdinalsOfCell;
            aMesh->get_sideset_cells_and_ords(aSideSetName,tSideSetCells,tSideOrdinalsOfCell);

            Matrix<IndexMat> tLinearHexMap = get_linear_hex_side_vertex_map();

            Matrix<IndexMat> tNonUniqueVerticesOnSideSet(1,tLinearHexMap.n_cols()*tSideSetCells.size());

            uint tCount = 0;
            // iterate through cells in side set
            for(moris::uint i =0; i < tSideSetCells.size(); i++)
            {
            	// get pointer to current cell/element
            	mtk::Cell* tCell = tSideSetCells(i);

            	// side ordinal of current cell/element
            	moris_index tSideOrdinal = tSideOrdinalsOfCell(i);

            	// vertex indices attached to the cell
            	Matrix<IndexMat> tCellVertices = tCell->get_vertex_inds();

            	for(moris::uint j = 0; j < tLinearHexMap.n_cols(); j++)
            	{
            		moris_index tVertexOrd = tLinearHexMap(tSideOrdinal,j);

            		tNonUniqueVerticesOnSideSet(tCount) = tCellVertices(tVertexOrd);
            		tCount++;
            	}
            }

            // make unique

            Matrix<IndexMat> tVerticesOnSideSet(1,1);
            unique( tNonUniqueVerticesOnSideSet,tVerticesOnSideSet );

            //moris::print(tVerticesOnSideSet,"tVerticesOnSideSet");
            //moris::print(tNonUniqueVerticesOnSideSet,"tNonUniqueVerticesOnSideSet");

            return tVerticesOnSideSet;
		}

        Matrix<DDRMat>
        compute_neumann_values(uint aNumNodes,
        					   Matrix<IndexMat> const & aVertexIndices,
        				       Matrix<DDRMat>   const & aNodeCoords)
		{
        	// allocate output
        	Matrix<DDRMat> tNeumannValues(aNumNodes,1);

        	// loop over nodes
        	for(moris::uint i=0; i<aVertexIndices.numel(); i++)
        	{
        		moris_index tVertexIndex = aVertexIndices(i);
        		tNeumannValues(tVertexIndex) = heat_flux_function(aNodeCoords.get_row(tVertexIndex));
        	}

        	return tNeumannValues;
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//2//2//2//3//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//2//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

        TEST_CASE( "Element_Diffusion_Bar", "[moris],[fem],[ElemDiff_Bar]" )
        {
            if(par_size() == 1 )
            {
                // Create a 3D mesh of HEX8 using MTK ------------------------------------------
                std::cout<<" Create a 3D mesh of five HEX8 using MTK "<<std::endl;
                //------------------------------------------------------------------------------
                uint aNumElemTypes = 1; // only 1 element type ( quad )
                uint aNumDim = 3;       // number of spatial dimensions

                // element connectivity
                Matrix< IdMat > aElementConnQuad = {{ 1, 2,  8,  7, 13, 14, 20, 19 },
                                                    { 2, 3,  9,  8, 14, 15, 21, 20 },
                                                    { 3, 4, 10,  9, 15, 16, 22, 21 },
                                                    { 4, 5, 11, 10, 16, 17, 23, 22 },
                                                    { 5, 6, 12, 11, 17, 18, 24, 23 }};


                // local to global element map for quads
                Matrix< IdMat > aElemLocalToGlobalQuad = { { 1, 2, 3, 4, 5} };

                // node coordinates
                Matrix< DDRMat > aCoords = {{ 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 },
                                            { 3.0, 0.0, 0.0 }, { 4.0, 0.0, 0.0 }, { 5.0, 0.0, 0.0 },
											{ 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 2.0, 1.0, 0.0 },
											{ 3.0, 1.0, 0.0 }, { 4.0, 1.0, 0.0 }, { 5.0, 1.0, 0.0 },
											{ 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 2.0, 0.0, 1.0 },
											{ 3.0, 0.0, 1.0 }, { 4.0, 0.0, 1.0 }, { 5.0, 0.0, 1.0 },
											{ 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 2.0, 1.0, 1.0 },
											{ 3.0, 1.0, 1.0 }, { 4.0, 1.0, 1.0 }, { 5.0, 1.0, 1.0 },};

                // specify the local to global map
                Matrix< IdMat > aNodeLocalToGlobal = {{ 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                       19, 20, 21, 22, 23, 24}};

                // Add sidesets to the mesh
                moris::mtk::MtkSetsInfo tMtkMeshSets;

                // Fixed bc
                moris::mtk::MtkSideSetInfo tFixed;
                Matrix<IndexMat> tFixedElementsAndOrds = {{1, 3}};
                tFixed.mElemIdsAndSideOrds = &tFixedElementsAndOrds;
                tFixed.mSideSetName        = "fixed" ;
                tMtkMeshSets.add_side_set(&tFixed);

                // Fixed bc
                moris::mtk::MtkSideSetInfo tHeatFluxSS;
                Matrix<IndexMat> tFluxCellAndOrds = {{5,1}};
                tHeatFluxSS.mElemIdsAndSideOrds = &tFluxCellAndOrds;
                tHeatFluxSS.mSideSetName        = "heat_flux" ;
                tMtkMeshSets.add_side_set(&tHeatFluxSS);

                // create mesh MTK database
                mtk::MtkMeshData tMeshData( aNumElemTypes );
                tMeshData.CreateAllEdgesAndFaces  = true;
                tMeshData.SpatialDim              = & aNumDim;
                tMeshData.ElemConn( 0 )           = & aElementConnQuad;
                tMeshData.NodeCoords              = & aCoords;
                tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
                tMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;
//                tMeshData.SetsInfo                = &tMtkMeshSets;

                mtk::Mesh* tMesh = create_mesh( MeshType::STK, tMeshData );


                std::string tMeshOutputFile = "./diffuse_test.exo";
                tMesh->create_output_mesh(tMeshOutputFile);

                //1) Create the fem nodes ------------------------------------------------------
                //std::cout<<" Create the fem nodes "<<std::endl;
                //------------------------------------------------------------------------------

                // number of mesh nodes
                uint tNumOfNodes = tMesh->get_num_nodes();
                //std::cout<<"Number of Nodes:"<< tNumOfNodes <<std::endl;

                //create a cell of fem nodes
                moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

                // loop over the mesh nodes
                for( uint k = 0; k < tNumOfNodes; k++ )
                {
                    // create a fem node for each mesh node
                    tNodes( k ) = new fem::Node( & tMesh->get_mtk_vertex( k ) );
                }

                //2) Create the IWGs -----------------------------------------------------------
                std::cout<<" 2) Create the IWGs "<<std::endl;
                //------------------------------------------------------------------------------

                // input a cell of IWG types to be created
                Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
                                                       fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                       fem::IWG_Type::SPATIALDIFF_NEUMANN,
  													   fem::IWG_Type::SPATIALDIFF_GHOST_LEFT,
													   fem::IWG_Type::SPATIALDIFF_GHOST_RIGHT};

                // number of IWGs to be created
                uint tNumOfIWGs = tIWGTypeList.size();

                std::cout<<"Length of IWG list:"<<tIWGTypeList.size()<<std::endl;

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

                //3) Create element blocks -------------------------------------------------------
                std::cout<<" 3) Create element blocks "<<std::endl;
                //------------------------------------------------------------------------------
                // get the number of elements on proc from the mesh
                luint tNumOfElements = tMesh->get_num_elems();

                // create equation objects
                Cell< MSI::Equation_Object* > tElements;
                tElements.reserve( tNumOfElements + 1 + 2 + 1 );

                // get the block names from the mesh
                moris::Cell<std::string> tBlockSetsNames = tMesh->get_set_names( EntityRank::ELEMENT);

                // Cell containing the block mesh cell ( a cell of mesh cells )
                moris::Cell<mtk::Cell*> tBlockSetElement( tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( 0 ) ).numel(), nullptr );

                // loop on the blocks
                for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
                {
                    Matrix< IndexMat > tBlockSetElementInd = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( Ik ) );

                    // loop on the elements in a block
                    for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                    {
                        // Cell containing the block mesh cell with mesh cells
                        tBlockSetElement( k ) = & tMesh->get_mtk_cell( k );
                    }
                }

                // create a fem element block
                moris::Cell< fem::Element_Block * > tElementBlocks( 5, nullptr );
                moris::Cell< IWG* > tIWGBulk = { tIWGs( 0 ) };
                tElementBlocks( 0 ) = new fem::Element_Block( tBlockSetElement,
                                                              fem::Element_Type::BULK,
                                                              tIWGBulk,
                                                              tNodes );

                // put the equation object of block 0 in the global list of equation objects
                tElements.append( tElementBlocks( 0 )->get_equation_object_list() );

                //4) Create Dirichlet side blocks -------------------------------------------------------
                std::cout<<" 4) Create Dirichlet side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Dirichlet sideset
                moris::Cell<mtk::Cell*> tSideSetDirichletElement( 1, nullptr );
                tSideSetDirichletElement( 0 ) = & tMesh->get_mtk_cell( 0 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGDirichlet = { tIWGs( 1 ) };
                tElementBlocks( 1 ) = new fem::Element_Block( tSideSetDirichletElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGDirichlet,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 1 )->get_equation_object_list() );

                //5) Create Neumann side blocks -------------------------------------------------------
                std::cout<<" 5) Create Neumann side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Neumann sideset
                moris::Cell<mtk::Cell*> tSideSetNeumannElement( 1, nullptr );
                tSideSetNeumannElement( 0 ) = & tMesh->get_mtk_cell( 4 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGNeumann = { tIWGs( 2 ) };
                tElementBlocks( 2 ) = new fem::Element_Block( tSideSetNeumannElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGNeumann,
                                                              tNodes );

                // put the equation object of block 2 in the global list of equation objects
                tElements.append( tElementBlocks( 2 )->get_equation_object_list() );

                //6) Create left Ghost side blocks -------------------------------------------------------
                std::cout<<" 6) Create left Ghost side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the left Ghost sideset
                moris::Cell<mtk::Cell*> tSideSetLeftGhostElement( 1, nullptr );
                tSideSetLeftGhostElement( 0 ) = & tMesh->get_mtk_cell( 3 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGLeftGhost = { tIWGs( 3 ) };
                tElementBlocks( 3 ) = new fem::Element_Block( tSideSetLeftGhostElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGLeftGhost,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 3 )->get_equation_object_list() );

                //7) Create right Ghost side blocks -------------------------------------------------------
                std::cout<<" 7) Create right Ghost side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the right Ghost sideset
                moris::Cell<mtk::Cell*> tSideSetRightGhostElement( 1, nullptr );
                tSideSetRightGhostElement( 0 ) = & tMesh->get_mtk_cell( 4 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGRightGhost = { tIWGs( 4 ) };
                tElementBlocks( 4 ) = new fem::Element_Block( tSideSetRightGhostElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGRightGhost,
                                                              tNodes );

                // put the equation object of block 2 in the global list of equation objects
                tElements.append( tElementBlocks( 4 )->get_equation_object_list() );

                //8) Impose boundary conditions -----------------------------------------------
                std::cout<<" 8) Impose boundary conditions "<<std::endl;
                //-----------------------------------------------------------------------------
                // nodal weak bc for Dirichlet
                real tTempValue = 0.0;
                Cell< moris_index > tDirichletElements = { 0 };
                for( uint iDirichlet = 0; iDirichlet < 1; iDirichlet++ )
                {
                     // set the list of face ordinals
                     tElements( tNumOfElements + iDirichlet )->set_list_of_side_ordinals( {{ 3 }} );

                     //get the nodal weak bcs of the element
                     Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumOfElements + iDirichlet )->get_weak_bcs();

                     // get the element number of nodes
                     uint tNumberOfNodes = tElements( tNumOfElements + iDirichlet )->get_num_nodes();

                     // set size of the element nodal weak bc
                     tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                     // loop over the element nodes
                     Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tDirichletElements( iDirichlet ) ).get_vertex_inds();

                     for( uint l = 0; l < tNumberOfNodes; l++ )
                     {
                         // copy weak bc into element
                         tNodalWeakBCs( l ) = tTempValue;
                     }
                }

                // nodal weak bc for Neumann
                real tFluxValue = 1.0;
                Cell< moris_index > tNeumannElements = { 4 };
                for( uint iNeumann = 0; iNeumann < 1; iNeumann++ )
                {
                     // set the list of face ordinals
                     tElements( tNumOfElements + tDirichletElements.size() + iNeumann )->set_list_of_side_ordinals( {{ 1 }} );

                     //get the nodal weak bcs of the element
                     Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumOfElements + tNeumannElements.size() + iNeumann )->get_weak_bcs();

                     // get the element number of nodes
                     uint tNumberOfNodes = tElements( tNumOfElements + tNeumannElements.size() + iNeumann )->get_num_nodes();

                     // set size of the element nodal weak bc
                     tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                     // loop over the element nodes
                     Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tNeumannElements( iNeumann ) ).get_vertex_inds();

                     for( uint l = 0; l < tNumberOfNodes; l++ )
                     {
                         // copy weak bc into element
                         tNodalWeakBCs( l ) = tFluxValue;
                     }
                }

                //  left Ghost
                     // set the list of face ordinals
                     tElements( 7 )->set_list_of_side_ordinals( {{ 1 }} );

                //  right Ghost
                     // set the list of face ordinals
                     tElements( 8 )->set_list_of_side_ordinals( {{ 3 }} );



                //9) Create the model solver interface -----------------------------------------
                std::cout<<" 9) Create the model solver interface "<<std::endl;
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

                tElementBlocks( 0 )->finalize( tModelSolverInterface );
                tElementBlocks( 1 )->finalize( tModelSolverInterface );
                tElementBlocks( 2 )->finalize( tModelSolverInterface );
                tElementBlocks( 3 )->finalize( tModelSolverInterface );
                tElementBlocks( 4 )->finalize( tModelSolverInterface );

                tModelSolverInterface->finalize();
                // calculate AdofMap
                Matrix< DDUMat > tAdofMap = tModelSolverInterface->get_dof_manager()->get_adof_ind_map();

                //10) Create solver interface ---------------------------------------------------
                std::cout<<" 10) Create solver interface "<<std::endl;
                //------------------------------------------------------------------------------

                MSI::MSI_Solver_Interface * tSolverInterface
                    = new moris::MSI::MSI_Solver_Interface( tModelSolverInterface );

                // 11) Create Nonlinear Problem -------------------------------------------------
                std::cout<<" 11) Create Nonlinear Problem "<<std::endl;
                //------------------------------------------------------------------------------

                NLA::Nonlinear_Problem* tNonlinearProblem
                    = new NLA::Nonlinear_Problem( tSolverInterface );

                // 12) Create Solvers and solver manager ----------------------------------------
                std::cout<<" 12) Create Solvers and solver manager "<<std::endl;
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

                // 13) Solve --------------------------------------------------------------------
                std::cout<<" 13) Solve "<<std::endl;
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

                // 14) postprocessing------------------------------------------------------------
                std::cout<<" 14) Postprocessing "<<std::endl;
                //------------------------------------------------------------------------------

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

                // checking the solution--------------------------------------------------------
                //------------------------------------------------------------------------------
                // Expected solution
                Matrix< DDRMat > tExpectedSolution = {{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }};

                print(tTempSolutionField,"Solution:");

                // define an epsilon for the deviation of the numerical solution
                double tEpsilon = 1E-04;
                double tEpsilonNodal = 1E-02;

                // define a bool for solution check
                bool tCheckNodalSolution = true;
                bool tCheckNormSolution  = true;

                // initialize L2-norms
                real tExpectedL2Norm = 0.0;
                real tSolutionL2Norm = 0.0;
                real tSolutionL2Error = 0.0;

                // loop over the nodes to compute norms and check node solutions
                for ( uint i = 0; i < 5; i++ )
                {
                	tExpectedL2Norm = tExpectedL2Norm + std::pow( (tExpectedSolution(i) + tExpectedSolution(i+1))/2 ,2);
                	tSolutionL2Norm = tSolutionL2Norm + std::pow( (tTempSolutionField(i) + tTempSolutionField(i+1))/2 ,2);
                	tSolutionL2Error = tSolutionL2Error + std::pow( (tTempSolutionField(i) + tTempSolutionField(i+1))/2
                		                                            -(tExpectedSolution(i) + tExpectedSolution(i+1))/2 ,2);

                	// check nodal solution
                    tCheckNodalSolution = tCheckNodalSolution
                                          && ( std::abs( tTempSolutionField( i ) - tExpectedSolution( i ) ) < tEpsilonNodal );
                }

                //L2-Norm, relative L2-error should be smaller than Epsilon
                tSolutionL2Norm = std::sqrt(tSolutionL2Norm);
                tExpectedL2Norm = std::sqrt(tExpectedL2Norm);
                tSolutionL2Error = std::sqrt(tSolutionL2Error);
                real tNormalizedL2Error = tSolutionL2Error / tExpectedL2Norm;

                std::cout<<"Expected L2-Norm: "<<tExpectedL2Norm<<std::endl;
                std::cout<<"Solution L2-Norm: "<<tSolutionL2Norm<<std::endl;
                std::cout<<"Normalized L2-Error: "<<tNormalizedL2Error<<std::endl;

                // check L2 Error
                tCheckNormSolution = tCheckNormSolution
                                      && ( std::abs( tNormalizedL2Error) < tEpsilon );

                // check if test has passed
                CHECK( tCheckNormSolution );
                CHECK( tCheckNodalSolution );

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

                // 10) Clean up -----------------------------------------------------------------
                //std::cout<<" 10) Clean up "<<std::endl;
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


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//3//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

        TEST_CASE( "Element_Diffusion_Bar_quadratic", "[moris],[fem],[ElemDiff_Bar_quad]" )
        {
            if(par_size() == 1 )
            {
                // Create a 3D mesh of HEX8 using MTK ------------------------------------------
                std::cout<<" Create a 3D mesh of five HEX8 using MTK "<<std::endl;
                //------------------------------------------------------------------------------
                uint aNumElemTypes = 1; // only 1 element type ( quad )
                uint aNumDim = 3;       // number of spatial dimensions

                // element connectivity - check
                Matrix< IdMat > aElementConnQuad = {{ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
                		                                  15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27},
                                                    { 2, 28, 29,  3,  6, 30, 31,  7, 32, 33, 34, 10, 14, 35,
                                                    	  36, 15, 37, 38, 39, 18, 40, 41, 42, 25, 43, 44, 45}};


                // local to global element map for quads - check
                Matrix< IdMat > aElemLocalToGlobalQuad = { { 1, 2} };

                // node coordinates - check
                Matrix< DDRMat > aCoords = {{  0,  0,  0}, {  1,  0,  0}, {  1,  1,  0}, {  0,  1,  0}, // 1-4
                		                    {  0,  0,  1}, {  1,  0,  1}, {  1,  1,  1}, {  0,  1,  1}, // 5-8
                                            {0.5,  0,  0}, {  1,0.5,  0}, {0.5,  1,  0}, {  0,0.5,  0}, // 9-12
											{  0,  0,0.5}, {  1,  0,0.5}, {  1,  1,0.5}, {  0,  1,0.5}, // 13-16
											{0.5,  0,  1}, {  1,0.5,  1}, {0.5,  1,  1}, {  0,0.5,  1}, // 17-20
											{0.5,0.5,0.5}, {0.5,0.5,  0}, {0.5,0.5,  1}, {  0,0.5,0.5}, // 21-24
											{  1,0.5,0.5}, {0.5,  0,0.5}, {0.5,  1,0.5}, {  2,  0,  0}, // 25-28
											{  2,  1,  0}, {  2,  0,  1}, {  2,  1,  1}, {1.5,  0,  0}, // 29-32
											{  2,0.5,  0}, {1.5,  1,  0}, {  2,  0,0.5}, {  2,  1,0.5}, // 33-36
											{1.5,  0,  1}, {  2,0.5,  1}, {1.5,  1,  1}, {1.5,0.5,0.5}, // 37-40
											{1.5,0.5,  0}, {1.5,0.5,  1}, {  2,0.5,0.5}, {1.5,  0,0.5}, // 41-44
											{1.5,  1,0.5}                                             };// 45

                // specify the local to global map - check
                Matrix< IdMat > aNodeLocalToGlobal = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                		                              12, 13, 14, 15, 16, 17, 18, 19, 20,
													  21, 22, 23, 24, 25, 26, 27, 28, 29,
													  30, 31, 32, 33, 34, 35, 36, 37, 38,
													          39, 40, 41, 42, 43, 44, 45}};

                // Add sidesets to the mesh
                moris::mtk::MtkSetsInfo tMtkMeshSets;

                // Fixed bc
                moris::mtk::MtkSideSetInfo tFixed;
                Matrix<IndexMat> tFixedElementsAndOrds = {{1, 3}};
                tFixed.mElemIdsAndSideOrds = &tFixedElementsAndOrds;
                tFixed.mSideSetName        = "fixed" ;
                tMtkMeshSets.add_side_set(&tFixed);

                // Fixed bc
                moris::mtk::MtkSideSetInfo tHeatFluxSS;
                Matrix<IndexMat> tFluxCellAndOrds = {{5,1}};
                tHeatFluxSS.mElemIdsAndSideOrds = &tFluxCellAndOrds;
                tHeatFluxSS.mSideSetName        = "heat_flux" ;
                tMtkMeshSets.add_side_set(&tHeatFluxSS);

                // create mesh MTK database
                mtk::MtkMeshData tMeshData( aNumElemTypes );
                tMeshData.CreateAllEdgesAndFaces  = true;
                tMeshData.SpatialDim              = & aNumDim;
                tMeshData.ElemConn( 0 )           = & aElementConnQuad;
                tMeshData.NodeCoords              = & aCoords;
                tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
                tMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal;
//                tMeshData.SetsInfo                = &tMtkMeshSets;

                mtk::Mesh* tMesh = create_mesh( MeshType::STK, tMeshData );


                std::string tMeshOutputFile = "./diffuse_test.exo";
                tMesh->create_output_mesh(tMeshOutputFile);

                //1) Create the fem nodes ------------------------------------------------------
                //std::cout<<" Create the fem nodes "<<std::endl;
                //------------------------------------------------------------------------------

                // number of mesh nodes
                uint tNumOfNodes = tMesh->get_num_nodes();
                //std::cout<<"Number of Nodes:"<< tNumOfNodes <<std::endl;

                //create a cell of fem nodes
                moris::Cell< Node_Base* > tNodes( tNumOfNodes, nullptr );

                // loop over the mesh nodes
                for( uint k = 0; k < tNumOfNodes; k++ )
                {
                    // create a fem node for each mesh node
                    tNodes( k ) = new fem::Node( & tMesh->get_mtk_vertex( k ) );
                }

                //2) Create the IWGs -----------------------------------------------------------
                std::cout<<" 2) Create the IWGs "<<std::endl;
                //------------------------------------------------------------------------------

                // input a cell of IWG types to be created
                Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
                                                       fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                       fem::IWG_Type::SPATIALDIFF_NEUMANN,
  													   fem::IWG_Type::SPATIALDIFF_GHOST_LEFT,
													   fem::IWG_Type::SPATIALDIFF_GHOST_RIGHT};

                // number of IWGs to be created
                uint tNumOfIWGs = tIWGTypeList.size();

                std::cout<<"Length of IWG list:"<<tIWGTypeList.size()<<std::endl;

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

                //3) Create element blocks -------------------------------------------------------
                std::cout<<" 3) Create element blocks "<<std::endl;
                //------------------------------------------------------------------------------
                // get the number of elements on proc from the mesh
                luint tNumOfElements = tMesh->get_num_elems();

                // create equation objects (bulk of each element + Dirichlet + Neumann + Ghosts)
                Cell< MSI::Equation_Object* > tElements;
                tElements.reserve( tNumOfElements + 1 + 1 + 2 );

                // get the block names from the mesh
                moris::Cell<std::string> tBlockSetsNames = tMesh->get_set_names( EntityRank::ELEMENT);

                // Cell containing the block mesh cell ( a cell of mesh cells )
                moris::Cell<mtk::Cell*> tBlockSetElement( tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT,
                		                                  tBlockSetsNames( 0 ) ).numel(), nullptr );

                // loop on the blocks
                for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
                {
                    Matrix< IndexMat > tBlockSetElementInd = tMesh->get_set_entity_loc_inds( EntityRank::ELEMENT,
                    		                                                                 tBlockSetsNames( Ik ) );

                    // loop on the elements in a block
                    for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                    {
                        // Cell containing the block mesh cell with mesh cells
                        tBlockSetElement( k ) = & tMesh->get_mtk_cell( k );
                    }
                }

                // create a fem element block
                moris::Cell< fem::Element_Block * > tElementBlocks( 5, nullptr );
                moris::Cell< IWG* > tIWGBulk = { tIWGs( 0 ) };
                tElementBlocks( 0 ) = new fem::Element_Block( tBlockSetElement,
                                                              fem::Element_Type::BULK,
                                                              tIWGBulk,
                                                              tNodes );

                // put the equation object of block 0 in the global list of equation objects
                tElements.append( tElementBlocks( 0 )->get_equation_object_list() );

                //4) Create Dirichlet side blocks -------------------------------------------------------
                std::cout<<" 4) Create Dirichlet side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Dirichlet sideset
                moris::Cell<mtk::Cell*> tSideSetDirichletElement( 1, nullptr );
                tSideSetDirichletElement( 0 ) = & tMesh->get_mtk_cell( 0 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGDirichlet = { tIWGs( 1 ) };
                tElementBlocks( 1 ) = new fem::Element_Block( tSideSetDirichletElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGDirichlet,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 1 )->get_equation_object_list() );

                //5) Create Neumann side blocks -------------------------------------------------------
                std::cout<<" 5) Create Neumann side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Neumann sideset
                moris::Cell<mtk::Cell*> tSideSetNeumannElement( 1, nullptr );
                tSideSetNeumannElement( 0 ) = & tMesh->get_mtk_cell( 1 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGNeumann = { tIWGs( 2 ) };
                tElementBlocks( 2 ) = new fem::Element_Block( tSideSetNeumannElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGNeumann,
                                                              tNodes );

                // put the equation object of block 2 in the global list of equation objects
                tElements.append( tElementBlocks( 2 )->get_equation_object_list() );

                //6) Create left Ghost side blocks -------------------------------------------------------
                std::cout<<" 6) Create left Ghost side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the left Ghost sideset
                moris::Cell<mtk::Cell*> tSideSetLeftGhostElement( 1, nullptr );
                tSideSetLeftGhostElement( 0 ) = & tMesh->get_mtk_cell( 0 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGLeftGhost = { tIWGs( 3 ) };
                tElementBlocks( 3 ) = new fem::Element_Block( tSideSetLeftGhostElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGLeftGhost,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 3 )->get_equation_object_list() );

                //7) Create right Ghost side blocks -------------------------------------------------------
                std::cout<<" 7) Create right Ghost side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the right Ghost sideset
                moris::Cell<mtk::Cell*> tSideSetRightGhostElement( 1, nullptr );
                tSideSetRightGhostElement( 0 ) = & tMesh->get_mtk_cell( 1 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGRightGhost = { tIWGs( 4 ) };
                tElementBlocks( 4 ) = new fem::Element_Block( tSideSetRightGhostElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGRightGhost,
                                                              tNodes );

                // put the equation object of block 2 in the global list of equation objects
                tElements.append( tElementBlocks( 4 )->get_equation_object_list() );

                //8) Impose boundary conditions -----------------------------------------------
                std::cout<<" 8) Impose boundary conditions "<<std::endl;
                //-----------------------------------------------------------------------------
                // nodal weak bc for Dirichlet
                real tTempValue = 0.0;
                Cell< moris_index > tDirichletElements = { 0 };
                for( uint iDirichlet = 0; iDirichlet < 1; iDirichlet++ )
                {
                     // set the list of face ordinals
                     tElements( tNumOfElements + iDirichlet )->set_list_of_side_ordinals( {{ 3 }} );

                     //get the nodal weak bcs of the element
                     Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumOfElements + iDirichlet )->get_weak_bcs();

                     // get the element number of nodes
                     uint tNumberOfNodes = tElements( tNumOfElements + iDirichlet )->get_num_nodes();

                     // set size of the element nodal weak bc
                     tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                     // loop over the element nodes
                     Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tDirichletElements( iDirichlet ) ).get_vertex_inds();

                     for( uint l = 0; l < tNumberOfNodes; l++ )
                     {
                         // copy weak bc into element
                         tNodalWeakBCs( l ) = tTempValue;
                     }
                }

                // nodal weak bc for Neumann
                real tFluxValue = 1.0;
                Cell< moris_index > tNeumannElements = { 1 };
                for( uint iNeumann = 0; iNeumann < 1; iNeumann++ )
                {
                     // set the list of face ordinals
                     tElements( tNumOfElements + tDirichletElements.size() + iNeumann )->set_list_of_side_ordinals( {{ 1 }} );

                     //get the nodal weak bcs of the element
                     Matrix< DDRMat > & tNodalWeakBCs = tElements( tNumOfElements + tNeumannElements.size() + iNeumann )->get_weak_bcs();

                     // get the element number of nodes
                     uint tNumberOfNodes = tElements( tNumOfElements + tNeumannElements.size() + iNeumann )->get_num_nodes();

                     // set size of the element nodal weak bc
                     tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                     // loop over the element nodes
                     Matrix< IndexMat > tNodeIndices = tMesh->get_mtk_cell( tNeumannElements( iNeumann ) ).get_vertex_inds();

                     for( uint l = 0; l < tNumberOfNodes; l++ )
                     {
                         // copy weak bc into element
                         tNodalWeakBCs( l ) = tFluxValue;
                     }
                }

                //  left Ghost
                     // set the list of face ordinals
                     tElements( tNumOfElements + 2 )->set_list_of_side_ordinals( {{ 3 }} );


                //  right Ghost
                     // set the list of face ordinals
                     tElements( tNumOfElements + 3 )->set_list_of_side_ordinals( {{ 1 }} );


                //9) Create the model solver interface -----------------------------------------
                std::cout<<" 9) Create the model solver interface "<<std::endl;
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

                tElementBlocks( 0 )->finalize( tModelSolverInterface );
                tElementBlocks( 1 )->finalize( tModelSolverInterface );
                tElementBlocks( 2 )->finalize( tModelSolverInterface );
                tElementBlocks( 3 )->finalize( tModelSolverInterface );
                tElementBlocks( 4 )->finalize( tModelSolverInterface );

                tModelSolverInterface->finalize();
                // calculate AdofMap
                Matrix< DDUMat > tAdofMap = tModelSolverInterface->get_dof_manager()->get_adof_ind_map();

                //10) Create solver interface ---------------------------------------------------
                std::cout<<" 10) Create solver interface "<<std::endl;
                //------------------------------------------------------------------------------

                MSI::MSI_Solver_Interface * tSolverInterface
                    = new moris::MSI::MSI_Solver_Interface( tModelSolverInterface );

                // 11) Create Nonlinear Problem -------------------------------------------------
                std::cout<<" 11) Create Nonlinear Problem "<<std::endl;
                //------------------------------------------------------------------------------

                NLA::Nonlinear_Problem* tNonlinearProblem
                    = new NLA::Nonlinear_Problem( tSolverInterface );

                // 12) Create Solvers and solver manager ----------------------------------------
                std::cout<<" 12) Create Solvers and solver manager "<<std::endl;
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

                // 13) Solve --------------------------------------------------------------------
                std::cout<<" 13) Solve "<<std::endl;
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

                // 14) postprocessing------------------------------------------------------------
                std::cout<<" 14) Postprocessing "<<std::endl;
                //------------------------------------------------------------------------------

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
                            = tElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i, tDofTypeList);
                        // add up the contribution of each element to the node value
                        tNodeVal = tNodeVal + tElemVal;
                    }
                    // fill the solution matrix with the node value
                    tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
                }

                // checking the solution--------------------------------------------------------
                //------------------------------------------------------------------------------
                // Expected solution (solution computed with same ghost penalty in 1D)
                Matrix< DDRMat > tExpectedSolution = {{0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                		                               0.5, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0,
													   0.5, 1.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0,
													   1.0, 0.5, 0.5, 2.0, 2.0, 2.0, 2.0, 1.5,
													   2.0, 1.5, 2.0, 2.0, 1.5, 2.0, 1.5, 1.5,
													   1.5, 1.5, 2.0, 1.5, 1.5}};


                //print(tTempSolutionField,"Solution:");
                //print(tExpectedSolution,"Expected Solution:");

                // define an epsilon for the deviation of the solution
                //double tEpsilon = 1E-04;
                double tEpsilonNodal = 2E-04;

                // define a bool for solution check
                bool tCheckNodalSolution = true;
                //bool tCheckNormSolution  = true;

                // initialize L2-norms
                //real tExpectedL2Norm = 0.0;
                //real tSolutionL2Norm = 0.0;
                //real tSolutionL2Error = 0.0;

                // loop over the nodes to compute norms and check node solutions
                for ( uint i = 0; i < 45; i++ )
                {
                	//tExpectedL2Norm = tExpectedL2Norm + std::pow( (tExpectedSolution(i) + tExpectedSolution(i+1))/2 ,2);
                	//tSolutionL2Norm = tSolutionL2Norm + std::pow( (tTempSolutionField(i) + tTempSolutionField(i+1))/2 ,2);
                	//tSolutionL2Error = tSolutionL2Error + std::pow( (tTempSolutionField(i) + tTempSolutionField(i+1))/2
                	//	                                            -(tExpectedSolution(i) + tExpectedSolution(i+1))/2 ,2);

                	//std::cout<<" value number: "<< i <<std::endl;
                	//std::cout<<" deviation from expected: "<< std::abs(tTempSolutionField(i)-tExpectedSolution(i)) <<std::endl;

                	// check nodal solution
                    tCheckNodalSolution = tCheckNodalSolution
                                          && ( std::abs( tTempSolutionField( i ) - tExpectedSolution( i ) ) < tEpsilonNodal );
                }

                //L2-Norm, relative L2-error should be smaller than Epsilon
                //tSolutionL2Norm = std::sqrt(tSolutionL2Norm);
                //tExpectedL2Norm = std::sqrt(tExpectedL2Norm);
                //tSolutionL2Error = std::sqrt(tSolutionL2Error);
                //real tNormalizedL2Error = tSolutionL2Error / tExpectedL2Norm;

                //std::cout<<"Expected L2-Norm: "<<tExpectedL2Norm<<std::endl;
                //std::cout<<"Solution L2-Norm: "<<tSolutionL2Norm<<std::endl;
                //std::cout<<"Normalized L2-Error: "<<tNormalizedL2Error<<std::endl;

                // check L2 Error
                //tCheckNormSolution = tCheckNormSolution
                //                      && ( std::abs( tNormalizedL2Error) < tEpsilon );

                // check if test has passed
                //REQUIRE( tCheckNormSolution );
                REQUIRE( tCheckNodalSolution );

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

                // 10) Clean up -----------------------------------------------------------------
                //std::cout<<" 10) Clean up "<<std::endl;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//4//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Test for 1x1x4 Bar with cut last element

        TEST_CASE( "Element_Diffusion_Bar_cut", "[moris],[fem],[ElemDiff_Bar_cut]" )
            {
        //namespace mtk {


            // Define the Interpolation Mesh ---------------------------------------------------------
            std::string tInterpString = "Bar_1x1x4";
            //moris::mtk::Mesh* tInterpMesh = create_mesh( MeshType::STK, tInterpString, NULL );

            uint aNumElemTypes = 1; // only 1 element type ( quad )
            uint aNumDim = 3;       // number of spatial dimensions

            // element connectivity
            Matrix< IdMat > aElementConnQuad_InterpolMesh = {{  1,  2,  3,  4,  5,  6,  7,  8 },
                                                             {  5,  6,  7,  8,  9, 10, 11, 12 },
                                                             {  9, 10, 11, 12, 13, 14, 15, 16 },
                                                             { 13, 14, 15, 16, 17, 18, 19, 20 } };


            // local to global element map for quads
            Matrix< IdMat > aElemLocalToGlobalQuad_InterpolMesh = { { 1, 2, 3, 4} };

            // node coordinates
            Matrix< DDRMat > aCoords_InterpolMesh = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
            		                                 {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
													 {0,0,2}, {1,0,2}, {1,1,2}, {0,1,2},
													 {0,0,3}, {1,0,3}, {1,1,3}, {0,1,3},
													 {0,0,4}, {1,0,4}, {1,1,4}, {0,1,4}};

            // specify the local to global map
            Matrix< IdMat > aNodeLocalToGlobal_InterpolMesh = {{ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
            		                                            11, 12, 13, 14, 15, 16, 17, 18, 19, 20}};

            // Add sidesets to the mesh
            moris::mtk::MtkSetsInfo tMtkMeshSets_InterpolMesh;

            // create mesh MTK database
            mtk::MtkMeshData tInterpolMeshData( aNumElemTypes );
            tInterpolMeshData.CreateAllEdgesAndFaces  = true;
            tInterpolMeshData.SpatialDim              = & aNumDim;
            tInterpolMeshData.ElemConn( 0 )           = & aElementConnQuad_InterpolMesh;
            tInterpolMeshData.NodeCoords              = & aCoords_InterpolMesh;
            tInterpolMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad_InterpolMesh;
            tInterpolMeshData.LocaltoGlobalNodeMap    = & aNodeLocalToGlobal_InterpolMesh;

            mtk::Mesh* tInterpolMesh = create_mesh( MeshType::STK, tInterpolMeshData );


            std::string tInterpolMeshOutputFile = "./DiffusionBar_Interpolation.exo";
            tInterpolMesh->create_output_mesh(tInterpolMeshOutputFile);

            // Define the Integration Mesh -----------------------------------------------------------
            Matrix<DDRMat> tIntegNodeCoordinates ={{0, 0, 0},
                                              {1, 0, 0},
                                              {0, 1, 0},
                                              {1, 1, 0},
                                              {0, 0, 1},
                                              {1, 0, 1},
                                              {0, 1, 1},
                                              {1, 1, 1},
                                              {0, 0, 2},
                                              {1, 0, 2},
                                              {0, 1, 2},
                                              {1, 1, 2},
                                              {0, 0, 3},
                                              {1, 0, 3},
                                              {0, 1, 3},
                                              {1, 1, 3},
                                              {0, 0, 4},
                                              {1, 0, 4},
                                              {0, 1, 4},
                                              {1, 1, 4},
                                              {0.5, 0, 3.5},
                                              {1, 0.5, 3.5},
                                              {0.5, 1, 3.5},
                                              {0, 0.5, 3.5},
                                              {0.5, 0.5, 3},
                                              {0.5, 0.5, 4},
                                               {0.5, 0.5, 3.5},
                                              {0.1, 0, 3.1},
                                              {0.9, 0, 3.1},
                                              {0.1, 0.1, 3.1},
                                              {0.9, 0.1, 3.1},
                                              {1, 0, 3.1},
                                              {0, 0, 3.1},
                                              {1, 0.1, 3.1},
                                              {1, 0.9, 3.1},
                                              {0.9, 0.9, 3.1},
                                              {1, 1, 3.1},
                                              {0.9, 1, 3.1},
                                              {0.1, 1, 3.1},
                                              {0.1, 0.9, 3.1},
                                              {0, 1, 3.1},
                                              {0, 0.9, 3.1},
                                              {0, 0.1, 3.1},
                                              {0.5, 0.5, 3.1}};


            Matrix<IndexMat> tIntegLocalToGlobalNodeMap = {{ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
            		                                   12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
													   23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
													   34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};

            Matrix<IndexMat> tElemsAsCellIds     = {{1,2,3,4}};
            Matrix<IndexMat> tElemsAsCellToNodes = {{1, 2, 4, 3, 5, 6, 8, 7},
                                                               {5, 6, 8, 7, 9, 10, 12, 11},
                                                               {9, 10, 12, 11, 13, 14, 16, 15},
                                                               {13, 14, 16, 15, 17, 18, 20, 19}};

            // Tetrahedral cells in material phase 0 (Material) ???????
            CellTopology     tPhase0ChildTopo  = CellTopology::TET4;
            Matrix<IndexMat> tCellIdsPhase0    = {{ 6,  8, 10, 12, 14, 16, 17, 18, 20, 31, 32,
            		                               33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74,
												   75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};

            Matrix<IndexMat> tCellToNodePhase0 = {{29, 14, 31, 32},
                                                  {13, 28, 30, 33},
                                                  {35, 16, 36, 37},
                                                  {14, 31, 32, 34},
                                                  {16, 36, 37, 38},
                                                  {39, 15, 40, 41},
                                                  {13, 30, 40, 43},
                                                  {15, 40, 41, 42},
                                                  {30, 13, 33, 43},
                                                  {13, 14, 31, 29},
                                                  {28, 13, 31, 29},
                                                  {13, 28, 31, 30},
                                                  {16, 14, 35, 36},
                                                  {14, 31, 34, 36},
                                                  {14, 34, 35, 36},
                                                  {15, 16, 39, 40},
                                                  {16, 36, 38, 40},
                                                  {16, 38, 39, 40},
                                                  {15, 13, 40, 43},
                                                  {15, 40, 42, 43},
                                                  {13, 30, 31, 44},
                                                  {14, 13, 31, 44},
                                                  {25, 13, 14, 44},
                                                  {14, 31, 36, 44},
                                                  {16, 14, 36, 44},
                                                  {25, 14, 16, 44},
                                                  {16, 36, 40, 44},
                                                  {15, 16, 40, 44},
                                                  {25, 16, 15, 44},
                                                  {30, 13, 40, 44},
                                                  {13, 15, 40, 44},
                                                  {13, 25, 15, 44}};


            // Tetrahedral cells in material phase 1 (Void) ???????
            CellTopology tPhase1ChildTopo = CellTopology::TET4;
            Matrix<IndexMat> tCellToNodeGhost0 = {{21, 27, 31, 30},
                                                        {17, 18, 21, 27},
                                                        {31, 27, 34, 36},
                                                        {18, 20, 22, 27},
                                                        {36, 27, 38, 40},
                                                        {20, 19, 23, 27},
                                                        {17, 24, 19, 27},
                                                        {30, 27, 31, 44},
                                                        {31, 27, 36, 44},
                                                        {36, 27, 40, 44},
                                                        {27, 30, 40, 44},
                                                        {17, 26, 18, 27},
                                                        {18, 26, 20, 27},
                                                        {20, 26, 19, 27},
                                                        {17, 19, 26, 27},
                                                        {21, 28, 30, 31},
                                                        {21, 28, 31, 29},
                                                        {21, 29, 31, 32},
                                                        {27, 21, 31, 32},
                                                        {18, 21, 27, 32},
                                                        {28, 21, 30, 33},
                                                        {21, 27, 30, 33},
                                                        {21, 17, 27, 33},
                                                        {27, 22, 34, 36},
                                                        {34, 22, 35, 36},
                                                        {22, 35, 36, 37},
                                                        {27, 22, 36, 37},
                                                        {20, 22, 27, 37},
                                                        {31, 27, 32, 34},
                                                        {27, 18, 32, 34},
                                                        {27, 22, 18, 34},
                                                        {27, 23, 38, 40},
                                                        {38, 23, 39, 40},
                                                        {36, 27, 37, 38},
                                                        {27, 20, 37, 38},
                                                        {27, 23, 20, 38},
                                                        {23, 39, 40, 41},
                                                        {27, 23, 40, 41},
                                                        {19, 23, 27, 41},
                                                        {27, 24, 42, 43},
                                                        {30, 27, 40, 43},
                                                        {40, 27, 42, 43},
                                                        {40, 27, 41, 42},
                                                        {27, 19, 41, 42},
                                                        {27, 24, 19, 42},
                                                        {27, 30, 33, 43},
                                                        {17, 27, 33, 43},
                                                        {24, 27, 17, 43}};

            // ???????
            Matrix<IndexMat> tCellIdsGhost0 = {{ 5,  7,  9, 11, 13, 15, 19, 21, 22, 23, 24,
            		                            25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38,
												39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52,
												56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68,
												69, 70, 71, 72}};


            // ???????
            moris::mtk::MtkSetsInfo tMtkMeshSets;
            // Define side sets on the integration mesh (i.e. fixed bc, interface and ghost)

            // interface side set
            moris::mtk::MtkSideSetInfo tInterfaceSideSet;
            Matrix<IndexMat> tInterfaceElemIdandSideOrd = {{6, 2},
                                                           {8, 1},
                                                           {9, 2},
                                                           {10, 2},
                                                           {12, 1},
                                                           {13, 2},
                                                           {14, 1},
                                                           {16, 2},
                                                           {17, 1},
                                                           {18, 1},
                                                           {20, 2},
                                                           {21, 2},
                                                           {22, 2},
                                                           {23, 2},
                                                           {24, 1},
                                                           {29, 1},
                                                           {30, 1},
                                                           {32, 2},
                                                           {33, 1},
                                                           {34, 1},
                                                           {37, 2},
                                                           {41, 2},
                                                           {43, 1},
                                                           {44, 1},
                                                           {45, 1},
                                                           {48, 2},
                                                           {52, 2},
                                                           {54, 1},
                                                           {55, 1},
                                                           {56, 2},
                                                           {59, 1},
                                                           {63, 1},
                                                           {65, 2},
                                                           {66, 2},
                                                           {67, 2},
                                                           {70, 1},
                                                           {73, 1},
                                                           {76, 1},
                                                           {79, 1},
                                                           {82, 2}};

            tInterfaceSideSet.mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd;
            tInterfaceSideSet.mSideSetName        = "iside" ;
            tMtkMeshSets.add_side_set(&tInterfaceSideSet);

            // Fixed bc
            moris::mtk::MtkSideSetInfo tFixed;
            Matrix<IndexMat> tFixedElementsAndOrds = {{1, 4}};
            tFixed.mElemIdsAndSideOrds = &tFixedElementsAndOrds;
            tFixed.mSideSetName        = "fixed" ;
            tMtkMeshSets.add_side_set(&tFixed);

            // Fixed bc
            moris::mtk::MtkSideSetInfo tGhost;
            Matrix<IndexMat> tGhostCellAndOrds = {{3,5},{4,4}};
            tGhost.mElemIdsAndSideOrds = &tGhostCellAndOrds;
            tGhost.mSideSetName        = "ghost_facets" ;
            tMtkMeshSets.add_side_set(&tGhost);

            // add block sets
            // Tet Cells in Omega 0
            moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
            tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
            tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
            tOmega0BlockSetTet.mBlockSetTopo = CellTopology::TET4;
            tMtkMeshSets.add_block_set(&tOmega0BlockSetTet);

            // Hex Cells in Omega 0
            Matrix<IdMat> tOmega0HexCellIds = {{1,2,3}};
            moris::mtk::MtkBlockSetInfo tOmega0BlockSetHex;
            tOmega0BlockSetHex.mCellIdsInSet = &tOmega0HexCellIds;
            tOmega0BlockSetHex.mBlockSetName = "Omega_0_hex";
            tOmega0BlockSetHex.mBlockSetTopo = CellTopology::HEX8;
            tMtkMeshSets.add_block_set(&tOmega0BlockSetHex);

            // Cells in the ghost domain of omega 1
            moris::mtk::MtkBlockSetInfo tOmega0GhostBlockSetTet;
            tOmega0GhostBlockSetTet.mCellIdsInSet = &tCellIdsGhost0;
            tOmega0GhostBlockSetTet.mBlockSetName = "Omega_0_Ghost";
            tOmega0GhostBlockSetTet.mBlockSetTopo = CellTopology::TET4;
            tMtkMeshSets.add_block_set(&tOmega0GhostBlockSetTet);

            // Integration Cells for Ghost penalization only
            Matrix<IdMat> tGhostCellIds = {{3,4}};
            moris::mtk::MtkBlockSetInfo tCellsForGhost;
            tCellsForGhost.mCellIdsInSet = &tGhostCellIds;
            tCellsForGhost.mBlockSetName = "Ghost_Cells_0";
            tCellsForGhost.mBlockSetTopo = CellTopology::HEX8;
            tMtkMeshSets.add_block_set(&tCellsForGhost);


            // Mesh data input structure
            moris::mtk::MtkMeshData tIntegMeshData(3);

            moris::uint tSpatialDim   = 3;
            moris::uint tNumElemTypes = 3;
            Matrix<IdMat> tNodeOwner(1,tIntegNodeCoordinates.n_rows(),moris::par_rank());
            tIntegMeshData.ElemConn(0) = &tElemsAsCellToNodes;
            tIntegMeshData.ElemConn(1) = &tCellToNodePhase0;
            tIntegMeshData.ElemConn(2) = &tCellToNodeGhost0;
            tIntegMeshData.LocaltoGlobalElemMap(0) = (&tElemsAsCellIds);
            tIntegMeshData.LocaltoGlobalElemMap(1) = (&tCellIdsPhase0);
            tIntegMeshData.LocaltoGlobalElemMap(2) = (&tCellIdsGhost0);
            tIntegMeshData.CreateAllEdgesAndFaces  = false;
            tIntegMeshData.Verbose                 = false;
            tIntegMeshData.SpatialDim              = &tSpatialDim;
            tIntegMeshData.NodeCoords              = &tIntegNodeCoordinates;
            tIntegMeshData.NodeProcOwner           = &tNodeOwner;
            tIntegMeshData.LocaltoGlobalNodeMap    = &tIntegLocalToGlobalNodeMap;
            tIntegMeshData.SetsInfo                = &tMtkMeshSets;
            tIntegMeshData.MarkNoBlockForIO        = false;

            moris::mtk::Mesh* tIntegMesh = moris::mtk::create_mesh( MeshType::STK, tIntegMeshData );

            // write to a file
            std::string tIntegFileOutput = "./DiffusionBar_Integration.exo";
            tIntegMesh->create_output_mesh(tIntegFileOutput,"tIntegFileOutput");

            // Define the relationship between integration mesh and interpolation mesh
            moris::moris_index tInterpCellIndex = 3;
            Matrix<IndexMat> tVertexIDs = {{13, 14, 16, 15, 17, 18, 20, 19, 21, 22, 23, 24, 25, 26, 27, 28,
            		                        29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};

            // element level parameter
            Matrix<DDRMat> tLocalCoordinatesWrtElem4 = {{-1, -1, -1},
                                                        { 1, -1, -1},
                                                        { 1, 1, -1},
                                                        {-1, 1, -1},
                                                        {-1, -1, 1},
                                                        { 1, -1, 1},
                                                        { 1, 1, 1},
                                                        {-1, 1, 1},
                                                        { 0, -1, 0},
                                                        { 1, 0, 0},
                                                        { 0, 1, 0},
                                                        { -1, 0, 0},
                                                        {0, 0, -1},
                                                        {0, 0, 1},
                                                        {0, 0, 0},
                                                        {-0.8, -1, -0.8},
                                                        {0.8, -1, -0.8},
                                                        {-0.8, -0.8, -0.8},
                                                        {0.8, -0.8, -0.8},
                                                        {1, -1, -0.8},
                                                        {-1, -1, -0.8},
                                                        {1, -0.8, -0.8},
                                                        {1, 0.8, -0.8},
                                                        {0.8, 0.8, -0.8},
                                                        {1, 1, -0.8},
                                                        {0.8, 1, -0.8},
                                                        {-0.8, 1, -0.8},
                                                        {-0.8, 0.8, -0.8},
                                                        {-1, 1, -0.8},
                                                        {-1, 0.8, -0.8},
                                                        {-1, -0.8, -0.8},
                                                        {0, 0, -0.8}};



            //delete tIntegMesh Data;
            //delete tIntegMeshData;


            //1) Create the fem nodes on Interpolation Mesh --------------------------------
            //std::cout<<" Create the fem nodes "<<std::endl;
            //------------------------------------------------------------------------------

            // number of mesh nodes
            uint tInterpolNumOfNodes = tInterpolMesh->get_num_nodes();
            //std::cout<<"Number of Nodes:"<< tNumOfNodes <<std::endl;

            //create a cell of fem nodes
            moris::Cell< Node_Base* > tInterpolNodes( tInterpolNumOfNodes, nullptr );

            // loop over the mesh nodes
            for( uint k = 0; k < tInterpolNumOfNodes; k++ )
            {
                // create a fem node for each mesh node
                tInterpolNodes( k ) = new fem::Node( & tInterpolMesh->get_mtk_vertex( k ) );
            }

            //2) Create the IWGs -----------------------------------------------------------
            std::cout<<" 2) Create the IWGs "<<std::endl;
            //------------------------------------------------------------------------------

            // input a cell of IWG types to be created
            Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::SPATIALDIFF_BULK,
                                                   fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                   fem::IWG_Type::SPATIALDIFF_NEUMANN,
                                                   fem::IWG_Type::SPATIALDIFF_GHOST_LEFT,
                                                   fem::IWG_Type::SPATIALDIFF_GHOST_RIGHT};

            // number of IWGs to be created
            uint tNumOfIWGs = tIWGTypeList.size();

            std::cout<<"Length of IWG list:"<<tIWGTypeList.size()<<std::endl;

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


            //3) Create element blocks for both Meshes. ------------------------------------
            std::cout<<" 3) Create element blocks "<<std::endl;
            //------------------------------------------------------------------------------
            // get the number of elements on proc from the interpolation mesh
            luint tNumOfInterpolElements = tInterpolMesh->get_num_elems();

            // get the number of elements on proc from the integration mesh
            uint tNumOfIntegElements = tIntegMesh->get_num_elems();

            // create equation objects for Interpolation Mesh
            Cell< MSI::Equation_Object* > tInterpolElements;
            tInterpolElements.reserve( tNumOfInterpolElements + 1 + 2 + 1 );
                                                           // + 1 Dirichlet + 2 Ghost + 1 Neumann elements

            // create equation objects for Integration Mesh
            Cell< MSI::Equation_Object* > tIntegElements;
            tIntegElements.reserve( tNumOfIntegElements + 1 + 2 + 1 );

            // get the block names from the mesh
            moris::Cell<std::string> tInterpolBlockSetsNames = tInterpolMesh->get_set_names( EntityRank::ELEMENT);

            // Cell containing the block mesh cell ( a cell of mesh cells )
            moris::Cell<mtk::Cell*> tInterpolBlockSetElement( tInterpolMesh->get_set_entity_loc_inds( EntityRank::ELEMENT,
            		                                          tInterpolBlockSetsNames( 0 ) ).numel(), nullptr );

            // loop on the blocks
            for( luint Ik=0; Ik < tInterpolBlockSetsNames.size(); ++Ik )
            {
                Matrix< IndexMat > tInterpolBlockSetElementInd = tInterpolMesh->get_set_entity_loc_inds( EntityRank::ELEMENT,
                		                                                                                 tInterpolBlockSetsNames( Ik ) );

                // loop on the elements in a block
                for( luint k=0; k < tInterpolBlockSetElementInd.numel(); ++k )
                {
                    // Cell containing the block mesh cell with mesh cells
                    tInterpolBlockSetElement( k ) = & tInterpolMesh->get_mtk_cell( k );
                }
            }


            // create a fem element block
            moris::Cell< fem::Element_Block * > tInterpolElementBlocks( 5, nullptr );
            moris::Cell< IWG* > tIWGBulk = { tIWGs( 0 ) };
            tInterpolElementBlocks( 0 ) = new fem::Element_Block( tInterpolBlockSetElement,
                                                                  fem::Element_Type::BULK,
                                                                  tIWGBulk,
                                                                  tInterpolNodes );

            // put the equation object of block 0 in the global list of equation objects
            tInterpolElements.append( tInterpolElementBlocks( 0 )->get_equation_object_list() );



            //4) Create Dirichlet side blocks ---------------------------------------------
            std::cout<<" 4) Create Dirichlet side block "<<std::endl;
            //-----------------------------------------------------------------------------
            // elements included in the Dirichlet sideset
            moris::Cell<mtk::Cell*> tSideSetDirichletElement( 1, nullptr );
            tSideSetDirichletElement( 0 ) = & tInterpolMesh->get_mtk_cell( 0 );

            // create a fem element sideblock
            moris::Cell< IWG* > tIWGDirichlet = { tIWGs( 1 ) };
            tInterpolElementBlocks( 1 ) = new fem::Element_Block( tSideSetDirichletElement,
                                                                  fem::Element_Type::SIDESET,
                                                                  tIWGDirichlet,
                                                                  tInterpolNodes );

            // put the equation object of block 1 in the global list of equation objects
            tInterpolElements.append( tInterpolElementBlocks( 1 )->get_equation_object_list() );

            //5) Create Neumann side blocks -------------------------------------------------------
            std::cout<<" 5) Create Neumann side block "<<std::endl;
            //-----------------------------------------------------------------------------
            // elements included in the Neumann sideset
            moris::Cell<mtk::Cell*> tSideSetNeumannElement( 1, nullptr );
            tSideSetNeumannElement( 0 ) = & tInterpolMesh->get_mtk_cell( 3 );

            // create a fem element sideblock
            moris::Cell< IWG* > tIWGNeumann = { tIWGs( 2 ) };
            tInterpolElementBlocks( 2 ) = new fem::Element_Block( tSideSetNeumannElement,
                                                                  fem::Element_Type::SIDESET,
                                                                  tIWGNeumann,
                                                                  tInterpolNodes );

            // put the equation object of block 2 in the global list of equation objects
            tInterpolElements.append( tInterpolElementBlocks( 2 )->get_equation_object_list() );

            //6) Create left Ghost side blocks -------------------------------------------------------
            std::cout<<" 6) Create left Ghost side block "<<std::endl;
            //-----------------------------------------------------------------------------
            // elements included in the left Ghost sideset
            moris::Cell<mtk::Cell*> tSideSetLeftGhostElement( 1, nullptr );
            tSideSetLeftGhostElement( 0 ) = & tInterpolMesh->get_mtk_cell( 2 );

            // create a fem element sideblock
            moris::Cell< IWG* > tIWGLeftGhost = { tIWGs( 3 ) };
            tInterpolElementBlocks( 3 ) = new fem::Element_Block( tSideSetLeftGhostElement,
                                                                  fem::Element_Type::SIDESET,
                                                                  tIWGLeftGhost,
                                                                  tInterpolNodes );

            // put the equation object of block 1 in the global list of equation objects
            tInterpolElements.append( tInterpolElementBlocks( 3 )->get_equation_object_list() );

            //7) Create right Ghost side blocks -------------------------------------------
            std::cout<<" 7) Create right Ghost side block "<<std::endl;
            //-----------------------------------------------------------------------------
            // elements included in the right Ghost sideset
            moris::Cell<mtk::Cell*> tSideSetRightGhostElement( 1, nullptr );
            tSideSetRightGhostElement( 0 ) = & tInterpolMesh->get_mtk_cell( 3 );

            // create a fem element sideblock
            moris::Cell< IWG* > tIWGRightGhost = { tIWGs( 4 ) };
            tInterpolElementBlocks( 4 ) = new fem::Element_Block( tSideSetRightGhostElement,
                                                                  fem::Element_Type::SIDESET,
                                                                  tIWGRightGhost,
                                                                  tInterpolNodes );

            // put the equation object of block 2 in the global list of equation objects
            tInterpolElements.append( tInterpolElementBlocks( 4 )->get_equation_object_list() );

            //8) Impose boundary conditions -----------------------------------------------
            std::cout<<" 8) Impose boundary conditions "<<std::endl;
            //-----------------------------------------------------------------------------
            // nodal weak bc for Dirichlet
            real tTempValue = 0.0;
            Cell< moris_index > tDirichletElements = { 0 };
            for( uint iDirichlet = 0; iDirichlet < 1; iDirichlet++ )
            {
                 // set the list of face ordinals
                 tInterpolElements( tNumOfInterpolElements + iDirichlet )->set_list_of_side_ordinals( {{ 4 }} );

                 //get the nodal weak bcs of the element
                 Matrix< DDRMat > & tNodalWeakBCs = tInterpolElements( tNumOfInterpolElements + iDirichlet )->get_weak_bcs();

                 // get the element number of nodes
                 uint tNumberOfNodes = tInterpolElements( tNumOfInterpolElements + iDirichlet )->get_num_nodes();

                 // set size of the element nodal weak bc
                 tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                 // loop over the element nodes
                 Matrix< IndexMat > tNodeIndices = tInterpolMesh->get_mtk_cell( tDirichletElements( iDirichlet ) ).get_vertex_inds();

                 for( uint l = 0; l < tNumberOfNodes; l++ )
                 {
                     // copy weak bc into element
                     tNodalWeakBCs( l ) = tTempValue;
                 }
            }

            // nodal weak bc for Neumann
            real tFluxValue = 1.0;
            Cell< moris_index > tNeumannElements = { 3 };
            for( uint iNeumann = 0; iNeumann < 1; iNeumann++ )
            {
                 // set the list of face ordinals
                 tInterpolElements( tNumOfInterpolElements + tDirichletElements.size() + iNeumann )->set_list_of_side_ordinals( {{ 5 }} );

                 //get the nodal weak bcs of the element
                 Matrix< DDRMat > & tNodalWeakBCs = tInterpolElements( tNumOfInterpolElements + tNeumannElements.size() + iNeumann )->get_weak_bcs();

                 // get the element number of nodes
                 uint tNumberOfNodes = tInterpolElements( tNumOfInterpolElements + tNeumannElements.size() + iNeumann )->get_num_nodes();

                 // set size of the element nodal weak bc
                 tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                 // loop over the element nodes
                 Matrix< IndexMat > tNodeIndices = tInterpolMesh->get_mtk_cell( tNeumannElements( iNeumann ) ).get_vertex_inds();

                 for( uint l = 0; l < tNumberOfNodes; l++ )
                 {
                     // copy weak bc into element
                     tNodalWeakBCs( l ) = tFluxValue;
                 }
            }

            //  left Ghost
                 // set the list of face ordinals
                 tInterpolElements( 6 )->set_list_of_side_ordinals( {{ 5 }} );

            //  right Ghost
                 // set the list of face ordinals
                 tInterpolElements( 7 )->set_list_of_side_ordinals( {{ 4 }} );


            //9) Create the model solver interface -----------------------------------------
            std::cout<<" 9) Create the model solver interface "<<std::endl;
            //------------------------------------------------------------------------------

            //FIXME force the communication table
            Matrix< IdMat > tCommunicationTable( 1, 1, 0 );

            // FIXME: get map from mesh
            uint tDofOrder = 1;
            map< moris_id, moris_index > tCoefficientsMap;
            //tMesh->get_adof_map( tDofOrder, tCoefficientsMap );

            uint tNumInterpolCoeff = 20;
            //= tMesh->get_num_coeffs( 1 )

            moris::MSI::Model_Solver_Interface* tModelSolverInterface
                = new moris::MSI::Model_Solver_Interface( tInterpolElements,
                                                          tCommunicationTable,
                                                          tCoefficientsMap,
                                                          tNumInterpolCoeff,
                                                          tInterpolMesh );

            tModelSolverInterface->set_param( "TEMP" )  = (sint)tDofOrder;

            tInterpolElementBlocks( 0 )->finalize( tModelSolverInterface );
            tInterpolElementBlocks( 1 )->finalize( tModelSolverInterface );
            tInterpolElementBlocks( 2 )->finalize( tModelSolverInterface );
            tInterpolElementBlocks( 3 )->finalize( tModelSolverInterface );
            tInterpolElementBlocks( 4 )->finalize( tModelSolverInterface );

            tModelSolverInterface->finalize();
            // calculate AdofMap
            Matrix< DDUMat > tAdofMap = tModelSolverInterface->get_dof_manager()->get_adof_ind_map();


            //10) Create solver interface ---------------------------------------------------
            std::cout<<" 10) Create solver interface "<<std::endl;
            //------------------------------------------------------------------------------

            MSI::MSI_Solver_Interface * tSolverInterface
                = new moris::MSI::MSI_Solver_Interface( tModelSolverInterface );

            // 11) Create Nonlinear Problem -------------------------------------------------
            std::cout<<" 11) Create Nonlinear Problem "<<std::endl;
            //------------------------------------------------------------------------------

            NLA::Nonlinear_Problem* tNonlinearProblem
                = new NLA::Nonlinear_Problem( tSolverInterface );

            // 12) Create Solvers and solver manager ----------------------------------------
            std::cout<<" 12) Create Solvers and solver manager "<<std::endl;
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


            // 13) Solve --------------------------------------------------------------------
            std::cout<<" 13) Solve "<<std::endl;
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

            // 14) postprocessing------------------------------------------------------------
            std::cout<<" 14) Postprocessing "<<std::endl;
            //------------------------------------------------------------------------------

            // dof type list for the solution to write on the mesh
            moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

            // create a matrix to be filled  with the solution
            Matrix< DDRMat > tTempSolutionField( tInterpolNumOfNodes, 1 );

            // loop over the nodes
            for( uint i = 0; i < tInterpolNumOfNodes; i++ )
            {
                // get a list of elements connected to the ith node
                Matrix<IndexMat> tConnectedElements =
                    tInterpolMesh->get_entity_connected_to_entity_loc_inds( static_cast< moris_index >( i ),
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
                        = tInterpolElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i,
                                                                                              tDofTypeList);
                    // add up the contribution of each element to the node value
                    tNodeVal = tNodeVal + tElemVal;
                }
                // fill the solution matrix with the node value
                tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
            }

            // checking the solution--------------------------------------------------------
            //------------------------------------------------------------------------------
            // Expected solution
            Matrix< DDRMat > tExpectedSolution = {{ 0.0, 1.0, 2.0, 3.0, 4.0}};
//            		                                0.0, 0.0, 0.0, 0.0,
//            		                                1.0, 1.0, 1.0, 1.0,
//													2.0, 2.0, 2.0, 2.0,
//													3.0, 3.0, 3.0, 3.0,
//													4.0, 4.0, 4.0, 4.0}};
            Matrix< DDUMat > tLookUp = {{ 0, 4, 8, 12, 16}};

            print(tTempSolutionField,"Solution:");

            // define an epsilon for the deviation of the numerical solution
            double tEpsilon = 1E-04;
            double tEpsilonNodal = 1E-02;

            // define a bool for solution check
            bool tCheckNodalSolution = true;
            bool tCheckNormSolution  = true;

            // initialize L2-norms
            real tExpectedL2Norm = 0.0;
            real tSolutionL2Norm = 0.0;
            real tSolutionL2Error = 0.0;

            // loop over the nodes to compute norms and check node solutions
            for ( uint i = 0; i < 4; i++ )
            {
            	tExpectedL2Norm = tExpectedL2Norm + std::pow( (tExpectedSolution(i) + tExpectedSolution(i+1))/2 , 2 );
            	tSolutionL2Norm = tSolutionL2Norm + std::pow( (tTempSolutionField(tLookUp(i)) + tTempSolutionField(tLookUp(i+1)))/2 , 2 );
            	tSolutionL2Error = tSolutionL2Error + std::pow( (tTempSolutionField(tLookUp(i)) + tTempSolutionField(tLookUp(i+1)))/2
            		                                            -(tExpectedSolution(i) + tExpectedSolution(i+1))/2 , 2 );

            	// check nodal solution
                tCheckNodalSolution = tCheckNodalSolution
                                      && ( std::abs( tTempSolutionField( tLookUp(i) ) - tExpectedSolution( i ) ) < tEpsilonNodal );
            }

                tCheckNodalSolution = tCheckNodalSolution
                                      && ( std::abs( tTempSolutionField( tLookUp(4) ) - tExpectedSolution( 4 ) ) < tEpsilonNodal );

            //L2-Norm, relative L2-error should be smaller than Epsilon
            tSolutionL2Norm = std::sqrt(tSolutionL2Norm);
            tExpectedL2Norm = std::sqrt(tExpectedL2Norm);
            tSolutionL2Error = std::sqrt(tSolutionL2Error);
            real tNormalizedL2Error = tSolutionL2Error / tExpectedL2Norm;

            std::cout<<"Expected L2-Norm: "<<tExpectedL2Norm<<std::endl;
            std::cout<<"Solution L2-Norm: "<<tSolutionL2Norm<<std::endl;
            std::cout<<"Normalized L2-Error: "<<tNormalizedL2Error<<std::endl;

            // check L2 Error
            tCheckNormSolution = tCheckNormSolution
                                  && ( std::abs( tNormalizedL2Error) < tEpsilon );

            // check if test has passed
            CHECK( tCheckNormSolution );
            CHECK( tCheckNodalSolution );

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
            tInterpolMeshData.FieldsInfo = &tFieldsInfo;

            // create the mesh
            mtk::Mesh* tMeshForOutput = create_mesh( MeshType::STK, tInterpolMeshData );

            // add field to the mesh
            tMeshForOutput->add_mesh_field_real_scalar_data_loc_inds( tTempFieldName,
                                                                      EntityRank::NODE,
                                                                      tTempSolutionField );

            // create output mesh
            std::string tOutputFile = "./int_ElemDiff_test.exo";
            tMeshForOutput->create_output_mesh( tOutputFile );


            // 10) Clean up -----------------------------------------------------------------
            //std::cout<<" 10) Clean up "<<std::endl;
            //------------------------------------------------------------------------------
            delete tInterpolMesh;
            delete tIntegMesh;
            delete tMeshForOutput;

            for( uint i = 0; i < tInterpolNumOfNodes; i++ )
            {
                delete tInterpolNodes( i );
            }

            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                delete tIWGs( i );
            }

            for( uint i = 0; i < tNumOfInterpolElements; i++ )
            {
                delete tInterpolElements( i );
            }

            delete tModelSolverInterface;
            delete tSolverInterface;
            delete tNonlinearProblem;
            delete tLinSolver;
            delete tNonlinearSolver;


            //}/*namespace mtk*/
            }/*TEST_CASE*/




    }/* namespace fem */
}/* namespace moris */
