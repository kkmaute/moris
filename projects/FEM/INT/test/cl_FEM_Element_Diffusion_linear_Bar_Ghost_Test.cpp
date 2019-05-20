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
#include "cl_FEM_Set.hpp"                      //FEM/INT/src

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

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

namespace moris
{
    namespace fem
    {

// test in this .cpp file:
// "ElemDiff_lin_Bar" - 5x1x1 'Bar' with linear elements and Ghost applied between elements 4 and 5
    /*
     *  Created on: Apr 30, 2019
     *      Author: wunsch
     */


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
            moris::Cell<mtk::Cell const*> tSideSetCells;
            Matrix<IndexMat> tSideOrdinalsOfCell;
            aMesh->get_sideset_cells_and_ords(aSideSetName,tSideSetCells,tSideOrdinalsOfCell);

            Matrix<IndexMat> tLinearHexMap = get_linear_hex_side_vertex_map();

            Matrix<IndexMat> tNonUniqueVerticesOnSideSet(1,tLinearHexMap.n_cols()*tSideSetCells.size());

            uint tCount = 0;
            // iterate through cells in side set
            for(moris::uint i =0; i < tSideSetCells.size(); i++)
            {
            	// get pointer to current cell/element
            	mtk::Cell const * tCell = tSideSetCells(i);

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
//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//TEST//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

        TEST_CASE( "Element_Diffusion_linear_Bar_Ghost", "[moris],[fem],[ElemDiff_lin_Bar]" )
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
                std::cout<<" 3) Create element blocks / sets "<<std::endl;
                //------------------------------------------------------------------------------
                // get the number of elements on proc from the mesh
                luint tNumOfElements = tMesh->get_num_elems();

                // create equation objects
                Cell< MSI::Equation_Object* > tElements;
                tElements.reserve( tNumOfElements + 1 + 2 + 1 );

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
                moris::Cell< MSI::Equation_Set * > tElementBlocks( 5, nullptr );
                moris::Cell< IWG* > tIWGBulk = { tIWGs( 0 ) };
                tElementBlocks( 0 ) = new fem::Set( tBlockSetElement,
                                                              fem::Element_Type::BULK,
                                                              tIWGBulk,
                                                              tNodes );

                // put the equation object of block 0 in the global list of equation objects
                tElements.append( tElementBlocks( 0 )->get_equation_object_list() );

                //4) Create Dirichlet side blocks -------------------------------------------------------
                std::cout<<" 4) Create Dirichlet side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Dirichlet sideset
                moris::Cell<mtk::Cell const *> tSideSetDirichletElement( 1, nullptr );
                tSideSetDirichletElement( 0 ) = & tMesh->get_mtk_cell( 0 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGDirichlet = { tIWGs( 1 ) };
                tElementBlocks( 1 ) = new fem::Set( tSideSetDirichletElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGDirichlet,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 1 )->get_equation_object_list() );

                //5) Create Neumann side blocks -------------------------------------------------------
                std::cout<<" 5) Create Neumann side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the Neumann sideset
                moris::Cell<mtk::Cell const *> tSideSetNeumannElement( 1, nullptr );
                tSideSetNeumannElement( 0 ) = & tMesh->get_mtk_cell( 4 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGNeumann = { tIWGs( 2 ) };
                tElementBlocks( 2 ) = new fem::Set( tSideSetNeumannElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGNeumann,
                                                              tNodes );

                // put the equation object of block 2 in the global list of equation objects
                tElements.append( tElementBlocks( 2 )->get_equation_object_list() );

                //6) Create left Ghost side blocks -------------------------------------------------------
                std::cout<<" 6) Create left Ghost side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the left Ghost sideset
                moris::Cell<mtk::Cell const *> tSideSetLeftGhostElement( 1, nullptr );
                tSideSetLeftGhostElement( 0 ) = & tMesh->get_mtk_cell( 3 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGLeftGhost = { tIWGs( 3 ) };
                tElementBlocks( 3 ) = new fem::Set( tSideSetLeftGhostElement,
                                                              fem::Element_Type::SIDESET,
                                                              tIWGLeftGhost,
                                                              tNodes );

                // put the equation object of block 1 in the global list of equation objects
                tElements.append( tElementBlocks( 3 )->get_equation_object_list() );

                //7) Create right Ghost side blocks -------------------------------------------------------
                std::cout<<" 7) Create right Ghost side block "<<std::endl;
                //-----------------------------------------------------------------------------
                // elements included in the right Ghost sideset
                moris::Cell<mtk::Cell const *> tSideSetRightGhostElement( 1, nullptr );
                tSideSetRightGhostElement( 0 ) = & tMesh->get_mtk_cell( 4 );

                // create a fem element sideblock
                moris::Cell< IWG* > tIWGRightGhost = { tIWGs( 4 ) };
                tElementBlocks( 4 ) = new fem::Set( tSideSetRightGhostElement,
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
                    = new moris::MSI::Model_Solver_Interface( tElementBlocks,
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


    }/* namespace fem */
}/* namespace moris */
