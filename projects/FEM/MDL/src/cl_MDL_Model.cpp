
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"

#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"               //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "../../INT/src/cl_FEM_Element_Bulk.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_Element_Block.hpp"

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

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_print.hpp" // for check

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        Model::Model(       mtk::Mesh *                   aMesh,
                      const uint                          aBSplineOrder,
                            Cell< Cell< fem::IWG_Type > > aIWGTypeList,
                            Cell< moris_index >           aSidesetList,
                            Cell< fem::BC_Type >          aSidesetBCTypeList ) : mMesh( aMesh )


//    Model::Model(       mtk::Mesh *                   aMesh,
//                  const uint                          aBSplineOrder,
//                        Cell< Cell< fem::IWG_Type > > aIWGTypeList ) : mMesh( aMesh )
        {
            // start timer
            tic tTimer1;

            mDofOrder = aBSplineOrder;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // ask mesh about number of nodes on proc
            luint tNumberOfNodes = aMesh->get_num_nodes();

            std::cout<<tNumberOfNodes<<" Number of nodes"<<std::endl;

            // create node objects
            mNodes.resize(  tNumberOfNodes, nullptr );

            for( luint k = 0; k<tNumberOfNodes; ++k )
            {
                mNodes( k ) = new fem::Node( &aMesh->get_mtk_vertex( k ) );
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created %u FEM nodes in %5.3f seconds.\n\n",
                        ( unsigned int ) tNumberOfNodes,
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1.5: create IWGs
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // number of IWGs to be created
            uint tNumOfIWGs = aIWGTypeList.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs for the problem considered
            mIWGs.resize( tNumOfIWGs );

            // loop over the IWG types
            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                mIWGs( i ).resize( aIWGTypeList( i ).size(), nullptr );

                for( uint Ki = 0; Ki < aIWGTypeList( i ).size(); Ki++)
                {
                    // create an IWG with the factory for the ith IWG type
                    mIWGs( i )( Ki ) = tIWGFactory.create_IWGs( aIWGTypeList( i )( Ki ) );
                }
            }

//            if ( tNumOfIWGs > 1 )
//            {
//                mIWGs1.resize( 2, nullptr );
//                mIWGs1( 0 ) = mIWGs( 1 )( 0 );
//                mIWGs1( 1 ) = mIWGs( 2 )( 0 );
//
//            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer2;

            // a factory to create the elements
            fem::Element_Factory tElementFactory;

            // get the number of element to create
            luint tNumberOfEquationObjects = aMesh->get_num_elems()
                                           + aMesh->get_sidesets_num_faces( aSidesetList );

            uint tNumberElementBlocks = 1 + aSidesetList.size();

            // create equation objects
            mElements.reserve( tNumberOfEquationObjects );

            mElementBlocks.resize( tNumberElementBlocks, nullptr );

            //  Create Blockset Elements ----------------------------------------
            std::cout<<" Create Blockset Elements "<<std::endl;
            //------------------------------------------------------------------------------
            // ask mesh about number of elements on proc
            moris::Cell<std::string> tBlockSetsNames = aMesh->get_set_names( EntityRank::ELEMENT);

            moris::Cell<mtk::Cell*> tBlockSetElement( aMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( 0 ) ).numel(), nullptr );

            moris::uint tEquationObjectCounter = 0;
            moris::uint tElementBlockCounter = 0;
            for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
            {
                Matrix< IndexMat > tBlockSetElementInd
                    = aMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( Ik ) );

                for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                {
                    tEquationObjectCounter++;

                    tBlockSetElement( k ) = & aMesh->get_mtk_cell( k );
                }

            }
            mElementBlocks( tElementBlockCounter ) = new fem::Element_Block( tBlockSetElement, fem::Element_Type::BULK, mIWGs( 0 ), mNodes );

            mElementBlocks( tElementBlockCounter )->finalize();

            mElements.append( mElementBlocks( tElementBlockCounter++ )->get_equation_object_list() );

            //  Create Blockset Elements ----------------------------------------
            std::cout<<" Create Sideset Elements "<<std::endl;
            //------------------------------------------------------------------------------

            moris::Cell<std::string> tSideSetsNames = aMesh->get_set_names( EntityRank::FACE );

            for( luint Ik = 0; Ik < aSidesetList.size(); ++Ik )
            {
                // get the treated sideset name
                std::string tTreatedSidesetName = tSideSetsNames( aSidesetList( Ik ) );

                // create a cell of sideset element
                moris::Cell< mtk::Cell * > tSideSetElement;

                // create a list of sideset ordinal
                Matrix< IndexMat >  aSidesetOrdinals;

                // get the treated sideset elements and ordinals
                aMesh->get_sideset_cells_and_ords( tTreatedSidesetName, tSideSetElement, aSidesetOrdinals );

                mElementBlocks( tElementBlockCounter ) = new fem::Element_Block( tSideSetElement, fem::Element_Type::SIDESET, mIWGs( Ik + 1 ), mNodes );

                mElementBlocks( tElementBlockCounter )->finalize();

                mElements.append( mElementBlocks( tElementBlockCounter++ )->get_equation_object_list() );

                for( luint k = 0; k < tSideSetElement.size(); ++k )
                {
                    mElements( tEquationObjectCounter )->set_list_of_side_ordinals( {{aSidesetOrdinals( k )}} ); //FIXME

                    // get the nodal weak bcs of the element
                    Matrix< DDRMat > & tNodalWeakBCs = mElements( tEquationObjectCounter )->get_weak_bcs();

                    // get the element number of nodes
                    uint tNumberOfNodes = mElements( tEquationObjectCounter++ )->get_num_nodes();

                    // set size of the element nodal weak bc
                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                    // loop over the element nodes
                    Matrix< IndexMat > tNodeIndices = tSideSetElement( k )->get_vertex_inds();

                    //--------------------------------------------------------------------------------------------
                    for( uint l = 0; l < tNumberOfNodes; l++ )
                    {
                        if ( aSidesetBCTypeList( Ik ) == fem::BC_Type::DIRICHLET )
                        {
                        // copy weak bc into element
                        tNodalWeakBCs( l ) = 5.0;
                        }
                        else if ( aSidesetBCTypeList( Ik ) == fem::BC_Type::NEUMANN )
                        {
                            // copy weak bc into element
                            tNodalWeakBCs( l ) = 20.0;
                        }
                    }
                }
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created %u FEM elements in %5.3f seconds.\n\n",
                        ( unsigned int ) tNumberOfEquationObjects,
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create Model Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer3;

            //--------------------------FIXME------------------------------------
            // This part should not be needed anymore when MTK has all the functionalities
            Matrix< IdMat > tCommTable;
            moris::map< moris::moris_id, moris::moris_index > tIdToIndMap;
            moris::uint tMaxNumAdofs;

            if ( mMesh->get_mesh_type() == MeshType::HMR )
            {
                if ( mDofOrder == 0 )
                {
                    mDofOrder  = this->get_lagrange_order_from_mesh();
                }

                // get map from mesh
                mMesh->get_adof_map( mDofOrder, mCoefficientsMap );

                tCommTable   = aMesh->get_communication_table();
                tIdToIndMap  = mCoefficientsMap;
                tMaxNumAdofs = aMesh->get_num_coeffs( mDofOrder );
            }
            else
            {
                tCommTable.set_size( 1, 1, 0 );
                tMaxNumAdofs = 1000000;
            }
            //--------------------------END FIXME--------------------------------

            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( mElements,
                                                                            tCommTable,
                                                                            tIdToIndMap,
                                                                            tMaxNumAdofs,
                                                                            mMesh );

            if ( mMesh->get_mesh_type() == MeshType::HMR )
            {
                mModelSolverInterface->set_param("L2")= (sint)mDofOrder;
                mModelSolverInterface->set_param("TEMP")= (sint)mDofOrder;
            }

            mModelSolverInterface->finalize();

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 5: create Nonlinear Problem
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mNonlinearProblem =  new NLA::Nonlinear_Problem( mSolverInterface );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 6: create Solvers and solver manager
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // create factory for nonlinear solver
            NLA::Nonlinear_Solver_Factory tNonlinFactory;

            // create nonlinear solver
            mNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            // create factory for linear solver
            dla::Solver_Factory  tSolFactory;

            // create linear solver
            mLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

            // set default parameters for linear solver
            mLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            mLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

            mLinearSolverAlgorithm->set_param("AZ_orthog") = AZ_modified;

            mLinearSolverAlgorithm->set_param("AZ_graph_fill") = 2;

            //mLinearSolverAlgorithm->set_param("AZ_keep_info") = 1;
//            mLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;

            // create solver manager
            mLinSolver = new dla::Linear_Solver();
            mNonlinearSolver = new NLA::Nonlinear_Solver();

            // set manager and settings
            mNonlinearSolverAlgorithm->set_linear_solver( mLinSolver );

            // set first solver
            mLinSolver->set_linear_algorithm( 0, mLinearSolverAlgorithm );

            mNonlinearSolver->set_nonlinear_algorithm( mNonlinearSolverAlgorithm, 0 );

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created Model-Solver Interface in %5.3f seconds.\n\n",
                        ( double ) tElapsedTime / 1000 );
            }

        }

//------------------------------------------------------------------------------

        Model::~Model()
        {

            // delete manager
            delete mLinSolver;

            // delete problem
            delete mNonlinearProblem;

            // delete NonLinSolverManager
            delete mNonlinearSolver;

            // delete SI
            delete mSolverInterface;

            // delete MSI
            delete mModelSolverInterface;

            // delete IWGs
            for( auto tIWG : mIWGs )
            {
                tIWG.clear();
            }

            // delete elements
            for( auto tElement : mElements )
            {
                delete tElement;
            }

            // delete nodes
            for( auto tNode : mNodes )
            {
                delete tNode;
            }
        }

//------------------------------------------------------------------------------

        void
        Model::set_weak_bcs( const Matrix<DDRMat> & aWeakBCs )
        {
            // set weak BCs
            for( auto tElement : mElements )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k < tNumberOfNodes; ++k )
                {
                    // copy weakbc into element
                    tNodalWeakBCs( k ) = aWeakBCs( tElement->get_node_index( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Model::set_weak_bcs_from_nodal_field( moris_index aFieldIndex )
        {
            for( auto tElement : mElements )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // copy weakbc into element
                    tNodalWeakBCs( k ) = mMesh->get_value_of_scalar_field( aFieldIndex,
                                                                           EntityRank::NODE,
                                                                           tElement->get_node_index( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Model::solve( Matrix<DDRMat> & aSolution )
        {

            // call solver
            mNonlinearSolver->solve( mNonlinearProblem );

            // temporary array for solver
            Matrix< DDRMat > tSolution;
            mNonlinearSolverAlgorithm->get_full_solution( tSolution );

            // get length of array
            uint tLength = tSolution.length();

            // make sure that length of vector is correct
//            MORIS_ERROR( tLength <= (uint)  mMesh->get_num_coeffs( mDofOrder ),
//                    "Number of ADOFs does not match" );

            // rearrange data into output
            aSolution.set_size( tLength, 1 );

            for( uint k=0; k<tLength; ++k )
            {
                aSolution( k ) = tSolution( mAdofMap( k ) );
            }
        }

//------------------------------------------------------------------------------

        real
        Model::compute_integration_error(
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
        {
            real aError = 0.0;
            for( auto tElement : mElements )
            {
                aError += tElement->compute_integration_error( aFunction );
            }
            return aError;
        }

//------------------------------------------------------------------------------

        uint
        Model::get_lagrange_order_from_mesh()
        {

            // set order of this model according to Lagrange order
            // of first element on mesh
           return mtk::interpolation_order_to_uint(
            mMesh->get_mtk_cell( 0 ).get_interpolation_order() );
        }

//------------------------------------------------------------------------------

        void
        Model::set_dof_order( const uint aOrder )
        {
            mDofOrder = aOrder;
//            MORIS_ASSERT( aOrder == mDofOrder,
//                    "Model: the functionality to change the order of the model has nor been implemented yet" );
        }

//------------------------------------------------------------------------------
        real
        Model::compute_element_average( const uint aElementIndex )
        {
            return mElements( aElementIndex )->compute_element_average_of_scalar_field();
        }

//------------------------------------------------------------------------------

        void
        Model::output_solution( const std::string & aFilePath )
        {
            if ( mMesh->get_mesh_type() == MeshType::HMR )
            {
                mSolHMR.set_size(mMesh->get_num_nodes(),1,-1.0);

                moris::Cell<std::string> tBlockSetsNames = mMesh->get_set_names( EntityRank::ELEMENT);

                for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
                {
                    Matrix< IndexMat > tBlockSetElementInd
                        = mMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( Ik ) );

                    for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                    {
                       uint tNumVert = mMesh->get_mtk_cell( k ).get_number_of_vertices();

                       //print( mElements(k)->get_pdof_values(), "Element");

                       for( luint Jk=0; Jk < tNumVert; ++Jk )
                       {
                           moris_index tID= mMesh->get_mtk_cell( k ).get_vertex_pointers()( Jk) ->get_index();

                           mSolHMR(tID) = mElements(k)->get_pdof_values()(Jk);
                       }
                    }
                }
            }
            else
            {
            // 8) Postprocessing
            // dof type list for the solution to write on the mesh
            moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

            uint tNumOfNodes = mNodes.size();

            // create a matrix to be filled  with the solution
            Matrix< DDRMat > tTempSolutionField( tNumOfNodes, 1 );

            // loop over the nodes
            for( uint i = 0; i < tNumOfNodes; i++ )
            {
                // get a list of elements connected to the ith node
                Matrix<IndexMat> tConnectedElements =
                    mMesh->get_entity_connected_to_entity_loc_inds( static_cast< moris_index >( i ),
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
                    real tElemVal = mElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i, tDofTypeList);
                    // add up the contribution of each element to the node value
                    tNodeVal = tNodeVal + tElemVal;
                }
                // fill the solution matrix with the node value
                tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
            }
            //print( tTempSolutionField, "tTempSolutionField" );

            // add field to the mesh
            mMesh->add_mesh_field_real_scalar_data_loc_inds( aFilePath,
                                                             EntityRank::NODE,
                                                             tTempSolutionField );

            // create output mesh
            std::string tOutputFile = "./int_ElemDiff_test_11.exo";
            mMesh->create_output_mesh( tOutputFile );
            }
        }

    } /* namespace mdl */
} /* namespace moris */
