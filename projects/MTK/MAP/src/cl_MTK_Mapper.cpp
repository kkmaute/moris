
#include "assert.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Mapper.hpp"

#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Set_User_Info.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

#include "cl_MSI_Solver_Interface.hpp"

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

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "op_elemwise_mult.hpp"
#include "op_div.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"

#include "cl_MDL_Model.hpp"

namespace moris
{
    namespace mapper
    {

//------------------------------------------------------------------------------

        Mapper::Mapper( mtk::Mesh_Manager* aMesh,
                        const moris_index  aMeshPairIndex,
                        const uint         aBSplineOrder ) : mSourceMeshPairIndex( aMeshPairIndex ),
                                                             mTargetMeshPairIndex( aMeshPairIndex ),
                                                             mMeshManager( aMesh ),
                                                             mBSplineOrder( aBSplineOrder )
        {
            // Retrieve source mesh pair
            mMeshManager->get_mesh_pair(aMeshPairIndex,mSourceInterpMesh,mSourceIntegMesh);

            // Retrieve target mesh pair
            mMeshManager->get_mesh_pair(aMeshPairIndex,mTargetInterpMesh,mTargetIntegMesh);
        }

//------------------------------------------------------------------------------

        Mapper::~Mapper()
        {
            // test if model and IWG have been created
            if( mHaveIwgAndModel )
            {
                // delete the fem model
                delete mModel;
            }
        }

//------------------------------------------------------------------------------
        void Mapper::create_iwg_and_model( const real aAlpha )
        {
            if( ! mHaveIwgAndModel )
            {
                // create a L2 IWG
                //FIXME should be provided to the function
                fem::IWG_Factory tIWGFactory;
                std::shared_ptr< fem::IWG > tIWGL2 = tIWGFactory.create_IWG( fem::IWG_Type::L2 );
                tIWGL2->set_residual_dof_type( { MSI::Dof_Type::L2 } );
                tIWGL2->set_dof_type_list( {{ MSI::Dof_Type::L2 }}, mtk::Master_Slave::MASTER );

                // define set info
                //FIXME should be provided to the function
                moris::Cell< fem::Set_User_Info > tSetInfo( 1 );
                tSetInfo( 0 ).set_mesh_index( 0 );
                tSetInfo( 0 ).set_IWGs( { tIWGL2 } );

                // create model
                mModel = new mdl::Model( mMeshManager,
                                         mBSplineOrder,
                                         tSetInfo,
                                         mTargetMeshPairIndex );

                // set bool for building IWG and model to true
                mHaveIwgAndModel = true;
            }
        }

//-----------------------------------------------------------------------------

        void Mapper::set_l2_alpha( const real & aAlpha )
        {
            // remove model
            if( ! mHaveIwgAndModel )
            {
                this->create_iwg_and_model( aAlpha );
            }
            else
            {
                MORIS_ERROR(false, "Model does exist set alpha there");
//                mIWG->set_alpha( aAlpha );
            }
        }

//-----------------------------------------------------------------------------

        void Mapper::perform_mapping( const std::string      & aSourceLabel,
                                      const enum EntityRank    aSourceEntityRank,
                                      const std::string      & aTargetLabel,
                                      const enum EntityRank    aTargetEntityRank )
        {
            // get index of source
            moris_index tSourceIndex = mSourceInterpMesh->get_field_ind( aSourceLabel,
                                                                         aSourceEntityRank );

            MORIS_ERROR( tSourceIndex != gNoIndex, "perform_mapping() Source Field not found");

            // get target index
            moris_index tTargetIndex = mTargetInterpMesh->get_field_ind( aTargetLabel,
                                                                         aTargetEntityRank );

            // test if output field has to be initialized
            if( tTargetIndex == gNoIndex )
            {
                // create target field
                tTargetIndex = mTargetInterpMesh->create_scalar_field( aTargetLabel,
                                                                       aTargetEntityRank );
            }

            switch( aSourceEntityRank )
            {
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                case( EntityRank::NODE ) :
                {
                   switch( aTargetEntityRank )
                   {
                       case( EntityRank::BSPLINE ) :
                       case( EntityRank::BSPLINE_2 ) :
                       case( EntityRank::BSPLINE_3 ) :
                       {
                          this->map_node_to_bspline_same_mesh( tSourceIndex,
                                                               tTargetIndex,
                                                               aTargetEntityRank );
                           break;
                       }
                       default :
                       {
                           MORIS_ERROR( false, "perform_mapping(): aTargetEntityRank not supported.");
                           break;
                       }
                   }
                   break;
                }
                case( EntityRank::BSPLINE ) :
                case( EntityRank::BSPLINE_2 ) :
                case( EntityRank::BSPLINE_3 ) :
                {
                    switch( aTargetEntityRank )
                    {
                        case( EntityRank::NODE ) :
                        {
                            this->map_bspline_to_node_same_mesh( tSourceIndex,
                                                                 aSourceEntityRank,
                                                                 tTargetIndex );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false, "perform_mapping(): aTargetEntityRank not supported.");
                            break;
                        }
                    }
                    break;
                }
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                default :
                {
                    MORIS_ERROR( false, "perform_mapping(): aSourceEntityRank not supported.");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void Mapper::map_node_to_bspline_same_mesh( const moris_index     aSourceIndex,
                                                    const moris_index     aTargetIndex,
                                                    const enum EntityRank aBSplineRank )
        {
            // create the model if it has not been created yet
            this->create_iwg_and_model();

            // set order to mBSplineOrder
            // mModel->set_dof_order( mBSplineOrder );

            // set weak bcs from field
            mModel->set_weak_bcs_from_nodal_field( aSourceIndex );

            moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::L2 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algortihm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

//            tLinearSolverAlgorithm->set_param("AZ_keep_info") = 1;
            //tLinearSolverAlgorithm->set_param("AZ_pre_calc") = AZ_reuse;
//            tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 5;

//            tLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;
            //tLinearSolverAlgorithm->set_param("ML_reuse") = true;

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algortihm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
            tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
            tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
            tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;

            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create time Solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tsa::Time_Solver_Factory tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            NLA::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface(mModel->get_solver_interface());

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( tDofTypes1 );
            tTimeSolver.set_dof_type_list( tDofTypes1 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tTimeSolver.solve();

            // test if output mesh is HMR
            if( mTargetInterpMesh->get_mesh_type() == MeshType::HMR )
            {
                // perform L2 projection
                tTimeSolver.get_full_solution( mTargetInterpMesh->get_field( aTargetIndex, aBSplineRank ) );

                Matrix< DDUMat > tAdofMap = mModel->get_adof_map();

                //print(tAdofMap,"tAdofMap");

                // temporary array for solver
                Matrix< DDRMat > tSolution = mTargetInterpMesh->get_field( aTargetIndex, aBSplineRank );

                uint tLength = tSolution.length();

                // rearrange data into output
                mTargetInterpMesh->get_field( aTargetIndex, aBSplineRank ).set_size( tLength, 1 );

                for( uint k=0; k<tLength; ++k )
                {
                    mTargetInterpMesh->get_field( aTargetIndex, aBSplineRank )( k ) = tSolution( tAdofMap( k ) );
                }
            }
            else
            {
                moris::Matrix< DDRMat > tSolution;
                tTimeSolver.get_full_solution( tSolution );

                // get number of coeffs of
                uint tNumberOfCoeffs = mTargetInterpMesh->get_num_coeffs( mBSplineOrder );

                // make sure that solution is correct
                MORIS_ERROR( tNumberOfCoeffs == tSolution.length(),
                                  "perform_mapping() number of coeffs does not match" );

                // copy solution into target
                for( uint k=0; k<tNumberOfCoeffs; ++k )
                {
                    // get ref to value
                    mTargetInterpMesh->get_value_of_scalar_field(
                            aTargetIndex,
                            aBSplineRank,
                            k ) = tSolution( k );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::map_bspline_to_node_same_mesh(
                         const moris_index     aSourceIndex,
                         const enum EntityRank aBSplineRank,
                         const moris_index     aTargetIndex )
        {

            // get number of nodes
            moris_index tNumberOfNodes = mTargetInterpMesh->get_num_nodes();

            // loop over all nodes
            for( moris_index k=0;  k<tNumberOfNodes; ++k )
            {
                // get weights
                const Matrix< DDRMat > & tTMatrix = mTargetInterpMesh->
                        get_t_matrix_of_node_loc_ind( k, aBSplineRank );

                // get indices
                Matrix< IndexMat > tBSplines = mTargetInterpMesh->
                        get_bspline_inds_of_node_loc_ind( k, aBSplineRank );

                // get number of coefficients
                uint tNumberOfCoeffs = tTMatrix.length();

                // value of node
                real & tValue = mTargetInterpMesh->get_value_of_scalar_field(
                        aTargetIndex,
                        EntityRank::NODE,
                        k );

                // reset value
                tValue = 0.0;

                for( uint i=0; i<tNumberOfCoeffs; ++i )
                {
                    tValue +=  mSourceInterpMesh->get_value_of_scalar_field(
                            aSourceIndex,
                            aBSplineRank,
                            tBSplines( i ) ) * tTMatrix( i );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Mapper::map_node_to_element_same_mesh(
                         const moris_index   aSourceIndex,
                         const moris_index   aTargetIndex )
        {
            // create the model if it has not been created yet
            this->create_iwg_and_model();

            // set weak bcs from field
            mModel->set_weak_bcs_from_nodal_field( aSourceIndex );

            // get number of elements
            uint tNumberOfElements = mTargetInterpMesh->get_num_elems();

            // loop over all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // get ref to entry in database
                real & tValue = mTargetInterpMesh->get_value_of_scalar_field(
                        aTargetIndex,
                        EntityRank::ELEMENT,
                        e );

                // calculate value
                tValue = mModel->compute_element_average( e );
            }
        }

 //------------------------------------------------------------------------------

//        void
//        Mapper::create_nodes_for_filter()
//        {
//            if( ! mHaveNodes )
//            {
//                // get number of nodes from mesh
//                uint tNumberOfNodes = mSourceInterpMesh->get_num_nodes();
//
//                // reserve node container
//                mNodes.resize( tNumberOfNodes, nullptr );
//
//
//                // populate container
//                for( uint k=0; k<tNumberOfNodes; ++k )
//                {
//                    mNodes( k ) = new Node( &mSourceInterpMesh->get_mtk_vertex( k ) );
//                }
//
//                // link to neighbors
//                /*for( uint k=0; k<tNumberOfNodes; ++k )
//                {
//                    Matrix< IndexMat > tNodeIndices =
//                            mSourceInterpMesh->get_entity_connected_to_entity_loc_inds(
//                                    k,
//                                    EntityRank::NODE,
//                                    EntityRank::NODE );
//
//                    uint tNumberOfConnectedNodes = tNodeIndices.length();
//                    mNodes( k )->init_neighbor_container( tNumberOfConnectedNodes );
//
//                    for( uint i=0; i<tNumberOfConnectedNodes; ++i )
//                    {
//                        mNodes( k )->insert_neighbor( mNodes( tNodeIndices( i ) ) );
//                    }
//
//                } */
//
//                // set node flag
//                mHaveNodes = true;
//            }
//        }
//------------------------------------------------------------------------------

//        void
//        Mapper::perform_filter(
//                        const std::string & aSourceLabel,
//                        const real        & aFilterRadius,
//                        Matrix< DDRMat >  & aValues )
//        {
//
//            MORIS_ERROR( par_size() == 1,
//                    "The filter is not written for parallel. In order do use it, mtk::Mapper needs access to node information from the aura.");
//
//            // fixme: the following two lines only work for HMR
//            moris_index tFieldIndex = mSourceInterpMesh->get_field_ind( aSourceLabel,
//                                                                  EntityRank::NODE );
//
//            const Matrix< DDRMat > & tSourceField = mSourceInterpMesh->get_field( tFieldIndex, EntityRank::NODE );
//
//            // calculate weights if this was not done already
//            this->calculate_filter_weights( aFilterRadius );
//
//            // get number of nodes on target
//            uint tNumberOfNodes = mNodes.size();
//
//            aValues.set_size( tNumberOfNodes, 1 );
//
//            for( uint k=0; k<tNumberOfNodes; ++k )
//            {
//
//                Matrix< IndexMat > & tIndices = mNodes( k )->get_node_indices();
//
//                uint tNumberOfIndices = tIndices.length();
//
//                Matrix< DDRMat > tValues( tNumberOfIndices , 1 );
//
//                for( uint i=0; i<tNumberOfIndices; ++i )
//                {
//                    tValues( i ) = tSourceField( tIndices( i ) );
//                }
//
//                // fill vector with values
//                aValues( k ) = dot ( mNodes( k )->get_weights(), tValues );
//            }
//        }

//------------------------------------------------------------------------------

//        void
//        Mapper::calculate_filter_weights( const real & aFilterRadius )
//        {
//            if( mFilterRadius != aFilterRadius )
//            {
//                // remember radius
//                mFilterRadius = aFilterRadius;
//
//                // create nodes for the filter
//                this->create_nodes_for_filter();
//
//                for( Node * tNode : mNodes )
//                {
//
//                    // flag myself
//                    tNode->flag();
//
//                    // cell containing neighbors
//                    Cell< Node * > tNeighbors;
//
//                    tNode->get_nodes_in_proximity( tNode->get_coords(), aFilterRadius, tNeighbors );
//
//                    uint tNumberOfNeighbors = tNeighbors.size();
//
//                    Matrix< DDRMat > & tWeights = tNode->get_weights();
//                    tWeights.set_size( tNumberOfNeighbors, 1 );
//
//                    Matrix< IndexMat > & tIndices = tNode->get_node_indices();
//                    tIndices.set_size( tNumberOfNeighbors, 1 );
//
//                    real tMyLevel = tNode->get_level();
//
//                    uint tCount = 0;
//                    for( Node * tNeighbor : tNeighbors )
//                    {
//                        // Kurt's formula with level based average
//                        tWeights( tCount )   =
//                                ( aFilterRadius - tNeighbor->get_distance() )
//                               *  ( tMyLevel + 1.0 ) / ( ( real ) tNeighbor->get_level() + 1.0);
//
//                        // Simple Weight by distance
//                        //tWeights( tCount )   =
//                        //        ( aFilterRadius - tNeighbor->get_distance() );
//
//                        // save index
//                        tIndices( tCount++ ) = tNeighbor->get_index();
//
//                        // unflag neighbors
//                        tNeighbor->unflag();
//                    }
//
//                    tWeights = tWeights / sum( tWeights );
//
//                    // unflag this node
//                    tNode->unflag();
//                }
//            }
//
//        }
    } /* namespace mtk */
} /* namespace moris */
