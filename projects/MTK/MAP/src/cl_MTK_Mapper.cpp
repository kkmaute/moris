#include "assert.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Mapper.hpp"

#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Set_User_Info.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mapper_Node.hpp"
#include "cl_MTK_Field.hpp"

#include "cl_HMR_Database.hpp"     //HMR/src
#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR_Mesh.hpp" //HMR/src
#include "cl_HMR_T_Matrix.hpp" //HMR/src

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_SOL_Warehouse.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "op_elemwise_mult.hpp"
#include "op_div.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_eye.hpp"

#include "cl_MDL_Model.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        Mapper::Mapper()
        {
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

        void Mapper::map_input_field_to_output_field(
                mtk::Field * aFieldSource,
                mtk::Field * aFieldTarget )
        {
            mFieldIn = aFieldSource;
            mFieldOut = aFieldTarget;

            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairIn = mFieldIn->get_mesh_pair();
            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairOut = mFieldOut->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn .second->get_interpolation_mesh( tMeshPairIn .first );
            moris::mtk::Mesh * tTargetMesh = tMeshPairOut.second->get_interpolation_mesh( tMeshPairOut.first );

            MORIS_ERROR( tSourceMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Source mesh is not and HMR mesh" );
            MORIS_ERROR( tTargetMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Target mesh is not and HMR mesh" );
            MORIS_ERROR( tSourceMesh->get_HMR_database() == tTargetMesh->get_HMR_database(),
                    "Mapper::map_input_field_to_output_field() Source and target mesh are based on two different instances of HMR" );

            std::shared_ptr< hmr::Database > tHMRDatabase = tSourceMesh->get_HMR_database();

            // grab orders of meshes
            uint tSourceLagrangeOrder = tSourceMesh->get_order();
            uint tTargetLagrangeOrder = tTargetMesh->get_order();

            // get order of Union Mesh
            uint tLagrangeOrder = std::max( tSourceLagrangeOrder, tTargetLagrangeOrder );

            uint tSourcePattern = tSourceMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            uint tTargetPattern = tTargetMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            uint tUnionPattern  = tHMRDatabase->get_parameters()->get_union_pattern();

            uint tUnionDescritizationOrder = mFieldOut->get_discretization_order();

            // create union pattern
            tHMRDatabase->create_union_pattern(
                    tSourcePattern,
                    tTargetPattern,
                    tUnionPattern );

            // create union mesh
            hmr::Interpolation_Mesh_HMR * tUnionInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                    tHMRDatabase,
                    tLagrangeOrder,
                    tUnionPattern,
                    tUnionDescritizationOrder,
                    tTargetPattern); // order, Lagrange pattern, bspline pattern

            //construct union integration mesh (note: this is not ever used but is needed for mesh manager)
            hmr::Integration_Mesh_HMR* tIntegrationUnionMesh = new hmr::Integration_Mesh_HMR(
                    tLagrangeOrder,
                    tUnionPattern,
                    tUnionInterpolationMesh);

            // Create mesh manager
            std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();

            // Register mesh pair
            uint tMeshIndexUnion = tMeshManager->register_mesh_pair( tUnionInterpolationMesh, tIntegrationUnionMesh );

            mtk::Field tFieldUnion( tMeshManager, tMeshIndexUnion );

            // map source Lagrange field to target Lagrange field
            if( tSourceLagrangeOrder >= tTargetLagrangeOrder )
            {
                // interpolate field onto union mesh
                this->interpolate_field(
                        mFieldIn,
                        &tFieldUnion );
            }
            else
            {
                // create union mesh. Bspline order will not be used
                hmr::Interpolation_Mesh_HMR * tHigherOrderInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                        tHMRDatabase,
                        tLagrangeOrder,
                        tSourcePattern,
                        tLagrangeOrder,
                        tSourcePattern); // order, Lagrange pattern, bspline order, bspline pattern

                //construct union integration mesh (note: this is not ever used but is needed for mesh manager)
                hmr::Integration_Mesh_HMR* tHigherOrderIntegrationMesh = new hmr::Integration_Mesh_HMR(
                        tLagrangeOrder,
                        tSourcePattern,
                        tHigherOrderInterpolationMesh);

                uint tMeshIndexUnion = tMeshManager->register_mesh_pair( tHigherOrderInterpolationMesh, tHigherOrderIntegrationMesh );

                mtk::Field tFieldHigerOrder( tMeshManager, tMeshIndexUnion );

                this->change_field_order( aFieldSource, &tFieldHigerOrder );

                // interpolate field onto union mesh
                this->interpolate_field(
                        &tFieldHigerOrder,
                        &tFieldUnion );
            }

            // project field to union
            this->perform_mapping(
                    &tFieldUnion,
                    EntityRank::NODE,
                    EntityRank::BSPLINE);

            // move coefficients to output field
            mFieldOut->get_coefficients() = std::move( tFieldUnion.get_coefficients() );

            this->perform_mapping(
                    mFieldOut,
                    EntityRank::BSPLINE,
                    EntityRank::NODE);
        }

        // -----------------------------------------------------------------------------

        // interpolate field values from source Lagrange to target Lagrange mesh
        void Mapper::interpolate_field(
                mtk::Field * aFieldSource,
                mtk::Field * aFieldTarget )
        {
            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairIn = aFieldSource->get_mesh_pair();
            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairOut = aFieldTarget->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn .second->get_interpolation_mesh( tMeshPairIn .first );
            moris::mtk::Mesh * tTargetMesh = tMeshPairOut.second->get_interpolation_mesh( tMeshPairOut.first );

            MORIS_ERROR( tSourceMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::interpolate_field() Source mesh is not and HMR mesh" );
            MORIS_ERROR( tTargetMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::interpolate_field() Target mesh is not and HMR mesh" );
            MORIS_ERROR( tSourceMesh->get_HMR_database() == tTargetMesh->get_HMR_database(),
                    "Mapper::interpolate_field() Source and target mesh are based on two different instances of HMR" );

            // grab orders of meshes
            uint tSourceLagrangeOrder = tSourceMesh->get_order();
            uint tTargetLagrangeOrder = tTargetMesh->get_order();

            MORIS_ERROR( tSourceLagrangeOrder == tTargetLagrangeOrder,
                    "Mapper::interpolate_field() Source and target mesh order are not the same" );

            MORIS_ERROR( aFieldTarget->get_number_of_dimensions() == aFieldSource->get_number_of_dimensions(),
                    "Mapper::interpolate_field() Source and target field number of dimmensions differ" );

            // pointer to mesh that is linked to input field
            hmr::Lagrange_Mesh_Base * tSourceLagrangeMesh = tSourceMesh->get_HMR_lagrange_mesh();
            hmr::Lagrange_Mesh_Base * tTargetLagrangeMesh = tTargetMesh->get_HMR_lagrange_mesh();

            uint tSourcePattern = tSourceLagrangeMesh->get_activation_pattern();
            //uint tTargetPattern = tTargetLagrangeMesh->get_activation_pattern();

            tTargetLagrangeMesh->select_activation_pattern();

            // unflag nodes on target
            tTargetLagrangeMesh->unflag_all_basis();

            // number of elements on target mesh
            auto tNumberOfElements = tTargetLagrangeMesh->get_number_of_elements();

            // number of nodes per element
            auto tNumberOfNodesPerElement = tTargetLagrangeMesh->get_number_of_basis_per_element();

            // create unity matrix
            Matrix< DDRMat > tEye;
            eye( tNumberOfNodesPerElement, tNumberOfNodesPerElement, tEye );

            // get values of source field
            const Matrix< DDRMat > & tSourceData = aFieldSource->get_node_values();

            // get target data
            Matrix< DDRMat > & tTargetData = aFieldTarget->get_node_values();

            // allocate value matrix
            tTargetData.set_size( tTargetLagrangeMesh->get_number_of_all_basis_on_proc(), aFieldTarget->get_number_of_dimensions() );

            // containers for source and target data
            Matrix< DDRMat > tElementSourceData( tNumberOfNodesPerElement, aFieldSource->get_number_of_dimensions() );

            hmr::T_Matrix * tTMatrix = new hmr::T_Matrix(
                    tSourceMesh->get_HMR_database()->get_parameters(),
                    tTargetLagrangeMesh );

            // loop over all elements
            for( luint Ie = 0; Ie < tNumberOfElements; ++Ie )
            {
                // get pointer to target element
                auto tTargetElement = tTargetLagrangeMesh->get_element( Ie );

                // get background element
                auto tBackgroundElement = tTargetElement->get_background_element();

                // initialize refinement Matrix
                Matrix< DDRMat > tR( tEye );

                while( ! tBackgroundElement->is_active( tSourcePattern ) )
                {
                    // right multiply refinement matrix
                    tR = tR * tTMatrix->get_refinement_matrix( tBackgroundElement->get_child_index() );

                    // jump to parent
                    tBackgroundElement = tBackgroundElement->get_parent();
                }

                // get pointer to source element
                auto tSourceElement = tSourceLagrangeMesh->get_element_by_memory_index( tBackgroundElement->get_memory_index() );

                // fill source data vector
                for( uint Ik = 0; Ik < tNumberOfNodesPerElement; ++Ik )
                {
                    // get pointer to source node
                    auto tNode  = tSourceElement->get_basis( Ik );
                    auto tIndex = tNode->get_index();

                    // copy data from source mesh
                    tElementSourceData.set_row( Ik, tSourceData.get_row( tIndex ) );
                }

                // copy target data to target mesh
                for( uint Ik = 0; Ik < tNumberOfNodesPerElement; ++Ik )
                {
                    // get pointer to target node
                    auto tNode = tTargetElement->get_basis( Ik );

                    // test if data has already been written to target
                    if ( ! tNode->is_flagged() )
                    {
                        // get node index
                        auto tIndex = tNode->get_index();

                        tTargetData.set_row( tIndex, tR.get_row( Ik ) * tElementSourceData );

                        // flag this node
                        tNode->flag();
                    }
                }
            }

            delete( tTMatrix );
        }

        // -----------------------------------------------------------------------------

        void Mapper::change_field_order(
                mtk::Field * aFieldSource,
                mtk::Field * aFieldTarget )
        {
            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairIn = aFieldSource->get_mesh_pair();
            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairOut = aFieldTarget->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn .second->get_interpolation_mesh( tMeshPairIn .first );
            moris::mtk::Mesh * tTargetMesh = tMeshPairOut.second->get_interpolation_mesh( tMeshPairOut.first );

            // pointer to mesh that is linked to input field
            hmr::Lagrange_Mesh_Base * tSourceLagrangeMesh = tSourceMesh->get_HMR_lagrange_mesh();
            hmr::Lagrange_Mesh_Base * tTargetLagrangeMesh = tTargetMesh->get_HMR_lagrange_mesh();

            // make sure that meshes are compatible
            MORIS_ASSERT( tSourceLagrangeMesh->get_activation_pattern() == tTargetLagrangeMesh->get_activation_pattern(),
                    "Mapper::change_field_order() Source and target field have to be build based on the same pattern" );

            tTargetLagrangeMesh->select_activation_pattern();

            // unflag all nodes on out mesh
            tTargetLagrangeMesh->unflag_all_basis();

            // source values
            const Matrix< DDRMat > & tSourceValues = aFieldSource->get_node_values();

            // target values
            Matrix< DDRMat > & tTargetValues = aFieldTarget->get_node_values();

            // allocate output memory
            tTargetValues.set_size( tTargetLagrangeMesh->get_number_of_nodes_on_proc(), 1 );

            // get number of elements
            uint tNumberOfElements = tSourceLagrangeMesh->get_number_of_elements();

            hmr::T_Matrix * tTMatrix = new hmr::T_Matrix(
                    tSourceMesh->get_HMR_database()->get_parameters(),
                    tSourceLagrangeMesh );

            uint tTargetMeshOrder = tTargetMesh->get_order();
            Matrix< DDRMat > tT = tTMatrix->get_change_order_matrix( tTargetMeshOrder );

            // delete t-Matrix object
            delete tTMatrix;

            uint tNumberOfNodesPerSourceElement = tSourceLagrangeMesh->get_number_of_basis_per_element();
            uint tNumberOfNodesPerTargetElement = tTargetLagrangeMesh->get_number_of_basis_per_element();

            Matrix< DDRMat > tLocalSourceValues( tNumberOfNodesPerSourceElement, 1 );
            Matrix< DDRMat > tN( 1, tNumberOfNodesPerSourceElement );

            // loop over all elements
            for( uint Ie = 0; Ie < tNumberOfElements; ++Ie )
            {
                // get pointer to source element
                const hmr::Element * tSourceElement = tSourceLagrangeMesh->get_element( Ie );

                for( uint Ii = 0; Ii < tNumberOfNodesPerSourceElement; ++Ii )
                {
                    tLocalSourceValues( Ii ) = tSourceValues( tSourceElement->get_basis( Ii )->get_index() );
                }

                // get pointer to target element
                hmr::Element * tTargetElement = tTargetLagrangeMesh->get_element( Ie );

                // loop over all nodes on target
                for( uint Ik = 0; Ik < tNumberOfNodesPerTargetElement; ++Ik )
                {
                    // get basis
                    hmr::Basis * tNode = tTargetElement->get_basis( Ik );

                    if( ! tNode->is_flagged() )
                    {
                        // copy row from T-Matrix
                        tT.get_row( Ik, tN );

                        // interpolate values
                        tTargetValues( tNode->get_index() ) = dot( tN, tLocalSourceValues );

                        // flag node
                        tNode->flag();
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        void Mapper::create_iwg_and_model(
                mtk::Field * aField,
                const real aAlpha )
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

                std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPairIn = aField->get_mesh_pair();

                // create model
                mModel = new mdl::Model(
                        tMeshPairIn.second.get(),
                        aField->get_discretization_mesh_index(),
                        tSetInfo,
                        tMeshPairIn.first );

                // set bool for building IWG and model to true
                mHaveIwgAndModel = true;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mapper::perform_mapping(
                mtk::Field * aField,
                const enum EntityRank aSourceEntityRank,
                const enum EntityRank aTargetEntityRank )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map");

            switch( aSourceEntityRank )
            {
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                case EntityRank::NODE :
                {
                    switch( aTargetEntityRank )
                    {
                        case EntityRank::BSPLINE:
                        {
                            this->map_node_to_bspline_from_field( aField );
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
                case EntityRank::BSPLINE:
                {
                    switch( aTargetEntityRank )
                    {
                        case EntityRank::NODE:
                        {
                            this->map_bspline_to_node_same_mesh( aField );
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

        void Mapper::map_node_to_bspline( Matrix<DDRMat> & aSolution )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map Node-to-Bspline");

            moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::L2 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algortihm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm =
                    tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm =
                    tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")                = 2;
            tNonlinearSolverAlgorithm->set_param("NLA_hard_break")              = false;
            tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
            tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian")        = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;

            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create time Solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tsa::Time_Solver_Factory tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm =
                    tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface(mModel->get_solver_interface());

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( tDofTypes1 );
            tTimeSolver.set_dof_type_list( tDofTypes1 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tTimeSolver.solve();
            tTimeSolver.get_full_solution( aSolution );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mapper::map_node_to_bspline_from_field( mtk::Field * aField )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map Node-to-Bspline");

            // create the model if it has not been created yet
            this->create_iwg_and_model( aField );

            // set weak bcs from field
            mModel->set_weak_bcs( aField->get_node_values() );

            this->map_node_to_bspline( aField->get_coefficients() );
        }

        //------------------------------------------------------------------------------

        void
        Mapper::map_bspline_to_node_same_mesh( mtk::Field * aField )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map Bspline-to-Node");

            std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > tMeshPair = aField->get_mesh_pair();

            moris::mtk::Mesh * tInterpolationMesh = tMeshPair.second->get_interpolation_mesh( tMeshPair.first );

            moris_index tDescritizationIndex= aField->get_discretization_mesh_index();

            // get number of nodes on block
            uint tNumberOfNodes= tInterpolationMesh->get_num_nodes();

            Matrix< DDRMat > & tNodalValues = aField->get_node_values();

            tNodalValues.set_size( tNumberOfNodes, 1 );

            const Matrix< DDRMat > & tCoefficients = aField->get_coefficients();

            for( uint Ik = 0; Ik < tNumberOfNodes; ++Ik )
            {
                // get pointer to node
                auto tNode = &tInterpolationMesh->get_mtk_vertex( Ik );

                // get PDOFs from node
                auto tBSplines = tNode->
                        get_interpolation( tDescritizationIndex )->
                        get_coefficients();

                // get T-Matrix
                const Matrix< DDRMat > & tTMatrix = *tNode->
                        get_interpolation( tDescritizationIndex )->
                        get_weights();

                // get number of coefficients
                uint tNumberOfCoeffs = tTMatrix.length();

                MORIS_ASSERT( tNumberOfCoeffs > 0, "No coefficients defined for node" ) ;

                // fill coeffs vector
                Matrix< DDRMat > tCoeffs( tNumberOfCoeffs, 1 );
                for( uint Ii = 0; Ii < tNumberOfCoeffs; ++Ii )
                {
                    tCoeffs( Ii ) = tCoefficients( tBSplines( Ii )->get_index() );
                }

                // write value into solution
                tNodalValues( Ik ) = moris::dot( tTMatrix, tCoeffs );
            }
        }

        // FIXME do not delete for future use
        ////------------------------------------------------------------------------------
        //
        //        void
        //        Mapper::map_node_to_element_same_mesh(
        //                         const moris_index   aSourceIndex,
        //                         const moris_index   aTargetIndex )
        //        {
        //            // create the model if it has not been created yet
        //            this->create_iwg_and_model();
        //
        //            // set weak bcs from field
        //            mModel->set_weak_bcs_from_nodal_field( aSourceIndex );
        //
        //            // get number of elements
        //            uint tNumberOfElements = mTargetMesh->get_num_elems();
        //
        //            // loop over all elements
        //            for( uint e=0; e<tNumberOfElements; ++e )
        //            {
        //                // get ref to entry in database
        //                real & tValue = mTargetMesh->get_value_of_scalar_field(
        //                        aTargetIndex,
        //                        EntityRank::ELEMENT,
        //                        e );
        //
        //                // calculate value
        //                tValue = mModel->compute_element_average( e );
        //            }
        //        }

        //------------------------------------------------------------------------------

        //        void
        //        Mapper::create_nodes_for_filter()
        //        {
        //            if( ! mHaveNodes )
        //            {
        //                // get number of nodes from mesh
        //                uint tNumberOfNodes = mSourceMesh->get_num_nodes();
        //
        //                // reserve node container
        //                mNodes.resize( tNumberOfNodes, nullptr );
        //
        //
        //                // populate container
        //                for( uint k=0; k<tNumberOfNodes; ++k )
        //                {
        //                    mNodes( k ) = new Node( &mSourceMesh->get_mtk_vertex( k ) );
        //                }
        //
        //                // link to neighbors
        //                /*for( uint k=0; k<tNumberOfNodes; ++k )
        //                {
        //                    Matrix< IndexMat > tNodeIndices =
        //                            mSourceMesh->get_entity_connected_to_entity_loc_inds(
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
        //            moris_index tFieldIndex = mSourceMesh->get_field_ind( aSourceLabel,
        //                                                                  EntityRank::NODE );
        //
        //            const Matrix< DDRMat > & tSourceField = mSourceMesh->get_field( tFieldIndex, EntityRank::NODE );
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
