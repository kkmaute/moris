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
#include "cl_MTK_Field_Discrete.hpp"

#include "cl_HMR_Database.hpp"     //HMR/src
#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR_Mesh.hpp" //HMR/src
#include "cl_HMR_T_Matrix.hpp" //HMR/src

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_SOL_Warehouse.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

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
            // cast output field to discrete field
            mtk::Field_Discrete * tDiscreteFieldTarget = dynamic_cast<mtk::Field_Discrete*> (aFieldTarget);

            // check that dynamic cast was successful
            MORIS_ERROR( tDiscreteFieldTarget != nullptr,
                    "Mapper::map_input_field_to_output_field - target field need to be discrete field.\n");

            mtk::Mesh_Pair tMeshPairIn  = aFieldSource->get_mesh_pair();
            mtk::Mesh_Pair tMeshPairOut = tDiscreteFieldTarget->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn.get_interpolation_mesh();
            moris::mtk::Mesh * tTargetMesh = tMeshPairOut.get_interpolation_mesh();

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

            uint tUnionDescritizationOrder = tDiscreteFieldTarget->get_discretization_order();

            // create union pattern
            tHMRDatabase->create_union_pattern(
                    tSourcePattern,
                    tTargetPattern,
                    tUnionPattern );

            uint tNumAdditionalMeshes = tHMRDatabase->get_number_of_additional_lagrange_meshes( tSourcePattern );

            for( uint Ik = 0; Ik < tNumAdditionalMeshes; Ik++ )
            {
                tHMRDatabase->get_additional_lagrange_mesh_by_index( Ik, tSourcePattern )->update_mesh();
            }

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

            mtk::Mesh_Pair tMeshPairUnion(tUnionInterpolationMesh, tIntegrationUnionMesh, true);

            mtk::Field_Discrete tFieldUnion( tMeshPairUnion, 0 );

            // map source Lagrange field to target Lagrange field
            if( tSourceLagrangeOrder >= tTargetLagrangeOrder )
            {
                // interpolate field onto union mesh
                this->interpolate_field(
                        aFieldSource,
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

                mtk::Mesh_Pair tMeshPairHigherOrder(tHigherOrderInterpolationMesh, tHigherOrderIntegrationMesh, true);

                mtk::Field_Discrete tFieldHigerOrder( tMeshPairHigherOrder, 0 );

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
            tDiscreteFieldTarget->unlock_field();
            tDiscreteFieldTarget->set_coefficients( tFieldUnion.get_coefficients() );

            this->perform_mapping(
                    tDiscreteFieldTarget,
                    EntityRank::BSPLINE,
                    EntityRank::NODE);
        }

        //------------------------------------------------------------------------------

        void Mapper::map_input_field_to_output_field_2( mtk::Field * aFieldSource )
        {
            // cast aFieldSource field to discrete field
            mtk::Field_Discrete * tDiscreteFieldSource = dynamic_cast<mtk::Field_Discrete*> (aFieldSource);

            // check that dynamic cast was successful
            MORIS_ERROR( tDiscreteFieldSource != nullptr,
                    "Mapper::map_input_field_to_output_field_2 - source field need to be discrete field.\n");

            mtk::Mesh_Pair tMeshPairIn = tDiscreteFieldSource->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn.get_interpolation_mesh();

            MORIS_ERROR( tSourceMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Source mesh is not and HMR mesh" );

            std::shared_ptr< hmr::Database > tHMRDatabase = tSourceMesh->get_HMR_database();

            // grab orders of meshes
            uint tSourceLagrangeOrder = tSourceMesh->get_order();
            uint tTargetOrder = tDiscreteFieldSource->get_discretization_order();

            // get order of Union Mesh
            uint tLagrangeOrder = std::max( tSourceLagrangeOrder, tTargetOrder );

            uint tSourcePattern = tSourceMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            uint tTargetPattern = tSourceMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            uint tUnionPattern  = tHMRDatabase->get_parameters()->get_union_pattern();

            uint tTargetBSPattern = tSourceMesh->
                    get_HMR_lagrange_mesh()->
                    get_bspline_pattern( tDiscreteFieldSource->get_discretization_mesh_index() );

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
                    tTargetOrder,
                    tTargetBSPattern); // order, Lagrange pattern, bspline pattern

            //construct union integration mesh (note: this is not ever used but is needed for mesh manager)
            hmr::Integration_Mesh_HMR* tIntegrationUnionMesh = new hmr::Integration_Mesh_HMR(
                    tLagrangeOrder,
                    tUnionPattern,
                    tUnionInterpolationMesh);

            mtk::Mesh_Pair tMeshPairUnion(tUnionInterpolationMesh, tIntegrationUnionMesh, true);

            mtk::Field_Discrete tFieldUnion( tMeshPairUnion, 0 );

            // map source Lagrange field to target Lagrange field
            if( tSourceLagrangeOrder >= tTargetOrder )
            {
                // interpolate field onto union mesh
                this->interpolate_field(
                        tDiscreteFieldSource,
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

                mtk::Mesh_Pair tMeshPairHigherOrder(tHigherOrderInterpolationMesh, nullptr, false);

                mtk::Field_Discrete tFieldHigerOrder( tMeshPairHigherOrder, 0 );

                this->change_field_order( tDiscreteFieldSource, &tFieldHigerOrder );

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
            tDiscreteFieldSource->unlock_field();
            tDiscreteFieldSource->set_coefficients( tFieldUnion.get_coefficients() );
        }

        // -----------------------------------------------------------------------------

        // interpolate field values from source Lagrange to target Lagrange mesh
        void Mapper::interpolate_field(
                mtk::Field * aFieldSource,
                mtk::Field * aFieldTarget )
        {
            mtk::Mesh_Pair tMeshPairIn  = aFieldSource->get_mesh_pair();
            mtk::Mesh_Pair tMeshPairOut = aFieldTarget->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn.get_interpolation_mesh();
            moris::mtk::Mesh * tTargetMesh = tMeshPairOut.get_interpolation_mesh();

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

            MORIS_ERROR( aFieldTarget->get_number_of_fields() == aFieldSource->get_number_of_fields(),
                    "Mapper::interpolate_field() - Number of fields differ in ource and target field.\n" );

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
            const Matrix< DDRMat > & tSourceData = aFieldSource->get_nodal_values();

            // get target data; note: the size of this vector is number of nodes times number of fields;
            // the number of nodes and number of basis differ as the number of basis includes the aura
            Matrix< DDRMat > tTargetData(
                    tTargetLagrangeMesh->get_number_of_nodes_on_proc(),
                    aFieldTarget->get_number_of_fields(),
                    MORIS_REAL_MAX);

            // containers for source and target data
            Matrix< DDRMat > tElementSourceData( tNumberOfNodesPerElement, aFieldSource->get_number_of_fields() );

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

            // check that target data is valid
            MORIS_ERROR( tTargetData.max() < MORIS_REAL_MAX,
                    "Mapper::interpolate_field - target data not valid.\n");

            // copy target data onto target field
            aFieldTarget->unlock_field();
            aFieldTarget->set_nodal_values( tTargetData );

            delete( tTMatrix );
        }

        // -----------------------------------------------------------------------------

        void Mapper::change_field_order(
                mtk::Field * aFieldSource,
                mtk::Field * aFieldTarget )
        {
            mtk::Mesh_Pair tMeshPairIn = aFieldSource->get_mesh_pair();
            mtk::Mesh_Pair tMeshPairOut = aFieldTarget->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn.get_interpolation_mesh();
            moris::mtk::Mesh * tTargetMesh = tMeshPairOut.get_interpolation_mesh();

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
            const Matrix< DDRMat > & tSourceValues = aFieldSource->get_nodal_values();

            // target values
            Matrix< DDRMat > tTargetValues( tTargetLagrangeMesh->get_number_of_nodes_on_proc(), 1 );

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

            aFieldTarget->unlock_field();
            aFieldTarget->set_nodal_values( tTargetValues );
        }

        //------------------------------------------------------------------------------
        void Mapper::create_iwg_and_model(
                mtk::Field * aField,
                const real aAlpha )
        {
            // cast aFieldSource field to discrete field
             mtk::Field_Discrete * tDiscreteField = dynamic_cast<mtk::Field_Discrete*> (aField);

             // check that dynamic cast was successful
             MORIS_ERROR( tDiscreteField != nullptr,
                     "Mapper::map_input_field_to_output_field_2 - source field need to be discrete field.\n");

            mtk::Mesh_Pair tMeshPair = tDiscreteField->get_mesh_pair();

            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();

            uint MeshPairIndex = tMeshManager->register_mesh_pair( tMeshPair );

            if( ! mHaveIwgAndModel )
            {
                // create a L2 IWG
                //FIXME should be provided to the function
                fem::IWG_Factory tIWGFactory;
                std::shared_ptr< fem::IWG > tIWGL2 = tIWGFactory.create_IWG( fem::IWG_Type::L2 );
                tIWGL2->set_residual_dof_type( { { MSI::Dof_Type::L2 } } );
                tIWGL2->set_dof_type_list( {{ MSI::Dof_Type::L2 }}, mtk::Master_Slave::MASTER );

                // define set info
                //FIXME should be provided to the function
                moris::Cell< fem::Set_User_Info > tSetInfo( 1 );
                tSetInfo( 0 ).set_mesh_index( 0 );
                tSetInfo( 0 ).set_IWGs( { tIWGL2 } );

                // create model
                mModel = new mdl::Model(
                        tMeshManager,
                        tDiscreteField->get_discretization_mesh_index(),
                        tSetInfo,
                        MeshPairIndex );

                // set bool for building IWG and model to true
                mHaveIwgAndModel = true;
            }
            // set weak bcs from field
            mModel->set_weak_bcs( tDiscreteField->get_nodal_values() );

            this->map_node_to_bspline( tDiscreteField );
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

        void Mapper::map_node_to_bspline(  mtk::Field * aField )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map Node-to-Bspline");

            // define time, nonlinear and linear solver
            sol::SOL_Warehouse tSolverWarehouse( mModel->get_solver_interface() );

            moris::Cell< moris::Cell< moris::ParameterList > > tParameterlist( 7 );
            for( uint Ik = 0; Ik < 7; Ik ++)
            {
                tParameterlist( Ik ).resize(1);
            }

            // choose solver type based on problem size
            // FIXME: solver should be received from solver warehouse
            uint tNumberOfCoefficients = aField->get_number_of_coefficients();

            if ( tNumberOfCoefficients*par_size() < 100000 )
            {
                tParameterlist( 0 )(0) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
            }
            else
            {
                tParameterlist( 0 )(0) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
                tParameterlist( 0 )( 0 ).set( "ifpack_prec_type", "ILU");
                tParameterlist( 0 )( 0 ).set( "fact: level-of-fill", 1);
            }

            tParameterlist( 1 )(0) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 2 )(0) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )(0).set( "NLA_max_iter", 2 );

            tParameterlist( 3 )(0) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )(0).set("NLA_DofTypes"      , "L2" );

            tParameterlist( 4 )(0) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 5 )(0) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )(0).set("TSA_DofTypes"      , "L2" );

            tParameterlist( 6 )(0) = moris::prm::create_solver_warehouse_parameterlist();
            tParameterlist( 6 )(0).set("SOL_TPL_Type"               , static_cast< uint >( sol::MapType::Epetra ) );

            tSolverWarehouse.set_parameterlist( tParameterlist );

            tSolverWarehouse.initialize();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tsa::Time_Solver * tTimeSolver = tSolverWarehouse.get_main_time_solver();

            tTimeSolver->solve();
            
            Matrix<DDRMat> tSolution;
            tTimeSolver->get_full_solution( tSolution );

            aField->unlock_field();
            aField->set_coefficients( tSolution );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Mapper::map_node_to_bspline_from_field( mtk::Field * aField )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map Node-to-Bspline");

            // create the model if it has not been created yet
            this->create_iwg_and_model( aField );
        }

        //------------------------------------------------------------------------------

        void
        Mapper::map_bspline_to_node_same_mesh( mtk::Field * aField )
        {
            // Tracer
            Tracer tTracer("MTK", "Mapper", "Map Bspline-to-Node");

            // cast aFieldSource field to discrete field
             mtk::Field_Discrete * tDiscreteField = dynamic_cast<mtk::Field_Discrete*> (aField);

             // check that dynamic cast was successful
             MORIS_ERROR( tDiscreteField != nullptr,
                     "Mapper::map_input_field_to_output_field_2 - source field need to be discrete field.\n");

            mtk::Mesh_Pair tMeshPair = tDiscreteField->get_mesh_pair();

            moris::mtk::Mesh * tInterpolationMesh = tMeshPair.get_interpolation_mesh();

            moris_index tDescritizationIndex= tDiscreteField->get_discretization_mesh_index();

            // get number of nodes on block
            uint tNumberOfNodes= tInterpolationMesh->get_num_nodes();

            Matrix< DDRMat > tNodalValues( tNumberOfNodes, 1, MORIS_REAL_MAX );

            const Matrix< DDRMat > & tCoefficients = tDiscreteField->get_coefficients();

            //-------------------------------------

            MORIS_ERROR( tInterpolationMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Source mesh is not and HMR mesh" );

            std::shared_ptr< hmr::Database > tHMRDatabase = tInterpolationMesh->get_HMR_database();

            // grab orders of meshes
            uint tSourceLagrangeOrder = tInterpolationMesh->get_order();
            uint tOrder = tDiscreteField->get_discretization_order();

            // get order of Union Mesh
            //uint tLagrangeOrder = std::max( tSourceLagrangeOrder, tOrder );

            //uint tSourcePattern = tInterpolationMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            //uint tTargetPattern = tSourceMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            //uint tPattern = 0;

            // map source Lagrange field to target Lagrange field
            if( tSourceLagrangeOrder >= tOrder )
            {

            }
//            else
//            {
//                // create union mesh. Bspline order will not be used
//                hmr::Interpolation_Mesh_HMR * tHigherOrderInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
//                        tHMRDatabase,
//                        tOrder,
//                        tSourcePattern,
//                        tOrder,
//                        tPattern); // order, Lagrange pattern, bspline order, bspline pattern
//
//                mtk::Mesh_Pair tMeshPairHigherOrder(tHigherOrderInterpolationMesh, nullptr, true);
//
//                mtk::Field_Discrete tFieldHigerOrder( tMeshPairHigherOrder, 0 );
//
//                //--------------------------------------------------
//
//                uint tUnionPattern  = tHMRDatabase->get_parameters()->get_union_pattern();
//
//                // create union pattern
//                tHMRDatabase->create_union_pattern(
//                        tSourcePattern,
//                        tSourcePattern,
//                        tUnionPattern );
//
//                // create union mesh
//                hmr::Interpolation_Mesh_HMR * tUnionInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
//                        tHMRDatabase,
//                        tOrder,
//                        tUnionPattern,
//                        tOrder,
//                        tPattern); // order, Lagrange pattern, bspline pattern
//
//                //construct union integration mesh (note: this is not ever used but is needed for mesh manager)
//                hmr::Integration_Mesh_HMR* tIntegrationUnionMesh = new hmr::Integration_Mesh_HMR(
//                        tOrder,
//                        tUnionPattern,
//                        tUnionInterpolationMesh);
//
//                mtk::Mesh_Pair tMeshPairUnion(tUnionInterpolationMesh, tIntegrationUnionMesh, true);
//
//                mtk::Field_Discrete tFieldUnion( tMeshPairUnion, 0 );
//
//                tFieldUnion.unlock_field();
//                tFieldUnion.set_coefficients( tDiscreteField->get_coefficients() );
//
//                // project field to union
//                this->perform_mapping(
//                        &tFieldUnion,
//                        EntityRank::BSPLINE,
//                        EntityRank::NODE);
//
//                //tFieldUnion.save_field_to_exodus( "Field_after1.exo");
//
//                tFieldHigerOrder.unlock_field();
//                tFieldHigerOrder.set_nodal_values( tFieldUnion.get_nodal_values() );
//
//                this->change_field_order( &tFieldHigerOrder, tDiscreteField );
//
//                return;
//            }

            //---------------------------------------------------------

            for( uint Ik = 0; Ik < tNumberOfNodes; ++Ik )
            {
                // get pointer to node
                auto tNode = &tInterpolationMesh->get_mtk_vertex( Ik );

                if ((uint) par_rank() == tInterpolationMesh->get_entity_owner(Ik, EntityRank::NODE ) )
                {
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

            tDiscreteField->unlock_field();
            tDiscreteField->set_nodal_values( tNodalValues );
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
