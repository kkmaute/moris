/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Solver_Interface.hpp
 *
 */
#pragma once

#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_FEM_Enums.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
        class SOL_Warehouse;
    }    // namespace sol

    namespace mtk
    {
        class Mesh;
    }
    namespace mdl
    {
        class Model;
    }

    namespace MSI
    {
        class MSI_Solver_Interface : public moris::Solver_Interface
        {
          private:
            moris::MSI::Model_Solver_Interface* mMSI    = nullptr;
            moris::MSI::Dof_Manager*            mDofMgn = nullptr;

            // pointer to solver warehouse
            std::shared_ptr< sol::SOL_Warehouse > mSolverWarehouse = nullptr;

            sol::Dist_Vector* mSolutionVector                = nullptr;
            sol::Dist_Vector* mPrevSolutionVector            = nullptr;
            sol::Dist_Vector* mAdjointSolutionVector         = nullptr;
            sol::Dist_Vector* mPreviousAdjointSolutionVector = nullptr;
            sol::Dist_Vector* mEigenSolutionVector           = nullptr;

            Matrix< DDRMat > mTime;
            Matrix< DDRMat > mPrevTime;

            Vector< enum MSI::Dof_Type > mListOfDofTypes;
            Vector< enum MSI::Dof_Type > mListOfSecondaryDofTypes;

            mdl::Model* mModel = nullptr;

            std::shared_ptr< Vector< real > > mEigenValues;

          public:
            MSI_Solver_Interface()
            {
                mTime     = { { 0.0 }, { 1.0 } };
                mPrevTime = { { 0.0 }, { 0.0 } };
            };

            //------------------------------------------------------------------------------
            MSI_Solver_Interface( moris::MSI::Model_Solver_Interface* aMSI )
                    : mMSI( aMSI )
                    , mDofMgn( mMSI->get_dof_manager() )
            {
                mTime     = { { 0.0 }, { 1.0 } };
                mPrevTime = { { 0.0 }, { 0.0 } };

                mMSI->set_solver_interface( this );
            };

            //------------------------------------------------------------------------------

            ~MSI_Solver_Interface() override {};

            //------------------------------------------------------------------------------

            void set_solution_vector( sol::Dist_Vector* aSolutionVector ) override;

            //------------------------------------------------------------------------------

            void set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector ) override;

            //------------------------------------------------------------------------------

            void set_eigen_solution_vector( sol::Dist_Vector* aSolutionVector ) override;

            //--------------------------------------------------------------------------------

            void set_eigen_values( const std::shared_ptr< Vector< real > >& aEigenValues ) override;

            //------------------------------------------------------------------------------

            void set_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector ) override;

            //------------------------------------------------------------------------------

            void set_previous_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector ) override;

            //------------------------------------------------------------------------------

            /**
             * @brief get the solution vector for a certain dof type
             *
             * @param aListOfDofTypes
             * @return sol::Dist_Vector*
             */

            sol::Dist_Vector*
            get_solution_vector(
                    const Vector< enum MSI::Dof_Type >& aListOfDofTypes,
                    Vector< moris_index > const &       aLocalCoefficientsIndices );

            //------------------------------------------------------------------------------

            sol::Dist_Vector*
            get_solution_vector_prev_time_step() override;

            //------------------------------------------------------------------------------
            /**
             * get previous solution vector
             * @param[ out ] aSolutionVector previous distributed solution vector
             */
            sol::Dist_Vector*
            get_eigen_solution_vector() override
            {
                return mEigenSolutionVector;
            }

            //------------------------------------------------------------------------------

            /**
             * get previous solution vector
             * @param[ out ] aSolutionVector previous distributed solution vector
             */
            std::shared_ptr< Vector< real > >&
            get_eigen_values() override
            {
                return mEigenValues;
            }

            //------------------------------------------------------------------------------

            void
            set_model( mdl::Model* aModel )
            {
                mModel = aModel;
            }

            void postmultiply_implicit_dQds() override;

            void compute_IQI() override;

            //------------------------------------------------------------------------------

            void get_exact_solution_from_hdf5_and_calculate_error( const char* aFilename );

            //------------------------------------------------------------------------------

            void get_residual_vector_for_output( const char* aFilename );

            //------------------------------------------------------------------------------

            void write_solution_to_hdf5_file( const char* aFilename );

            //------------------------------------------------------------------------------

            void set_time( const Matrix< DDRMat >& aTime ) override;

            //------------------------------------------------------------------------------

            Matrix< DDRMat > get_time() override;

            //------------------------------------------------------------------------------

            void set_previous_time( const Matrix< DDRMat >& aTime ) override;

            //------------------------------------------------------------------------------

            Matrix< DDRMat > get_previous_time() override;

            //------------------------------------------------------------------------------

            void free_block_memory( const uint aMyEquSetInd ) override;

            //------------------------------------------------------------------------------

            void initialize_set(
                    const uint                             aMyEquSetInd,
                    const bool                             aIsStaggered,
                    const moris::fem::Time_Continuity_Flag aTimeContinuityOnlyFlag               = moris::fem::Time_Continuity_Flag::DEFAULT,
                    const bool                             aIsAdjointOffDiagonalTimeContribution = false ) override;

            //------------------------------------------------------------------------------

            void update_model() override;

            //------------------------------------------------------------------------------

            void report_beginning_of_assembly() override;

            //------------------------------------------------------------------------------

            void report_end_of_assembly() override;

            //------------------------------------------------------------------------------

            void initiate_output(
                    const uint aOutputIndex,
                    const real aTime,
                    const bool aEndOfTimeIteration ) override;

            /**
             * Updates the underlying problem via the model.
             */
            void update_problem() override;

            //------------------------------------------------------------------------------

            void
            set_requested_dof_types( const Vector< enum MSI::Dof_Type >& aListOfDofTypes ) override
            {
                mListOfDofTypes = aListOfDofTypes;
            };

            //------------------------------------------------------------------------------

            void
            set_secondary_dof_types( const Vector< enum MSI::Dof_Type >& aListOfDofTypes ) override
            {
                mListOfSecondaryDofTypes = aListOfDofTypes;
            };

            //------------------------------------------------------------------------------

            const Vector< enum MSI::Dof_Type >&
            get_requested_dof_types() override
            {
                return mListOfDofTypes;
            };

            //------------------------------------------------------------------------------

            const Vector< enum MSI::Dof_Type >&
            get_secondary_dof_types()
            {
                return mListOfSecondaryDofTypes;
            };

            //------------------------------------------------------------------------------

            void
            set_time_levels_for_type(
                    const enum Dof_Type aDofType,
                    const moris::uint   aNumTimeLevels ) override
            {
                mDofMgn->set_time_levels_for_type( aDofType, aNumTimeLevels );
            };

            //------------------------------------------------------------------------------

            // number of elements blocks on proc
            moris::uint
            get_num_sets() override
            {
                return mMSI->get_num_equation_sets();
            };

            //------------------------------------------------------------------------------

            // number of elements on proc
            moris::uint
            get_num_my_elements() override
            {
                return mMSI->get_num_eqn_objs();
            };

            //------------------------------------------------------------------------------

            moris::uint
            get_num_equation_objects_on_set( uint aMyEquSetInd ) override
            {
                return mMSI->get_equation_set( aMyEquSetInd )->get_num_equation_objects();
            };

            //------------------------------------------------------------------------------

            fem::Element_Type
            get_set_type( uint aMyEquSetInd ) override
            {
                return mMSI->get_equation_set( aMyEquSetInd )->get_element_type();
            };

            //------------------------------------------------------------------------------

            moris::uint
            get_num_my_dofs() override
            {
                return mDofMgn->get_num_owned_adofs();
            };

            //------------------------------------------------------------------------------

            moris::uint get_num_rhs() override;

            //------------------------------------------------------------------------------

            moris::uint
            get_num_eigen_vectors() override
            {
                return mMSI->get_num_eigen_vectors();
            }

            //------------------------------------------------------------------------------

            moris::uint
            get_max_num_global_dofs() override
            {
                moris::uint tNumMyDofs = mDofMgn->get_num_owned_adofs();

                // sum up all distributed dofs
                moris::uint tMaxNumGlobalDofs = sum_all( tNumMyDofs );

                return tMaxNumGlobalDofs;
            };

            //------------------------------------------------------------------------------
            // local-to-global map
            moris::Matrix< DDSMat >
            get_my_local_global_map() override
            {
                return mDofMgn->get_local_adof_ids();
            };

            //------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_my_local_global_map( const Vector< enum Dof_Type >& aListOfDofTypes ) override
            {
                return mDofMgn->get_local_adof_ids( aListOfDofTypes );
            };

            //------------------------------------------------------------------------------
            Matrix< DDSMat >
            get_my_local_global_overlapping_map() override
            {
                return mDofMgn->get_local_overlapping_adof_ids();
            };

            //------------------------------------------------------------------------------

            Matrix< DDSMat >
            get_my_local_global_overlapping_map( const Vector< enum Dof_Type >& aListOfDofTypes ) override
            {
                return mDofMgn->get_local_overlapping_adof_ids( aListOfDofTypes );
            };

            //------------------------------------------------------------------------------
            void
            get_element_topology(
                    const moris::uint& aMyElementInd,
                    Matrix< DDSMat >&  aElementTopology ) override
            {
                mMSI->get_eqn_obj( aMyElementInd )
                        ->get_equation_obj_dof_ids( aElementTopology );
            };

            //------------------------------------------------------------------------------
            void
            get_element_topology(
                    const moris::uint& aMyEquSetInd,
                    const moris::uint& aMyElementInd,
                    Matrix< DDSMat >&  aElementTopology ) override
            {
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_equation_obj_dof_ids( aElementTopology );
            };

            //------------------------------------------------------------------------------

            Matrix< DDUMat >
            get_constrained_Ids() override
            {
                // Matrix< DDUMat > tLocalConstrIds;// = mDofMgn->get_full_to_free_constraints();
                return Matrix< DDUMat >( 0, 0 );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_operator(
                    const moris::uint& aMyElementInd,
                    Matrix< DDRMat >&  aElementMatrix ) override
            {
                // mMSI->get_eqn_obj( aMyElementInd )->set_time( mTime );
                mMSI->get_eqn_obj( aMyElementInd )
                        ->get_egn_obj_jacobian( aElementMatrix );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_operator(
                    const moris::uint& aMyEquSetInd,
                    const moris::uint& aMyElementInd,
                    Matrix< DDRMat >&  aElementMatrix ) override
            {
                // mMSI->get_equation_set( aMyEquSetInd )->get_equation_object_list()( aMyElementInd )->set_time( mTime );
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_egn_obj_jacobian( aElementMatrix );
            };
            //------------------------------------------------------------------------------

            void
            get_equation_object_rhs(
                    const moris::uint&          aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override
            {
                mMSI->get_eqn_obj( aMyElementInd )
                        ->get_equation_obj_residual( aElementRHS );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override
            {
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_equation_obj_residual( aElementRHS );
            };

            void
            get_equation_object_staggered_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override
            {
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_staggered_equation_obj_residual( aElementRHS );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_off_diag_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override
            {
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_equation_obj_off_diagonal_residual( aElementRHS );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_operator_and_rhs(
                    const moris::uint&          aMyElementInd,
                    Matrix< DDRMat >&           aElementMatrix,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override
            {
                // mMSI->get_eqn_obj( aMyElementInd )->set_time( mTime );
                mMSI->get_eqn_obj( aMyElementInd )
                        ->get_egn_obj_jacobian_and_residual( aElementMatrix, aElementRHS );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_operator_and_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Matrix< DDRMat >&           aElementMatrix,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override
            {
                // mMSI->get_equation_set( aMyEquSetInd )->get_equation_object_list()( aMyElementInd )->set_time( mTime );
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_egn_obj_jacobian_and_residual( aElementMatrix, aElementRHS );
            };

            //------------------------------------------------------------------------------

            mtk::Mesh*
            get_mesh_pointer_for_multigrid() override
            {
                return mMSI->get_mesh_pointer_for_multigrid();
            };

            //------------------------------------------------------------------------------

            void
            read_multigrid_maps(
                    const moris::uint              aLevel,
                    const moris::Matrix< DDSMat >& aExtFineIndices,
                    const moris::sint              aTypeTimeIdentifier,
                    moris::Matrix< DDSMat >&       aInternalFineIndices ) override
            {
                mMSI->get_msi_multigrid_pointer()->read_multigrid_maps( aLevel, aExtFineIndices, aTypeTimeIdentifier, aInternalFineIndices );
            };

            //------------------------------------------------------------------------------

            const Vector< Matrix< DDSMat > >&
            get_lists_of_multigrid_identifiers() override
            {
                return mMSI->get_msi_multigrid_pointer()->get_lists_of_multigrid_identifiers();
            };

            //------------------------------------------------------------------------------

            const Vector< Matrix< DDUMat > >&
            get_lists_of_ext_index_multigrid() override
            {
                return mMSI->get_msi_multigrid_pointer()->get_lists_of_ext_index_multigrid();
            };

            //------------------------------------------------------------------------------

            const Vector< Vector< Matrix< DDSMat > > >&
            get_multigrid_map() override
            {
                return mMSI->get_msi_multigrid_pointer()->get_multigrid_map();
            };

            //------------------------------------------------------------------------------

            const moris::Matrix< DDUMat >&
            get_number_remaining_dofs() override
            {
                return mMSI->get_msi_multigrid_pointer()->get_number_remaining_dofs();
            };

            //------------------------------------------------------------------------------

            const Matrix< DDSMat >&
            get_type_time_identifier_to_type_map() override
            {
                return mDofMgn->get_typetime_identifier_to_type_map();
            };

            //------------------------------------------------------------------------------

            moris::sint
            get_adof_index_for_type( moris::uint aDofType ) override
            {
                return mMSI->get_adof_index_for_type( aDofType );
            };

            //------------------------------------------------------------------------------

            void
            calculate_criteria(
                    const moris::uint& aMySetInd,
                    const moris::uint& aMyEquationObjectInd ) override
            {
                mMSI->get_equation_set( aMySetInd )->get_equation_object_list()( aMyEquationObjectInd )->compute_QI();
            };

            //------------------------------------------------------------------------------

            const Vector< moris::Matrix< DDRMat > >& get_criteria( const moris::uint& aMySetInd ) override;

            //------------------------------------------------------------------------------

            void set_requested_IQI_names( const Vector< std::string >& aIQINames ) override;

            //------------------------------------------------------------------------------

            /**
             * @brief Set the solver warehouse object store a pointer of the warehouse
             *
             * @param aSolverWarehouse
             */
            void
            set_solver_warehouse( const std::shared_ptr< sol::SOL_Warehouse >& aSolverWarehouse ) override;

            //------------------------------------------------------------------------------

            /**
             * @brief compute number of zeros id diogonal and off disgonal part
             * of the matrix based on adof ids
             *
             */

            void compute_sparsity_pattern() override;

            //------------------------------------------------------------------------------

            /**
             * @brief initial estimate size of the stencil based on the spatial dimension and
             *
             * @return uint
             */

            virtual uint estimate_number_of_nonzero_columns();

            //-------------------------------------------------------------------------------

            /**
             * @brief This function receives the shared adof connectivity
             *
             * @param [in] aSharedAdofConn outer cell: neighbor processor index, inner cell: shared aof index , most inner cell: list of adofs connected to the adof
             * @param [out] aAdofConnectivityReceive outer cell: received from neighbor processor index, inner cells: list of all adofs that are connoted to shared adofs consecutively
             * @param [out] aAdofConnectivityOffsetReceive outer cell: received from neighbor processor index, inner cells: offset indicating position of the connected adofs
             */

            void communicate_shared_adof_connectivity(
                    Vector< Vector< Vector< uint > > > const & aSharedAdofConn,
                    Vector< Vector< uint > >&                  aAdofConnectivityReceive,
                    Vector< Vector< uint > >&                  aAdofConnectivityOffsetReceive,
                    Vector< moris_index > const &              aCommCell );
        };
    }    // namespace MSI
}    // namespace moris