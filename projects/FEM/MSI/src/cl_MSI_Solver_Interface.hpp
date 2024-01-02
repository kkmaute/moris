/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Solver_Interface.hpp
 *
 */
#ifndef SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_
#define SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_

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

            virtual ~MSI_Solver_Interface(){};

            //------------------------------------------------------------------------------

            void set_solution_vector( sol::Dist_Vector* aSolutionVector );

            //------------------------------------------------------------------------------

            void set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector );

            //------------------------------------------------------------------------------

            void set_eigen_solution_vector( sol::Dist_Vector* aSolutionVector );

            //------------------------------------------------------------------------------

            void set_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector );

            //------------------------------------------------------------------------------

            void set_previous_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector );

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
            get_solution_vector_prev_time_step();

            //------------------------------------------------------------------------------
            /**
             * get previous solution vector
             * @param[ out ] aSolutionVector previous distributed solution vector
             */
            sol::Dist_Vector*
            get_eigen_solution_vector()
            {
                return mEigenSolutionVector;
            }

            //------------------------------------------------------------------------------

            void
            set_model( mdl::Model* aModel )
            {
                mModel = aModel;
            }

            void postmultiply_implicit_dQds();

            void compute_IQI();

            //------------------------------------------------------------------------------

            void get_exact_solution_from_hdf5_and_calculate_error( const char* aFilename );

            //------------------------------------------------------------------------------

            void get_residual_vector_for_output( const char* aFilename );

            //------------------------------------------------------------------------------

            void write_solution_to_hdf5_file( const char* aFilename );

            //------------------------------------------------------------------------------

            void set_time( const Matrix< DDRMat >& aTime );

            //------------------------------------------------------------------------------

            Matrix< DDRMat > get_time();

            //------------------------------------------------------------------------------

            void set_previous_time( const Matrix< DDRMat >& aTime );

            //------------------------------------------------------------------------------

            Matrix< DDRMat > get_previous_time();

            //------------------------------------------------------------------------------

            void free_block_memory( const uint aMyEquSetInd );

            //------------------------------------------------------------------------------

            void initialize_set(
                    const uint                             aMyEquSetInd,
                    const bool                             aIsStaggered,
                    const moris::fem::Time_Continuity_Flag aTimeContinuityOnlyFlag               = moris::fem::Time_Continuity_Flag::DEFAULT,
                    const bool                             aIsAdjointOffDiagonalTimeContribution = false );

            //------------------------------------------------------------------------------

            void report_beginning_of_assembly();

            //------------------------------------------------------------------------------

            void report_end_of_assembly();

            //------------------------------------------------------------------------------

            void initiate_output(
                    const uint aOutputIndex,
                    const real aTime,
                    const bool aEndOfTimeIteration );

            //------------------------------------------------------------------------------

            void
            set_requested_dof_types( const Vector< enum MSI::Dof_Type > aListOfDofTypes )
            {
                mListOfDofTypes = aListOfDofTypes;
            };

            //------------------------------------------------------------------------------

            void
            set_secondary_dof_types( const Vector< enum MSI::Dof_Type > aListOfDofTypes )
            {
                mListOfSecondaryDofTypes = aListOfDofTypes;
            };

            //------------------------------------------------------------------------------

            const Vector< enum MSI::Dof_Type >&
            get_requested_dof_types()
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
                    const moris::uint   aNumTimeLevels )
            {
                mDofMgn->set_time_levels_for_type( aDofType, aNumTimeLevels );
            };

            //------------------------------------------------------------------------------

            // number of elements blocks on proc
            moris::uint
            get_num_my_blocks()
            {
                return mMSI->get_num_eqn_blocks();
            };

            //------------------------------------------------------------------------------

            // number of elements on proc
            moris::uint
            get_num_my_elements()
            {
                return mMSI->get_num_eqn_objs();
            };

            //------------------------------------------------------------------------------

            moris::uint
            get_num_equation_objects_on_set( uint aMyEquSetInd )
            {
                return mMSI->get_equation_set( aMyEquSetInd )->get_num_equation_objects();
            };

            //------------------------------------------------------------------------------

            enum fem::Element_Type
            get_set_type( uint aMyEquSetInd )
            {
                return mMSI->get_equation_set( aMyEquSetInd )->get_element_type();
            };

            //------------------------------------------------------------------------------

            moris::uint
            get_num_my_dofs()
            {
                return mDofMgn->get_num_owned_adofs();
            };

            //------------------------------------------------------------------------------

            moris::uint get_num_rhs();

            //------------------------------------------------------------------------------

            moris::uint
            get_num_eigen_vectors()
            {
                return mMSI->get_num_eigen_vectors();
            }

            //------------------------------------------------------------------------------

            moris::uint
            get_max_num_global_dofs()
            {
                moris::uint tNumMyDofs = mDofMgn->get_num_owned_adofs();

                // sum up all distributed dofs
                moris::uint tMaxNumGlobalDofs = sum_all( tNumMyDofs );

                return tMaxNumGlobalDofs;
            };

            //------------------------------------------------------------------------------
            // local-to-global map
            moris::Matrix< DDSMat >
            get_my_local_global_map()
            {
                return mDofMgn->get_local_adof_ids();
            };

            //------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_my_local_global_map( const Vector< enum Dof_Type >& aListOfDofTypes )
            {
                return mDofMgn->get_local_adof_ids( aListOfDofTypes );
            };

            //------------------------------------------------------------------------------
            Matrix< DDSMat >
            get_my_local_global_overlapping_map()
            {
                return mDofMgn->get_local_overlapping_adof_ids();
            };

            //------------------------------------------------------------------------------

            Matrix< DDSMat >
            get_my_local_global_overlapping_map( const Vector< enum Dof_Type >& aListOfDofTypes )
            {
                return mDofMgn->get_local_overlapping_adof_ids( aListOfDofTypes );
            };

            //------------------------------------------------------------------------------
            void
            get_element_topology(
                    const moris::uint& aMyElementInd,
                    Matrix< DDSMat >&  aElementTopology )
            {
                mMSI->get_eqn_obj( aMyElementInd )
                        ->get_equation_obj_dof_ids( aElementTopology );
            };

            //------------------------------------------------------------------------------
            void
            get_element_topology(
                    const moris::uint& aMyEquSetInd,
                    const moris::uint& aMyElementInd,
                    Matrix< DDSMat >&  aElementTopology )
            {
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_equation_obj_dof_ids( aElementTopology );
            };

            //------------------------------------------------------------------------------

            Matrix< DDUMat >
            get_constrained_Ids()
            {
                // Matrix< DDUMat > tLocalConstrIds;// = mDofMgn->get_full_to_free_constraints();
                return Matrix< DDUMat >( 0, 0 );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_operator(
                    const moris::uint& aMyElementInd,
                    Matrix< DDRMat >&  aElementMatrix )
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
                    Matrix< DDRMat >&  aElementMatrix )
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
                    Vector< Matrix< DDRMat > >& aElementRHS )
            {
                mMSI->get_eqn_obj( aMyElementInd )
                        ->get_equation_obj_residual( aElementRHS );
            };

            //------------------------------------------------------------------------------

            void
            get_equation_object_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS )
            {
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_equation_obj_residual( aElementRHS );
            };

            void
            get_equation_object_staggered_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS )
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
                    Vector< Matrix< DDRMat > >& aElementRHS )
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
                    Vector< Matrix< DDRMat > >& aElementRHS )
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
                    Vector< Matrix< DDRMat > >& aElementRHS )
            {
                // mMSI->get_equation_set( aMyEquSetInd )->get_equation_object_list()( aMyElementInd )->set_time( mTime );
                mMSI->get_equation_set( aMyEquSetInd )
                        ->get_equation_object_list()( aMyElementInd )
                        ->get_egn_obj_jacobian_and_residual( aElementMatrix, aElementRHS );
            };

            //------------------------------------------------------------------------------

            mtk::Mesh*
            get_mesh_pointer_for_multigrid()
            {
                return mMSI->get_mesh_pointer_for_multigrid();
            };

            //------------------------------------------------------------------------------

            void
            read_multigrid_maps(
                    const moris::uint              aLevel,
                    const moris::Matrix< DDSMat >& aExtFineIndices,
                    const moris::sint              aTypeTimeIdentifier,
                    moris::Matrix< DDSMat >&       aInternalFineIndices )
            {
                mMSI->get_msi_multigrid_pointer()->read_multigrid_maps( aLevel, aExtFineIndices, aTypeTimeIdentifier, aInternalFineIndices );
            };

            //------------------------------------------------------------------------------

            const Vector< Matrix< DDSMat > >&
            get_lists_of_multigrid_identifiers()
            {
                return mMSI->get_msi_multigrid_pointer()->get_lists_of_multigrid_identifiers();
            };

            //------------------------------------------------------------------------------

            const Vector< Matrix< DDUMat > >&
            get_lists_of_ext_index_multigrid()
            {
                return mMSI->get_msi_multigrid_pointer()->get_lists_of_ext_index_multigrid();
            };

            //------------------------------------------------------------------------------

            const Vector< Vector< Matrix< DDSMat > > >&
            get_multigrid_map()
            {
                return mMSI->get_msi_multigrid_pointer()->get_multigrid_map();
            };

            //------------------------------------------------------------------------------

            const moris::Matrix< DDUMat >&
            get_number_remaining_dofs()
            {
                return mMSI->get_msi_multigrid_pointer()->get_number_remaining_dofs();
            };

            //------------------------------------------------------------------------------

            const Matrix< DDSMat >&
            get_type_time_identifier_to_type_map()
            {
                return mDofMgn->get_typetime_identifier_to_type_map();
            };

            //------------------------------------------------------------------------------

            moris::sint
            get_adof_index_for_type( moris::uint aDofType )
            {
                return mMSI->get_adof_index_for_type( aDofType );
            };

            //------------------------------------------------------------------------------

            void
            calculate_criteria(
                    const moris::uint& aMySetInd,
                    const moris::uint& aMyEquationObjectInd )
            {
                mMSI->get_equation_set( aMySetInd )->get_equation_object_list()( aMyEquationObjectInd )->compute_QI();
            };

            //------------------------------------------------------------------------------

            const Vector< moris::Matrix< DDRMat > >& get_criteria( const moris::uint& aMySetInd );

            //------------------------------------------------------------------------------

            void set_requested_IQI_names( const Vector< std::string >& aIQINames );

            //------------------------------------------------------------------------------

            /**
             * @brief Set the solver warehouse object store a pointer of the warehouse
             *
             * @param aSolverWarehouse
             */
            void
            set_solver_warehouse( std::shared_ptr< sol::SOL_Warehouse > aSolverWarehouse );


            //------------------------------------------------------------------------------

            /**
             * @brief compute number of zeros id diogonal and off disgonal part
             * of the matrix based on adof ids
             *
             */

            virtual void compute_sparsity_pattern() override;

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

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
