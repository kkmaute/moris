/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Solver_Interface.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_
#define SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_

#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Enums.hpp"

#include "cl_DLA_Geometric_Multigrid.hpp"

#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
        class Dist_Matrix;
        class SOL_Warehouse;
    }    // namespace sol

    namespace mtk
    {
        class Mesh;
    }

    namespace MSI
    {
        enum class Dof_Type;
    }

    class Solver_Interface
    {
      private:
        // Dummy member variable

        bool mIsForwardAnalysis            = true;
        bool mIsAdjointSensitivityAnalysis = true;

      protected:
        Vector< moris_id > mNonZeroDigonal;
        Vector< moris_id > mNonZeroOffDigonal;

      public:
        /** Destructor */
        virtual ~Solver_Interface() {};

        //------------------------------------------------------------------------------
        /**
         * indicated that this equation model is used for the sensitivity analysis
         */
        void
        set_sensitivity_analysis_type( bool tIsAdjointSensitivityAnalysis )
        {
            mIsForwardAnalysis            = false;
            mIsAdjointSensitivityAnalysis = tIsAdjointSensitivityAnalysis;
        };

        //------------------------------------------------------------------------------
        /**
         * indicated that this equation model is used for the forward analysis
         */
        void
        set_is_forward_analysis()
        {
            mIsForwardAnalysis = true;
        };

        //------------------------------------------------------------------------------
        /**
         * returns if this is the a forward analysis
         * @param[ out ] mIsForwardAnalysis
         */
        bool
        is_forward_analysis()
        {
            return mIsForwardAnalysis;
        };

        //------------------------------------------------------------------------------

        virtual void
        postmultiply_implicit_dQds()
        {
            MORIS_ERROR( false, "Solver_Interface::postmultiply_implicit_dQds: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        compute_IQI()
        {
            MORIS_ERROR( false, "Solver_Interface::compute_IQI: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            MORIS_ERROR( false, "Solver_Interface::set_solution_vector: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_eigen_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            MORIS_ERROR( false, "Solver_Interface::set_eigen_solution_vector: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_eigen_values( const std::shared_ptr< Vector< real > >& aEigenValues )
        {
            MORIS_ERROR( false, "Solver_Interface::set_eigen_values: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            MORIS_ERROR( false, "Solver_Interface::set_adjoint_solution_vector: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_previous_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            MORIS_ERROR( false, "Solver_Interface::set_adjoint_solution_vector: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_time_levels_for_type(
                const enum MSI::Dof_Type aDofType,
                const moris::uint        aNumTimeLevels )
        {
            MORIS_ERROR( false, "Solver_Interface::set_time_levels_for_type: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector )
        {
            MORIS_ERROR( false, "Solver_Interface::set_solution_vector_prev_time_step: not set." );
        };

        //------------------------------------------------------------------------------

        virtual sol::Dist_Vector*
        get_solution_vector_prev_time_step()
        {
            MORIS_ERROR( false, "Solver_Interface::get_solution_vector_prev_time_step: not set." );
            return nullptr;
        }

        //------------------------------------------------------------------------------

        virtual sol::Dist_Vector*
        get_eigen_solution_vector()
        {
            MORIS_ERROR( false, "Solver_Interface::get_eigen_solution_vector: not set." );
            return nullptr;
        }

        //------------------------------------------------------------------------------

        virtual std::shared_ptr< Vector< real > >&
        get_eigen_values()
        {
            MORIS_ERROR( false, "Solver_Interface::get_eigen_solution_vector: not set." );
            return *( new std::shared_ptr< Vector< real > >() );
        }

        //------------------------------------------------------------------------------

        virtual void
        set_time( const Matrix< DDRMat >& aTime )
        {
            MORIS_ERROR( false, "Solver_Interface::set_time: not set." );
        };

        //------------------------------------------------------------------------------

        virtual Matrix< DDRMat >
        get_time()
        {
            MORIS_ERROR( false, "Solver_Interface::get_time: not set." );
            return Matrix< DDRMat >( 0, 0 );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_previous_time( const Matrix< DDRMat >& aTime )
        {
            MORIS_ERROR( false, "Solver_Interface::set_previous_time: not set." );
        };

        //------------------------------------------------------------------------------

        virtual Matrix< DDRMat >
        get_previous_time()
        {
            MORIS_ERROR( false, "Solver_Interface::get_previous_time: not set." );
            return Matrix< DDRMat >( 0, 0 );
        };

        //------------------------------------------------------------------------------

        virtual void set_residual_norm( const real& aResNorm ) {
            // MORIS_ERROR( false, "Solver_Interface::set_residual_norm: not set.");
        };

        //------------------------------------------------------------------------------

        virtual void free_block_memory( const uint aBlockInd ) = 0;

        //------------------------------------------------------------------------------

        virtual void initialize_set(
                const uint                             aBlockInd,
                const bool                             aIsStaggered                          = false,
                const moris::fem::Time_Continuity_Flag aTimeContinuityOnlyFlag               = moris::fem::Time_Continuity_Flag::DEFAULT,
                const bool                             aIsAdjointOffDiagonalTimeContribution = false ) {};

        //------------------------------------------------------------------------------

        virtual void update_model() {};

        //------------------------------------------------------------------------------

        virtual void report_beginning_of_assembly() {};

        //------------------------------------------------------------------------------

        virtual void report_end_of_assembly() {};

        //------------------------------------------------------------------------------

        virtual void
        set_requested_dof_types( const Vector< enum MSI::Dof_Type >& aListOfDofTypes )
        {
            MORIS_ERROR( false, "Solver_Interface::set_requested_dof_types: not set." );
        };

        //------------------------------------------------------------------------------

        virtual void
        set_secondary_dof_types( const Vector< enum MSI::Dof_Type >& aListOfDofTypes )
        {
            MORIS_ERROR( false, "Solver_Interface::set_secondary_dof_types: not set." );
        };

        //------------------------------------------------------------------------------

        virtual const Vector< enum MSI::Dof_Type >&
        get_requested_dof_types()
        {
            MORIS_ERROR( false, "Solver_Interface::get_requested_dof_types: not set." );
            return *( new Vector< enum MSI::Dof_Type >() );
        };

        //------------------------------------------------------------------------------

        virtual void
        initiate_output(
                const uint aOutputIndex,
                const real aTime,
                const bool aEndOfTimeIteration )
        {
            MORIS_ERROR( false, "Solver_Interface::initiate_output: not set." );
        };

        //------------------------------------------------------------------------------

        // local dimension of the problem
        virtual moris::uint get_num_my_dofs() = 0;

        //------------------------------------------------------------------------------

        virtual moris::uint get_num_rhs() = 0;

        //------------------------------------------------------------------------------

        virtual moris::uint get_num_eigen_vectors() = 0;

        //------------------------------------------------------------------------------

        virtual uint get_max_num_global_dofs() = 0;

        //------------------------------------------------------------------------------

        // number local elements blocks
        virtual moris::uint get_num_sets() = 0;

        //------------------------------------------------------------------------------

        virtual moris::uint get_num_equation_objects_on_set( uint aBlockInd ) = 0;

        //------------------------------------------------------------------------------

        virtual fem::Element_Type
        get_set_type( uint aMyEquSetInd )
        {
            MORIS_ERROR( false, "Solver_Interface::initiate_output: not set." );
            return fem::Element_Type::UNDEFINED;
        };

        //------------------------------------------------------------------------------

        // number local elements
        virtual moris::uint get_num_my_elements() = 0;

        //------------------------------------------------------------------------------

        // local-to-global map // FIXME pass return value in as reference
        virtual moris::Matrix< DDSMat > get_my_local_global_map() = 0;

        //------------------------------------------------------------------------------

        // FIXME pass return value in as reference
        virtual moris::Matrix< DDSMat >
        get_my_local_global_map( const Vector< enum MSI::Dof_Type >& aListOfDofTypes )
        {
            MORIS_ERROR( false, "Solver_Interface::get_my_local_global_map: not set." );
            return Matrix< DDSMat >( 0, 0 );
        }

        //------------------------------------------------------------------------------

        // FIXME pass return value in as reference
        virtual moris::Matrix< DDSMat >
        get_my_local_global_overlapping_map()
        {
            // MORIS_ERROR( false, "Solver_Interface::get_my_local_global_overlapping_map(): Virtual class not overwritten" );
            return Matrix< DDSMat >( 0, 0 );
        };

        //------------------------------------------------------------------------------

        /**
         * @brief Get the my local global overlapping map of the adofs
         *
         * @param aListOfDofTypes
         * @return moris::Matrix< DDSMat >
         */

        virtual moris::Matrix< DDSMat >
        get_my_local_global_overlapping_map( const Vector< enum MSI::Dof_Type >& aListOfDofTypes )
        {
            MORIS_ERROR( false, "Solver_Interface::get_my_local_global_overlapping_map(): Virtual class not overwritten" );
            return Matrix< DDSMat >( 0, 0 );
        };

        //------------------------------------------------------------------------------

        virtual moris::Matrix< DDUMat > get_constrained_Ids() = 0;

        //------------------------------------------------------------------------------

        virtual void get_equation_object_operator(
                const moris::uint&       aMyElementInd,
                moris::Matrix< DDRMat >& aElementMatrix ) = 0;

        //------------------------------------------------------------------------------

        virtual void get_equation_object_operator(
                const moris::uint&       aMyBlockInd,
                const moris::uint&       aMyElementInd,
                moris::Matrix< DDRMat >& aElementMatrix ) = 0;

        //------------------------------------------------------------------------------

        virtual void get_element_topology(
                const moris::uint&       aMyElementInd,
                moris::Matrix< DDSMat >& aElementTopology ) = 0;

        //------------------------------------------------------------------------------

        virtual void get_element_topology(
                const moris::uint&       aMyBlockInd,
                const moris::uint&       aMyElementInd,
                moris::Matrix< DDSMat >& aElementTopology ) = 0;

        //------------------------------------------------------------------------------

        virtual void get_equation_object_rhs(
                const moris::uint&          aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS ) = 0;

        //------------------------------------------------------------------------------

        virtual void get_equation_object_rhs(
                const moris::uint&          aMyBlockInd,
                const moris::uint&          aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS ) = 0;

        //------------------------------------------------------------------------------

        // fixme move that into get_equation_object_rhs()
        virtual void
        get_equation_object_off_diag_rhs(
                const moris::uint&          aMyBlockInd,
                const moris::uint&          aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS )
        {
            MORIS_ERROR( false, "not implemented" );
        };

        //------------------------------------------------------------------------------

        virtual void
        get_equation_object_staggered_rhs(
                const moris::uint&          aMyEquSetInd,
                const moris::uint&          aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS )
        {
            MORIS_ERROR( false, "not implemented" );
        };

        //------------------------------------------------------------------------------

        virtual void get_equation_object_operator_and_rhs(
                const moris::uint&          aMyElementInd,
                moris::Matrix< DDRMat >&    aElementMatrix,
                Vector< Matrix< DDRMat > >& aElementRHS ) = 0;

        //------------------------------------------------------------------------------

        virtual void get_equation_object_operator_and_rhs(
                const moris::uint&          aMyBlockInd,
                const moris::uint&          aMyElementInd,
                moris::Matrix< DDRMat >&    aElementMatrix,
                Vector< Matrix< DDRMat > >& aElementRHS ) = 0;

        //------------------------------------------------------------

        virtual void
        set_time_value(
                const moris::real& aLambda,
                moris::uint        aPos )
        {
            MORIS_ERROR( false, "Solver_Interface::set_time_value: not set." );
        }

        //------------------------------------------------------------

        virtual void
        use_matrix_market_files()
        {
            MORIS_ERROR( false, "error in use_matrix_market_files" );
        }

        //------------------------------------------------------------------------------

        virtual const char*
        get_matrix_market_path()
        {
            // assert(0);
            return nullptr;
        }

        //------------------------------------------------------------------------------

        virtual mtk::Mesh*
        get_mesh_pointer_for_multigrid()
        {
            MORIS_ERROR( false, "Solver_Interface::get_mesh_pointer_for_multigrid, Only works with MSI and multigrid" );
            return nullptr;
        };

        //------------------------------------------------------------------------------

        virtual void
        read_multigrid_maps(
                const moris::uint              aLevel,
                const moris::Matrix< DDSMat >& aExtFineIndices,
                const moris::sint              aTypeTimeIdentifier,
                moris::Matrix< DDSMat >&       aInternalFineIndices )
        {
            MORIS_ERROR( false, "Solver_Interface::read_multigrid_maps, Only works with MSI and multigrid" );
        };

        //------------------------------------------------------------------------------

        virtual const Vector< Matrix< DDUMat > >&
        get_lists_of_ext_index_multigrid()
        {
            MORIS_ERROR( false, "Solver_Interface::get_lists_of_ext_index_multigrid, Only works with MSI and multigrid" );
            return *( new Vector< Matrix< DDUMat > >() );
        };

        //------------------------------------------------------------------------------

        virtual const Vector< Matrix< DDSMat > >&
        get_lists_of_multigrid_identifiers()
        {
            MORIS_ERROR( false, "Solver_Interface::get_lists_of_ext_index_multigrid, Only works with MSI and multigrid" );
            return *( new Vector< Matrix< DDSMat > >() );
        };

        //------------------------------------------------------------------------------

        virtual const Vector< Vector< Matrix< DDSMat > > >&
        get_multigrid_map()
        {
            MORIS_ERROR( false, "Solver_Interface::get_multigrid_map, Only works with MSI and multigrid" );
            return *( new Vector< Vector< Matrix< DDSMat > > >() );
        };

        //------------------------------------------------------------------------------

        virtual const moris::Matrix< DDUMat >&
        get_number_remaining_dofs()
        {
            MORIS_ERROR( false, "Solver_Interface::get_number_remaining_dofs, Only works with MSI and multigrid" );
            return *( new Matrix< DDUMat >() );
        };

        //------------------------------------------------------------------------------

        virtual const Matrix< DDSMat >&
        get_type_time_identifier_to_type_map()
        {
            MORIS_ERROR( false, "Solver_Interface::get_type_time_identifier_to_type_map, Only works with MSI and multigrid" );
            return *( new Matrix< DDSMat >() );
        };

        //------------------------------------------------------------------------------

        virtual moris::sint
        get_adof_index_for_type( moris::uint aDofType )
        {
            MORIS_ERROR( false, "Solver_Interface::get_adof_index_for_type, Only works with MSI and multigrid" );
            return 0;
        };

        //---------------------------------------------------------------------------------------------------------

        /**
         * @brief builds thr graph based on the precomputed sparsity pattern
         * This is for Petsc right now
         *
         * @param aMat
         * @param aUseSparsityPattern
         */
        void build_graph( moris::sol::Dist_Matrix* aMat, bool aUseSparsityPattern = false );

        //---------------------------------------------------------------------------------------------------------

        void fill_matrix_and_RHS(
                moris::sol::Dist_Matrix* aMat,
                moris::sol::Dist_Vector* aVectorRHS );

        //---------------------------------------------------------------------------------------------------------

        void fill_matrix_and_RHS(
                moris::sol::Dist_Matrix* aMat,
                moris::sol::Dist_Vector* aVectorRHS,
                moris::sol::Dist_Vector* aFullSolutionVector );

        //---------------------------------------------------------------------------------------------------------

        void assemble_jacobian(
                moris::sol::Dist_Matrix* aMat,
                const fem::Time_Continuity_Flag = fem::Time_Continuity_Flag::DEFAULT );

        //---------------------------------------------------------------------------------------------------------

        void assemble_RHS(
                moris::sol::Dist_Vector* aVectorRHS,
                const fem::Time_Continuity_Flag = fem::Time_Continuity_Flag::DEFAULT );

        //------------------------------------------------------------------------------

        void assemble_staggered_RHS_contribution( moris::sol::Dist_Vector* aVectorRHS );

        //------------------------------------------------------------------------------

        // FIXME first draft. change this after diagonal jac gets moved into solver
        void assemble_additional_DqDs_RHS_contribution( moris::sol::Dist_Vector* aVectorRHS );

        //---------------------------------------------------------------------------------------------------------

        void get_adof_ids_based_on_criteria(
                Vector< moris::Matrix< IdMat > >& aCriteriaIds,
                const moris::real                 aThreshold );

        //---------------------------------------------------------------------------------------------------------

        virtual void
        calculate_criteria(
                const moris::uint& aMySetInd,
                const moris::uint& aMyElementInd )
        {
            MORIS_ERROR( false, "Solver_Interface::calculate_criteria(), not implemented for base class" );
        };

        //---------------------------------------------------------------------------------------------------------

        virtual const Vector< moris::Matrix< DDRMat > >&
        get_criteria( const moris::uint& aMySetInd )
        {
            MORIS_ERROR( false, "Solver_Interface::get_criteria(), not implemented for base class" );
            return *( new Vector< moris::Matrix< DDRMat > >() );
        };

        //---------------------------------------------------------------------------------------------------------

        virtual void
        set_requested_IQI_names( const Vector< std::string >& aIQINames )
        {
            MORIS_ERROR( false, "Solver_Interface::set_requested_IQI_type(), not implemented for base class" );
        };

        //------------------------------------------------------------------------------

        /**
         * @brief Set the solver warehouse object
         *
         * @param aSolverWarehouse
         */

        virtual void
        set_solver_warehouse( const std::shared_ptr< sol::SOL_Warehouse >& aSolverWarehouse );

        //------------------------------------------------------------------------------

        /**
         * @brief virtual method to be overloaded by MSI child class
         *
         */
        virtual void compute_sparsity_pattern()
        {
            MORIS_ERROR( false, "Solver_Interface::compute_sparsity_pattern(), not implemented for base class" );
        };

        /**
         * Updates the underlying problem being solved.
         */
        virtual void update_problem() {
        };
    };
}    // namespace moris

#endif /* SRC_DISTLINALG_CL_SOLVER_INPUT_HPP_ */
