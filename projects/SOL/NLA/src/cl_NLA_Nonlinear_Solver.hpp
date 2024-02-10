/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Solver.hpp
 *
 */
#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_HPP_

// MORIS header files.
#include "moris_typedefs.hpp"    // CON/src
#include "cl_Vector.hpp"
#include <memory>
#include "cl_Param_List.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"

namespace moris
{
    class Solver_Interface;
    namespace sol
    {
        class SOL_Warehouse;
    }
    namespace tsa
    {
        class Time_Solver_Algorithm;
    }
    namespace NLA
    {
        class Nonlinear_Problem;
        class Nonlinear_Algorithm;

        class Nonlinear_Solver
        {
          private:
            //! List of list of dof types
            Vector< Vector< enum MSI::Dof_Type > > mStaggeredDofTypeList;

            //! List of secondary dependencies
            Vector< Vector< enum MSI::Dof_Type > > mSecondaryDofTypeList;

            //! List with nonlinear solvers
            Vector< std::shared_ptr< Nonlinear_Algorithm > > mNonlinearSolverAlgorithmList;

            //! List with nonlinear solvers
            Vector< Nonlinear_Solver* > mNonLinearSubSolverList;

            //! Pointer to solver database
            sol::SOL_Warehouse* mSolverWarehouse = nullptr;

            //! Pointer to nonlinear problem
            Nonlinear_Problem* mNonlinearProblem = nullptr;

            //! Pointer to solver interface
            Solver_Interface* mSolverInput = nullptr;

            // Nonlinear solver manager index
            sint mNonlinearSolverManagerIndex = -1;

            //! Reference norm
            real mRefNorm = -1.0;

            //! Current residual norm
            real mResidualNorm = -1.0;

            //! Static reference norm
            real mStaticRefNorm = -1.0;

            //! Current residual norm
            real mStaticResidualNorm = -1.0;

            //! Flag for computing static residual
            bool mComputeStaticResidual = false;

            //! Number of iterations of algorithm needed to converge / stop relative to maximum
            real mRelNumIterations = 0;

            ParameterList mParameterListNonLinearSolver;

            enum NonlinearSolverType mNonLinSolverType = NonlinearSolverType::END_ENUM;

            uint mCallCounter                = 0;
            uint mCallCounterNonlinearSolver = 0;

            sint mLevel = 0;

            sint mTimeIter = 0;

            tsa::Time_Solver_Algorithm* mMyTimeSolverAlgorithm = nullptr;
            //        enum NonlinearSolverType mMyNonLinSolverType = NonlinearSolverType::UNDEFINED;

          protected:

          public:
            /**
             * @brief Constructor. Creates a default nonlinear solver
             *
             * @param[in] aNonLinSolverType Nonlinear solver type. Default is Newton
             */
            Nonlinear_Solver( const enum NonlinearSolverType aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

            Nonlinear_Solver( const enum NonlinearSolverType aNonLinSolverType, const ParameterList aParameterlist );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Constructor. Sets given List of nonlinear solvers to this nonlinear solver manager
             *
             * @param[in] aNonlinearSolverList List of nonlinear solvers.
             * @param[in] aNonLinSolverType Nonlinear solver type. Default is Newton
             */
            Nonlinear_Solver(
                    Vector< std::shared_ptr< Nonlinear_Algorithm > >& aNonlinearSolverList,
                    const enum NonlinearSolverType                  aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

            //--------------------------------------------------------------------------------------------------

            ~Nonlinear_Solver();

            //--------------------------------------------------------------------------------------------------

            void free_memory();

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Sets one of the lists this nonlinear solver manager is operating on. Should be called multiple times for black solvers
             *
             * @param[in] aStaggeredDofTypeList List of dof types.
             * @param[in] aLevel                Solver level in the block structure. Default is 0
             */
            void set_dof_type_list(
                    const Vector< enum MSI::Dof_Type > aStaggeredDofTypeList,
                    const sint                       aLevel = 0 );

            //--------------------------------------------------------------------------------------------------

            void set_secondary_dof_type_list( const Vector< enum MSI::Dof_Type > aStaggeredDofTypeList );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Sets sub-non-linear solver this nonlinear solver is operating on
             *
             * @param[in] aNonlinearSolver Pointer to nonlinear solver
             */
            void set_sub_nonlinear_solver( NLA::Nonlinear_Solver* aNonlinearSolver );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Sets sub-non-linear solver this nonlinear solver is operating on
             *
             * @param[in] aNonlinearSolver Pointer to nonlinear solver
             * @param[in] aListEntry       List entry
             */
            void set_sub_nonlinear_solver(
                    NLA::Nonlinear_Solver* aNonlinearSolver,
                    const uint             aListEntry );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Sets solver interface
             *
             * @param[in] aSolverInterface Pointer to solver interface
             */
            void
            set_solver_interface( Solver_Interface* aSolverInterface )
            {
                mSolverInput = aSolverInterface;
            };

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Gets solver interface
             *
             * @param[out] Pointer to solver interface
             */
            Solver_Interface*
            get_solver_interface()
            {
                return mSolverInput;
            };

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Sets sub-non-linear solver this nonlinear solver is operating on
             *
             * @param[in] aNonlinearSolver Pointer to nonlinear solver
             * @param[in] aListEntry       List entry
             */
            Nonlinear_Solver*
            get_sub_nonlinear_solver( const uint aListEntry )
            {
                return mNonLinearSubSolverList( aListEntry );
            };

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set nonlinear solver. Uses push back to add the given nonlinear solver to the list.
             *
             * @param[in] aNonLinSolver Pointer to nonlinear solver.
             */
            void set_nonlinear_algorithm( std::shared_ptr< Nonlinear_Algorithm > aNonLinSolver );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set nonlinear solver on position in list
             *
             * @param[in] aNonLinSolver Pointer to nonlinear solver.
             * @param[in] aListEntry Pointer to nonlinear solver.
             */
            void set_nonlinear_algorithm(
                    std::shared_ptr< Nonlinear_Algorithm > aLinSolver,
                    const uint                             aListEntry );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Get function for list of list of dof types
             *
             * @param[out] rListOfListsOfDofTypes Returns the nonlinear solver managers list of list of dof types
             */
            Vector< Vector< enum MSI::Dof_Type > >
            get_dof_type_list()
            {
                return mStaggeredDofTypeList;
            };

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Get the union of this nonlinear solver managers dof types
             *
             * @param[out] rUnionListOfDofTypes Returns the union list of this nonlinear solver managers dof types
             */
            Vector< enum MSI::Dof_Type > get_dof_type_union();

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Get secondary dof type list
             *
             * @param[out] rUnionListOfDofTypes Returns the union list of this nonlinear solver managers dof types
             */
            Vector< Vector< enum MSI::Dof_Type > >
            get_secondary_dof_type_list()
            {
                return mSecondaryDofTypeList;
            };

            //--------------------------------------------------------------------------------------------------

            Vector< enum MSI::Dof_Type > get_sec_dof_type_union();

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Sets the index of this nonlinear solver manager
             *
             * @param[in] aNonlinearSolverManagerIndex Nonlinear solver manager index
             */
            void set_nonlinear_solver_manager_index( const sint aNonlinearSolverManagerIndex );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Get the nonlinear solver manager index
             *
             * @param[out] rNonlinearSolverManagerIndex Returns the nonlinear solver manager index
             */
            sint get_nonlinear_solver_manager_index();

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Returns pointer to the solver database
             *
             * @param[out] rSolverDatabase Returns the pointer to the solver database
             */
            sol::SOL_Warehouse*
            get_solver_warehouse()
            {
                return mSolverWarehouse;
            };

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set pointer to the solver database
             *
             * @param[in] rSolverDatabase Pointer to the solver database
             */
            void set_solver_warehouse( sol::SOL_Warehouse* aSolverWarehouse );

            //--------------------------------------------------------------------------------------------------
            /**
             * @brief sets the time iteration
             *
             * @param[in] aTimeIter time iteration number
             */
            void set_time_step_iter( const sint aTimeIter );

            //--------------------------------------------------------------------------------------------------
            /**
             * @brief returns the time iteration value
             */
            sint get_time_step_iter();

            //--------------------------------------------------------------------------------------------------

            void solve( sol::Dist_Vector* aFullVector );

            void solve( Nonlinear_Problem* aNonlinearProblem );

            //--------------------------------------------------------------------------------------------------

            void get_full_solution( Matrix< DDRMat >& LHSValues );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Returns reference residual norm
             *
             * @return  reference residual norm
             */
            real
            get_ref_norm()
            {
                MORIS_ASSERT( mRefNorm >= 0.0,
                        "Nonlinear_Solver::get_ref_norm - invalid reference residual norm." );

                return mRefNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief set reference residual norm
             *
             * @param[in] aRefNorm  reference norm
             */
            void
            set_ref_norm( const real& aRefNorm )
            {
                mRefNorm = aRefNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Returns current residual norm
             *
             * @return residual norm
             */
            real
            get_residual_norm()
            {
                MORIS_ASSERT( mResidualNorm >= 0.0,
                        "Nonlinear_Solver::get_residual_norm - invalid residual norm." );

                return mResidualNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set current residual norm
             *
             * @param[in] aResidualNorm  current residual norm
             */
            void
            set_residual_norm( const real& aResidualNorm )
            {
                mResidualNorm = aResidualNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Returns static residual reference  norm
             *
             * @return  reference residual norm
             */
            real
            get_static_ref_norm()
            {
                MORIS_ASSERT( mStaticRefNorm >= 0.0,
                        "Nonlinear_Solver::get_static_ref_norm - invalid static residual reference norm." );

                return mStaticRefNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief set static residual reference norm
             *
             * @param[in] aRefNorm  reference norm
             */
            void
            set_static_ref_norm( const real& aStaticRefNorm )
            {
                mStaticRefNorm = aStaticRefNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Returns current residual norm
             *
             * @return static residual norm
             */
            real
            get_static_residual_norm()
            {
                MORIS_ASSERT( mStaticResidualNorm >= 0.0,
                        "Nonlinear_Solver::get_static_residual_norm - invalid static residual norm." );

                return mStaticResidualNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set current static residual norm
             *
             * @param[in] aResidualNorm  current static residual norm
             */
            void
            set_static_residual_norm( const real& aStaticResidualNorm )
            {
                mStaticResidualNorm = aStaticResidualNorm;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Returns number of iterations need to converge relative to maximum number of iterations
             *
             * @return relative number of iterations
             */
            real
            get_relative_number_iterations()
            {
                return mRelNumIterations;
            }

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief set number of iterations need to converge relative to maximum number of iterations
             *
             * @param[in] aRelNumIterations  relative number of iterations
             */
            void
            set_relative_number_iterations( const real& aRelNumIterations )
            {
                mRelNumIterations = aRelNumIterations;
            }

            //--------------------------------------------------------------------------------------------------

            void set_nonlinear_solver_manager_parameters();

            //--------------------------------------------------------------------------------------------------

            void set_compute_static_residual_flag( bool aFlag );

            //--------------------------------------------------------------------------------------------------

            bool get_compute_static_residual_flag();

            //--------------------------------------------------------------------------------------------------

            Nonlinear_Problem* get_my_nonlin_problem();

            void set_nonlin_solver_type( enum NonlinearSolverType aNonLinSolverType );

            enum NonlinearSolverType get_nonlin_solver_type();

            void set_time_solver_type( tsa::Time_Solver_Algorithm* aTimeSolverAlgorithm );

            //--------------------------------------------------------------------------------------------------

            ParameterListTypes&
            set_param( char const * aKey )
            {
                return mParameterListNonLinearSolver( aKey );
            }
        };
    }    // namespace NLA
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_HPP_ */
