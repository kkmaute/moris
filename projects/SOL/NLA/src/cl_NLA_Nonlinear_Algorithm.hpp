/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Algorithm.hpp
 *
 */

#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_ALGORITHM_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_ALGORITHM_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;
        class Dist_Vector;
    }    // namespace sol
    class Solver_Interface;
    namespace tsa
    {
        class Time_Solver_Algorithm;
    }
    namespace dla
    {
        class Linear_Solver_Algorithm;
        class Linear_Solver;
    }    // namespace dla

    namespace NLA
    {
        class Nonlinear_Solver;
        class Nonlinear_Problem;
        class Nonlinear_Algorithm
        {
          private:

          protected:
            //! Pointer to my nonlinear solver manager
            Nonlinear_Solver* mMyNonLinSolverManager = nullptr;

            //! Pointer to the linear solver manager
            dla::Linear_Solver* mLinSolverManager = nullptr;

            //! Pointer to the linear solver manager for adjoint solve. Can be the same than mLinSolverManager
            dla::Linear_Solver* mLinSolverManagerForAdjoint = nullptr;

            //! pointer to the nonlinear problem
            Nonlinear_Problem* mNonlinearProblem = nullptr;

            //! Parameterlist for this nonlinear solver
            moris::Parameter_List mParameterListNonlinearSolver;

            bool mLinSolverOwned = false;

            friend class Convergence;

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Member function which keeps track of used time for a particular purpose.
             */
            moris::real calculate_time_needed( const clock_t aTime );

          public:

            //--------------------------------------------------------------------------------------------------

            Nonlinear_Algorithm( const Parameter_List& aParameterlist = prm::create_nonlinear_algorithm_parameter_list() )
                    : mParameterListNonlinearSolver( aParameterlist ){};

            //--------------------------------------------------------------------------------------------------

            virtual ~Nonlinear_Algorithm(){};

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set the linear solver
             *
             * @param[in] aLinSolverManager Linear solver manager
             */
            void set_linear_solver( dla::Linear_Solver* aLinSolver );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set the linear solver for adjoint solve
             *
             * @param[in] aLinSolverManager Linear solver manager
             */
            void set_linear_solver_for_adjoint_solve( dla::Linear_Solver* aLinSolver );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Call to solve the nonlinear system
             *
             * @param[in] aNonlinearProblem Nonlinear problem
             */
            virtual void solver_nonlinear_system( Nonlinear_Problem* aNonlinearProblem ) = 0;

            //--------------------------------------------------------------------------------------------------

            virtual void get_full_solution( moris::Matrix< DDRMat >& LHSValues ) = 0;

            //--------------------------------------------------------------------------------------------------

            virtual void extract_my_values(
                    const moris::uint&                      aNumIndices,
                    const moris::Matrix< DDSMat >&          aGlobalBlockRows,
                    const moris::uint&                      aBlockRowOffsets,
                    Vector< moris::Matrix< DDRMat > >& LHSValues ) = 0;

            //--------------------------------------------------------------------------------------------------

            void set_nonlinear_solver_manager( Nonlinear_Solver* aNonlinSolverManager );

            //--------------------------------------------------------------------------------------------------

            Nonlinear_Solver* get_my_nonlin_solver();

            virtual void
            set_my_time_solver_algorithm( std::shared_ptr< tsa::Time_Solver_Algorithm > aMyTimeSolverAlgorithm )
            {
                MORIS_ASSERT( false, "set_my_time_solver_algorithm(): function not implemented" );
            }

            virtual void
            initialize_variables( Nonlinear_Problem* aNonlinearProblem )
            {
                MORIS_ASSERT( false, "initialize_variables(): this function is to be used for the arc length algorithm" );
            }
        };
    }    // namespace NLA
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_ALGORITHM_HPP_ */
