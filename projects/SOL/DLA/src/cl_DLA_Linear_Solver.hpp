/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver.hpp
 *
 */

#ifndef MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_
#define MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_

// MORIS header files.
#include "typedefs.hpp"    // CON/src
#include "cl_Cell.hpp"
#include <memory>
#include "cl_Param_List.hpp"

namespace moris
{
    namespace dla
    {
        class Linear_Solver_Algorithm;
        class Linear_Problem;
        class Preconditioner_Trilinos;
        class Linear_Solver
        {
          private:
            //! Linear solver list
            moris::Cell< std::shared_ptr< Linear_Solver_Algorithm > > mLinearSolverList;

            moris::uint mCallCounter = 0;

            std::string mLhsOutputFileName;

          protected:
            moris::ParameterList mParameterListLinearSolver;

          public:
            //--------------------------------------------------------------------------------------------------
            /**
             * @brief Constructor. Creates a default linear solver.
             */
            Linear_Solver();

            Linear_Solver( const moris::ParameterList aParameterlist );

            //--------------------------------------------------------------------------------------------------

            virtual ~Linear_Solver(){};

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set linear solver. Uses push back to add the given linear solver to the list.
             *
             * @param[in] aNonLinSolver Pointer to nonlinear solver.
             */
            void set_linear_algorithm( std::shared_ptr< Linear_Solver_Algorithm > aLinSolverAlgorithm );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Set linear solver on position in list
             *
             * @param[in] aLinSolver Pointer to nonlinear solver.
             * @param[in] aListEntry Pointer to nonlinear solver.
             */
            void set_linear_algorithm( const moris::uint       aListEntry,
                    std::shared_ptr< Linear_Solver_Algorithm > aLinSolverAlgorithm );

            //--------------------------------------------------------------------------------------------------

            /**
             * @brief Solve linear system
             *
             * @param[in] aLinearProblem Pointer to linear problem.
             * @param[in] aLinearProblem Iteration number.
             */
            void solver_linear_system( dla::Linear_Problem* aLinearProblem,
                    const moris::sint                       aIter );

            //--------------------------------------------------------------------------------------------------

            void set_linear_solver_manager_parameters();

            //--------------------------------------------------------------------------------------------------

            void
            set_LHS_output_filename( std::string aLhsOutputFileName )
            {
                mLhsOutputFileName = aLhsOutputFileName;
            }

            const std::string&
            get_LHS_output_filename()
            {
                return mLhsOutputFileName;
            }

            //--------------------------------------------------------------------------------------------------

            ParameterListTypes&
            set_param( char const * aKey )
            {
                return mParameterListLinearSolver( aKey );
            }

            moris::ParameterList const &
            get_param_list()
            {
                return mParameterListLinearSolver;
            }

            //--------------------------------------------------------------------------------------------------
        };
    }    // namespace dla
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_ */
