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
#include "moris_typedefs.hpp"    // CON/src
#include "cl_Vector.hpp"
#include <memory>
#include <utility>
#include "fn_PRM_SOL_Parameters.hpp"

namespace moris::dla
{
    class Linear_Solver_Algorithm;
    class Linear_Problem;
    class Preconditioner_Trilinos;
    class Linear_Solver
    {
      private:
        //! Linear solver list
        Vector< std::shared_ptr< Linear_Solver_Algorithm > > mLinearSolverList;

        moris::uint mCallCounter = 0;

        std::string mLhsOutputFileName;

      protected:
        moris::Parameter_List mParameterListLinearSolver;

      public:
        Linear_Solver( const moris::Parameter_List& aParameterList = prm::create_linear_solver_parameter_list() );

        //--------------------------------------------------------------------------------------------------

        virtual ~Linear_Solver(){};

        //--------------------------------------------------------------------------------------------------

        /**
         * @brief Set linear solver. Uses push back to add the given linear solver to the list.
         *
         * @param[in] aNonLinSolver Pointer to nonlinear solver.
         */
        void set_linear_algorithm( const std::shared_ptr< Linear_Solver_Algorithm >& aLinSolverAlgorithm );

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
            mLhsOutputFileName = std::move( aLhsOutputFileName );
        }

        const std::string&
        get_LHS_output_filename()
        {
            return mLhsOutputFileName;
        }

        //--------------------------------------------------------------------------------------------------
    };
}    // namespace moris::dla

#endif /* MORIS_DISTLINALG_CL_DLA_Linear_Solver_MANAGER_HPP_ */
