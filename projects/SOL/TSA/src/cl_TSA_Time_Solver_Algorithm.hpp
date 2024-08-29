/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Algorithm.hpp
 *
 */

#pragma once

#include <iostream>
#include <utility>

// MORIS header files.
#include "cl_Vector.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;
        class Dist_Vector;
    }    // namespace sol
    class Solver_Interface;
    namespace NLA
    {
        class Nonlinear_Solver;
    }
    namespace sol
    {
        class SOL_Warehouse;
    }

    namespace tsa
    {
        class Time_Solver;
        class Time_Solver_Algorithm
        {
          private:

          protected:

            uint mTimeSteps; // Number of time steps
            real mTimeIncrements; // Increments in time

            //! Pointer to database
            sol::SOL_Warehouse* mSolverWarehouse = nullptr;

            //! Pointer to solver interface
            Solver_Interface* mSolverInterface = nullptr;

            Vector< Time_Solver* > mTimeSolverList;

            //!  Pointer to time solver
            Time_Solver* mMyTimeSolver;

            //!  Nonlinear Solver
            NLA::Nonlinear_Solver* mNonlinearSolver = nullptr;

            //!  Nonlinear Solver for adjoint solve. Can be the same than mNonlinearSolver
            NLA::Nonlinear_Solver* mNonlinearSolverForAdjoint = nullptr;

            sol::Dist_Map* mFullMap = nullptr;

            //! name for output file for solution vector if set by user
            std::string mOutputSolVecFileName = "";

            /**
             * @brief Member function which keeps track of used time for a particular purpose.
             */
            static real calculate_time_needed( clock_t aTime );

          public:

            /**
             * @brief Constructor using a given parameter list
             *
             * @param[in] aParameterList     User defined parameter list
             */
            explicit Time_Solver_Algorithm(
                    const Parameter_List& aParameterList = prm::create_time_solver_algorithm_parameter_list() );

            //-------------------------------------------------------------------------------
            /**
             * @brief Destructor         *
             */
            virtual ~Time_Solver_Algorithm();

            //-------------------------------------------------------------------------------

            void delete_pointers();

            //-------------------------------------------------------------------------------
            /**
             * @brief Solve call using a given solution vector
             *
             * @param[in] aFullVector     Solution Vector
             */
            virtual void solve( Vector< sol::Dist_Vector* >& aFullVector ){};

            //-------------------------------------------------------------------------------
            /**
             * @brief finalize call for algorithm
             *
             */
            void finalize();

            //-------------------------------------------------------------------------------
            /**
             * @brief set warehouse to algorithm
             *
             * @param[in] aFullVector     Solution Vector
             */
            void
            set_solver_warehouse( sol::SOL_Warehouse* aSolverWarehouse )
            {
                mSolverWarehouse = aSolverWarehouse;
            };

            //-------------------------------------------------------------------------------
            /**
             * @brief set nonlinear solver
             *
             * @param[in] aNonlinearSolver     Nonlinear solver
             */
            void
            set_nonlinear_solver( NLA::Nonlinear_Solver* aNonlinearSolver )
            {
                mNonlinearSolver = aNonlinearSolver;
            };

            //-------------------------------------------------------------------------------
            /**
             * @brief set nonlinear solver for adjoint solve
             *
             * @param[in] aNonlinearSolver     Nonlinear solver
             */
            void
            set_nonlinear_solver_for_adjoint_solve( NLA::Nonlinear_Solver* aNonlinearSolver )
            {
                mNonlinearSolverForAdjoint = aNonlinearSolver;
            };

            //-------------------------------------------------------------------------------
            /**
             * @brief set time solver pointer
             *
             * @param[in] aTimeSolver     Time solver
             */
            void
            set_time_solver( Time_Solver* aTimeSolver )
            {
                mMyTimeSolver = aTimeSolver;
            };

            //-------------------------------------------------------------------------------
            /**
             * @brief set file name for outputting solution vectors
             *
             * @param[in] aOutputSolVecFileName   Name for output files containing solution vector
             */
            void
            set_output_filename( std::string aOutputSolVecFileName )
            {
                mOutputSolVecFileName = std::move( aOutputSolVecFileName );
            };

            //-------------------------------------------------------------------------------
            // arc-length functions
            virtual void
            set_lambda_increment( moris::real aLambdaInc )
            {
                MORIS_ASSERT( false, "Time_Solver_Algorithm:set_lambda_increment(): arc-length uses the monolithic time solver" );
            };

            //-------------------------------------------------------------------------------

            virtual moris::real
            get_new_lambda()
            {
                MORIS_ASSERT( false, "Time_Solver_Algorithm:get_new_lambda(): arc-length uses the monolithic time solver" );
                return 0;
            };
        };
    }    // namespace tsa
}    // namespace moris
