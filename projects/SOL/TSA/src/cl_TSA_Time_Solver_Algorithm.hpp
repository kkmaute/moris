/*
 * cl_TSA_Time_Solver_Algorithm.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_ALGORITHM_HPP_
#define MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_ALGORITHM_HPP_

#include <iostream>

// MORIS header files.
#include "cl_Cell.hpp"
#include "cl_Param_List.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
class Dist_Map;
class Dist_Vector;
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
        //! Pointer to database
        sol::SOL_Warehouse * mSolverWarehouse = nullptr;

        //! Pointer to solver interface
        Solver_Interface * mSolverInterface = nullptr;

        moris::Cell< Time_Solver * > mTimeSolverList;

        Time_Solver * mMyTimeSolver;

        //! Full Vector
        Dist_Vector * mFullVector = nullptr;

        //! Full Vector
        Dist_Vector * mPrevFullVector = nullptr;

        moris::uint mCallCounter = 0;

        NLA::Nonlinear_Solver * mNonlinearSolver = nullptr;

        //!  flag indicating if this is the master time solver
        bool mIsMasterTimeSolver = false;

        Dist_Map * mFullMap = nullptr;

        //! Parameterlist for this nonlinear solver
        moris::ParameterList mParameterListTimeSolver;

        /**
         * @brief Member function which keeps track of used time for a particular purpose.
         */
        moris::real calculate_time_needed( const clock_t aTime );

    public:
        //-------------------------------------------------------------------------------
        /**
         * @brief default constructor
         *
         * @param[in] rSolverDatabase Poiner to the solver database
         */
        Time_Solver_Algorithm( const enum sol::MapType aMapType = sol::MapType::Epetra );

        //-------------------------------------------------------------------------------
        /**
         * @brief Constructor using a given parameter list
         *
         * @param[in] aParameterlist     User defined parameter list
         */
        Time_Solver_Algorithm( const ParameterList aParameterlist,
                               const enum sol::MapType aMapType = sol::MapType::Epetra );

        //-------------------------------------------------------------------------------
        /**
         * @brief Destructor         *
         */
        virtual ~Time_Solver_Algorithm();

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call
         *
         * @param[in] aFullVector
         */
        virtual void solve(){};

        //-------------------------------------------------------------------------------
        /**
         * @brief Solve call using a given soltion vector
         *
         * @param[in] aFullVector     Solution Vector
         */
        virtual void solve( Dist_Vector * aFullVector ){};

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
        void set_solver_warehouse( sol::SOL_Warehouse * aSolverWarehouse )
        {
            mSolverWarehouse = aSolverWarehouse;
        };

        //-------------------------------------------------------------------------------
        /**
         * @brief Get solution vector from algorithm
         *
         * @param[in] aLHSValues     Solution Vector
         */
        void get_full_solution( moris::Matrix< DDRMat > & aLHSValues );

        //-------------------------------------------------------------------------------
        /**
          * @brief set nonlinear solver
          *
          * @param[in] aNonlinearSolver     Nonlinear solver
          */
        void set_nonlinear_solver( NLA::Nonlinear_Solver * aNonlinearSolver )
        {
            mNonlinearSolver = aNonlinearSolver;
        };

        //-------------------------------------------------------------------------------
        /**
          * @brief set time solver pointer
          *
          * @param[in] aTimeSolver     Time solver
          */
        void set_time_solver( Time_Solver * aTimeSolver )
        {
            mMyTimeSolver = aTimeSolver;
        };

        //-------------------------------------------------------------------------------

        void set_time_solver_parameters();

        ParameterListTypes&  set_param( char const* aKey )
        {
            return mParameterListTimeSolver( aKey );
        }

        //-------------------------------------------------------------------------------
        // arc-length functions
        virtual void set_lambda_increment( moris::real aLambdaInc )
        {
            MORIS_ASSERT(false, "Time_Solver_Algorithm:set_lambda_increment(): arc-length uses the monolithic time solver");
        };

        //-------------------------------------------------------------------------------

        virtual moris::real get_new_lambda()
        {
            MORIS_ASSERT(false, "Time_Solver_Algorithm:get_new_lambda(): arc-length uses the monolithic time solver");
            return 0;
        };
    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_ALGORITHM_HPP_ */
