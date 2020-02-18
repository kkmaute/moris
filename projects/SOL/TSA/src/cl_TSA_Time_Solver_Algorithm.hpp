/*
 * cl_TSA_Time_Solver_Algorithm.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_ALGORITHM_HPP_
#define MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_ALGORITHM_HPP_

#include <iostream>

#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_SOL_Warehouse.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_SOL_Enums.hpp"

// MORIS header files.

namespace moris
{
class Dist_Map;
class Dist_Vector;
class Solver_Interface;

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

        //moris::Cell< moris::Cell< enum MSI::Dof_Type > > mMyDofTypeList;

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

        Time_Solver_Algorithm( const enum sol::MapType aMapType = sol::MapType::Epetra );

        Time_Solver_Algorithm( const ParameterList aParameterlist,
                               const enum sol::MapType aMapType = sol::MapType::Epetra );

        //-------------------------------------------------------------------------------

        ~Time_Solver_Algorithm();

        //-------------------------------------------------------------------------------

        virtual void solve(){};

        //-------------------------------------------------------------------------------

        virtual void solve( Dist_Vector * aFullVector ){};

        //-------------------------------------------------------------------------------

        void finalize();

        //-------------------------------------------------------------------------------

        void set_solver_warehouse( sol::SOL_Warehouse * aSolverWarehouse )
        {
            mSolverWarehouse = aSolverWarehouse;
        };

        void get_full_solution( moris::Matrix< DDRMat > & LHSValues );

        /**
         * @brief Sets one of the lists this nonlinear solver manager is operating on. Should be called multiple times for black solvers
         *
         * @param[in] aStaggeredDofTypeList List of dof types.
         * @param[in] aLevel                Solver level in the block structure. Default is 0
         */
//        void set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
//                                const moris::sint                       aLevel =  0);

        void set_nonlinear_solver( NLA::Nonlinear_Solver * aNonlinearSolver )
        {
            mNonlinearSolver = aNonlinearSolver;
        };

        void set_time_solver( Time_Solver * aTimeSolver )
        {
            mMyTimeSolver = aTimeSolver;
        };

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

        virtual moris::real get_new_lambda()
        {
            MORIS_ASSERT(false, "Time_Solver_Algorithm:get_new_lambda(): arc-length uses the monolithic time solver");
            return 0;
        };
    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_ALGORITHM_HPP_ */
