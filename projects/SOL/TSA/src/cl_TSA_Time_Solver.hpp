/*
 * cl_TSA_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_HPP_

#include <iostream>

#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Database.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_DLA_Enums.hpp"

// MORIS header files.

namespace moris
{
class Map_Class;
class Dist_Vector;
class Solver_Interface;

namespace tsa
{
    class Time_Solver
    {
    private:

    protected:
        //! Pointer to database
        NLA::Nonlinear_Database * mDatabase = nullptr;

        //! Pointer to solver interface
        Solver_Interface * mSolverInterface = nullptr;

        moris::Cell< Time_Solver * > mTimeSolverList;

        //! Full Vector
        Dist_Vector * mFullVector = nullptr;

        //! Full Vector
        Dist_Vector * mPrevFullVector = nullptr;

        moris::Cell< moris::Cell< enum MSI::Dof_Type > > mMyDofTypeList;

        moris::uint mCallCounter = 0;

        NLA::Nonlinear_Solver * mNonlinearSolver = nullptr;

        //!  flag indicating if this is the master time solver
        bool mIsMasterTimeSolver = false;

        Map_Class       * mFullMap;

        /**
         * @brief Member function which keeps track of used time for a particular purpose.
         */
        moris::real calculate_time_needed( const clock_t aTime );

    public:
        //-------------------------------------------------------------------------------

        Time_Solver( const enum MapType aMapType = MapType::Epetra );

        //-------------------------------------------------------------------------------

        ~Time_Solver();

        //-------------------------------------------------------------------------------

        virtual void solve(){};

        //-------------------------------------------------------------------------------

        virtual void solve( Dist_Vector * aFullVector ){};

        //-------------------------------------------------------------------------------

        void finalize();

        //-------------------------------------------------------------------------------

        void set_database( NLA::Nonlinear_Database * aDatabase )
        {
            mDatabase = aDatabase;
        };

        void get_full_solution( moris::Matrix< DDRMat > & LHSValues );

        /**
         * @brief Sets one of the lists this nonlinear solver manager is operating on. Should be called multiple times for black solvers
         *
         * @param[in] aStaggeredDofTypeList List of dof types.
         * @param[in] aLevel                Solver level in the block structure. Default is 0
         */
        void set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
                                const moris::sint                       aLevel =  0);

        void set_nonlinear_solver( NLA::Nonlinear_Solver * aNonlinearSolver )
        {
            mNonlinearSolver = aNonlinearSolver;
        };

        void set_time_solver(       Time_Solver * aTimeSolver )
        {
            if( mCallCounter == 0 )
            {
                // removes all elements from the Cell and destroy them
                mTimeSolverList.clear();

                // Resize the Cell to size = 1
                mTimeSolverList.resize( 1, nullptr );

                // Set nonlinear solver on first entry
                mTimeSolverList( 0 ) =  aTimeSolver;
            }
            else
            {
                // set nonlinear solver on next entry
                mTimeSolverList.push_back( aTimeSolver );
            }

            mCallCounter = mCallCounter + 1;
        }
    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_TIME_SOLVER_HPP_ */
