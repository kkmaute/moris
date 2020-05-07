/*
 * cl_TSA_Time_Solver_Algorithm.cpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#include <ctime>

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver.hpp"

#include "cl_SOL_Dist_Map.hpp"

extern moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace tsa;

//-------------------------------------------------------------------------------

Time_Solver_Algorithm::Time_Solver_Algorithm( const enum sol::MapType aMapType )
{
    this->set_time_solver_parameters();
}

Time_Solver_Algorithm::Time_Solver_Algorithm( const ParameterList aParameterlist,
                                              const enum sol::MapType aMapType ) : mParameterListTimeSolver( aParameterlist )
{}

//-------------------------------------------------------------------------------

Time_Solver_Algorithm::~Time_Solver_Algorithm()
{
    if ( mIsMasterTimeSolver )
    {
        delete( mFullVector );
//        delete( mPrevFullVector );
    }
    delete( mPrevFullVector );                 // FIXME There's a delete somewhere in HMR which need this memory leak. has to be fixed
    delete( mFullMap );
}
//-------------------------------------------------------------------------------

void Time_Solver_Algorithm::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
{
    mFullVector->extract_copy( LHSValues );
}

//-------------------------------------------------------------------------------

moris::real Time_Solver_Algorithm::calculate_time_needed( const clock_t aTime )
{
    moris::real tDeltaTime = (moris::real) ( clock() - aTime ) / CLOCKS_PER_SEC;

    moris::real tDeltaTimeMax   = tDeltaTime;

    max_all( tDeltaTime, tDeltaTimeMax );

    return tDeltaTimeMax;
}

//-------------------------------------------------------------------------------

//void Time_Solver_Algorithm::set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
//                                               const moris::sint                       aLevel )
//{
//    mMyDofTypeList.push_back( aStaggeredDofTypeList );
//}

//-------------------------------------------------------------------------------

void Time_Solver_Algorithm::finalize()
{
    if ( mIsMasterTimeSolver )
    {
        // create map object
        Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

        mSolverInterface = mSolverWarehouse->get_solver_interface();

        uint tNumRHMS = mSolverInterface->get_num_rhs();

        MORIS_LOG_INFO( "Creating main time solver system with %-5i dofs.", mSolverInterface->get_my_local_global_overlapping_map().numel() );

        mFullMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map() );

        // full vector and prev full vector
        mFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );

        mSolverInterface->set_solution_vector( mFullVector );

        mFullVector->vec_put_scalar( 0.0 );

        mPrevFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );
    }
    else
    {
        // create map object
        Matrix_Vector_Factory tMatFactory( mMyTimeSolver->get_solver_warehouse()->get_tpl_type() );

        mSolverInterface = mMyTimeSolver->get_solver_interface();
//        mSolverInterface = mMyTimeSolver->get_solver_warehouse()->get_solver_interface();

        uint tNumRHMS = mSolverInterface->get_num_rhs();

        mFullMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map() );

        mPrevFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, tNumRHMS );
    }

    mPrevFullVector->vec_put_scalar( 0.0 );
}

//-------------------------------------------------------------------------------

void Time_Solver_Algorithm::set_time_solver_parameters()
{
    // Number of time steps
    mParameterListTimeSolver.insert( "TSA_Num_Time_Steps", 1 );

    // Time Frame
    mParameterListTimeSolver.insert( "TSA_Time_Frame", 1.0 );
}
