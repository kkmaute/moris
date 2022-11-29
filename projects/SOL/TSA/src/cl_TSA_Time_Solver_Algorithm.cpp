/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Algorithm.cpp
 *
 */

#include <ctime>

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_TSA_Time_Solver_Algorithm.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
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

//-------------------------------------------------------------------------------

Time_Solver_Algorithm::Time_Solver_Algorithm(
        const ParameterList     aParameterlist,
        const enum sol::MapType aMapType )
        : mParameterListTimeSolver( aParameterlist )
{
}

//-------------------------------------------------------------------------------

Time_Solver_Algorithm::~Time_Solver_Algorithm()
{
    this->delete_pointers();
}

//-------------------------------------------------------------------------------

void
Time_Solver_Algorithm::delete_pointers()
{
    delete mFullMap;
    mFullMap = nullptr;
}

//-------------------------------------------------------------------------------

moris::real
Time_Solver_Algorithm::calculate_time_needed( const clock_t aTime )
{
    moris::real tDeltaTime = ( moris::real )( clock() - aTime ) / CLOCKS_PER_SEC;

    moris::real tDeltaTimeMax = max_all( tDeltaTime );

    return tDeltaTimeMax;
}

//-------------------------------------------------------------------------------

void
Time_Solver_Algorithm::finalize()
{
    this->delete_pointers();

    // create map object
    sol::Matrix_Vector_Factory tMatFactory( mMyTimeSolver->get_solver_warehouse()->get_tpl_type() );

    mSolverInterface = mMyTimeSolver->get_solver_interface();

    // build map for full vector
    mFullMap = tMatFactory.create_full_map(
            mSolverInterface->get_my_local_global_map(),
            mSolverInterface->get_my_local_global_overlapping_map() );
}

//-------------------------------------------------------------------------------

void
Time_Solver_Algorithm::set_time_solver_parameters()
{
    // Number of time steps
    mParameterListTimeSolver.insert( "TSA_Num_Time_Steps", 1 );

    // Time Frame
    mParameterListTimeSolver.insert( "TSA_Time_Frame", 1.0 );
}
