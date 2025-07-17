/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Solver.cpp
 *
 */
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver.hpp"

#include <utility>
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"

#include "cl_SOL_Warehouse.hpp"

#include "cl_SOL_Dist_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

// Detailed Logging package
#include "cl_Tracer.hpp"

using namespace moris;
using namespace NLA;

//--------------------------------------------------------------------------------------------------------------------------

Nonlinear_Solver::Nonlinear_Solver(
        const Parameter_List& aSolverParameterList,
        const Parameter_List& aAlgorithmParameterList )
        : mSecondaryDofTypeList( Vector< Vector< enum MSI::Dof_Type > >( 0 ) )
        , mParameterListNonLinearSolver( aSolverParameterList )
{
    // create solver object
    Nonlinear_Solver_Factory               tSolFactory;
    std::shared_ptr< Nonlinear_Algorithm > tNonLinSolver = tSolFactory.create_nonlinear_solver( aAlgorithmParameterList );
    mNonlinearSolverAlgorithmList.resize( 1, nullptr );
    mNonlinearSolverAlgorithmList( 0 ) = tNonLinSolver;

    mStaggeredDofTypeList.resize( 0 );
}

//--------------------------------------------------------------------------------------------------

Nonlinear_Solver::Nonlinear_Solver(
        Vector< std::shared_ptr< Nonlinear_Algorithm > >& aNonlinearSolverList,
        const Parameter_List&                             aSolverParameterList )
        : mSecondaryDofTypeList( Vector< Vector< enum MSI::Dof_Type > >( 0 ) )
        , mParameterListNonLinearSolver( aSolverParameterList )
{
    mNonlinearSolverAlgorithmList = aNonlinearSolverList;
}

//--------------------------------------------------------------------------------------------------

Nonlinear_Solver::~Nonlinear_Solver()
{
    this->free_memory();
}

//--------------------------------------------------------------------------------------------------

void Nonlinear_Solver::free_memory()
{
    if ( mNonlinearProblem != nullptr )
    {
        delete ( mNonlinearProblem );

        mNonlinearProblem = nullptr;
    }
}

//--------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_dof_type_list(
        const Vector< enum MSI::Dof_Type >& aStaggeredDofTypeList,
        const moris::sint                   aLevel )
{
    mStaggeredDofTypeList.push_back( aStaggeredDofTypeList );
}

//--------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_secondary_dof_type_list( const Vector< enum MSI::Dof_Type >& aStaggeredDofTypeList )
{
    if ( mSecondaryDofTypeList.size() == 0 )
    {
        mSecondaryDofTypeList.clear();
    }

    mSecondaryDofTypeList.push_back( aStaggeredDofTypeList );
}

//--------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_nonlinear_algorithm( const std::shared_ptr< Nonlinear_Algorithm >& aNonLinSolver )
{
    if ( mCallCounter == 0 )
    {
        // removes all elements from the Cell and destroy them
        mNonlinearSolverAlgorithmList.clear();

        // Resize the Cell to size = 1
        mNonlinearSolverAlgorithmList.resize( 1, nullptr );

        // Set nonlinear solver on first entry
        mNonlinearSolverAlgorithmList( 0 ) = aNonLinSolver;
    }
    else
    {
        // set nonlinear solver on next entry
        mNonlinearSolverAlgorithmList.push_back( aNonLinSolver );
    }

    mCallCounter = mCallCounter + 1;
}

//-------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_nonlinear_algorithm(
        std::shared_ptr< Nonlinear_Algorithm > aNonLinSolver,
        const moris::uint                      aListEntry )
{
    // Check if list is smaller than given entry
    if ( aListEntry >= mNonlinearSolverAlgorithmList.size() )
    {
        // Resize to new entry value and set nullptr on new entries
        mNonlinearSolverAlgorithmList.resize( aListEntry + 1, nullptr );
    }
    // Set nonlinear solver on entry
    mNonlinearSolverAlgorithmList( aListEntry ) = std::move( aNonLinSolver );
}

//-------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_sub_nonlinear_solver( Nonlinear_Solver* aNonLinSolver )
{
    if ( mCallCounterNonlinearSolver == 0 )
    {
        // removes all elements from the Cell and destroy them
        mNonLinearSubSolverList.clear();

        // Resize the Cell to size = 1
        mNonLinearSubSolverList.resize( 1, nullptr );

        // Set nonlinear solver on first entry
        mNonLinearSubSolverList( 0 ) = aNonLinSolver;
    }
    else
    {
        // set nonlinear solver on next entry
        mNonLinearSubSolverList.push_back( aNonLinSolver );
    }

    mCallCounterNonlinearSolver = mCallCounterNonlinearSolver + 1;
}

//-------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_sub_nonlinear_solver(
        Nonlinear_Solver* aNonLinSolver,
        const moris::uint aListEntry )
{
    // Check if list is smaller than given entry
    if ( aListEntry >= mNonLinearSubSolverList.size() )
    {
        // Resize to new entry value and set nullptr on new entries
        mNonLinearSubSolverList.resize( aListEntry + 1, nullptr );
    }
    // Set nonlinear solver on entry
    mNonLinearSubSolverList( aListEntry ) = aNonLinSolver;
}

//-------------------------------------------------------------------------------------------------------

Vector< enum MSI::Dof_Type >
Nonlinear_Solver::get_dof_type_union()
{
    moris::sint tCounter = 0;

    // Loop over all dof type lists to determine the total number of dof types
    for ( moris::uint Ik = 0; Ik < mStaggeredDofTypeList.size(); ++Ik )
    {
        tCounter = tCounter + mStaggeredDofTypeList( Ik ).size();
    }

    // Create list of dof types with earlier determines size
    Vector< enum MSI::Dof_Type > tUnionEnumList( tCounter );
    tCounter = 0;

    // Loop over all dof types. Add them to union list
    for ( moris::uint Ik = 0; Ik < mStaggeredDofTypeList.size(); ++Ik )
    {
        for ( moris::uint Ii = 0; Ii < mStaggeredDofTypeList( Ik ).size(); ++Ii )
        {
            tUnionEnumList( tCounter++ ) = mStaggeredDofTypeList( Ik )( Ii );
        }
    }

    return tUnionEnumList;
}

//-------------------------------------------------------------------------------------------------------

Vector< enum MSI::Dof_Type >
Nonlinear_Solver::get_sec_dof_type_union()
{
    moris::sint tCounter = 0;

    // Loop over all dof type lists to determine the total number of dof types
    for ( moris::uint Ik = 0; Ik < mSecondaryDofTypeList.size(); ++Ik )
    {
        tCounter = tCounter + mSecondaryDofTypeList( Ik ).size();
    }

    // Create list of dof types with earlier determines size
    Vector< enum MSI::Dof_Type > tUnionEnumList( tCounter );
    tCounter = 0;

    // Loop over all dof types. Add them to union list
    for ( moris::uint Ik = 0; Ik < mSecondaryDofTypeList.size(); ++Ik )
    {
        for ( moris::uint Ii = 0; Ii < mSecondaryDofTypeList( Ik ).size(); ++Ii )
        {
            tUnionEnumList( tCounter++ ) = mSecondaryDofTypeList( Ik )( Ii );
        }
    }

    return tUnionEnumList;
}

//--------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_solver_warehouse( sol::SOL_Warehouse* aSolverWarehouse )
{
    mSolverWarehouse = aSolverWarehouse;

    mSolverInput = mSolverWarehouse->get_solver_interface();
}

//-------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::solve( sol::Dist_Vector* aFullVector )
{
    Tracer tTracer( "NonLinearSolver", "NonlinearSolver", "Solve" );

    mSolverInput = mSolverWarehouse->get_solver_interface();

    Vector< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

    mSolverInput->set_requested_dof_types( tDofTypeUnion );
    mSolverInput->set_secondary_dof_types( tDofTypeUnion );

    if ( mNonLinSolverType == NonlinearSolverType::NLBGS_SOLVER )
    {
        this->free_memory();

        mNonlinearProblem = new Nonlinear_Problem(
                mSolverWarehouse,
                mSolverInput,
                aFullVector,
                mNonlinearSolverManagerIndex,
                false,
                mSolverWarehouse->get_tpl_type() );
    }
    else
    {
        this->free_memory();

        mNonlinearProblem = new Nonlinear_Problem(
                mSolverWarehouse,
                mSolverInput,
                aFullVector,
                mNonlinearSolverManagerIndex,
                true,
                mSolverWarehouse->get_tpl_type() );
    }

    mNonlinearProblem->set_nonlinear_solver( this );

    map< enum MSI::Dof_Type, std::string > tDofTypeToNameMap = MSI::get_dof_type_name_map();

    std::string tDofTypeNames = "";
    for ( uint Ik = 0; Ik < tDofTypeUnion.size(); Ik++ )
    {
        tDofTypeNames = tDofTypeNames + tDofTypeToNameMap.find( tDofTypeUnion( Ik ) ) + ", ";
    }

    tDofTypeNames.pop_back();
    tDofTypeNames.pop_back();

    MORIS_LOG_SPEC( "Nonlinear solver operates on DOF types: ", tDofTypeNames );

    for ( const auto& iNonLinearAlgorithm : mNonlinearSolverAlgorithmList )
    {
        iNonLinearAlgorithm->set_nonlinear_solver_manager( this );
        iNonLinearAlgorithm->solver_nonlinear_system( mNonlinearProblem );
    }

    this->free_memory();
}

//-------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::solve( Nonlinear_Problem* aNonlinearProblem )
{
    moris::sint tErrorStatus        = 0;
    moris::sint tMaxNumLinRestarts  = mParameterListNonLinearSolver.get< moris::sint >( "NLA_max_non_lin_solver_restarts" );
    moris::sint tTryRestartOnFailIt = 1;

    MORIS_ERROR( mNonLinSolverType != NonlinearSolverType::NLBGS_SOLVER, "Nonlinear_Solver::solve(); Nonlinear Solver is NLBGS" );

    mNonlinearSolverAlgorithmList( 0 )->set_nonlinear_solver_manager( this );

    mNonlinearSolverAlgorithmList( 0 )->solver_nonlinear_system( aNonlinearProblem );

    // Restart the nonlinear solver
    while ( tErrorStatus != 0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && (moris::sint)mNonlinearSolverAlgorithmList.size() <= tMaxNumLinRestarts )
    {
        if ( par_rank() == 0 )
        {
            MORIS_LOG( " ... Previous nonlinear solve failed. Trying restart %i of %i", tTryRestartOnFailIt, tMaxNumLinRestarts );
        }

        // Re-solve scaled linear system with current solution as an initial guess
        // tErrorStatus = mNonlinearSolverAlgorithmList( tTryRestartOnFailIt )->solver_nonlinear_system( aNonlinearProblem );
        mNonlinearSolverAlgorithmList( tTryRestartOnFailIt )->solver_nonlinear_system( aNonlinearProblem );

        // Iterate TryRestartOnFailIt counter
        tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
    }
}

//--------------------------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::get_full_solution( moris::Matrix< DDRMat >& LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_time_step_iter( const sint aTimeIter )
{
    mTimeIter = aTimeIter;
}

//--------------------------------------------------------------------------------------------------------------------------

moris::sint
Nonlinear_Solver::get_time_step_iter()
{
    return mTimeIter;
}
//--------------------------------------------------------------------------------------------------------------------------
Nonlinear_Problem*
Nonlinear_Solver::get_my_nonlin_problem()
{
    return mNonlinearProblem;
}

//--------------------------------------------------------------------------------------------------------------------------

void Nonlinear_Solver::set_compute_static_residual_flag( bool aComputeFlag, bool mUseStaticResidual )
{
    mComputeStaticResidual = aComputeFlag;
    mUseStaticResidual     = mUseStaticResidual;

    // set flag for all sub-solvers
    for ( uint Ik = 0; Ik < mNonLinearSubSolverList.size(); ++Ik )
    {
        // get sub-solver
        Nonlinear_Solver* tSubSolver = mNonLinearSubSolverList( Ik );

        // check that sub-solver exists
        MORIS_ERROR( tSubSolver != nullptr,
                "Nonlinear_Solver::set_compute_static_residual_flag - sub-solver %d not defined.",
                Ik );

        // set compute static residual flag
        tSubSolver->set_compute_static_residual_flag( mComputeStaticResidual, mUseStaticResidual );
    }
}

//--------------------------------------------------------------------------------------------------------------------------

bool Nonlinear_Solver::get_compute_static_residual_flag()
{
    return mComputeStaticResidual;
}

//--------------------------------------------------------------------------------------------------------------------------

bool Nonlinear_Solver::get_use_static_residual_flag()
{
    return mUseStaticResidual;
}
