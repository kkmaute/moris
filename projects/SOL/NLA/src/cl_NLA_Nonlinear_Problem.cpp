/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Problem.cpp
 *
 */
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include <ctime>

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_Tracer.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

//-----------------------------------------------------------------------------

Nonlinear_Problem::Nonlinear_Problem(
        sol::SOL_Warehouse*     aSolverWarehouse,
        Solver_Interface*       aSolverInterface,
        sol::Dist_Vector*       aFullVector,
        const moris::sint       aNonlinearSolverManagerIndex,
        const bool              aBuildLinerSystemFlag,
        const enum sol::MapType aMapType )
        : mFullVector( aFullVector )
        , mBuildLinerSystemFlag( aBuildLinerSystemFlag )
        , mMapType( aMapType )
        , mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
    mSolverWarehouse = aSolverWarehouse;

    mSolverInterface = aSolverInterface;

    //    if( mMapType == sol::MapType::Petsc )
    //    {
    //        // Initialize petsc solvers
    //        PetscInitializeNoArguments();
    //    }

    const Vector< enum MSI::Dof_Type >& tRequestedDofTypes = mSolverInterface->get_requested_dof_types();

    // delete pointers if they already exist
    this->delete_pointers();

    // Build Matrix vector factory
    sol::Matrix_Vector_Factory tMatFactory( mMapType );

    MORIS_LOG_SPEC( "Total number of DOFs for non-linear system",
            sum_all( aSolverInterface->get_my_local_global_map( tRequestedDofTypes ).numel() ) );

    // create map object FIXME ask linear problem for map
    mMap = tMatFactory.create_map( aSolverInterface->get_my_local_global_map( tRequestedDofTypes ),
            aSolverInterface->get_my_local_global_overlapping_map( tRequestedDofTypes ) );

    // create map object FIXME ask linear problem for map
    mMapFull = tMatFactory.create_full_map(
            aSolverInterface->get_my_local_global_map(),
            aSolverInterface->get_my_local_global_overlapping_map() );

    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        // create solver factory
        Solver_Factory tSolFactory;

        mLinearProblem = tSolFactory.create_linear_system(
                aSolverInterface,
                mSolverWarehouse,
                mMap,
                mMapFull,
                mMapType );
    }
}

//-----------------------------------------------------------------------------

Nonlinear_Problem::Nonlinear_Problem(
        Solver_Interface*       aSolverInterface,
        const moris::sint       aNonlinearSolverManagerIndex,
        const bool              aBuildLinerSystemFlag,
        const enum sol::MapType aMapType )
        : mBuildLinerSystemFlag( aBuildLinerSystemFlag )
        , mMapType( aMapType )
        , mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
    mSolverInterface = aSolverInterface;

    if ( mMapType == sol::MapType::Petsc )
    {
#ifdef MORIS_HAVE_PETSC
        // Initialize petsc solvers
        PetscInitializeNoArguments();
#else
        MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
    }

    // Build Matrix vector factory
    sol::Matrix_Vector_Factory tMatFactory( mMapType );

    mMap = tMatFactory.create_full_map(
            aSolverInterface->get_my_local_global_map(),
            aSolverInterface->get_my_local_global_overlapping_map() );

    uint tNumRHMS = aSolverInterface->get_num_rhs();

    // full vector
    mFullVector = tMatFactory.create_vector(
            aSolverInterface,
            mMap,
            tNumRHMS );

    mSolverInterface->set_solution_vector( mFullVector );

    mFullVector->vec_put_scalar( 0.0 );

    // delete pointers if they already exist
    this->delete_pointers();

    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        MORIS_LOG_INFO( "Build linear problem with index %-5i", mNonlinearSolverManagerIndex );

        // create solver factory
        Solver_Factory tSolFactory;

        mLinearProblem = tSolFactory.create_linear_system( aSolverInterface, mMapType );
    }

    // set flag that interface has been set
    mIsLeaderSystem = true;
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::set_interface( Solver_Interface* aSolverInterface )
{
}

//-----------------------------------------------------------------------------

Nonlinear_Problem::~Nonlinear_Problem()
{
    this->delete_pointers();

    delete mMap;
    mMap = nullptr;

    delete mMapFull;
    mMapFull = nullptr;

    if ( mIsLeaderSystem )
    {
        delete mFullVector;
    }

    if ( mIsLeaderSystem )
    {
        if ( mMapType == sol::MapType::Petsc )
        {
#ifdef MORIS_HAVE_PETSC
            PetscFinalize();
#else
            MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
        }
    }
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::delete_pointers()
{
    delete ( mLinearProblem );

    mLinearProblem = nullptr;
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::update_fem_model()
{
    mSolverInterface->update_model();
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::build_linearized_problem(
        const bool& aRebuildJacobian,
        const bool& aCombinedResJacAssembly,
        const sint  aNonLinearIt )
{
    Tracer tTracer( "NonLinearProblem", "Build" );

    if ( aCombinedResJacAssembly )
    {
        // build residual and jacobian
        mLinearProblem->assemble_residual_and_jacobian();
    }
    else
    {
        // build jacobian only if requested
        if ( aRebuildJacobian )
        {
            mLinearProblem->assemble_jacobian();
        }

        // build residual
        mLinearProblem->assemble_residual();
    }

    // in case of sensitivity analysis and staggered solver: assemble contribution contributions of
    // of secondary dof types to residual
    if ( !mSolverInterface->is_forward_analysis() )
    {
        //        std::cout << "need fix in Nonlinear_Problem::build_linearized_problem \n";

        Vector< enum MSI::Dof_Type > tSecDofTypes = mMyNonLinSolver->get_sec_dof_type_union();

        // in case of secondary dof types
        if ( tSecDofTypes.size() != 0 )
        {
            Vector< enum MSI::Dof_Type > tDofTypeUnion = mMyNonLinSolver->get_dof_type_union();

            mSolverInterface->set_requested_dof_types( tSecDofTypes );
            mSolverInterface->set_secondary_dof_types( tDofTypeUnion );

            mLinearProblem->assemble_staggered_residual_contribution();

            mSolverInterface->set_requested_dof_types( tDofTypeUnion );
            mSolverInterface->set_secondary_dof_types( tSecDofTypes );
        }
    }
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::build_linearized_problem(
        const bool& aRebuildJacobian,
        const sint  aNonLinearIt,
        const sint  aRestart )
{
    Tracer tTracer( "NonLinearProblem", "Build" );

    this->restart_from_sol_vec( aRestart );

    if ( aRebuildJacobian )
    {
        mLinearProblem->assemble_jacobian();
    }

    mLinearProblem->assemble_residual();
}

//-----------------------------------------------------------------------------

sol::Dist_Vector*
Nonlinear_Problem::get_full_vector()
{
    return mFullVector;
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::extract_my_values(
        const moris::uint&                 aNumIndices,
        const moris::Matrix< DDSMat >&     aGlobalBlockRows,
        const moris::uint&                 aBlockRowOffsets,
        Vector< moris::Matrix< DDRMat > >& LHSValues )
{
    mFullVector->extract_my_values(
            aNumIndices,
            aGlobalBlockRows,
            aBlockRowOffsets,
            LHSValues );
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::print_sol_vec( const sint aNonLinearIt )
{
    char NonLinNum[ 100 ];
    std::sprintf( NonLinNum, "NonLIt.%04u", aNonLinearIt );

    char SolVector[ 100 ];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector, ".h5\0" );

    mFullVector->save_vector_to_HDF5( SolVector );
}

//-----------------------------------------------------------------------------

void Nonlinear_Problem::restart_from_sol_vec( const sint aRestart )
{
    char NonLinNum[ 100 ];
    std::sprintf( NonLinNum, "NonLIt.%04u", aRestart );

    char SolVector[ 100 ];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector, ".h5\0" );

    mFullVector->read_vector_from_HDF5( SolVector );
}

//--------------------------------------------------------------------------------------------------

void Nonlinear_Problem::set_time_value(
        const moris::real& aLambda,
        moris::uint        aPos )
{
    mSolverInterface->set_time_value( aLambda, aPos );
}

//--------------------------------------------------------------------------------------------------

real Nonlinear_Problem::get_static_residual_norm()
{
    return mLinearProblem->compute_static_residual_norm();
}
