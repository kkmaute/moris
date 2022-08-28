/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_Trilinos.cpp
 *
 */

#include "cl_DLA_Preconditioner_Trilinos.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"

#include "cl_Stopwatch.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include "Teuchos_ParameterList.hpp"

#include "Ifpack.h"
#include "Ifpack_ValidParameters.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_LocalFilter.h"

#include "Ifpack_AdditiveSchwarz.h"

using namespace moris;
using namespace dla;

//-------------------------------------------------------------------------------

Preconditioner_Trilinos::Preconditioner_Trilinos()
{

}

//-------------------------------------------------------------------------------

Preconditioner_Trilinos::Preconditioner_Trilinos(
        const moris::ParameterList  aParameterlist,
        Linear_Problem            * aLinearSystem )
{
    this->initialize(
            aParameterlist,
            aLinearSystem );
}

//-------------------------------------------------------------------------------

Preconditioner_Trilinos::~Preconditioner_Trilinos()
{

}

//-------------------------------------------------------------------------------

void Preconditioner_Trilinos::initialize(
        const moris::ParameterList  aParameterlist,
        Linear_Problem            * aLinearSystem )
{
    // set parameter list
    mParameterList = aParameterlist;

    // check whether one preconditioner has been defined
    bool tIsIfpack = ! mParameterList.get< std::string >( "ifpack_prec_type" ).empty();
    bool tIsMl     = ! mParameterList.get< std::string >( "ml_prec_type" ).empty();

    if ( !tIsIfpack && !tIsMl )
    {
        mIsInitialized = false;
        return;
    }
    else
    {
        mIsInitialized = true;
    }

    // check that only one preconditioner is defined
    MORIS_ERROR( ( tIsIfpack && !tIsMl ) || ( !tIsIfpack && tIsMl ),
            "Preconditioner_Trilinos::initialize - One and only one preconditioner must be specified.\n");

    // store linear system
    mLinearSystem = aLinearSystem;
}

//-------------------------------------------------------------------------------

void Preconditioner_Trilinos::build(const sint & aIter)
{
    // check if preconditioner is initialized; if not do nothing
    if (!mIsInitialized)
    {
        return;
    }

    // check that linear system is defined
    MORIS_ERROR( mLinearSystem,
            "Preconditioner_Trilinos::build - Preconditioner has not been initialized.\n");

    // build Ifpack preconditioner
    if( ! mParameterList.get< std::string >( "ifpack_prec_type" ).empty() )
    {
        // initialize and build, and compute preconditioner in first iteration
        // or if preconditioner should not be reused
        if ( aIter == 1 || mParameterList.get< bool >( "prec_reuse" ) == false )
        {
            this->build_ifpack_preconditioner();
            this->compute_ifpack_preconditioner();
        }
        else
        {
            // recompute preconditioner if reused
            this->compute_ifpack_preconditioner();
        }
    }
    // build ml preconditioner
    else if( ! mParameterList.get< std::string >( "ml_prec_type" ).empty() )
    {
        // initialize and build, and compute preconditioner in first iteration
        // or if preconditioner should not be reused
        if ( aIter == 1 || mParameterList.get< bool >( "prec_reuse" ) == false )
        {
            this->build_ml_preconditioner();
            this->compute_ml_preconditioner();
        }
        else
        {
            // recompute preconditioner if reused
            this->compute_ml_preconditioner( true );
        }
    }
}

//-------------------------------------------------------------------------------

bool Preconditioner_Trilinos::exists()
{
    if (!mIsInitialized)
    {
        return false;
    }
    else
    {
        return !mIfPackPrec.is_null() || !mMlPrec.is_null();
    }
}

//-------------------------------------------------------------------------------

Teuchos::RCP< Epetra_Operator > Preconditioner_Trilinos::get_operator()
{
    MORIS_ERROR( this->exists(),
            "Preconditioner_Trilinos::get_prec - no preconditioner has been built.\n");

    // return pointer to preconditioner
    if ( !mIfPackPrec.is_null() )
    {
        return mIfPackPrec;
    }
    else
    {
        return mMlPrec;
    }
}

//-------------------------------------------------------------------------------

moris::sint Preconditioner_Trilinos::build_ifpack_preconditioner()
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;

    ParameterList tIfpackParameterlist;

    // To get all valid IFPACK parameters:
    // tIfpackParameterlist = Ifpack_GetValidParameters(); std::cout << tIfpackParameterlist;

    // Allocate an IFPACK factory.
    Ifpack Factory;

    // check that linear system is set
    MORIS_ERROR(mLinearSystem,
            "Preconditioner_Trilinos::build_ifpack_preconditioner - linear system not set.\n" );

    // Get pointer to operator
    Epetra_RowMatrix * tOperator = mLinearSystem->get_matrix()->get_matrix();

    // Get preconditioner type
    std::string PrecType = mParameterList.get< std::string >( "ifpack_prec_type" );

    // Get overlap across processors
    int OverlapLevel = mParameterList.get< moris::sint >( "overlap-level" ) ;

    // Create the preconditioner.
    mIfPackPrec = rcp ( Factory.Create ( PrecType, tOperator, OverlapLevel ) );

    // Specify local solver specific parameters
    if ( PrecType == "ILU" )
    {
        tIfpackParameterlist.set ( "fact: level-of-fill"     , mParameterList.get< moris::sint >( "fact: level-of-fill" ) );
        tIfpackParameterlist.set ( "fact: absolute threshold", mParameterList.get< moris::real >( "fact: absolute threshold" ) );
        tIfpackParameterlist.set ( "fact: relative threshold", mParameterList.get< moris::real >( "fact: relative threshold" ) );
        tIfpackParameterlist.set ( "fact: relax value"       , mParameterList.get< moris::real >( "fact: relax value" ) );
    }
    else if ( PrecType == "ILUT" )
    {
        tIfpackParameterlist.set ( "fact: ilut level-of-fill", mParameterList.get< moris::real >( "fact: ilut level-of-fill" ) );
        tIfpackParameterlist.set ( "fact: absolute threshold", mParameterList.get< moris::real >( "fact: absolute threshold" ) );
        tIfpackParameterlist.set ( "fact: relative threshold", mParameterList.get< moris::real >( "fact: relative threshold" ) );
        tIfpackParameterlist.set ( "fact: drop tolerance"    , mParameterList.get< moris::real >( "fact: drop tolerance" ) );
        tIfpackParameterlist.set ( "fact: relax value"       , mParameterList.get< moris::real >( "fact: relax value" )   );
    }
    else if ( PrecType == "IC" )
    {
        tIfpackParameterlist.set ( "fact: level-of-fill",      mParameterList.get< moris::sint >( "fact: level-of-fill" ) );
        tIfpackParameterlist.set ( "fact: absolute threshold", mParameterList.get< moris::real >( "fact: absolute threshold" ) );
        tIfpackParameterlist.set ( "fact: relative threshold", mParameterList.get< moris::real >( "fact: relative threshold" ) );
        tIfpackParameterlist.set ( "fact: drop tolerance"    , mParameterList.get< moris::real >( "fact: drop tolerance" ) );
    }
    else if ( PrecType == "ICT" )
    {
        tIfpackParameterlist.set ( "fact: ict level-of-fill",  mParameterList.get< moris::real >( "fact: ict level-of-fill" ) );
        tIfpackParameterlist.set ( "fact: absolute threshold", mParameterList.get< moris::real >( "fact: absolute threshold" ) );
        tIfpackParameterlist.set ( "fact: relative threshold", mParameterList.get< moris::real >( "fact: relative threshold" ) );
        tIfpackParameterlist.set ( "fact: drop tolerance"    , mParameterList.get< moris::real >( "fact: drop tolerance" ) );
        tIfpackParameterlist.set ( "fact: relax value"       , mParameterList.get< moris::real >( "fact: relax value" ) );
    }
    else if ( PrecType == "Amesos" )
    {
        tIfpackParameterlist.set ( "amesos: solver type", mParameterList.get< std::string > ( "amesos: solver type" ) );
    }
    else if ( PrecType == "point relaxation" )
    {
        tIfpackParameterlist.set ( "relaxation: type"                  , mParameterList.get< std::string >( "relaxation: type" ) );
        tIfpackParameterlist.set ( "relaxation: sweeps"                , mParameterList.get< moris::sint >( "relaxation: sweeps" ) );
        tIfpackParameterlist.set ( "relaxation: damping factor"        , mParameterList.get< moris::real >( "relaxation: damping factor" ) );
        tIfpackParameterlist.set ( "relaxation: zero starting solution", mParameterList.get< bool        >( "relaxation: zero starting solution" ) );

        tIfpackParameterlist.set ( "relaxation: min diagonal value"    , mParameterList.get< std::string >( "relaxation: min diagonal value" ) );
        tIfpackParameterlist.set ( "relaxation: backward mode"         , mParameterList.get< bool        >( "relaxation: backward mode" ) );
        tIfpackParameterlist.set ( "relaxation: use l1"                , mParameterList.get< bool        >( "relaxation: use l1" ) );
        tIfpackParameterlist.set ( "relaxation: l1 eta"                , mParameterList.get< moris::real >( "relaxation: l1 eta") );
    }
    else if ( PrecType == "block relaxation" )
    {
        tIfpackParameterlist.set ( "relaxation: type"                  , mParameterList.get< std::string >( "relaxation: type" ) );
        tIfpackParameterlist.set ( "relaxation: sweeps"                , mParameterList.get< moris::sint >( "relaxation: sweeps" ) );
        tIfpackParameterlist.set ( "relaxation: damping factor"        , mParameterList.get< moris::real >( "relaxation: damping factor" ) );
        tIfpackParameterlist.set ( "relaxation: zero starting solution", mParameterList.get< bool        >( "relaxation: zero starting solution" ) );

        tIfpackParameterlist.set ( "partitioner: type"                 , mParameterList.get< std::string >( "partitioner: type" ) );
        tIfpackParameterlist.set ( "partitioner: overlap"              , mParameterList.get< moris::sint >( "partitioner: overlap" ) );
        tIfpackParameterlist.set ( "partitioner: local parts"          , mParameterList.get< moris::sint >( "partitioner: local parts" ) );
        tIfpackParameterlist.set ( "partitioner: print level"          , mParameterList.get< moris::sint >( "partitioner: print level" ) );
        tIfpackParameterlist.set ( "partitioner: use symmetric graph"  , mParameterList.get< bool        >( "partitioner: use symmetric graph" ) );
    }
    else if ( PrecType == "SPARSKIT" )
    {
        tIfpackParameterlist.set ( "fact: sparskit: lfil"    , mParameterList.get< moris::sint >( "fact: sparskit: lfil" ) );
        tIfpackParameterlist.set ( "fact: sparskit: tol"     , mParameterList.get< moris::real >( "fact: sparskit: tol" ) );
        tIfpackParameterlist.set ( "fact: sparskit: droptol" , mParameterList.get< moris::real >( "fact: sparskit: droptol" ) );
        tIfpackParameterlist.set ( "fact: sparskit: permtol" , mParameterList.get< moris::real >( "fact: sparskit: permtol" ) );
        tIfpackParameterlist.set ( "fact: sparskit: alph"    , mParameterList.get< moris::real >( "fact: sparskit: alph" ) );
        tIfpackParameterlist.set ( "fact: sparskit: mbloc"   , mParameterList.get< moris::sint >( "fact: sparskit: mbloc" ) );
        tIfpackParameterlist.set ( "fact: sparskit: type"    , mParameterList.get< std::string >( "fact: sparskit: type" ) );
    }
    else if ( PrecType == "Krylov" )
    {
        tIfpackParameterlist.set ( "krylov: iterations"             ,  mParameterList.get< moris::real >( "krylov: iterations" ) );
        tIfpackParameterlist.set ( "krylov: tolerance"              ,  mParameterList.get< std::string >( "krylov: tolerance" ) );
        tIfpackParameterlist.set ( "krylov: solver"                 ,  mParameterList.get< moris::sint >( "krylov: solver" ) );
        tIfpackParameterlist.set ( "krylov: preconditioner"         ,  mParameterList.get< moris::sint >( "krylov: preconditioner" ) );
        tIfpackParameterlist.set ( "krylov: number of sweeps"       ,  mParameterList.get< moris::sint >( "krylov: number of sweeps" ) );
        tIfpackParameterlist.set ( "krylov: block size"             ,  mParameterList.get< moris::sint >( "krylov: block size" ) );
        tIfpackParameterlist.set ( "krylov: damping parameter"      ,  mParameterList.get< moris::real >( "krylov: damping parameter" ) );
        tIfpackParameterlist.set ( "krylov: zero starting solution" ,  mParameterList.get< moris::sint >( "krylov: zero starting solution" ) );
    }
    else
    {
        MORIS_ERROR( false, "Ifpack preconditioner type %s not implemented", PrecType.c_str() );
    }

    tIfpackParameterlist.set ( "schwarz: combine mode"     , mParameterList.get< std::string >( "schwarz: combine mode" ) );
    tIfpackParameterlist.set ( "schwarz: compute condest"  , mParameterList.get< bool        >( "schwarz: compute condest" ) );
    tIfpackParameterlist.set ( "schwarz: filter singletons", mParameterList.get< bool        >( "schwarz: filter singletons" ) );
    tIfpackParameterlist.set ( "schwarz: reordering type"  , mParameterList.get< std::string >( "schwarz: reordering type" ) );

    // start timer
    tic tTimer;

    // Set the parameters.
    IFPACK_CHK_ERR( mIfPackPrec->SetParameters( tIfpackParameterlist ) );

    // Initialize the preconditioner.
    IFPACK_CHK_ERR( mIfPackPrec->Initialize() );

    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    moris::real tElapsedTimeMax = max_all( tElapsedTime );

    if ( par_rank() == 0 )
    {
        MORIS_LOG_INFO( "SOL: Total time to initialize ifpack preconditioner is %5.3f seconds.",
                ( double ) tElapsedTimeMax / 1000);
    }

    return 0;
}

//-------------------------------------------------------------------------------

moris::sint Preconditioner_Trilinos::compute_ifpack_preconditioner( bool tRecompute )
{
    // start timer
    tic tTimer;

    // Compute preconditioner for current matrix
    IFPACK_CHK_ERR( mIfPackPrec->Compute() );

    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    moris::real tElapsedTimeMax = max_all( tElapsedTime );

    if ( par_rank() == 0 )
    {
        MORIS_LOG_INFO( "SOL: Total time to compute ifpack preconditioner is %5.3f seconds.",
                ( double ) tElapsedTimeMax / 1000);
    }

    return 0;
}

//-------------------------------------------------------------------------------

moris::sint Preconditioner_Trilinos::build_ml_preconditioner()
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;

    // check that linear system is set
    MORIS_ERROR(mLinearSystem,
            "Preconditioner_Trilinos::build_ml_preconditioner - linear system not set.\n" );

    // ml parameter list
    ParameterList tMlParams;

    // Get pointer to operator
    Epetra_RowMatrix * tOperator = mLinearSystem->get_matrix()->get_matrix();

    // Get default setting
    ML_Epetra::SetDefaults (  mParameterList.get< std::string >( "ml_prec_type" ), tMlParams );

    // General Parameters
    tMlParams.set ( "PDE equations", mParameterList.get< moris::sint >( "PDE equations" ) );

    tMlParams.set ( "ML output"                 , mParameterList.get< moris::sint >( "ML output" ) );
    tMlParams.set ( "ML print initial list"     , mParameterList.get< moris::sint >( "ML print initial list" ) );
    tMlParams.set ( "ML print final list"       , mParameterList.get< moris::sint >( "ML print final list" ) );
    tMlParams.set ( "print unused"              , mParameterList.get< moris::sint >( "print unused" ) );
    tMlParams.set ( "eigen-analysis: type"      , mParameterList.get< std::string >( "eigen-analysis: type" ) );
    tMlParams.set ( "eigen-analysis: iterations", mParameterList.get< moris::sint >( "eigen-analysis: iterations" ) );

    // Multigrid Cycle Parameters
    tMlParams.set ( "cycle applications"      , mParameterList.get< moris::sint >( "cycle applications" ) );
    tMlParams.set ( "max levels"              , mParameterList.get< moris::sint >( "max levels" ) );
    tMlParams.set ( "increasing or decreasing", mParameterList.get< std::string >( "increasing or decreasing" ) );
    tMlParams.set ( "prec type"               , mParameterList.get< std::string >( "prec type" ) );

    // Aggregation and Prolongator Parameters
    tMlParams.set ( "aggregation: type"            , mParameterList.get< std::string >( "aggregation: type" ) );
    tMlParams.set ( "aggregation: threshold"       , mParameterList.get< moris::real >( "aggregation: threshold" ) );
    tMlParams.set ( "aggregation: damping factor"  , mParameterList.get< moris::real >( "aggregation: damping factor" ) );
    tMlParams.set ( "aggregation: smoothing sweeps", mParameterList.get< moris::sint >( "aggregation: smoothing sweeps" ) );

    tMlParams.set ( "aggregation: use tentative restriction", mParameterList.get< bool >       ( "aggregation: use tentative restriction" ) );
    tMlParams.set ( "aggregation: symmetrize"               , mParameterList.get< bool >       ( "aggregation: symmetrize" ) );
    tMlParams.set ( "aggregation: global aggregates"        , mParameterList.get< moris::sint >( "aggregation: global aggregates" ) );
    tMlParams.set ( "aggregation: local aggregates"         , mParameterList.get< moris::sint >( "aggregation: local aggregates" ) );
    tMlParams.set ( "aggregation: nodes per aggregate"      , mParameterList.get< moris::sint >( "aggregation: nodes per aggregate" ) );

    tMlParams.set ( "energy minimization: enable" , mParameterList.get< bool >       ( "energy minimization: enable" ) );
    tMlParams.set ( "energy minimization: type"   , mParameterList.get< moris::sint >( "energy minimization: type" ) );
    tMlParams.set ( "energy minimization: droptol", mParameterList.get< moris::real >( "energy minimization: droptol" ) );
    tMlParams.set ( "energy minimization: cheap"  , mParameterList.get< bool >       ( "energy minimization: cheap" ) );

    // Smoother Parameters
    tMlParams.set ( "smoother: type"          , mParameterList.get< std::string >( "smoother: type" ) );
    tMlParams.set ( "smoother: sweeps"        , mParameterList.get< moris::sint >( "smoother: sweeps" ) );
    tMlParams.set ( "smoother: damping factor", mParameterList.get< moris::real >( "smoother: damping factor" ) );
    tMlParams.set ( "smoother: pre or post"   , mParameterList.get< std::string >( "smoother: pre or post" ) );

    // Aztec smoothing
    tMlParams.set ( "smoother: Aztec as solver",false );

    // Ifpack smoothing
    tMlParams.set ( "smoother: ifpack level-of-fill"     , mParameterList.get< moris::real >( "smoother: ifpack level-of-fill" ) );
    tMlParams.set ( "smoother: ifpack overlap"           , mParameterList.get< moris::sint >( "smoother: ifpack overlap" ) );
    tMlParams.set ( "smoother: ifpack absolute threshold", mParameterList.get< moris::real >( "smoother: ifpack absolute threshold" ) );
    tMlParams.set ( "smoother: ifpack relative threshold", mParameterList.get< moris::real >( "smoother: ifpack relative threshold" ) );

    // Coarsest Grid Parameters
    tMlParams.set ( "coarse: type"          , mParameterList.get< std::string >( "coarse: type" ) );
    tMlParams.set ( "coarse: max size"      , mParameterList.get< moris::sint >( "coarse: max size" ) );
    tMlParams.set ( "coarse: pre or post"   , mParameterList.get< std::string >( "coarse: pre or post" ) );
    tMlParams.set ( "coarse: sweeps"        , mParameterList.get< moris::sint >( "coarse: sweeps" ) );
    tMlParams.set ( "coarse: damping factor", mParameterList.get< moris::real >( "coarse: damping factor" ) );

    // Load-Balancing Parameters
    tMlParams.set ( "repartition: enable"     , mParameterList.get< moris::sint >( "repartition: enable" ) );
    tMlParams.set ( "repartition: partitioner", mParameterList.get< std::string >( "repartition: partitioner" ) );

    // Analysis Parameters
    tMlParams.set ( "analyze memory"              , mParameterList.get< bool >       ( "analyze memory" ) );
    tMlParams.set ( "viz: enable"                 , mParameterList.get< bool >       ( "viz: enable" ) );
    tMlParams.set ( "viz: output format"          , mParameterList.get< std::string >( "viz: output format" ) );
    tMlParams.set ( "viz: print starting solution", mParameterList.get< bool >       ( "viz: print starting solution" ) );

    // Smoothed Aggregation and the Null Space
    tMlParams.set ( "null space: type"               , mParameterList.get< std::string >( "null space: type" ) );
    tMlParams.set ( "null space: dimension"          , mParameterList.get< moris::sint >( "null space: dimension" ) );
    tMlParams.set ( "null space: vectors to compute" , mParameterList.get< moris::sint >( "null space: vectors to compute" ) );
    tMlParams.set ( "null space: add default vectors", mParameterList.get< bool >       ( "null space: add default vectors" ) );

    // start timer
    tic tTimer;

    // create preconditioner
    mMlPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*tOperator, tMlParams));

    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    moris::real tElapsedTimeMax = max_all( tElapsedTime );

    if ( par_rank() == 0 )
    {
        MORIS_LOG_INFO( "SOL: Total time to initialize ml preconditioner is %5.3f seconds.",
                ( double ) tElapsedTimeMax / 1000);
    }

    return 0;
}

//-------------------------------------------------------------------------------

moris::sint Preconditioner_Trilinos::compute_ml_preconditioner( bool tRecompute )
{
    // start timer
    tic tTimer;

    // Compute or recompute preconditioner for current matrix
    if ( tRecompute )
    {
        mMlPrec->ReComputePreconditioner();
    }
    else
    {
        mMlPrec->ComputePreconditioner();
    }
    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    moris::real tElapsedTimeMax = max_all( tElapsedTime );

    if ( par_rank() == 0 )
    {
        MORIS_LOG_INFO( "SOL: Total time to compute ml preconditioner is %5.3f seconds.",
                ( double ) tElapsedTimeMax / 1000);
    }

    return 0;
}

//-------------------------------------------------------------------------------

