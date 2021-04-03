/*
 * cl_DLA_Preconditioner_Trilinos.cpp
 *
 *  Created on: Feb 06, 2020
 *      Author: schmidt
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
#include "Ifpack_Preconditioner.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_LocalFilter.h"

#include "Ifpack_AdditiveSchwarz.h"

using namespace moris;
using namespace dla;

Preconditioner_Trilinos::Preconditioner_Trilinos(
        const moris::ParameterList aParameterlist,
        Linear_Problem * aLinearSystem ) : mLinearSystem( aLinearSystem )
{
    mParameterList = aParameterlist;
}

//-------------------------------------------------------------------------------

Preconditioner_Trilinos::~Preconditioner_Trilinos()
{

}

//-------------------------------------------------------------------------------

moris::sint Preconditioner_Trilinos::build_ifpack_preconditioner()
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;

    ParameterList tIfpackParameterlist;

    // Allocate an IFPACK factory.  The object contains no data, only
    // the Create() method for creating preconditioners.
    Ifpack Factory;

    // Get pointer to operator
    Epetra_RowMatrix * tOperator = mLinearSystem->get_matrix()->get_matrix();

    // Create the preconditioner.  For the list of PrecType values that
    // Create() accepts, please check the IFPACK documentation.
    std::string PrecType = mParameterList.get< std::string >( "ifpack_prec_type" );

    int OverlapLevel = mParameterList.get< moris::sint >( "overlap-level" ) ;

    mPreconditioner = rcp ( Factory.Create ( PrecType, tOperator, OverlapLevel ) );

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
        tIfpackParameterlist.set ( "fact: relax value"       , mParameterList.get< moris::real >( "fact: relax value" )   );
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
        tIfpackParameterlist.set ( "partitioner: root node"            , mParameterList.get< moris::sint >( "partitioner: root node" ) );
        tIfpackParameterlist.set ( "partitioner: use symmetric graph"  , mParameterList.get< bool        >( "partitioner: use symmetric graph" ) );
    }
    else
    {
        MORIS_ERROR( false, "Ifpack preconditioner type not implemented" );
    }

    tIfpackParameterlist.set ( "schwarz: combine mode"     , mParameterList.get< std::string >( "schwarz: combine mode" ) );
    tIfpackParameterlist.set ( "schwarz: compute condest"  , mParameterList.get< bool        >( "schwarz: compute condest" ) );
    tIfpackParameterlist.set ( "schwarz: filter singletons", mParameterList.get< bool        >( "schwarz: filter singletons" ) );
    tIfpackParameterlist.set ( "schwarz: reordering type"  , mParameterList.get< std::string >( "schwarz: reordering type" ) );

    // start timer
    tic tTimer;

    // Set the parameters.
    IFPACK_CHK_ERR( mPreconditioner->SetParameters( tIfpackParameterlist ) );

    // Initialize the preconditioner.
    IFPACK_CHK_ERR( mPreconditioner->Initialize() );

    // Build the preconditioner, by looking at the values of the matrix.
    IFPACK_CHK_ERR( mPreconditioner->Compute() );

    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    moris::real tElapsedTimeMax = max_all( tElapsedTime );

    if ( par_rank() == 0 )
    {
        MORIS_LOG_INFO( "SOL: Total time to build ifpack preconditioner is %5.3f seconds.",
                ( double ) tElapsedTimeMax / 1000);
    }

    std::cout << *mPreconditioner;

    return 0;
}

//-------------------------------------------------------------------------------

moris::sint Preconditioner_Trilinos::build_ml_preconditioner()
{
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Get pointer to operator
    Epetra_RowMatrix * tOperator = mLinearSystem->get_matrix()->get_matrix();

    mlParams.set ( "PDE equations", mParameterList.get< moris::sint >( "PDE equations" ) );

    // General Parameters
    mlParams.set ( "ML output"                 , mParameterList.get< moris::sint >( "ML output" ) );
    mlParams.set ( "print unused"              , mParameterList.get< moris::sint >( "print unused" ) );
    mlParams.set ( "ML print initial list"     , mParameterList.get< moris::sint >( "ML print initial list" ) );
    mlParams.set ( "ML print final list"       , mParameterList.get< moris::sint >( "ML print final list" ) );
    mlParams.set ( "eigen-analysis: type"      , mParameterList.get< std::string >( "eigen-analysis: type" ) );
    mlParams.set ( "eigen-analysis: iterations", mParameterList.get< moris::sint >( "eigen-analysis: iterations" ) );

    // Multigrid Cycle Parameters
    mlParams.set ( "cycle applications"      , mParameterList.get< moris::sint >( "cycle applications" ) );
    mlParams.set ( "max levels"              , mParameterList.get< moris::sint >( "max levels" ) );
    mlParams.set ( "increasing or decreasing", mParameterList.get< std::string >( "increasing or decreasing" ) );
    mlParams.set ( "prec type"               , mParameterList.get< std::string >( "prec type" ) );

    // Aggregation and Prolongator Parameters
    mlParams.set ( "aggregation: type"            , mParameterList.get< std::string >( "aggregation: type" ) );
    mlParams.set ( "aggregation: threshold"       , mParameterList.get< moris::real >( "aggregation: threshold" ) );
    mlParams.set ( "aggregation: damping factor"  , mParameterList.get< moris::real >( "aggregation: damping factor" ) );
    mlParams.set ( "aggregation: smoothing sweeps", mParameterList.get< moris::sint >( "aggregation: smoothing sweeps" ) );

    mlParams.set ( "aggregation: use tentative restriction", mParameterList.get< bool >( "aggregation: use tentative restriction" ) );
    mlParams.set ( "aggregation: symmetrize"               , mParameterList.get< bool >( "aggregation: symmetrize" ) );
    mlParams.set ( "aggregation: global aggregates"        , mParameterList.get< moris::sint >( "aggregation: global aggregates" ) );
    mlParams.set ( "aggregation: local aggregates"         , mParameterList.get< moris::sint >( "aggregation: local aggregates" ) );
    mlParams.set ( "aggregation: nodes per aggregate"      , mParameterList.get< moris::sint >( "aggregation: nodes per aggregate" ) );

    mlParams.set ( "energy minimization: enable" , mParameterList.get< bool >( "energy minimization: enable" ) );
    mlParams.set ( "energy minimization: type"   , mParameterList.get< moris::sint >( "energy minimization: type" ) );
    mlParams.set ( "energy minimization: droptol", mParameterList.get< moris::real >( "energy minimization: droptol" ) );
    mlParams.set ( "energy minimization: cheap"  , mParameterList.get< bool >( "energy minimization: cheap" ) );

    // Smoother Parameters
    mlParams.set ( "smoother: type"          , mParameterList.get< std::string >( "smoother: type" ) );
    mlParams.set ( "smoother: sweeps"        , mParameterList.get< moris::sint >( "smoother: sweeps" ) );
    mlParams.set ( "smoother: damping factor", mParameterList.get< moris::real >( "smoother: damping factor" ) );
    mlParams.set ( "smoother: pre or post"   , mParameterList.get< std::string >( "smoother: pre or post" ) );

    // Aztec smoothing
    mlParams.set ( "smoother: Aztec as solver",false );

    // Ifpack smoothing
    mlParams.set ( "smoother: ifpack level-of-fill"     , mParameterList.get< moris::real >( "smoother: ifpack level-of-fill" ) );
    mlParams.set ( "smoother: ifpack overlap"           , mParameterList.get< moris::sint >( "smoother: ifpack overlap" ) );
    mlParams.set ( "smoother: ifpack absolute threshold", mParameterList.get< moris::real >( "smoother: ifpack absolute threshold" ) );
    mlParams.set ( "smoother: ifpack relative threshold", mParameterList.get< moris::real >( "smoother: ifpack relative threshold" ) );

    // Coarsest Grid Parameters
    mlParams.set ( "coarse: type"          , mParameterList.get< std::string >( "coarse: type" ) );
    mlParams.set ( "coarse: max size"      , mParameterList.get< moris::sint >( "coarse: max size" ) );
    mlParams.set ( "coarse: pre or post"   , mParameterList.get< std::string >( "coarse: pre or post" ) );
    mlParams.set ( "coarse: sweeps"        , mParameterList.get< moris::sint >( "coarse: sweeps" ) );
    mlParams.set ( "coarse: damping factor", mParameterList.get< moris::real >( "coarse: damping factor" ) );

    // Load-Balancing Parameters
    mlParams.set ( "repartition: enable"     , mParameterList.get< moris::sint >( "repartition: enable" ) );
    mlParams.set ( "repartition: partitioner", mParameterList.get< std::string >( "repartition: partitioner" ) );

    // Analysis Parameters
    mlParams.set ( "analyze memory"              , mParameterList.get< bool >( "analyze memory" ) );
    mlParams.set ( "viz: enable"                 , mParameterList.get< bool >( "viz: enable" ) );
    mlParams.set ( "viz: output format"          , mParameterList.get< std::string >( "viz: output format" ) );
    mlParams.set ( "viz: print starting solution", mParameterList.get< bool >( "viz: print starting solution" ) );

    // Smoothed Aggregation and the Null Space
    mlParams.set ( "null space: type"               , mParameterList.get< std::string >( "null space: type" ) );
    mlParams.set ( "null space: dimension"          , mParameterList.get< moris::sint >( "null space: dimension" ) );
    //mlParams.set ( "null space: vectors"            , mParameterList.get< moris::sint >( "null space: vectors" ) );
    mlParams.set ( "null space: vectors to compute" , mParameterList.get< moris::sint >( "null space: vectors to compute" ) );
    mlParams.set ( "null space: add default vectors", mParameterList.get< bool >( "null space: add default vectors" ) );

    // start timer
    tic tTimer;

    mMlPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*tOperator, mlParams));

    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
    moris::real tElapsedTimeMax = max_all( tElapsedTime );

    if ( par_rank() == 0 )
    {
        MORIS_LOG_INFO( "SOL: Total time to build ml preconditioner is %5.3f seconds.",
                ( double ) tElapsedTimeMax / 1000);
    }

    return 0;
}

//-------------------------------------------------------------------------------



