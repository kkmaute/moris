/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_query_test.cpp
 *
 */

#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

// query function definitions
#include <string>

// Tracer Enums to understand keywords
#include "cl_Tracer_Enums.hpp"

#include "cl_Query.hpp"

// ----------------------------------------------------------------------------

TEST_CASE( "TreeQuery", "[moris],[ios],[query],[tree],[TreeQuery]")
{
    // specify file to read
    std::string tFileNameRead = "gClock.log";

    // specify file to output to
    std::string tFileNameWrite = "gClock_Tree.log";

    // create query
    moris::ios::Query tQuery;
    tQuery.initialize(tFileNameRead);

    // execute tree query
    tQuery.tree_query(tFileNameWrite);

}

TEST_CASE( "TimeSolverQuery", "[moris],[ios],[query],[table],[TimeSolverQuery]")
{
    // specify file to read
    std::string tFileNameRead = "gClock.log";

    // specify file to output to
    std::string tFileNameWrite = "gClock_TimeSolver.log";

    // create query
    moris::ios::Query tQuery;
    tQuery.initialize(tFileNameRead);

    // execute table query for arbitrary time solvers,
    tQuery.table_query(tFileNameWrite,
                       moris::EntityBase::TimeSolver,
                       moris::EntityType::Arbitrary,
                       moris::EntityAction::Solve );
}

TEST_CASE( "NonLinearSolverQuery", "[moris],[ios],[query],[table],[NonLinearSolverQuery]")
{
    // specify file to read
    std::string tFileNameRead = "gClock.log";

    // specify file to output to
    std::string tFileNameWrite = "gClock_NonLinearSolver.log";

    // create query
    moris::ios::Query tQuery;
    tQuery.initialize(tFileNameRead);

    // execute table query for arbitrary non-linear solvers
    tQuery.table_query(tFileNameWrite,
                       moris::EntityBase::NonLinearSolver,
                       moris::EntityType::Arbitrary,
                       moris::EntityAction::Solve );
}

TEST_CASE( "NonLinearProblemQuery", "[moris],[ios],[query],[table],[NonLinearProblemQuery]")
{
    // specify file to read
    std::string tFileNameRead = "gClock.log";

    // specify file to output to
    std::string tFileNameWrite = "gClock_NonLinearProblem.log";

    // create query
    moris::ios::Query tQuery;
    tQuery.initialize(tFileNameRead);

    // execute table query for arbitrary non-linear solvers
    tQuery.table_query(tFileNameWrite,
                       moris::EntityBase::NonLinearProblem,
                       moris::EntityType::Arbitrary,
                       moris::EntityAction::Arbitrary );
}

