/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Logger.cpp
 *
 */

#include <mpi.h>
#include <iostream>

#include "cl_Logger.hpp"
#include "catch.hpp"

TEST_CASE("Test the logger","[LOGGER][!throws]")
{
    // This is not an exception just a informational log, therefore this should not throw.
    REQUIRE_NOTHROW(XTK_INFO<<"Test the logger");

    // std::cout does not raise an exception
    REQUIRE_NOTHROW(std::cout<<"Error warning");

    int tProcRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);

    // See if it works in parallel
    if (tProcRank == 0)
    {
        REQUIRE_NOTHROW(XTK_INFO<<"My rank is "<< tProcRank);
    }

    else
    {
        REQUIRE_NOTHROW(XTK_INFO<<"I am not processor 0");
    }

}

