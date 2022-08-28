/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Waitbar.cpp
 *
 */

#include <thread>

// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "chronos.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::chronos::Waitbar",
		"[moris],[chronos],[cl_Waitbar],[Waitbar]")
{
	#include "chronos/cl_Waitbar.inc"

	REQUIRE(wall_time_microseconds >= big_container.size());
}

