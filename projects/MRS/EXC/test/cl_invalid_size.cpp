/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_invalid_size.cpp
 *
 */

#include <sstream>

// Third-party header files.
#include <catch.hpp>

// moris header files.
#include "exceptions.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::exceptions::invalid_size",
		"[moris],[except],[cl_invalid_size],[invalid_size]")
{
	#include "exceptions/cl_invalid_size.inc"

	REQUIRE_THROWS_AS(throw error, moris::exceptions::invalid_size);
	REQUIRE(error_msg.str() == std::string(error.what()));
}

