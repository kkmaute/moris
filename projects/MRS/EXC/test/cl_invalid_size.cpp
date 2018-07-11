// C++ header files.
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
