// C++ header files.
#include <sstream>

// Third-party header files.
#include <catch.hpp>

// moris header files.
#include "exceptions.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::exceptions::tic_toc_error",
		"[moris],[exceptions],[cl_tic_toc_error],[tic_toc_error]")
{
	#include "exceptions/cl_tic_toc_error.inc"

	REQUIRE_THROWS_AS(throw error, moris::exceptions::tic_toc_error);
	REQUIRE(error_msg.str() == std::string(error.what()));
}
