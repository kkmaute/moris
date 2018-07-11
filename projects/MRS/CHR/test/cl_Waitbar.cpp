// C++ header files.
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
