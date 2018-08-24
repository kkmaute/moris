// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::ios::to_stdio",
        "[moris],[ios],[fn_to_stdio],[to_stdio]")
{
    #include "fn_to_stdio.inc" // snippets IOS

    REQUIRE(read_mode == std::string("r"));
    REQUIRE(write_mode == std::string("w"));
    REQUIRE(append_mode == std::string("a"));
    REQUIRE(readwrite_mode == std::string("r+"));
    REQUIRE(trunc_mode == std::string("w+"));
    REQUIRE(update_mode == std::string("a+"));
}
