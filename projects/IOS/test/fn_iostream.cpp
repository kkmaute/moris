// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::ios::iostream",
        "[moris],[ios],[fn_iostream],[iostream]")
{
    SECTION("moris::clog")
    {
        #include "fn_iostream/clog.inc" // snippets IOS
    }

    SECTION("moris::cout")
    {
        #include "fn_iostream/cout.inc" // snippets IOS
    }

    SECTION("moris::cerr")
    {
        #include "fn_iostream/cerr.inc" // snippets IOS
    }
}
