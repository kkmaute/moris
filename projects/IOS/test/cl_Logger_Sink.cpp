// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::ios::Logger_Sink",
        "[moris],[ios],[cl_Logger_Sink],[Logger_Sink]")
{
    #include "cl_Logger_Sink/clog.inc" // snippets IOS
    REQUIRE(clog_msg1 == 24);
    REQUIRE(clog_msg2 == 7);

    #include "cl_Logger_Sink/cout.inc" // snippets IOS
    REQUIRE(cout_msg1 == 33);
    REQUIRE(cout_msg2 == 16);

    #include "cl_Logger_Sink/cerr.inc" // snippets IOS
    REQUIRE(cerr_msg1 == 25);
    REQUIRE(cerr_msg2 == 8);
}
