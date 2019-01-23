// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "MORIS_LOG",
        "[moris],[ios],[cl_Logger],[MORIS_LOG]")
{
    SECTION("MORIS_LOG_FUNCTION")
    {
//        #include "cl_Logger/log_function.inc" // snippets IOS
    }

    SECTION("MORIS_LOG")
    {
        #include "cl_Logger/log.inc" // snippets IOS
    }

    SECTION("MORIS_LOG_TRACE")
    {
        #include "cl_Logger/log_trace.inc" // snippets IOS
    }

    SECTION("MORIS_LOG_DEBUG")
    {
        #include "cl_Logger/log_debug.inc" // snippets IOS
    }

    SECTION("MORIS_LOG_INFO")
    {
        #include "cl_Logger/log_info.inc" // snippets IOS
    }

    SECTION("MORIS_LOG_WARNING")
    {
        #include "cl_Logger/log_warning.inc" // snippets IOS
    }

    SECTION("MORIS_LOG_ERROR")
    {
        #include "cl_Logger/log_error.inc" // snippets IOS
    }

    SECTION("MORIS_LOG_FATAL")
    {
        #include "cl_Logger/log_fatal.inc" // snippets IOS
    }
}
