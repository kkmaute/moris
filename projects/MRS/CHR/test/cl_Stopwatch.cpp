/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Stopwatch.cpp
 *
 */

#include <thread>

// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "chronos.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::chronos::Stopwatch",
        "[moris],[chronos],[cl_Stopwatch],[Stopwatch]")
{
        #include "chronos/cl_Stopwatch/Stopwatch.inc"

    SECTION("moris::chronos::Stopwatch::reset")
    {
        #include "chronos/cl_Stopwatch/reset.inc"

        // Assert only 1 second has passed, but not 2, because we reset the timer.
        REQUIRE(wall_time_hours == 0);
        REQUIRE(wall_time_minutes == 0);
        REQUIRE((wall_time_seconds >= 1 && wall_time_seconds < 2));
        REQUIRE((wall_time_milliseconds >= 1.0e+03 && wall_time_milliseconds < 2.0e+03));
        REQUIRE((wall_time_microseconds >= 1.0e+06 && wall_time_microseconds < 2.0e+06));
        REQUIRE((wall_time_nanoseconds >= 1.0e+09 && wall_time_nanoseconds < 2.0e+09));
    }

    SECTION("moris::chronos::Stopwatch::stop")
    {
        #include "chronos/cl_Stopwatch/stop.inc"

        // Assert 2ms have passed, but not 3ms.
        REQUIRE(wall_time_hours == 0);
        REQUIRE(wall_time_minutes == 0);
        REQUIRE(wall_time_seconds == 0);
        REQUIRE((wall_time_milliseconds >= 1 && wall_time_milliseconds < 2));
        REQUIRE((wall_time_microseconds >= 1.0e+03 && wall_time_microseconds < 2.0e+03));
        REQUIRE((wall_time_nanoseconds >= 1.0e+06 && wall_time_nanoseconds < 2.0e+06));
    }

    SECTION("moris::chronos::Stopwatch::is_stopped")
    {
        #include "chronos/cl_Stopwatch/is_stopped.0.inc"

        REQUIRE(t.is_stopped() == true);

        #include "chronos/cl_Stopwatch/is_stopped.1.inc"

        REQUIRE(t.is_stopped() == false);
    }

    SECTION("moris::chronos::Stopwatch::resume")
    {
        #include "chronos/cl_Stopwatch/resume.inc"

        // Assert only 1ms has passed, but not 2ms, because we stopped the timer.
        REQUIRE(wall_time_hours == 0);
        REQUIRE(wall_time_minutes == 0);
        REQUIRE(wall_time_seconds == 0);
        REQUIRE((wall_time_milliseconds >= 1 && wall_time_milliseconds < 2));
        REQUIRE((wall_time_microseconds >= 1e+03 && wall_time_microseconds < 2e+03));
        REQUIRE((wall_time_nanoseconds >= 1e+06 && wall_time_nanoseconds < 2e+06));
    }

    SECTION("moris::chronos::Stopwatch::toc")
    {
        #include "chronos/cl_Stopwatch/toc.inc"

        // Assert 1 second has passed.
        REQUIRE(wall_time_hours == 0);
        REQUIRE(wall_time_minutes == 0);
        REQUIRE(wall_time_seconds >= 1);
        REQUIRE(wall_time_milliseconds >= 1.0e+03);
        REQUIRE(wall_time_microseconds >= 1.0e+06);
        REQUIRE(wall_time_nanoseconds >= 1.0e+09);
    }

    SECTION("moris::chronos::Stopwatch::elapsed")
    {
        #include "chronos/cl_Stopwatch/elapsed.inc"

        // Assert formatting is different for each case.
        REQUIRE(wall_time_1.find("wall") != std::string::npos);
        REQUIRE(wall_time_1.find("user") != std::string::npos);
        REQUIRE(wall_time_1.find("system") != std::string::npos);
        REQUIRE(wall_time_1.find("CPU") != std::string::npos);
        REQUIRE(wall_time_1.find("%") != std::string::npos);

        REQUIRE(wall_time_2.find("wall") != std::string::npos);
        REQUIRE(wall_time_2.find("user") == std::string::npos);
        REQUIRE(wall_time_2.find("system") == std::string::npos);
        REQUIRE(wall_time_2.find("CPU") == std::string::npos);
        REQUIRE(wall_time_2.find("%") == std::string::npos);

        REQUIRE(wall_time_3.find("wall") == std::string::npos);
        REQUIRE(wall_time_3.find("user") != std::string::npos);
        REQUIRE(wall_time_3.find("system") == std::string::npos);
        REQUIRE(wall_time_3.find("CPU") == std::string::npos);
        REQUIRE(wall_time_3.find("%") == std::string::npos);

        REQUIRE(wall_time_4.find("wall") == std::string::npos);
        REQUIRE(wall_time_4.find("user") == std::string::npos);
        REQUIRE(wall_time_4.find("system") != std::string::npos);
        REQUIRE(wall_time_4.find("CPU") == std::string::npos);
        REQUIRE(wall_time_4.find("%") == std::string::npos);

        REQUIRE(wall_time_5.find("wall") == std::string::npos);
        REQUIRE(wall_time_5.find("user") == std::string::npos);
        REQUIRE(wall_time_5.find("system") == std::string::npos);
        REQUIRE(wall_time_5.find("CPU") != std::string::npos);
        REQUIRE(wall_time_5.find("%") == std::string::npos);

        REQUIRE(wall_time_6.find("wall") == std::string::npos);
        REQUIRE(wall_time_6.find("user") == std::string::npos);
        REQUIRE(wall_time_6.find("system") == std::string::npos);
        REQUIRE(wall_time_6.find("CPU") == std::string::npos);
        REQUIRE(wall_time_6.find("%") != std::string::npos);
    }
}

