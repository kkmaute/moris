// Third-party header files.
#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

extern moris::Logger gLogger;

// ----------------------------------------------------------------------------

TEST_CASE(
        "MORIS_LOG_TEST",
        "[moris],[ios],[cl_Logger],[MORIS_LOG],[log_test]")
{
        gLogger.log_info( "Test output message for MORIS Logger  || %-5i   ||   %-5.2e \n", 124, 0.00456789 );
        gLogger.log_info( "Test output message for MORIS Logger  || %-5i   ||   %-5.2e \n", 124, 0.00456789 );

        MORIS_LOG_INFO("Test output message for MORIS Logger  || %-5i   ||   %-5.2e \n", 122, 0.005555555);

        MORIS_LOG("Test output message for MORIS Logger  || %-5i   ||   %-5.2e \n", 121, 0.005555555);
}
