// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "Opt_Input.hpp" // OPT/src

// ---------------------------------------------------------------------
using namespace moris;

TEST_CASE( "[optimization]" )
{
    moris::real tObjective = 0.0;
    Matrix< DDRMat >  tADVs;

    SECTION( "GCMMA" )
    {
        // objective  = 3*x - (x^3)
        // constraint = x - 2
        // Initial value of x = 1.95
        // One possible solution to this problem is x = 2
        // This is a problem with one adv and to change the setup, one
        // needs to perform changes in opt_input.cpp

        moris::opt_input::create_solve_opt_problem( 1, tObjective, tADVs );

        REQUIRE( std::abs( tObjective - ( -2.0 ) ) < 1.0e-02 ); // check value of objective
        REQUIRE( std::abs( tADVs(0,0) - ( +2.0 ) ) < 1.0e-03 ); // check value of design variable, x
    }

    SECTION( "SQP" )
    {
        // objective  = 3*x - (x^3)
        // constraint = x - 2
        // Initial value of x = 1.95
        // One possible solution to this problem is x = 2
        // This is a problem with one adv and to change the setup, one
        // needs to perform changes in opt_input.cpp

        moris::opt_input::create_solve_opt_problem( 2, tObjective, tADVs );

        REQUIRE( std::abs( tObjective - ( -2.0 ) ) < 1.0e-08 ); // check value of objective
        REQUIRE( std::abs( tADVs(0,0) - ( +2.0 ) ) < 1.0e-08 ); // check value of design variable, x
    }

    SECTION( "SQP_GCMMA" )
    {
        // objective  = 3*x - (x^3)
        // constraint = x - 2
        // Initial value of x = 1.95
        // One possible solution to this problem is x = 2
        // This is a problem with one adv and to change the setup, one
        // needs to perform changes in opt_input.cpp

        moris::opt_input::create_solve_opt_problem( 3, tObjective, tADVs );

        REQUIRE( std::abs( tObjective - ( -2.0 ) ) < 1.0e-08 ); // check value of objective
        REQUIRE( std::abs( tADVs(0,0) - ( +2.0 ) ) < 1.0e-08 ); // check value of design variable, x
    }

    SECTION( "SWEEP" )
    {
        // objective  = 3*x - (x^3)
        // This is a problem with one adv. You start with ADV value 1.9, and end
        // at adv value 2.1. To change the setup, one needs to perform changes
        // in opt_input.cpp

        moris::opt_input::create_solve_opt_problem( 5, tObjective, tADVs );

        REQUIRE( std::abs( tObjective - ( -2.961 ) ) < 1.0e-14 ); // check value of objective
        REQUIRE( std::abs( tADVs(0,0) - ( +2.100 ) ) < 1.0e-14 ); // check value of design variable, x
    }
}

