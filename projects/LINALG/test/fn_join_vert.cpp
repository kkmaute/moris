/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include "catch.hpp"
#include <utility>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_join_vert.hpp"

namespace moris
{
    TEST_CASE( "Vertical concatenation preserves rows and operand lifetimes", "[LINALG][join_vert]" )
    {
        // Contact geometry maps concatenate leader/follower rows. Check a
        // rectangular layout, an empty accumulator and overlapping assignment.
        Matrix< DDSMat > tLeader   = { { 1, 2 }, { 3, 4 } };
        Matrix< DDSMat > tFollower = { { 5, 6 } };
        auto             tDeferred = join_vert( tLeader, tFollower );
        Matrix< DDSMat > tJoined   = tDeferred;
        REQUIRE( tJoined.n_rows() == 3 );
        REQUIRE( tJoined.n_cols() == 2 );
        CHECK( tJoined( 0, 1 ) == 2 );
        CHECK( tJoined( 1, 0 ) == 3 );
        CHECK( tJoined( 2, 0 ) == 5 );

        Matrix< DDSMat > tEmpty;
        tEmpty = join_vert( tEmpty, tFollower );
        REQUIRE( tEmpty.n_rows() == 1 );
        CHECK( tEmpty( 0, 1 ) == 6 );
        tLeader = join_vert( tLeader, tFollower );
        REQUIRE( tLeader.n_rows() == 3 );
        CHECK( tLeader( 2, 1 ) == 6 );

        // Return concatenations of temporary wrappers from a function. Each
        // result must own its values after the temporary operands are gone.
        auto tBothTemporary = []() {
            return join_vert( Matrix< DDSMat >{ { 7, 8 } }, Matrix< DDSMat >{ { 9, 10 } } );
        }();
        CHECK( tBothTemporary( 1, 1 ) == 10 );
        auto tLeftTemporary = join_vert( Matrix< DDSMat >{ { 7, 8 } }, tFollower );
        CHECK( tLeftTemporary( 0, 0 ) == 7 );
        auto tRightTemporary = join_vert( tFollower, Matrix< DDSMat >{ { 9, 10 } } );
        CHECK( tRightTemporary( 1, 0 ) == 9 );

        // A fixed-size input cannot also be the return type: joining two
        // three-vectors produces an owning six-vector on either backend.
        Matrix< F31RMat > tFixed( 3, 1, 2.0 );
        auto              tFixedJoined = join_vert( std::move( tFixed ), Matrix< F31RMat >( 3, 1, 4.0 ) );
        REQUIRE( tFixedJoined.n_rows() == 6 );
        REQUIRE( tFixedJoined.n_cols() == 1 );
        CHECK( tFixedJoined( 0 ) == 2.0 );
        CHECK( tFixedJoined( 5 ) == 4.0 );
    }
}    // namespace moris
