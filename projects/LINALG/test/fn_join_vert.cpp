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
    TEST_CASE( "Vertical concatenation handles all const and temporary combinations", "[LINALG][join_vert]" )
    {
        // Larger operands also exercise heap-backed storage. Evaluate on the
        // next statement, after temporary wrappers returned by either factory
        // have expired; ASan must find no borrowed temporary in the result.
        Matrix< DDRMat > tLeft( 64, 2, 2.0 );
        Matrix< DDRMat > tRight( 64, 2, 4.0 );

        auto tCheck = []( auto aLeftFactory, auto aRightFactory ) {
            auto                   tJoined = join_vert( aLeftFactory(), aRightFactory() );
            const Matrix< DDRMat > tResult = tJoined;
            REQUIRE( tResult.n_rows() == 128 );
            REQUIRE( tResult.n_cols() == 2 );

            for ( uint iColumn = 0; iColumn < 2; ++iColumn )
            {
                for ( uint iRow = 0; iRow < 64; ++iRow )
                {
                    CHECK( tResult( iRow, iColumn ) == 2.0 );
                    CHECK( tResult( 64 + iRow, iColumn ) == 4.0 );
                }
            }
        };

        auto tCheckLeft = [ & ]( auto aLeftFactory ) {
            tCheck( aLeftFactory, [ & ]() -> Matrix< DDRMat >& { return tRight; } );
            tCheck( aLeftFactory, [ & ]() -> const Matrix< DDRMat >& { return tRight; } );
            tCheck( aLeftFactory, []() { return Matrix< DDRMat >( 64, 2, 4.0 ); } );
            tCheck( aLeftFactory, []() -> const Matrix< DDRMat > { return Matrix< DDRMat >( 64, 2, 4.0 ); } );
        };

        tCheckLeft( [ & ]() -> Matrix< DDRMat >& { return tLeft; } );
        tCheckLeft( [ & ]() -> const Matrix< DDRMat >& { return tLeft; } );
        tCheckLeft( []() { return Matrix< DDRMat >( 64, 2, 2.0 ); } );
        tCheckLeft( []() -> const Matrix< DDRMat > { return Matrix< DDRMat >( 64, 2, 2.0 ); } );
    }

    TEST_CASE( "Vertical concatenation owns const rvalue operands", "[LINALG][join_vert]" )
    {
        // Moving a const reference must still snapshot both operands while a
        // temporary participates. Mutations expose a borrowed result even when
        // its storage remains alive, without relying on an allocator reusing it.
        Matrix< DDRMat >        tLeft       = { { 1.0, 2.0 } };
        Matrix< DDRMat >        tRight      = { { 3.0, 4.0 } };
        const Matrix< DDRMat >& tConstLeft  = tLeft;
        const Matrix< DDRMat >& tConstRight = tRight;
        auto                    tLeftMoved  = join_vert( std::move( tConstLeft ), tRight );
        auto                    tRightMoved = join_vert( tLeft, std::move( tConstRight ) );
        auto                    tBothMoved  = join_vert( std::move( tConstLeft ), std::move( tConstRight ) );
        tLeft.fill( 9.0 );
        tRight.fill( 8.0 );

        for ( const auto* tJoined : { &tLeftMoved, &tRightMoved, &tBothMoved } )
        {
            CHECK( ( *tJoined )( 0, 0 ) == 1.0 );
            CHECK( ( *tJoined )( 0, 1 ) == 2.0 );
            CHECK( ( *tJoined )( 1, 0 ) == 3.0 );
            CHECK( ( *tJoined )( 1, 1 ) == 4.0 );
        }
    }

    TEST_CASE( "Vertical concatenation retains native lvalue evaluation", "[LINALG][join_vert]" )
    {
        // Armadillo safely defers a join of persistent wrappers. Eigen's native
        // join is owning; preserve each backend's existing expression behavior.
        Matrix< DDRMat > tLeft         = { { 1.0 } };
        Matrix< DDRMat > tRight        = { { 2.0 } };
        auto             tJoined       = join_vert( tLeft, tRight );
        tLeft( 0 )                     = 3.0;
        const Matrix< DDRMat > tResult = tJoined;
#ifdef MORIS_USE_ARMA
        CHECK( tResult( 0 ) == 3.0 );
#else
        CHECK( tResult( 0 ) == 1.0 );
#endif
        CHECK( tResult( 1 ) == 2.0 );
    }
}    // namespace moris
