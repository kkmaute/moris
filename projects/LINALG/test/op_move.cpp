/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_move.cpp
 *
 */

#include <catch.hpp>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_move.hpp"
#include "op_equal_equal.hpp"    //ALG
#include "fn_all_true.hpp"
#include "fn_isempty.hpp"

#include <type_traits>
#include <utility>
#include <vector>

namespace moris
{
    TEST_CASE( "moris::op_move", "[linalgebra],[op_move]" )
    {
        Matrix< DDRMat > Amatrix( 50, 50, 0.0 );
        Matrix< DDRMat > Bmatrix( 50, 50, 0.0 );
        Matrix< DDRMat > Cmatrix( 50, 50 );

        Amatrix( 0, 0 ) = 1.0;
        Amatrix( 0, 1 ) = 2.0;
        Amatrix( 0, 2 ) = 3.0;
        Amatrix( 1, 0 ) = 4.0;
        Amatrix( 1, 1 ) = 5.0;
        Amatrix( 1, 2 ) = 6.0;
        Amatrix( 2, 0 ) = 7.0;
        Amatrix( 2, 1 ) = 8.0;
        Amatrix( 2, 2 ) = 9.0;

        Bmatrix = Amatrix.copy();

        Cmatrix = move( Bmatrix );

        REQUIRE( all_true( Cmatrix == Amatrix ) );

        REQUIRE( isempty( Bmatrix ) );
    }

    TEST_CASE( "moris::Matrix native move constructor", "[linalgebra],[move_semantics]" )
    {
        STATIC_REQUIRE( std::is_nothrow_move_constructible< Matrix< DDRMat > >::value
                        == std::is_nothrow_move_constructible< DDRMat >::value );

        // Large owning matrices avoid the backend's inline storage optimization.
        Matrix< DDRMat > tSource( 50, 50, 1.5 );
        tSource( 0, 0 )   = 42.0;
        tSource( 49, 49 ) = 99.0;

        auto* tOriginalData = tSource.data();

        Matrix< DDRMat > tMoved( std::move( tSource ) );

        REQUIRE( tMoved.data() == tOriginalData );
        REQUIRE( tMoved.n_rows() == 50 );
        REQUIRE( tMoved.n_cols() == 50 );
        REQUIRE( tMoved( 0, 0 ) == 42.0 );
        REQUIRE( tMoved( 49, 49 ) == 99.0 );
        REQUIRE( tMoved( 5, 5 ) == 1.5 );

        // A moved-from matrix must remain reusable; its dimensions are unspecified.
        tSource = Matrix< DDRMat >( 3, 3, 7.0 );

        REQUIRE( tSource.n_rows() == 3 );
        REQUIRE( tSource.n_cols() == 3 );
        REQUIRE( tSource( 0, 0 ) == 7.0 );
        REQUIRE( tMoved( 0, 0 ) == 42.0 );
    }

    TEST_CASE( "moris::Matrix native move assignment", "[linalgebra],[move_semantics]" )
    {
#ifdef MORIS_USE_EIGEN
        // A strict borrowed destination can reject a size-changing assignment.
        STATIC_REQUIRE_FALSE( std::is_nothrow_move_assignable< Matrix< DDRMat > >::value );
#else
        STATIC_REQUIRE( std::is_nothrow_move_assignable< Matrix< DDRMat > >::value
                        == std::is_nothrow_move_assignable< DDRMat >::value );
#endif

        Matrix< DDRMat > tSource( 50, 50, 2.0 );
        tSource( 2, 2 ) = 7.0;

        auto* tOriginalData = tSource.data();

        // Exercise replacement of an existing allocation.
        Matrix< DDRMat > tTarget( 30, 30, -1.0 );
        tTarget = std::move( tSource );

        REQUIRE( tTarget.data() == tOriginalData );
        REQUIRE( tTarget.n_rows() == 50 );
        REQUIRE( tTarget.n_cols() == 50 );
        REQUIRE( tTarget( 2, 2 ) == 7.0 );
        REQUIRE( tTarget( 0, 0 ) == 2.0 );

        tSource = Matrix< DDRMat >( 3, 3, 9.0 );

        REQUIRE( tSource.n_rows() == 3 );
        REQUIRE( tSource.n_cols() == 3 );
        REQUIRE( tSource( 0, 0 ) == 9.0 );
        REQUIRE( tTarget( 2, 2 ) == 7.0 );
    }

    TEST_CASE( "moris::Matrix copy preserves independent storage", "[linalgebra],[move_semantics]" )
    {
        Matrix< DDRMat > tSource( 50, 50, 2.0 );
        Matrix< DDRMat > tConstructed( tSource );
        Matrix< DDRMat > tAssigned( 30, 30, -1.0 );
        tAssigned = tSource;

        REQUIRE( all_true( tConstructed == tSource ) );
        REQUIRE( all_true( tAssigned == tSource ) );
        REQUIRE( tConstructed.data() != tSource.data() );
        REQUIRE( tAssigned.data() != tSource.data() );
        REQUIRE( tConstructed.data() != tAssigned.data() );

        tConstructed( 0, 0 ) = 3.0;
        tAssigned( 0, 0 )    = 4.0;

        REQUIRE( tSource( 0, 0 ) == 2.0 );
        REQUIRE( tConstructed( 0, 0 ) == 3.0 );
        REQUIRE( tAssigned( 0, 0 ) == 4.0 );
    }

    TEST_CASE( "moris::Matrix vector push_back with move semantics", "[linalgebra],[move_semantics]" )
    {
        std::vector< Matrix< DDRMat > > tMatrices;
        tMatrices.reserve( 1 );
        tMatrices.emplace_back( 50, 50, 1.0 );

        auto*            tFirstData        = tMatrices.front().data();
        const auto       tOriginalCapacity = tMatrices.capacity();
        Matrix< DDRMat > tSource( 50, 50, 2.0 );
        auto*            tSourceData = tSource.data();

        // Fill the capacity so the following push_back must relocate existing entries.
        while ( tMatrices.size() < tOriginalCapacity )
        {
            tMatrices.emplace_back( 50, 50, 3.0 );
        }

        tMatrices.push_back( std::move( tSource ) );

        REQUIRE( tMatrices.capacity() > tOriginalCapacity );
        REQUIRE( tMatrices.back().data() == tSourceData );
        REQUIRE( tMatrices.front().n_rows() == 50 );
        REQUIRE( tMatrices.front().n_cols() == 50 );
        REQUIRE( tMatrices.front()( 0, 0 ) == 1.0 );
        REQUIRE( tMatrices.back()( 0, 0 ) == 2.0 );

        // Vector may copy existing entries when the backend move can throw.
        if constexpr ( std::is_nothrow_move_constructible< Matrix< DDRMat > >::value )
        {
            REQUIRE( tMatrices.front().data() == tFirstData );
        }
    }

}    // namespace moris

#ifdef MORIS_USE_EIGEN
TEST_CASE( "Eigen move assignment from a borrowed view retains its owner", "[linalgebra],[move]" )
{
    moris::Matrix< moris::DDRMat > tOwner( 64, 64, 3.0 );
    auto*                          tStorage = tOwner.data();
    moris::Matrix< moris::DDRMat > tView( tStorage, 64, 64, false, true );
    tOwner = std::move( tView );
    CHECK( tOwner.data() == tStorage );
    CHECK( tOwner( 63, 63 ) == 3.0 );
    tOwner( 0, 0 ) = 7.0;
    CHECK( tStorage[ 0 ] == 7.0 );
}
#endif
