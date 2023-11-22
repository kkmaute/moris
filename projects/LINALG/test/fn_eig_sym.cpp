/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eig_sym.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_eig_sym.hpp"

namespace moris
{
TEST_CASE( "moris::eig_sym", "[linalgebra],[eig_sym]" )
    {
    Matrix< DDRMat > A( 3, 3 );

    A( 0, 0 ) = -1.0; A( 0, 1 ) = 3.0; A( 0, 2 ) = 0.0;
    A( 1, 0 ) =  3.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 2.0;
    A( 2, 0 ) =  0.0; A( 2, 1 ) = 2.0; A( 2, 2 ) = 4.0;

    Matrix< DDRMat > eig_vec;
    Matrix< DDRMat > eig_val;

    moris::eig_sym( eig_val, eig_vec, A );

    // NOTE: This test evalues the values, not the signs. To test/return certain signs and get consistent results,
    // one could require that the first element in each eigenvector be of a certain sign. Then, if you get an
    // eigenvector that does not follow this rule, you multiply it by -1.

    REQUIRE( moris::equal_to( std::abs( eig_val( 0, 0 ) ), 2.341200659105694 ) );
    REQUIRE( moris::equal_to( std::abs( eig_val( 1, 0 ) ), 3.043564371011120 ) );
    REQUIRE( moris::equal_to( std::abs( eig_val( 2, 0 ) ), 7.297636288094574 ) );

    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 0 ) ), 0.905449906553528 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 1 ) ), 0.304846459001434 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 2 ) ), 0.295345734955654 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 0 ) ), 0.404796670485593 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 1 ) ), 0.410888760082367 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 2 ) ), 0.816890495967332 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 0 ) ), 0.127671932256022 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 1 ) ), 0.859208393390255 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 2 ) ), 0.495440021033577 ) );
    }
}

