/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eig_gen.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_eig_gen.hpp"

namespace moris
{
TEST_CASE( "moris::eig_gen", "[linalgebra],[eig_gen]" )
    {
    moris::Matrix< moris::DDCMat > A( 3, 3 );

    A( 0, 0 ) = -1.0; A( 0, 1 ) = 1.0; A( 0, 2 ) = 3.0;
    A( 1, 0 ) =  1.0; A( 1, 1 ) = 2.0; A( 1, 2 ) = 0.0;
    A( 2, 0 ) =  3.0; A( 2, 1 ) = 0.0; A( 2, 2 ) = 2.0;

    moris::Matrix< moris::DDCMat > eig_vec;
    moris::Matrix< moris::DDCMat > eig_val;

   moris::eig_gen( eig_val, eig_vec, A );
   // NOTE: This test evaluates the values, not the signs. To test/return certain signs and get consistent results,
   // one could require that the first element in each eigenvector be of a certain sign. Then, if you get an
   // eigenvector that does not follow this rule, you multiply it by -1.

   // NOTE: Eigen and Armadillo both return the eigenvalues and vectores in a different order

#ifdef MORIS_USE_EIGEN
    REQUIRE( moris::equal_to( std::abs( eig_val( 0, 0 )), 3.000000000000000 ) );
    REQUIRE( moris::equal_to( std::abs( eig_val( 1, 0 )), 2.000000000000000 ) );
    REQUIRE( moris::equal_to( std::abs( eig_val( 2, 0 )), 4.000000000000000 ) );

    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 0 ) ), 0.8451542547285165 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 1 ) ), 0.0000000000000000 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 2 ) ), 0.5345224838248487 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 0 ) ), 0.1690308509457034 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 1 ) ), 0.9486832980505139 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 2 ) ), 0.2672612419124244 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 0 ) ), 0.5070925528371100 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 1 ) ), 0.3162277660168378 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 2 ) ), 0.8017837257372730 ) );
#endif

#ifdef MORIS_USE_ARMA
    REQUIRE( moris::equal_to( std::abs( eig_val( 0, 0 )), 3.000000000000000 ) );
    REQUIRE( moris::equal_to( std::abs( eig_val( 1, 0 )), 4.000000000000000 ) );
    REQUIRE( moris::equal_to( std::abs( eig_val( 2, 0 )), 2.000000000000000 ) );

    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 0 ) ), 0.8451542547285165 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 2 ) ), 0.0000000000000000 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 1 ) ), 0.5345224838248487 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 0 ) ), 0.1690308509457034 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 2 ) ), 0.9486832980505139 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 1 ) ), 0.2672612419124244 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 0 ) ), 0.5070925528371100 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 2 ) ), 0.3162277660168378 ) );
    REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 1 ) ), 0.8017837257372730 ) );
#endif
    }
}

