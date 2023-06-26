/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_2D.cpp
 *
 */

#include "cl_HMR_T_Matrix.hpp"

namespace moris::hmr
{

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 2 >::evaluate_geometry_interpolation(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       aN )
    {
        // unpack xi and eta from input vector
        auto xi  = aXi( 0 );
        auto eta = aXi( 1 );

        // populate matrix with values
        aN( 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
        aN( 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
        aN( 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
        aN( 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
    }

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 2 >::evaluate_shape_function(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       aN )
    {
        // Initialize xi and eta matrices
        Matrix< DDRMat > tNxi( mLagrangeOrder + 1, 1 );
        Matrix< DDRMat > tNeta( mLagrangeOrder + 1, 1 );

        // evaluate contributions for xi and eta
        for ( uint i = 0; i <= mLagrangeOrder; i++ )
        {
            tNxi( i ) = this->lagrange_shape_1d( i, aXi( 0 ) );
        }
        for ( uint j = 0; j <= mLagrangeOrder; j++ )
        {
            tNeta( j ) = this->lagrange_shape_1d( j, aXi( 1 ) );
        }

        // create shape vector in correct order
        for ( uint iNodeIndex = 0; iNodeIndex < mNumberOfNodes; iNodeIndex++ )
        {
            aN( iNodeIndex ) = tNxi( mLagrangeIJK( 0, iNodeIndex ) ) * tNeta( mLagrangeIJK( 1, iNodeIndex ) );
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 2 >::get_child_corner_nodes(
            uint       aChildIndex,
            Matrix< DDRMat >& aXi )
    {
        switch ( aChildIndex )
        {
            case ( 0 ):
            {
                aXi( 0, 0 ) = -1.0;
                aXi( 1, 0 ) = 0.0;
                aXi( 2, 0 ) = 0.0;
                aXi( 3, 0 ) = -1.0;

                aXi( 0, 1 ) = -1.0;
                aXi( 1, 1 ) = -1.0;
                aXi( 2, 1 ) = 0.0;
                aXi( 3, 1 ) = 0.0;

                break;
            }
            case ( 1 ):
            {
                aXi( 0, 0 ) = 0.0;
                aXi( 1, 0 ) = 1.0;
                aXi( 2, 0 ) = 1.0;
                aXi( 3, 0 ) = 0.0;

                aXi( 0, 1 ) = -1.0;
                aXi( 1, 1 ) = -1.0;
                aXi( 2, 1 ) = 0.0;
                aXi( 3, 1 ) = 0.0;

                break;
            }

            case ( 2 ):
            {
                aXi( 0, 0 ) = -1.0;
                aXi( 1, 0 ) = 0.0;
                aXi( 2, 0 ) = 0.0;
                aXi( 3, 0 ) = -1.0;

                aXi( 0, 1 ) = 0.0;
                aXi( 1, 1 ) = 0.0;
                aXi( 2, 1 ) = 1.0;
                aXi( 3, 1 ) = 1.0;

                break;
            }

            case ( 3 ):
            {
                aXi( 0, 0 ) = 0.0;
                aXi( 1, 0 ) = 1.0;
                aXi( 2, 0 ) = 1.0;
                aXi( 3, 0 ) = 0.0;

                aXi( 0, 1 ) = 0.0;
                aXi( 1, 1 ) = 0.0;
                aXi( 2, 1 ) = 1.0;
                aXi( 3, 1 ) = 1.0;

                break;
            }
            default:
            {
                MORIS_ERROR( false, "invalid child index" );
            }
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

}
