/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_3D.cpp
 *
 */

#include "cl_HMR_T_Matrix.hpp"

namespace moris::hmr
{

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 3 >::evaluate_geometry_interpolation(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       aN )
    {
        // unpack xi and eta from input vector
        auto xi   = aXi( 0 );
        auto eta  = aXi( 1 );
        auto zeta = aXi( 2 );

        // populate output matrix
        aN( 0 ) = -( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aN( 1 ) = ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aN( 2 ) = -( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aN( 3 ) = ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 ) * 0.125;
        aN( 4 ) = ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
        aN( 5 ) = -( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
        aN( 6 ) = ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 ) * 0.125;
        aN( 7 ) = -( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 ) * 0.125;
    }

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 3 >::evaluate_shape_function(
            const Matrix< DDRMat >& aXi,
            Matrix< DDRMat >&       aN )
    {
        // evaluate contributions for xi and eta and zeta
        Matrix< DDRMat > tNxi( mLagrangeOrder + 1, 1 );
        Matrix< DDRMat > tNeta( mLagrangeOrder + 1, 1 );
        Matrix< DDRMat > tNzeta( mLagrangeOrder + 1, 1 );

        for ( uint i = 0; i <= mLagrangeOrder; i++ )
        {
            tNxi( i ) = this->lagrange_shape_1d( i, aXi( 0 ) );
        }

        for ( uint j = 0; j <= mLagrangeOrder; j++ )
        {
            tNeta( j ) = this->lagrange_shape_1d( j, aXi( 1 ) );
        }

        for ( uint k = 0; k <= mLagrangeOrder; k++ )
        {
            tNzeta( k ) = this->lagrange_shape_1d( k, aXi( 2 ) );
        }

        // create shape vector in correct order
        for ( uint iNodeIndex = 0; iNodeIndex < mNumberOfNodes; iNodeIndex++ )
        {
            aN( iNodeIndex ) =
                    tNxi( mLagrangeIJK( 0, iNodeIndex ) )
                            * tNeta( mLagrangeIJK( 1, iNodeIndex ) )
                            * tNzeta( mLagrangeIJK( 2, iNodeIndex ) );
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 3 >::get_child_corner_nodes(
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
                aXi( 4, 0 ) = -1.0;
                aXi( 5, 0 ) = 0.0;
                aXi( 6, 0 ) = 0.0;
                aXi( 7, 0 ) = -1.0;

                aXi( 0, 1 ) = -1.0;
                aXi( 1, 1 ) = -1.0;
                aXi( 2, 1 ) = 0.0;
                aXi( 3, 1 ) = 0.0;
                aXi( 4, 1 ) = -1.0;
                aXi( 5, 1 ) = -1.0;
                aXi( 6, 1 ) = 0.0;
                aXi( 7, 1 ) = 0.0;

                aXi( 0, 2 ) = -1.0;
                aXi( 1, 2 ) = -1.0;
                aXi( 2, 2 ) = -1.0;
                aXi( 3, 2 ) = -1.0;
                aXi( 4, 2 ) = 0.0;
                aXi( 5, 2 ) = 0.0;
                aXi( 6, 2 ) = 0.0;
                aXi( 7, 2 ) = 0.0;

                break;
            }
            case ( 1 ):
            {
                aXi( 0, 0 ) = 0.0;
                aXi( 1, 0 ) = 1.0;
                aXi( 2, 0 ) = 1.0;
                aXi( 3, 0 ) = 0.0;
                aXi( 4, 0 ) = 0.0;
                aXi( 5, 0 ) = 1.0;
                aXi( 6, 0 ) = 1.0;
                aXi( 7, 0 ) = 0.0;

                aXi( 0, 1 ) = -1.0;
                aXi( 1, 1 ) = -1.0;
                aXi( 2, 1 ) = 0.0;
                aXi( 3, 1 ) = 0.0;
                aXi( 4, 1 ) = -1.0;
                aXi( 5, 1 ) = -1.0;
                aXi( 6, 1 ) = 0.0;
                aXi( 7, 1 ) = 0.0;

                aXi( 0, 2 ) = -1.0;
                aXi( 1, 2 ) = -1.0;
                aXi( 2, 2 ) = -1.0;
                aXi( 3, 2 ) = -1.0;
                aXi( 4, 2 ) = 0.0;
                aXi( 5, 2 ) = 0.0;
                aXi( 6, 2 ) = 0.0;
                aXi( 7, 2 ) = 0.0;

                break;
            }
            case ( 2 ):
            {
                aXi( 0, 0 ) = -1.0;
                aXi( 1, 0 ) = 0.0;
                aXi( 2, 0 ) = 0.0;
                aXi( 3, 0 ) = -1.0;
                aXi( 4, 0 ) = -1.0;
                aXi( 5, 0 ) = 0.0;
                aXi( 6, 0 ) = 0.0;
                aXi( 7, 0 ) = -1.0;

                aXi( 0, 1 ) = 0.0;
                aXi( 1, 1 ) = 0.0;
                aXi( 2, 1 ) = 1.0;
                aXi( 3, 1 ) = 1.0;
                aXi( 4, 1 ) = 0.0;
                aXi( 5, 1 ) = 0.0;
                aXi( 6, 1 ) = 1.0;
                aXi( 7, 1 ) = 1.0;

                aXi( 0, 2 ) = -1.0;
                aXi( 1, 2 ) = -1.0;
                aXi( 2, 2 ) = -1.0;
                aXi( 3, 2 ) = -1.0;
                aXi( 4, 2 ) = 0.0;
                aXi( 5, 2 ) = 0.0;
                aXi( 6, 2 ) = 0.0;
                aXi( 7, 2 ) = 0.0;

                break;
            }
            case ( 3 ):
            {
                aXi( 0, 0 ) = 0.0;
                aXi( 1, 0 ) = 1.0;
                aXi( 2, 0 ) = 1.0;
                aXi( 3, 0 ) = 0.0;
                aXi( 4, 0 ) = 0.0;
                aXi( 5, 0 ) = 1.0;
                aXi( 6, 0 ) = 1.0;
                aXi( 7, 0 ) = 0.0;

                aXi( 0, 1 ) = 0.0;
                aXi( 1, 1 ) = 0.0;
                aXi( 2, 1 ) = 1.0;
                aXi( 3, 1 ) = 1.0;
                aXi( 4, 1 ) = 0.0;
                aXi( 5, 1 ) = 0.0;
                aXi( 6, 1 ) = 1.0;
                aXi( 7, 1 ) = 1.0;

                aXi( 0, 2 ) = -1.0;
                aXi( 1, 2 ) = -1.0;
                aXi( 2, 2 ) = -1.0;
                aXi( 3, 2 ) = -1.0;
                aXi( 4, 2 ) = 0.0;
                aXi( 5, 2 ) = 0.0;
                aXi( 6, 2 ) = 0.0;
                aXi( 7, 2 ) = 0.0;

                break;
            }
            case ( 4 ):
            {
                aXi( 0, 0 ) = -1.0;
                aXi( 1, 0 ) = 0.0;
                aXi( 2, 0 ) = 0.0;
                aXi( 3, 0 ) = -1.0;
                aXi( 4, 0 ) = -1.0;
                aXi( 5, 0 ) = 0.0;
                aXi( 6, 0 ) = 0.0;
                aXi( 7, 0 ) = -1.0;

                aXi( 0, 1 ) = -1.0;
                aXi( 1, 1 ) = -1.0;
                aXi( 2, 1 ) = 0.0;
                aXi( 3, 1 ) = 0.0;
                aXi( 4, 1 ) = -1.0;
                aXi( 5, 1 ) = -1.0;
                aXi( 6, 1 ) = 0.0;
                aXi( 7, 1 ) = 0.0;

                aXi( 0, 2 ) = 0.0;
                aXi( 1, 2 ) = 0.0;
                aXi( 2, 2 ) = 0.0;
                aXi( 3, 2 ) = 0.0;
                aXi( 4, 2 ) = 1.0;
                aXi( 5, 2 ) = 1.0;
                aXi( 6, 2 ) = 1.0;
                aXi( 7, 2 ) = 1.0;

                break;
            }
            case ( 5 ):
            {
                aXi( 0, 0 ) = 0.0;
                aXi( 1, 0 ) = 1.0;
                aXi( 2, 0 ) = 1.0;
                aXi( 3, 0 ) = 0.0;
                aXi( 4, 0 ) = 0.0;
                aXi( 5, 0 ) = 1.0;
                aXi( 6, 0 ) = 1.0;
                aXi( 7, 0 ) = 0.0;

                aXi( 0, 1 ) = -1.0;
                aXi( 1, 1 ) = -1.0;
                aXi( 2, 1 ) = 0.0;
                aXi( 3, 1 ) = 0.0;
                aXi( 4, 1 ) = -1.0;
                aXi( 5, 1 ) = -1.0;
                aXi( 6, 1 ) = 0.0;
                aXi( 7, 1 ) = 0.0;

                aXi( 0, 2 ) = 0.0;
                aXi( 1, 2 ) = 0.0;
                aXi( 2, 2 ) = 0.0;
                aXi( 3, 2 ) = 0.0;
                aXi( 4, 2 ) = 1.0;
                aXi( 5, 2 ) = 1.0;
                aXi( 6, 2 ) = 1.0;
                aXi( 7, 2 ) = 1.0;

                break;
            }
            case ( 6 ):
            {
                aXi( 0, 0 ) = -1.0;
                aXi( 1, 0 ) = 0.0;
                aXi( 2, 0 ) = 0.0;
                aXi( 3, 0 ) = -1.0;
                aXi( 4, 0 ) = -1.0;
                aXi( 5, 0 ) = 0.0;
                aXi( 6, 0 ) = 0.0;
                aXi( 7, 0 ) = -1.0;

                aXi( 0, 1 ) = 0.0;
                aXi( 1, 1 ) = 0.0;
                aXi( 2, 1 ) = 1.0;
                aXi( 3, 1 ) = 1.0;
                aXi( 4, 1 ) = 0.0;
                aXi( 5, 1 ) = 0.0;
                aXi( 6, 1 ) = 1.0;
                aXi( 7, 1 ) = 1.0;

                aXi( 0, 2 ) = 0.0;
                aXi( 1, 2 ) = 0.0;
                aXi( 2, 2 ) = 0.0;
                aXi( 3, 2 ) = 0.0;
                aXi( 4, 2 ) = 1.0;
                aXi( 5, 2 ) = 1.0;
                aXi( 6, 2 ) = 1.0;
                aXi( 7, 2 ) = 1.0;

                break;
            }
            case ( 7 ):
            {
                aXi( 0, 0 ) = 0.0;
                aXi( 1, 0 ) = 1.0;
                aXi( 2, 0 ) = 1.0;
                aXi( 3, 0 ) = 0.0;
                aXi( 4, 0 ) = 0.0;
                aXi( 5, 0 ) = 1.0;
                aXi( 6, 0 ) = 1.0;
                aXi( 7, 0 ) = 0.0;

                aXi( 0, 1 ) = 0.0;
                aXi( 1, 1 ) = 0.0;
                aXi( 2, 1 ) = 1.0;
                aXi( 3, 1 ) = 1.0;
                aXi( 4, 1 ) = 0.0;
                aXi( 5, 1 ) = 0.0;
                aXi( 6, 1 ) = 1.0;
                aXi( 7, 1 ) = 1.0;

                aXi( 0, 2 ) = 0.0;
                aXi( 1, 2 ) = 0.0;
                aXi( 2, 2 ) = 0.0;
                aXi( 3, 2 ) = 0.0;
                aXi( 4, 2 ) = 1.0;
                aXi( 5, 2 ) = 1.0;
                aXi( 6, 2 ) = 1.0;
                aXi( 7, 2 ) = 1.0;

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
