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

    template<>
    Matrix< DDRMat > T_Matrix< 3 >::get_supporting_points( const uint aOrder )
    {
        // the following lines are needed to get the interpolation points
        Matrix< DDRMat > aXihat( 3, std::pow( aOrder + 1, 3 ) );
        switch ( aOrder )
        {
            case 1:
            {
                // hex 8
                aXihat( 0, 0 ) = -1.000000;
                aXihat( 1, 0 ) = -1.000000;
                aXihat( 2, 0 ) = -1.000000;
                aXihat( 0, 1 ) = 1.000000;
                aXihat( 1, 1 ) = -1.000000;
                aXihat( 2, 1 ) = -1.000000;
                aXihat( 0, 2 ) = 1.000000;
                aXihat( 1, 2 ) = 1.000000;
                aXihat( 2, 2 ) = -1.000000;
                aXihat( 0, 3 ) = -1.000000;
                aXihat( 1, 3 ) = 1.000000;
                aXihat( 2, 3 ) = -1.000000;
                aXihat( 0, 4 ) = -1.000000;
                aXihat( 1, 4 ) = -1.000000;
                aXihat( 2, 4 ) = 1.000000;
                aXihat( 0, 5 ) = 1.000000;
                aXihat( 1, 5 ) = -1.000000;
                aXihat( 2, 5 ) = 1.000000;
                aXihat( 0, 6 ) = 1.000000;
                aXihat( 1, 6 ) = 1.000000;
                aXihat( 2, 6 ) = 1.000000;
                aXihat( 0, 7 ) = -1.000000;
                aXihat( 1, 7 ) = 1.000000;
                aXihat( 2, 7 ) = 1.000000;
                break;
            }
            case 2:
            {
                // hex 27
                aXihat( 0, 0 )  = -1.000000;
                aXihat( 1, 0 )  = -1.000000;
                aXihat( 2, 0 )  = -1.000000;
                aXihat( 0, 1 )  = 1.000000;
                aXihat( 1, 1 )  = -1.000000;
                aXihat( 2, 1 )  = -1.000000;
                aXihat( 0, 2 )  = 1.000000;
                aXihat( 1, 2 )  = 1.000000;
                aXihat( 2, 2 )  = -1.000000;
                aXihat( 0, 3 )  = -1.000000;
                aXihat( 1, 3 )  = 1.000000;
                aXihat( 2, 3 )  = -1.000000;
                aXihat( 0, 4 )  = -1.000000;
                aXihat( 1, 4 )  = -1.000000;
                aXihat( 2, 4 )  = 1.000000;
                aXihat( 0, 5 )  = 1.000000;
                aXihat( 1, 5 )  = -1.000000;
                aXihat( 2, 5 )  = 1.000000;
                aXihat( 0, 6 )  = 1.000000;
                aXihat( 1, 6 )  = 1.000000;
                aXihat( 2, 6 )  = 1.000000;
                aXihat( 0, 7 )  = -1.000000;
                aXihat( 1, 7 )  = 1.000000;
                aXihat( 2, 7 )  = 1.000000;
                aXihat( 0, 8 )  = 0.000000;
                aXihat( 1, 8 )  = -1.000000;
                aXihat( 2, 8 )  = -1.000000;
                aXihat( 0, 9 )  = 1.000000;
                aXihat( 1, 9 )  = 0.000000;
                aXihat( 2, 9 )  = -1.000000;
                aXihat( 0, 10 ) = 0.000000;
                aXihat( 1, 10 ) = 1.000000;
                aXihat( 2, 10 ) = -1.000000;
                aXihat( 0, 11 ) = -1.000000;
                aXihat( 1, 11 ) = 0.000000;
                aXihat( 2, 11 ) = -1.000000;
                aXihat( 0, 12 ) = -1.000000;
                aXihat( 1, 12 ) = -1.000000;
                aXihat( 2, 12 ) = 0.000000;
                aXihat( 0, 13 ) = 1.000000;
                aXihat( 1, 13 ) = -1.000000;
                aXihat( 2, 13 ) = 0.000000;
                aXihat( 0, 14 ) = 1.000000;
                aXihat( 1, 14 ) = 1.000000;
                aXihat( 2, 14 ) = 0.000000;
                aXihat( 0, 15 ) = -1.000000;
                aXihat( 1, 15 ) = 1.000000;
                aXihat( 2, 15 ) = 0.000000;
                aXihat( 0, 16 ) = 0.000000;
                aXihat( 1, 16 ) = -1.000000;
                aXihat( 2, 16 ) = 1.000000;
                aXihat( 0, 17 ) = 1.000000;
                aXihat( 1, 17 ) = 0.000000;
                aXihat( 2, 17 ) = 1.000000;
                aXihat( 0, 18 ) = 0.000000;
                aXihat( 1, 18 ) = 1.000000;
                aXihat( 2, 18 ) = 1.000000;
                aXihat( 0, 19 ) = -1.000000;
                aXihat( 1, 19 ) = 0.000000;
                aXihat( 2, 19 ) = 1.000000;
                aXihat( 0, 20 ) = 0.000000;
                aXihat( 1, 20 ) = 0.000000;
                aXihat( 2, 20 ) = 0.000000;
                aXihat( 0, 21 ) = 0.000000;
                aXihat( 1, 21 ) = 0.000000;
                aXihat( 2, 21 ) = -1.000000;
                aXihat( 0, 22 ) = 0.000000;
                aXihat( 1, 22 ) = 0.000000;
                aXihat( 2, 22 ) = 1.000000;
                aXihat( 0, 23 ) = -1.000000;
                aXihat( 1, 23 ) = 0.000000;
                aXihat( 2, 23 ) = 0.000000;
                aXihat( 0, 24 ) = 1.000000;
                aXihat( 1, 24 ) = 0.000000;
                aXihat( 2, 24 ) = 0.000000;
                aXihat( 0, 25 ) = 0.000000;
                aXihat( 1, 25 ) = -1.000000;
                aXihat( 2, 25 ) = 0.000000;
                aXihat( 0, 26 ) = 0.000000;
                aXihat( 1, 26 ) = 1.000000;
                aXihat( 2, 26 ) = 0.000000;
                break;
            }
            case 3:
            {
                // hex 64
                real c = 1.0 / 3.0;

                aXihat( 0, 0 )  = -1.0;
                aXihat( 1, 0 )  = -1.0;
                aXihat( 2, 0 )  = -1.0;
                aXihat( 0, 1 )  = 1.0;
                aXihat( 1, 1 )  = -1.0;
                aXihat( 2, 1 )  = -1.0;
                aXihat( 0, 2 )  = 1.0;
                aXihat( 1, 2 )  = 1.0;
                aXihat( 2, 2 )  = -1.0;
                aXihat( 0, 3 )  = -1.0;
                aXihat( 1, 3 )  = 1.0;
                aXihat( 2, 3 )  = -1.0;
                aXihat( 0, 4 )  = -1.0;
                aXihat( 1, 4 )  = -1.0;
                aXihat( 2, 4 )  = 1.0;
                aXihat( 0, 5 )  = 1.0;
                aXihat( 1, 5 )  = -1.0;
                aXihat( 2, 5 )  = 1.0;
                aXihat( 0, 6 )  = 1.0;
                aXihat( 1, 6 )  = 1.0;
                aXihat( 2, 6 )  = 1.0;
                aXihat( 0, 7 )  = -1.0;
                aXihat( 1, 7 )  = 1.0;
                aXihat( 2, 7 )  = 1.0;
                aXihat( 0, 8 )  = -c;
                aXihat( 1, 8 )  = -1.0;
                aXihat( 2, 8 )  = -1.0;
                aXihat( 0, 9 )  = c;
                aXihat( 1, 9 )  = -1.0;
                aXihat( 2, 9 )  = -1.0;
                aXihat( 0, 10 ) = -1.0;
                aXihat( 1, 10 ) = -c;
                aXihat( 2, 10 ) = -1.0;
                aXihat( 0, 11 ) = -1.0;
                aXihat( 1, 11 ) = c;
                aXihat( 2, 11 ) = -1.0;
                aXihat( 0, 12 ) = -1.0;
                aXihat( 1, 12 ) = -1.0;
                aXihat( 2, 12 ) = -c;
                aXihat( 0, 13 ) = -1.0;
                aXihat( 1, 13 ) = -1.0;
                aXihat( 2, 13 ) = c;
                aXihat( 0, 14 ) = 1.0;
                aXihat( 1, 14 ) = -c;
                aXihat( 2, 14 ) = -1.0;
                aXihat( 0, 15 ) = 1.0;
                aXihat( 1, 15 ) = c;
                aXihat( 2, 15 ) = -1.0;
                aXihat( 0, 16 ) = 1.0;
                aXihat( 1, 16 ) = -1.0;
                aXihat( 2, 16 ) = -c;
                aXihat( 0, 17 ) = 1.0;
                aXihat( 1, 17 ) = -1.0;
                aXihat( 2, 17 ) = c;
                aXihat( 0, 18 ) = c;
                aXihat( 1, 18 ) = 1.0;
                aXihat( 2, 18 ) = -1.0;
                aXihat( 0, 19 ) = -c;
                aXihat( 1, 19 ) = 1.0;
                aXihat( 2, 19 ) = -1.0;
                aXihat( 0, 20 ) = 1.0;
                aXihat( 1, 20 ) = 1.0;
                aXihat( 2, 20 ) = -c;
                aXihat( 0, 21 ) = 1.0;
                aXihat( 1, 21 ) = 1.0;
                aXihat( 2, 21 ) = c;
                aXihat( 0, 22 ) = -1.0;
                aXihat( 1, 22 ) = 1.0;
                aXihat( 2, 22 ) = -c;
                aXihat( 0, 23 ) = -1.0;
                aXihat( 1, 23 ) = 1.0;
                aXihat( 2, 23 ) = c;
                aXihat( 0, 24 ) = -c;
                aXihat( 1, 24 ) = -1.0;
                aXihat( 2, 24 ) = 1.0;
                aXihat( 0, 25 ) = c;
                aXihat( 1, 25 ) = -1.0;
                aXihat( 2, 25 ) = 1.0;
                aXihat( 0, 26 ) = -1.0;
                aXihat( 1, 26 ) = -c;
                aXihat( 2, 26 ) = 1.0;
                aXihat( 0, 27 ) = -1.0;
                aXihat( 1, 27 ) = c;
                aXihat( 2, 27 ) = 1.0;
                aXihat( 0, 28 ) = 1.0;
                aXihat( 1, 28 ) = -c;
                aXihat( 2, 28 ) = 1.0;
                aXihat( 0, 29 ) = 1.0;
                aXihat( 1, 29 ) = c;
                aXihat( 2, 29 ) = 1.0;
                aXihat( 0, 30 ) = c;
                aXihat( 1, 30 ) = 1.0;
                aXihat( 2, 30 ) = 1.0;
                aXihat( 0, 31 ) = -c;
                aXihat( 1, 31 ) = 1.0;
                aXihat( 2, 31 ) = 1.0;
                aXihat( 0, 32 ) = -c;
                aXihat( 1, 32 ) = -c;
                aXihat( 2, 32 ) = -1.0;
                aXihat( 0, 33 ) = -c;
                aXihat( 1, 33 ) = c;
                aXihat( 2, 33 ) = -1.0;
                aXihat( 0, 34 ) = c;
                aXihat( 1, 34 ) = c;
                aXihat( 2, 34 ) = -1.0;
                aXihat( 0, 35 ) = c;
                aXihat( 1, 35 ) = -c;
                aXihat( 2, 35 ) = -1.0;
                aXihat( 0, 36 ) = -c;
                aXihat( 1, 36 ) = -1.0;
                aXihat( 2, 36 ) = -c;
                aXihat( 0, 37 ) = c;
                aXihat( 1, 37 ) = -1.0;
                aXihat( 2, 37 ) = -c;
                aXihat( 0, 38 ) = c;
                aXihat( 1, 38 ) = -1.0;
                aXihat( 2, 38 ) = c;
                aXihat( 0, 39 ) = -c;
                aXihat( 1, 39 ) = -1.0;
                aXihat( 2, 39 ) = c;
                aXihat( 0, 40 ) = -1.0;
                aXihat( 1, 40 ) = -c;
                aXihat( 2, 40 ) = -c;
                aXihat( 0, 41 ) = -1.0;
                aXihat( 1, 41 ) = -c;
                aXihat( 2, 41 ) = c;
                aXihat( 0, 42 ) = -1.0;
                aXihat( 1, 42 ) = c;
                aXihat( 2, 42 ) = c;
                aXihat( 0, 43 ) = -1.0;
                aXihat( 1, 43 ) = c;
                aXihat( 2, 43 ) = -c;
                aXihat( 0, 44 ) = 1.0;
                aXihat( 1, 44 ) = -c;
                aXihat( 2, 44 ) = -c;
                aXihat( 0, 45 ) = 1.0;
                aXihat( 1, 45 ) = c;
                aXihat( 2, 45 ) = -c;
                aXihat( 0, 46 ) = 1.0;
                aXihat( 1, 46 ) = c;
                aXihat( 2, 46 ) = c;
                aXihat( 0, 47 ) = 1.0;
                aXihat( 1, 47 ) = -c;
                aXihat( 2, 47 ) = c;
                aXihat( 0, 48 ) = c;
                aXihat( 1, 48 ) = 1.0;
                aXihat( 2, 48 ) = -c;
                aXihat( 0, 49 ) = -c;
                aXihat( 1, 49 ) = 1.0;
                aXihat( 2, 49 ) = -c;
                aXihat( 0, 50 ) = -c;
                aXihat( 1, 50 ) = 1.0;
                aXihat( 2, 50 ) = c;
                aXihat( 0, 51 ) = c;
                aXihat( 1, 51 ) = 1.0;
                aXihat( 2, 51 ) = c;
                aXihat( 0, 52 ) = -c;
                aXihat( 1, 52 ) = -c;
                aXihat( 2, 52 ) = 1.0;
                aXihat( 0, 53 ) = c;
                aXihat( 1, 53 ) = -c;
                aXihat( 2, 53 ) = 1.0;
                aXihat( 0, 54 ) = c;
                aXihat( 1, 54 ) = c;
                aXihat( 2, 54 ) = 1.0;
                aXihat( 0, 55 ) = -c;
                aXihat( 1, 55 ) = c;
                aXihat( 2, 55 ) = 1.0;
                aXihat( 0, 56 ) = -c;
                aXihat( 1, 56 ) = -c;
                aXihat( 2, 56 ) = -c;
                aXihat( 0, 57 ) = c;
                aXihat( 1, 57 ) = -c;
                aXihat( 2, 57 ) = -c;
                aXihat( 0, 58 ) = c;
                aXihat( 1, 58 ) = c;
                aXihat( 2, 58 ) = -c;
                aXihat( 0, 59 ) = -c;
                aXihat( 1, 59 ) = c;
                aXihat( 2, 59 ) = -c;
                aXihat( 0, 60 ) = -c;
                aXihat( 1, 60 ) = -c;
                aXihat( 2, 60 ) = c;
                aXihat( 0, 61 ) = c;
                aXihat( 1, 61 ) = -c;
                aXihat( 2, 61 ) = c;
                aXihat( 0, 62 ) = c;
                aXihat( 1, 62 ) = c;
                aXihat( 2, 62 ) = c;
                aXihat( 0, 63 ) = -c;
                aXihat( 1, 63 ) = c;
                aXihat( 2, 63 ) = c;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "something went wrong while creating T-Matrices." );
                break;
            }
        }

        return aXihat;
    }

    // -----------------------------------------------------------------------------------------------------------------

    template<>
    void T_Matrix< 3 >::populate_child_matrices( const Matrix< DDRMat >& aTL, const Matrix< DDRMat >& aTR )
    {
        // Number of coefficients
        uint tNumCoefficients = aTL.n_rows();

        // Tensor product
        uint tChildCol = 0;
        for ( uint iTRow3 = 0; iTRow3 < tNumCoefficients; iTRow3++ )
        {
            for ( uint iTRow2 = 0; iTRow2 < tNumCoefficients; iTRow2++ )
            {
                for ( uint iTRow1 = 0; iTRow1 < tNumCoefficients; iTRow1++ )
                {
                    uint tChildRow = 0;
                    for ( uint iTCol3 = 0; iTCol3 < tNumCoefficients; iTCol3++ )
                    {
                        for ( uint iTCol2 = 0; iTCol2 < tNumCoefficients; iTCol2++ )
                        {
                            for ( uint iTCol1 = 0; iTCol1 < tNumCoefficients; iTCol1++ )
                            {
                                mChild( 0 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTL( iTRow1, iTCol1 ) * aTL( iTRow2, iTCol2 ) * aTL( iTRow3, iTCol3 );
                                mChild( 1 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTR( iTRow1, iTCol1 ) * aTL( iTRow2, iTCol2 ) * aTL( iTRow3, iTCol3 );
                                mChild( 2 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTL( iTRow1, iTCol1 ) * aTR( iTRow2, iTCol2 ) * aTL( iTRow3, iTCol3 );
                                mChild( 3 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTR( iTRow1, iTCol1 ) * aTR( iTRow2, iTCol2 ) * aTL( iTRow3, iTCol3 );
                                mChild( 4 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTL( iTRow1, iTCol1 ) * aTL( iTRow2, iTCol2 ) * aTR( iTRow3, iTCol3 );
                                mChild( 5 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTR( iTRow1, iTCol1 ) * aTL( iTRow2, iTCol2 ) * aTR( iTRow3, iTCol3 );
                                mChild( 6 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTL( iTRow1, iTCol1 ) * aTR( iTRow2, iTCol2 ) * aTR( iTRow3, iTCol3 );
                                mChild( 7 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTR( iTRow1, iTCol1 ) * aTR( iTRow2, iTCol2 ) * aTR( iTRow3, iTCol3 );
                                tChildRow++;
                            }
                        }
                    }
                    tChildCol++;
                }
            }
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

}
