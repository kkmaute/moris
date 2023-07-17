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

    template<>
    Matrix< DDRMat > T_Matrix< 2 >::get_supporting_points( const uint aOrder )
    {
        // the following lines are needed to get the interpolation points
        Matrix< DDRMat > aXihat( 2, std::pow( aOrder + 1, 2 ) );

        switch ( aOrder )
        {
            case 1:
            {
                // quad 4
                aXihat( 0, 0 ) = -1.000000;
                aXihat( 1, 0 ) = -1.000000;
                aXihat( 0, 1 ) = 1.000000;
                aXihat( 1, 1 ) = -1.000000;
                aXihat( 0, 2 ) = 1.000000;
                aXihat( 1, 2 ) = 1.000000;
                aXihat( 0, 3 ) = -1.000000;
                aXihat( 1, 3 ) = 1.000000;
                break;
            }
            case 2:
            {
                // quad 9
                aXihat( 0, 0 ) = -1.000000;
                aXihat( 1, 0 ) = -1.000000;
                aXihat( 0, 1 ) = 1.000000;
                aXihat( 1, 1 ) = -1.000000;
                aXihat( 0, 2 ) = 1.000000;
                aXihat( 1, 2 ) = 1.000000;
                aXihat( 0, 3 ) = -1.000000;
                aXihat( 1, 3 ) = 1.000000;
                aXihat( 0, 4 ) = 0.000000;
                aXihat( 1, 4 ) = -1.000000;
                aXihat( 0, 5 ) = 1.000000;
                aXihat( 1, 5 ) = 0.000000;
                aXihat( 0, 6 ) = 0.000000;
                aXihat( 1, 6 ) = 1.000000;
                aXihat( 0, 7 ) = -1.000000;
                aXihat( 1, 7 ) = 0.000000;
                aXihat( 0, 8 ) = 0.000000;
                aXihat( 1, 8 ) = 0.000000;
                break;
            }
            case 3:
            {
                // quad 16
                real c = 1.0 / 3.0;

                aXihat( 0, 0 )  = -1.000000;
                aXihat( 1, 0 )  = -1.000000;
                aXihat( 0, 1 )  = 1.000000;
                aXihat( 1, 1 )  = -1.000000;
                aXihat( 0, 2 )  = 1.000000;
                aXihat( 1, 2 )  = 1.000000;
                aXihat( 0, 3 )  = -1.000000;
                aXihat( 1, 3 )  = 1.000000;
                aXihat( 0, 4 )  = -c;
                aXihat( 1, 4 )  = -1.000000;
                aXihat( 0, 5 )  = c;
                aXihat( 1, 5 )  = -1.000000;
                aXihat( 0, 6 )  = 1.000000;
                aXihat( 1, 6 )  = -c;
                aXihat( 0, 7 )  = 1.000000;
                aXihat( 1, 7 )  = c;
                aXihat( 0, 8 )  = c;
                aXihat( 1, 8 )  = 1.000000;
                aXihat( 0, 9 )  = -c;
                aXihat( 1, 9 )  = 1.000000;
                aXihat( 0, 10 ) = -1.000000;
                aXihat( 1, 10 ) = c;
                aXihat( 0, 11 ) = -1.000000;
                aXihat( 1, 11 ) = -c;
                aXihat( 0, 12 ) = -c;
                aXihat( 1, 12 ) = -c;
                aXihat( 0, 13 ) = c;
                aXihat( 1, 13 ) = -c;
                aXihat( 0, 14 ) = c;
                aXihat( 1, 14 ) = c;
                aXihat( 0, 15 ) = -c;
                aXihat( 1, 15 ) = c;
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
    void T_Matrix< 2 >::populate_child_matrices( const Matrix< DDRMat >& aTL, const Matrix< DDRMat >& aTR )
    {
        // Number of coefficients
        uint tNumCoefficients = aTL.n_rows();

        // Tensor product
        uint tChildCol = 0;
        for ( uint iTRow2 = 0; iTRow2 < tNumCoefficients; iTRow2++ )
        {
            for ( uint iTRow1 = 0; iTRow1 < tNumCoefficients; iTRow1++ )
            {
                uint tChildRow = 0;
                for ( uint iTCol2 = 0; iTCol2 < tNumCoefficients; iTCol2++ )
                {
                    for ( uint iTCol1 = 0; iTCol1 < tNumCoefficients; iTCol1++ )
                    {
                        mChild( 0 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTL( iTRow1, iTCol1 ) * aTL( iTRow2, iTCol2 );
                        mChild( 1 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTR( iTRow1, iTCol1 ) * aTL( iTRow2, iTCol2 );
                        mChild( 2 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTL( iTRow1, iTCol1 ) * aTR( iTRow2, iTCol2 );
                        mChild( 3 )( mBasisIndex( tChildRow ), mBasisIndex( tChildCol ) ) = aTR( iTRow1, iTCol1 ) * aTR( iTRow2, iTCol2 );
                        tChildRow++;
                    }
                }
                tChildCol++;
            }
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

}
