/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator_Test_Polynomial.cpp
 *
 */

#include "fn_load_matrix_from_binary_file.hpp"      //LNA/src
#include "cl_MTK_Integrator_Test_Polynomial.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Integrator_Test_Polynomial::Integrator_Test_Polynomial(
                const std::string &aPath )
        {
            // load matrix from file ( use double, otherwise test would break in 64bit mode)

            Matrix< DDRMat > tData;
            load_matrix_from_binary_file( tData, aPath );

            // initialize counter
            uint tCount = 0;

            // read number of dimensions
            mNumberOfDimensions = tData( tCount++ );

            // read number of coeffs
            mNumberOfCoeffs = tData( tCount++ );

            // read solution
            mIntegral = tData( tCount++ );

            // initialize coefficient matrix
            mCoeffs.set_size( mNumberOfCoeffs, 1 );

            // read coefficients
            for ( uint k = 0; k < mNumberOfCoeffs; ++k )
            {
                mCoeffs( k ) = tData( tCount++ );
            }

            // initialize exponents
            mExponents.set_size( mNumberOfDimensions, mNumberOfCoeffs );

            // read exponents
            for ( uint k = 0; k < mNumberOfCoeffs; ++k )
            {
                for ( uint i = 0; i < mNumberOfDimensions; ++i )
                {
                    mExponents( i, k ) = tData( tCount++ );
                }
            }
        }

        //------------------------------------------------------------------------------

        real
        Integrator_Test_Polynomial::eval( const Matrix< DDRMat > &aXi )
        {
            real aValue = 0;
            for ( uint k = 0; k < mNumberOfCoeffs; ++k )
            {
                real tValue = mCoeffs( k );
                for ( uint i = 0; i < mNumberOfDimensions; ++i )
                {
                    tValue *= std::pow( aXi( i ), mExponents( i, k ) );
                }
                aValue += tValue;
            }
            return aValue;
        }

        //------------------------------------------------------------------------------

        real
        Integrator_Test_Polynomial::get_integral()
        {
            return mIntegral;
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
