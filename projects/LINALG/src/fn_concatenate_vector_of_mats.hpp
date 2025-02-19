/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_concatenate_vector_of_mats.hpp
 *
 */

#pragma once

// MORIS library header files.
#include "cl_Matrix.hpp"

namespace moris
{
    template< typename MatrixType >
    Matrix< MatrixType >
    concatenate_vector_of_mats( Vector< Matrix< MatrixType > > aMat,
            moris_index                                        aFixedDim )
    {
        moris_index tFixedDimSize = 0;
        moris_index tTotalSize    = 0;

        bool tFixedDimSizeSet = false;

        Vector< moris_index > tMatSize( 2 );

        MORIS_ASSERT( aFixedDim <= 2, "Fixed dim can be 0 for row or 1 for col." );

        for ( moris::uint i = 0; i < aMat.size(); i++ )
        {

            tMatSize( 0 ) = aMat( i ).n_rows();
            tMatSize( 1 ) = aMat( i ).n_cols();

            if ( tMatSize( aFixedDim ) == 0 )
            {
                continue;
            }

            if ( not tFixedDimSizeSet )
            {
                tFixedDimSize    = tMatSize( aFixedDim );
                tFixedDimSizeSet = true;
            }
            else
            {
                MORIS_ASSERT( tMatSize( aFixedDim ) == tFixedDimSize, "fixed dimension not consistent." );
            }

            tTotalSize = tTotalSize + aMat( i ).numel();
        }

        if ( tTotalSize == 0 )
        {
            return Matrix< MatrixType >();
        }

        // size the concatenated matrix
        Matrix< MatrixType > tConcatenatedMat;

        moris_index tStart = 0;
        moris_index tEnd   = 0;

        if ( aFixedDim == 0 )
        {
            tConcatenatedMat.resize( tFixedDimSize, tTotalSize / tFixedDimSize );

            for ( moris::uint i = 0; i < aMat.size(); i++ )
            {
                if ( aMat( i ).n_rows() == 0 )
                {
                    continue;
                }

                tEnd = tStart + aMat( i ).n_cols() - 1;

                tConcatenatedMat( { 0, tFixedDimSize - 1 }, { tStart, tEnd } ) = aMat( i ).matrix_data();
                tStart                                                         = tEnd + 1;
            }
        }
        else if ( aFixedDim == 1 )
        {
            tConcatenatedMat.resize( tTotalSize / tFixedDimSize, tFixedDimSize );
            for ( moris::uint i = 0; i < aMat.size(); i++ )
            {
                if ( aMat( i ).n_cols() == 0 )
                {
                    continue;
                }

                tEnd = tStart + aMat( i ).n_rows() - 1;

                tConcatenatedMat( { tStart, tEnd }, { 0, tFixedDimSize - 1 } ) = aMat( i ).matrix_data();
                tStart                                                         = tEnd + 1;
            }
        }

        return tConcatenatedMat;
    }
}    // namespace moris
