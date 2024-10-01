/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Domain.hpp
 *
 */

#pragma once

#include "HMR_Globals.hpp" //HMR/src
#include "moris_typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src

namespace moris::hmr
{
    /**
     * the Domain struct contains basic information needed for
     * mesh generation and calculation of IDs.
     */
    template < uint N >
    struct Domain
    {
        // span of calculation domain, global ijk
        luint mDomainIJK[ gMaxNumberOfLevels ][ N ][ 2 ];

        // span of calculation domain plus aura, global ijk
        luint mAuraIJK[ gMaxNumberOfLevels ][ N ][ 2 ];

        // this information is needed for ID calculation
        luint mNumberOfElementsPerDimension[ gMaxNumberOfLevels ][ N ];

        luint mNumberOfElementsOnLevelZero;

        // this information is needed for ID calculation
        luint mLevelOffset[ gMaxNumberOfLevels ];

        // frame of elements within domain, local ijk
        luint mFrameIJK[ gMaxNumberOfLevels ][ N ][ 2 ];

//--------------------------------------------------------------------------------
    public:
//--------------------------------------------------------------------------------
        // empty constructor
        Domain(){};
//--------------------------------------------------------------------------------
        Domain( const Matrix< DDLUMat > & aDomainIJK,
                const luint               aPaddingSize )
        {
            mNumberOfElementsOnLevelZero = 1;

            // loop over all dimensions
            for( uint k=0; k< N; ++k )
            {
                // copy data into DomainIJK
                mDomainIJK[ 0 ][ k ][ 0 ] = aDomainIJK( 0, k );
                mDomainIJK[ 0 ][ k ][ 1 ] = aDomainIJK( 1, k );

                // add padding for Aura
                mAuraIJK[ 0 ][ k ][ 0 ]   = aDomainIJK( 0, k ) - aPaddingSize;
                mAuraIJK[ 0 ][ k ][ 1 ]   = aDomainIJK( 1, k ) + aPaddingSize;

                // element per dimension counter for first level
                mNumberOfElementsPerDimension[ 0 ][ k ] = mAuraIJK[ 0 ][ k ][ 1 ] - mAuraIJK[ 0 ][ k ][ 0 ] + 1;

                // element counter for fist level
                mNumberOfElementsOnLevelZero *=  mNumberOfElementsPerDimension[ 0 ][ k ] ;

                // create frame
                mFrameIJK[ 0 ][ k ][ 0 ] = aPaddingSize;
                mFrameIJK[ 0 ][ k ][ 1 ] = mNumberOfElementsPerDimension[ 0 ][ k ] - aPaddingSize - 1;
            }

            // create element table for higher levels
            for ( luint l=1; l<gMaxNumberOfLevels; ++l )
            {
                for ( luint k=0; k<N; ++k )
                {
                    mNumberOfElementsPerDimension[ l ][ k ] = 2 * mNumberOfElementsPerDimension[ l-1 ][ k ];

                    mDomainIJK[ l ][ k ][ 0 ] = 2*mDomainIJK[ l-1 ][ k ][ 0 ];
                    mDomainIJK[ l ][ k ][ 1 ] = 2*mDomainIJK[ l-1 ][ k ][ 1 ] + 1;

                    mAuraIJK[ l ][ k ][ 0 ] = 2*mAuraIJK[ l-1 ][ k ][ 0 ];
                    mAuraIJK[ l ][ k ][ 1 ] = 2*mAuraIJK[ l-1 ][ k ][ 1 ] + 1;

                    mFrameIJK[ l ][ k ][ 0 ] = 2*mFrameIJK[ l-1 ][ k ][ 0 ];
                    mFrameIJK[ l ][ k ][ 1 ] = 2*mFrameIJK[ l-1 ][ k ][ 1 ] + 1;
                }
            }

            // element offset counter
            auto tCount = mNumberOfElementsOnLevelZero;

            // children per element, needer for offset increase
            luint tNumberOfChildrenPerElement = std::pow( 2, N );

            // first entry of offset
            mLevelOffset[ 0 ] = 0;

            // all higher levels
            for ( luint l=1; l<gMaxNumberOfLevels; ++l )
            {
                // add counter to offset
                mLevelOffset[ l ] = mLevelOffset[ l-1 ] + tCount;

                // multiply counter
                tCount *= tNumberOfChildrenPerElement;
            }
        }
    };
//--------------------------------------------------------------------------------

}
