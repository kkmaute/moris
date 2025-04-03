/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Mesh_1D.hpp
 *
 */

#pragma once

#include "cl_HMR_Background_Mesh.hpp"    //HMR/src

namespace moris::hmr
{
    //-------------------------------------------------------------------------------

    template<>
    void
    Background_Mesh< 1 >::create_coarsest_frame()
    {
        // calculate number of elements in frame
        luint tNumberOfElements = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ]
                                - mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + 1;

        // clear coarsest element list
        mCoarsestElements.clear();

        // assign memory
        mCoarsestElements.resize( tNumberOfElements, nullptr );

        // initialize counter
        luint tCount = 0;

        // loop over domain
        for ( luint i = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
                i <= mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
                ++i )
        {
            mCoarsestElements( tCount++ ) =
                    mCoarsestElementsIncludingAura( calc_subdomain_id_of_element( 0, i ) );
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline luint
    Background_Mesh< 1 >::calc_domain_id_of_element(
            uint  aLevel,
            luint aI ) const
    {
        if ( aLevel < gMaxNumberOfLevels )
        {
            luint tI = aI + mMySubDomain.mAuraIJK[ aLevel ][ 0 ][ 0 ];

            // test if input is valid
            if ( tI >= mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ] )
            {
                // return error value
                return gNoEntityID;
            }
            else
            {
                // calculate domain id
                return mDomain.mLevelOffset[ aLevel ] + tI;
            }
        }
        else
        {
            // return error value
            return gNoEntityID;
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline luint
    Background_Mesh< 1 >::calc_domain_id_of_element(
            uint  aLevel,
            luint aI ) const
    {
        if ( aLevel < gMaxNumberOfLevels )
        {
            luint tI = aI + mMySubDomain.mAuraIJK[ aLevel ][ 0 ][ 0 ];

            // test if input is valid
            if ( tI >= mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ] )
            {
                // return error value
                return gNoEntityID;
            }
            else
            {
                // calculate domain id
                return mDomain.mLevelOffset[ aLevel ] + tI;
            }
        }
        else
        {
            // return error value
            return gNoEntityID;
        }
    }

    //-------------------------------------------------------------------------------
    template<>
    inline luint
    Background_Mesh< 1 >::calc_subdomain_id_of_element(
            uint  aLevel,
            luint aI ) const
    {
        if ( aLevel < gMaxNumberOfLevels )
        {
            // test if input is valid
            if ( aI >= mMySubDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ] )
            {
                // return no value
                return gNoEntityID;
            }
            else
            {
                // calculate element ID
                return mMySubDomain.mLevelOffset[ aLevel ] + aI;
            }
        }
        else
        {
            // return no value
            return gNoEntityID;
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline luint
    Background_Mesh< 1 >::calc_subdomain_id_from_global_id(
            uint  aLevel,
            luint aID ) const
    {
        if ( aLevel < gMaxNumberOfLevels )
        {
            // calculate global i position
            luint tI = aID - mDomain.mLevelOffset[ aLevel ];

            // calculate local i
            tI -= mMySubDomain.mAuraIJK[ aLevel ][ 0 ][ 0 ];

            // calculate local ID
            return this->calc_subdomain_id_of_element( aLevel, tI );
        }
        else
        {
            // return no value
            return gNoEntityID;
        }
    }

    //--------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 1 >::calc_ijk_from_global_id(
            const uint&  aLevel,
            const luint& aID,
            luint*       aIJK ) const
    {
        if ( aLevel < gMaxNumberOfLevels )
        {
            // calculate global i position
            luint tI = aID - mDomain.mLevelOffset[ aLevel ];
        }
        else
        {
            // return no value
            return gNoEntityID;
        }
    }

    //-------------------------------------------------------------------------------
}    // namespace moris::hmr
