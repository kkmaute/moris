/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Mesh_3D.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BACKGROUND_MESH_3D_HPP_
#define SRC_HMR_CL_HMR_BACKGROUND_MESH_3D_HPP_

#include "cl_HMR_Background_Mesh.hpp"    //HMR/src

namespace moris::hmr
{
    //-------------------------------------------------------------------------------

    template<>
    inline luint
    Background_Mesh< 3 >::calc_domain_id_of_element(
            uint  aLevel,
            luint aI,
            luint aJ,
            luint aK ) const
    {
        MORIS_ASSERT( aLevel < gMaxNumberOfLevels, "calc_domain_id_of_element(), Requested refinement level larger than maximal refinement level" );

        luint tI = aI + mMySubDomain.mAuraIJK[ aLevel ][ 0 ][ 0 ];
        luint tJ = aJ + mMySubDomain.mAuraIJK[ aLevel ][ 1 ][ 0 ];
        luint tK = aK + mMySubDomain.mAuraIJK[ aLevel ][ 2 ][ 0 ];

        MORIS_ASSERT( ( tI < mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ] )
                              && ( tJ < mDomain.mNumberOfElementsPerDimension[ aLevel ][ 1 ] )
                              && ( tK < mDomain.mNumberOfElementsPerDimension[ aLevel ][ 2 ] ),
                "calc_domain_id_of_element(), I, J or K position of this element outside of the domain" );

        // calculate domain id
        return mDomain.mLevelOffset[ aLevel ] + tI
             + mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ]
                       * ( tJ + tK * mDomain.mNumberOfElementsPerDimension[ aLevel ][ 1 ] );
    }

    //-------------------------------------------------------------------------------

    template<>
    inline luint
    Background_Mesh< 3 >::calc_subdomain_id_of_element(
            uint  aLevel,
            luint aI,
            luint aJ,
            luint aK ) const
    {
        MORIS_ASSERT( aLevel < gMaxNumberOfLevels, "calc_subdomain_id_of_element(), Requested refinement level larger than maximal refinement level" );

        if ( ( aI >= mMySubDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ] )
                || ( aJ >= mMySubDomain.mNumberOfElementsPerDimension[ aLevel ][ 1 ] )
                || ( aK >= mMySubDomain.mNumberOfElementsPerDimension[ aLevel ][ 2 ] ) )
        {
            // return no value
            return gNoEntityID;
        }
        else
        {
            // calculate element ID
            return mMySubDomain.mLevelOffset[ aLevel ] + aI
                 + mMySubDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ] *    //
                           ( aJ + aK * mMySubDomain.mNumberOfElementsPerDimension[ aLevel ][ 1 ] );
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline luint
    Background_Mesh< 3 >::calc_subdomain_id_from_global_id(
            uint  aLevel,
            luint aID ) const
    {
        MORIS_ASSERT( aLevel < gMaxNumberOfLevels,
                "calc_subdomain_id_from_global_id(), Requested refinement level larger than maximal refinement level" );

        // subtract level offset from ID
        luint tID = aID - mDomain.mLevelOffset[ aLevel ];

        // help variable
        luint tPlane = mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ]
                     * mDomain.mNumberOfElementsPerDimension[ aLevel ][ 1 ];

        // calculate global ijk position
        luint tK = tID / tPlane;
        tID -= tK * tPlane;

        luint tJ = tID / mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ];
        luint tI = tID - tJ * mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ];

        // calculate local ijk position
        tI -= mMySubDomain.mAuraIJK[ aLevel ][ 0 ][ 0 ];
        tJ -= mMySubDomain.mAuraIJK[ aLevel ][ 1 ][ 0 ];
        tK -= mMySubDomain.mAuraIJK[ aLevel ][ 2 ][ 0 ];

        return this->calc_subdomain_id_of_element( aLevel, tI, tJ, tK );
    }

    //--------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::calc_ijk_from_global_id(
            const uint&  aLevel,
            const luint& aID,
            luint*       aIJK ) const
    {
        MORIS_ASSERT( aLevel < gMaxNumberOfLevels,
                "calc_ijk_from_global_id(), Requested refinement level larger than maximal refinement level" );

        // subtract level offset from ID
        luint tID = aID - mDomain.mLevelOffset[ aLevel ];

        // help variable
        luint tPlane = mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ]
                     * mDomain.mNumberOfElementsPerDimension[ aLevel ][ 1 ];

        // calculate global ijk position
        aIJK[ 2 ] = tID / tPlane;
        tID -= aIJK[ 2 ] * tPlane;

        aIJK[ 1 ] = tID / mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ];
        aIJK[ 0 ] = tID - aIJK[ 1 ] * mDomain.mNumberOfElementsPerDimension[ aLevel ][ 0 ];
    }

    //--------------------------------------------------------------------------------


    template<>
    inline void
    Background_Mesh< 3 >::calc_element_ids(
            uint                     aLevel,
            const Matrix< DDLUMat >& aIJK,
            Matrix< DDLUMat >&       aIDs ) const
    {
        // reserve memory for output
        aIDs.set_size( 8, 1 );

        // child 0
        aIDs( 0 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 0 ),
                aIJK( 1, 0 ),
                aIJK( 2, 0 ) );

        // child 1
        aIDs( 1 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 1 ),
                aIJK( 1, 1 ),
                aIJK( 2, 1 ) );

        // child 2
        aIDs( 2 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 2 ),
                aIJK( 1, 2 ),
                aIJK( 2, 2 ) );

        // child 3
        aIDs( 3 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 3 ),
                aIJK( 1, 3 ),
                aIJK( 2, 3 ) );

        // child 4
        aIDs( 4 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 4 ),
                aIJK( 1, 4 ),
                aIJK( 2, 4 ) );

        // child 5
        aIDs( 5 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 5 ),
                aIJK( 1, 5 ),
                aIJK( 2, 5 ) );

        // child 6
        aIDs( 6 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 6 ),
                aIJK( 1, 6 ),
                aIJK( 2, 6 ) );

        // child 7
        aIDs( 7 ) = calc_domain_id_of_element(
                aLevel,
                aIJK( 0, 7 ),
                aIJK( 1, 7 ),
                aIJK( 2, 7 ) );
    }

    //--------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::create_coarsest_frame()
    {
        // calculate number of elements in frame
        luint tNumberOfElements = ( mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ]
                                          - mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + 1 )
                                * ( mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ]
                                        - mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + 1 )
                                * ( mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ]
                                        - mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + 1 );

        // clear coarsest element list
        mCoarsestElements.clear();

        // assign memory
        mCoarsestElements.resize( tNumberOfElements, nullptr );

        // initialize counter
        luint tCount = 0;

        // loop over domain
        for ( luint k = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
                k <= mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];
                ++k )
        {
            for ( luint j = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
                    j <= mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
                    ++j )
            {
                for ( luint i = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
                        i <= mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
                        ++i )
                {
                    // add element from mCoarsestElementsIncludingAura
                    mCoarsestElements( tCount++ ) =
                            mCoarsestElementsIncludingAura(
                                    calc_subdomain_id_of_element( 0, i, j, k ) );
                }
            }
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::initialize_coarsest_elements()
    {
        // assign memory for coarsest elements
        mCoarsestElementsIncludingAura.resize( mMySubDomain.mNumberOfElementsOnLevelZero,
                nullptr );

        // calculate number of elements on level zero pre direction
        Matrix< DDLUMat > tNumberOfElements = get_number_of_subdomain_elements_per_direction_on_level_zero();

        luint* tIJK = new luint[ 3 ];

        // initialize counter
        luint tCount = 0;

        // add elements on level zero
        for ( luint k = 0; k < tNumberOfElements( 2 ); ++k )
        {
            tIJK[ 2 ] = k;
            for ( luint j = 0; j < tNumberOfElements( 1 ); ++j )
            {
                tIJK[ 1 ] = j;
                for ( luint i = 0; i < tNumberOfElements( 0 ); ++i )
                {
                    tIJK[ 0 ] = i;
                    this->insert_zero_level_element( tCount++,
                            new Background_Element< 3 >( (Background_Element_Base*)nullptr,
                                    mActivePattern,
                                    tIJK,
                                    this->calc_domain_id_of_element( 0, i, j, k ),
                                    (uint)0,
                                    (uint)0,
                                    (uint)gNoProcOwner ) );
                }
            }
        }
        delete[] tIJK;
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::finalize_coarsest_elements()
    {
        // set boundaries for loop
        luint tImin[ 27 ];
        luint tImax[ 27 ];
        luint tJmin[ 27 ];
        luint tJmax[ 27 ];
        luint tKmin[ 27 ];
        luint tKmax[ 27 ];

        // create quadrants for aura elements
        // quadrant 0
        tImin[ 0 ] = 0;
        tImax[ 0 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 0 ] = 0;
        tJmax[ 0 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 0 ] = 0;
        tKmax[ 0 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 1
        tImin[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 1 ] = 0;
        tJmax[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 1 ] = 0;
        tKmax[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 2
        tImin[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 2 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 2 ] = 0;
        tJmax[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 2 ] = 0;
        tKmax[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 3
        tImin[ 3 ] = 0;
        tImax[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 3 ] = 0;
        tKmax[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 4
        tImin[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 4 ] = 0;
        tKmax[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 5
        tImin[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 5 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 5 ] = 0;
        tKmax[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 6
        tImin[ 6 ] = 0;
        tImax[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 6 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 6 ] = 0;
        tKmax[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 7
        tImin[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 7 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 7 ] = 0;
        tKmax[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 8
        tImin[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 8 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 8 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 8 ] = 0;
        tKmax[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] - 1;

        // quadrant 9
        tImin[ 9 ] = 0;
        tImax[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 9 ] = 0;
        tJmax[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 10
        tImin[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 10 ] = 0;
        tJmax[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 11
        tImin[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 11 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 11 ] = 0;
        tJmax[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 12
        tImin[ 12 ] = 0;
        tImax[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 13
        tImin[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 14
        tImin[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 14 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 15
        tImin[ 15 ] = 0;
        tImax[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 15 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 16
        tImin[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 16 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 17
        tImin[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 17 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 17 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 18
        tImin[ 18 ] = 0;
        tImax[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 18 ] = 0;
        tJmax[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 18 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 19
        tImin[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 19 ] = 0;
        tJmax[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 19 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 20
        tImin[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 20 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 20 ] = 0;
        tJmax[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] - 1;
        tKmin[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 20 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 21
        tImin[ 21 ] = 0;
        tImax[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 21 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 22
        tImin[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 22 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 23
        tImin[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 23 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 23 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 24
        tImin[ 24 ] = 0;
        tImax[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] - 1;
        tJmin[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 24 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 24 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 25
        tImin[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 25 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 25 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // quadrant 26
        tImin[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] + 1;
        tImax[ 26 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ] - 1;
        tJmin[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] + 1;
        tJmax[ 26 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ] - 1;
        tKmin[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] + 1;
        tKmax[ 26 ] = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ] - 1;

        // loop over all quadrants
        for ( uint q = 0; q < 27; ++q )
        {
            // loop over my own elments
            if ( mMyProcNeighbors( q ) == par_rank() )
            {
                // this is an active element
                for ( auto k = tKmin[ q ]; k <= tKmax[ q ]; ++k )
                {
                    for ( auto j = tJmin[ q ]; j <= tJmax[ q ]; ++j )
                    {
                        for ( auto i = tImin[ q ]; i <= tImax[ q ]; ++i )
                        {
                            // get local id == local index for level 0
                            auto tIndex = calc_subdomain_id_of_element( 0, i, j, k );

                            // assume that all non owned elements are padding elements
                            mCoarsestElementsIncludingAura( tIndex )->set_active_flag( mActivePattern );

                            // set owner ( in this case: myself )
                            mCoarsestElementsIncludingAura( tIndex )->set_owner( mMyProcNeighbors( q ) );
                        }
                    }
                }
            }
            else if ( mMyProcNeighbors( q ) != gNoProcNeighbor )
            {
                // count elements for aura
                uint tCount = ( tImax[ q ] - tImin[ q ] + 1 )
                            * ( tJmax[ q ] - tJmin[ q ] + 1 )
                            * ( tKmax[ q ] - tKmin[ q ] + 1 );

                // set size for matrix

                mCoarsestAura( q ).set_size( tCount, 1 );

                // reset counter
                tCount = 0;

                // this is an active element of a neigbor proc
                for ( auto k = tKmin[ q ]; k <= tKmax[ q ]; ++k )
                {
                    for ( auto j = tJmin[ q ]; j <= tJmax[ q ]; ++j )
                    {
                        for ( auto i = tImin[ q ]; i <= tImax[ q ]; ++i )
                        {
                            // get local id == local index for level 0
                            auto tIndex = calc_subdomain_id_of_element( 0, i, j, k );

                            // assume that all non owned elements are padding elements
                            mCoarsestElementsIncludingAura( tIndex )->set_padding_flag();

                            // set owner
                            mCoarsestElementsIncludingAura( tIndex )->set_owner( mMyProcNeighbors( q ) );

                            // add element index to aura
                            mCoarsestAura( q )( tCount++ ) = tIndex;
                        }
                    }
                }
            }
            else
            {
                // this is definetly a padding element
                for ( auto k = tKmin[ q ]; k <= tKmax[ q ]; ++k )
                {
                    for ( auto j = tJmin[ q ]; j <= tJmax[ q ]; ++j )
                    {
                        for ( auto i = tImin[ q ]; i <= tImax[ q ]; ++i )
                        {
                            // get local id == local index for level 0
                            auto tIndex = calc_subdomain_id_of_element( 0, i, j, k );

                            // deactivate element
                            mCoarsestElementsIncludingAura( tIndex )->set_padding_flag();

                            // padding element does not belong to anybody
                            mCoarsestElementsIncludingAura( tIndex )->set_owner( gNoProcOwner );
                        }
                    }
                }
            }
        }

        // create quadrants for inverse aura

        // quadrant 0
        luint tDelta = mParameters->get_padding_size() - 1;
        tImin[ 0 ]   = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 0 ]   = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 0 ]   = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 0 ]   = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 0 ]   = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 0 ]   = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 1
        tImin[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 1 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 2
        tImin[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 2 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 3
        tImin[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 3 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 4
        tImin[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 4 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 5
        tImin[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 5 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 6
        tImin[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 6 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 7
        tImin[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 7 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 8
        tImin[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 8 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ] + tDelta;

        // quadrant 9
        tImin[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 9 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 10
        tImin[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 10 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 11
        tImin[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 11 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 12
        tImin[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 12 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 13
        tImin[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 13 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 14
        tImin[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 14 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 15
        tImin[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 15 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 16
        tImin[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 16 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 17
        tImin[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 0 ];
        tKmax[ 17 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 18
        tImin[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 18 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 19
        tImin[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 19 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 20
        tImin[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ] + tDelta;
        tKmin[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 20 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 21
        tImin[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 21 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 22
        tImin[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 22 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 23
        tImin[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 0 ];
        tJmax[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 23 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 24
        tImin[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ] + tDelta;
        tJmin[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 24 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 25
        tImin[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 0 ];
        tImax[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 25 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        // quadrant 26
        tImin[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ] - tDelta;
        tImax[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 0 ][ 1 ];
        tJmin[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ] - tDelta;
        tJmax[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 1 ][ 1 ];
        tKmin[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ] - tDelta;
        tKmax[ 26 ] = mMySubDomain.mFrameIJK[ 0 ][ 2 ][ 1 ];

        for ( uint q = 0; q < 27; ++q )
        {
            if ( mMyProcNeighbors( q ) != gNoProcNeighbor
                    && mMyProcNeighbors( q ) != par_rank() )
            {
                // count elements for inverse aura
                luint tCount = ( tImax[ q ] - tImin[ q ] + 1 )
                             * ( tJmax[ q ] - tJmin[ q ] + 1 )
                             * ( tKmax[ q ] - tKmin[ q ] + 1 );

                mCoarsestInverseAura( q ).set_size( tCount, 1 );

                tCount = 0;
                for ( auto k = tKmin[ q ]; k <= tKmax[ q ]; ++k )
                {
                    for ( auto j = tJmin[ q ]; j <= tJmax[ q ]; ++j )
                    {
                        for ( auto i = tImin[ q ]; i <= tImax[ q ]; ++i )
                        {
                            mCoarsestInverseAura( q )( tCount++ ) =
                                    calc_subdomain_id_of_element( 0, i, j, k );
                        }
                    }
                }
            }
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::refine_element(
            Background_Element_Base* aElement,
            const bool               aKeepState )
    {
        // only perform if element is not refined already
        // and element is below max defined level

        if ( aElement->get_level() < gMaxNumberOfLevels - 1 )
        {
            if ( !aElement->has_children() )
            {
                // get owner of element
                uint tOwner = aElement->get_owner();

                // get level of new element
                auto tLevel = aElement->get_level() + 1;

                // get ijk positions of children
                Matrix< DDLUMat > tIJK;
                aElement->get_ijk_of_children( tIJK );

                // ask background mesh for IDs of children
                Matrix< DDLUMat > tIDs;
                this->calc_element_ids(
                        tLevel,
                        tIJK,
                        tIDs );

                // tIDs.print("tIDs");
                //  temporary array for ijk position
                luint tCIJK[ 3 ];

                // child 0
                tCIJK[ 0 ] = tIJK( 0, 0 );
                tCIJK[ 1 ] = tIJK( 1, 0 );
                tCIJK[ 2 ] = tIJK( 2, 0 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 0 ),
                        tLevel,
                        (uint)0,
                        tOwner ) );

                // child 1
                tCIJK[ 0 ] = tIJK( 0, 1 );
                tCIJK[ 1 ] = tIJK( 1, 1 );
                tCIJK[ 2 ] = tIJK( 2, 1 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 1 ),
                        tLevel,
                        (uint)1,
                        tOwner ) );

                // child 2
                tCIJK[ 0 ] = tIJK( 0, 2 );
                tCIJK[ 1 ] = tIJK( 1, 2 );
                tCIJK[ 2 ] = tIJK( 2, 2 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 2 ),
                        tLevel,
                        (uint)2,
                        tOwner ) );

                // child 3
                tCIJK[ 0 ] = tIJK( 0, 3 );
                tCIJK[ 1 ] = tIJK( 1, 3 );
                tCIJK[ 2 ] = tIJK( 2, 3 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 3 ),
                        tLevel,
                        (uint)3,
                        tOwner ) );

                // child 4
                tCIJK[ 0 ] = tIJK( 0, 4 );
                tCIJK[ 1 ] = tIJK( 1, 4 );
                tCIJK[ 2 ] = tIJK( 2, 4 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 4 ),
                        tLevel,
                        (uint)4,
                        tOwner ) );

                // child 5
                tCIJK[ 0 ] = tIJK( 0, 5 );
                tCIJK[ 1 ] = tIJK( 1, 5 );
                tCIJK[ 2 ] = tIJK( 2, 5 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 5 ),
                        tLevel,
                        (uint)5,
                        tOwner ) );

                // child 6
                tCIJK[ 0 ] = tIJK( 0, 6 );
                tCIJK[ 1 ] = tIJK( 1, 6 );
                tCIJK[ 2 ] = tIJK( 2, 6 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 6 ),
                        tLevel,
                        (uint)6,
                        tOwner ) );

                // child 7
                tCIJK[ 0 ] = tIJK( 0, 7 );
                tCIJK[ 1 ] = tIJK( 1, 7 );
                tCIJK[ 2 ] = tIJK( 2, 7 );
                aElement->insert_child( new Background_Element< 3 >(
                        aElement,
                        mActivePattern,
                        tCIJK,
                        tIDs( 7 ),
                        tLevel,
                        (uint)7,
                        tOwner ) );

                // set refined switch
                if ( aKeepState )
                {
                    // loop over all children
                    for ( uint k = 0; k < 8; ++k )
                    {
                        // get pointer to child and deactivate element
                        aElement->get_child( k )->deactivate( mActivePattern );
                    }
                }
                else
                {
                    // set refined switch
                    aElement->set_refined_flag( mActivePattern );
                }

                // test if this is a padding element
                if ( aElement->is_padding() )
                {
                    // loop over all children
                    for ( uint k = 0; k < 8; ++k )
                    {
                        // get pointer to child and set refinement flag
                        aElement->get_child( k )->set_padding_flag();
                    }
                }
            }
            else    // element has children
            {
                // activate children if they are deactivated
                for ( uint k = 0; k < 8; ++k )
                {
                    // get child
                    auto tChild = aElement->get_child( k );

                    // test if child is deactivated
                    if ( tChild->is_neither_active_nor_refined( mActivePattern ) )
                    {
                        // activate child
                        tChild->set_active_flag( mActivePattern );
                    }
                }

                // refine element
                aElement->set_refined_flag( mActivePattern );
            }
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::collect_neighbors_on_level_zero()
    {

        //       k = -1                      k = 0                k = +1
        //   --------------              --------------       --------------
        //  | 21 |  8 | 20 |            | 13 |  2 | 12 |     | 25 | 16 | 24 |
        //   --------------   j          --------------       --------------
        //  |  9 |  4 |  7 |  ^         |  3 |    |  1 |     | 17 |  5 | 15 |
        //   --------------   |          --------------       --------------
        //  | 18 |  6 | 19 |  o-->i     | 10 |  0 | 11 |     | 22 | 14 | 23 |
        //   --------------              --------------       --------------
        // maximum i
        luint tIMax = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ];

        // maximum j
        luint tJMax = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ];

        // maximum k
        luint tKMax = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ];

        // matrix containing neighbor ids
        luint tNeighborIndex[ 26 ];

        // initialize element counter
        luint tCount = 0;

        // loop over all k
        for ( luint k = 0; k < tKMax; ++k )
        {
            // loop over all j
            for ( luint j = 0; j < tJMax; ++j )
            {
                // loop over all i
                for ( luint i = 0; i < tIMax; ++i )
                {

                    // calculate all 26 neighbor indices

                    tNeighborIndex[ 0 ] = this->calc_subdomain_id_of_element( 0, i, j - 1, k );

                    tNeighborIndex[ 1 ] = this->calc_subdomain_id_of_element( 0, i + 1, j, k );

                    tNeighborIndex[ 2 ] = this->calc_subdomain_id_of_element( 0, i, j + 1, k );

                    tNeighborIndex[ 3 ] = this->calc_subdomain_id_of_element( 0, i - 1, j, k );

                    tNeighborIndex[ 4 ] = this->calc_subdomain_id_of_element( 0, i, j, k - 1 );

                    tNeighborIndex[ 5 ] = this->calc_subdomain_id_of_element( 0, i, j, k + 1 );

                    tNeighborIndex[ 6 ] = this->calc_subdomain_id_of_element( 0, i, j - 1, k - 1 );

                    tNeighborIndex[ 7 ] = this->calc_subdomain_id_of_element( 0, i + 1, j, k - 1 );

                    tNeighborIndex[ 8 ] = this->calc_subdomain_id_of_element( 0, i, j + 1, k - 1 );

                    tNeighborIndex[ 9 ] = this->calc_subdomain_id_of_element( 0, i - 1, j, k - 1 );

                    tNeighborIndex[ 10 ] = this->calc_subdomain_id_of_element( 0, i - 1, j - 1, k );

                    tNeighborIndex[ 11 ] = this->calc_subdomain_id_of_element( 0, i + 1, j - 1, k );

                    tNeighborIndex[ 12 ] = this->calc_subdomain_id_of_element( 0, i + 1, j + 1, k );

                    tNeighborIndex[ 13 ] = this->calc_subdomain_id_of_element( 0, i - 1, j + 1, k );

                    tNeighborIndex[ 14 ] = this->calc_subdomain_id_of_element( 0, i, j - 1, k + 1 );

                    tNeighborIndex[ 15 ] = this->calc_subdomain_id_of_element( 0, i + 1, j, k + 1 );

                    tNeighborIndex[ 16 ] = this->calc_subdomain_id_of_element( 0, i, j + 1, k + 1 );

                    tNeighborIndex[ 17 ] = this->calc_subdomain_id_of_element( 0, i - 1, j, k + 1 );

                    tNeighborIndex[ 18 ] = this->calc_subdomain_id_of_element( 0, i - 1, j - 1, k - 1 );

                    tNeighborIndex[ 19 ] = this->calc_subdomain_id_of_element( 0, i + 1, j - 1, k - 1 );

                    tNeighborIndex[ 20 ] = this->calc_subdomain_id_of_element( 0, i + 1, j + 1, k - 1 );

                    tNeighborIndex[ 21 ] = this->calc_subdomain_id_of_element( 0, i - 1, j + 1, k - 1 );

                    tNeighborIndex[ 22 ] = this->calc_subdomain_id_of_element( 0, i - 1, j - 1, k + 1 );

                    tNeighborIndex[ 23 ] = this->calc_subdomain_id_of_element( 0, i + 1, j - 1, k + 1 );

                    tNeighborIndex[ 24 ] = this->calc_subdomain_id_of_element( 0, i + 1, j + 1, k + 1 );

                    tNeighborIndex[ 25 ] = this->calc_subdomain_id_of_element( 0, i - 1, j + 1, k + 1 );

                    // get element
                    Background_Element_Base* tElement = mCoarsestElementsIncludingAura( tCount++ );

                    // loop over all possible neighbors
                    for ( uint n = 0; n < 26; ++n )
                    {
                        // test if neighbor exists
                        if ( tNeighborIndex[ n ] != gNoEntityID )
                        {
                            // insert neighbor
                            tElement->insert_neighbor( n,
                                    mCoarsestElementsIncludingAura( tNeighborIndex[ n ] ) );
                        }
                    }
                }
            }
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::calc_corner_nodes_of_element(
            const Background_Element_Base* aElement,
            Matrix< DDRMat >&              aNodeCoords )
    {
        // get ijk position of element
        const luint* tIJK = aElement->get_ijk();

        // get level of element
        const uint tLevel = aElement->get_level();

        aNodeCoords.set_size( 3, 8 );

        // calculate node 0
        aNodeCoords( 0, 0 ) = ( tIJK[ 0 ]
                                      + mMySubDomain.mAuraIJK[ tLevel ][ 0 ][ 0 ] )
                                    * mElementLength[ tLevel ][ 0 ]
                            + mDomainOffset[ 0 ];

        aNodeCoords( 1, 0 ) = ( tIJK[ 1 ]
                                      + mMySubDomain.mAuraIJK[ tLevel ][ 1 ][ 0 ] )
                                    * mElementLength[ tLevel ][ 1 ]
                            + mDomainOffset[ 1 ];

        aNodeCoords( 2, 0 ) = ( tIJK[ 2 ]
                                      + mMySubDomain.mAuraIJK[ tLevel ][ 2 ][ 0 ] )
                                    * mElementLength[ tLevel ][ 2 ]
                            + mDomainOffset[ 2 ];

        // node 1
        aNodeCoords( 0, 1 ) = aNodeCoords( 0, 0 )
                            + mElementLength[ tLevel ][ 0 ];
        aNodeCoords( 1, 1 ) = aNodeCoords( 1, 0 );
        aNodeCoords( 2, 1 ) = aNodeCoords( 2, 0 );

        // node 2
        aNodeCoords( 0, 2 ) = aNodeCoords( 0, 1 );
        aNodeCoords( 1, 2 ) = aNodeCoords( 1, 1 )
                            + mElementLength[ tLevel ][ 1 ];
        aNodeCoords( 2, 2 ) = aNodeCoords( 2, 1 );

        // node 3
        aNodeCoords( 0, 3 ) = aNodeCoords( 0, 0 );
        aNodeCoords( 1, 3 ) = aNodeCoords( 1, 2 );
        aNodeCoords( 2, 3 ) = aNodeCoords( 2, 2 );

        // node 4
        aNodeCoords( 0, 4 ) = aNodeCoords( 0, 0 );
        aNodeCoords( 1, 4 ) = aNodeCoords( 1, 0 );
        aNodeCoords( 2, 4 ) = aNodeCoords( 2, 0 )
                            + mElementLength[ tLevel ][ 2 ];

        // node 5
        aNodeCoords( 0, 5 ) = aNodeCoords( 0, 1 );
        aNodeCoords( 1, 5 ) = aNodeCoords( 1, 1 );
        aNodeCoords( 2, 5 ) = aNodeCoords( 2, 4 );

        // node 6
        aNodeCoords( 0, 6 ) = aNodeCoords( 0, 2 );
        aNodeCoords( 1, 6 ) = aNodeCoords( 1, 2 );
        aNodeCoords( 2, 6 ) = aNodeCoords( 2, 4 );

        // node 7
        aNodeCoords( 0, 7 ) = aNodeCoords( 0, 3 );
        aNodeCoords( 1, 7 ) = aNodeCoords( 1, 3 );
        aNodeCoords( 2, 7 ) = aNodeCoords( 2, 4 );
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::calc_center_of_element(
            const Background_Element_Base* aElement,
            Matrix< DDRMat >&              aNodeCoords )
    {
        // get ijk position of element
        const luint* tIJK = aElement->get_ijk();

        // get level of element
        const uint tLevel = aElement->get_level();

        aNodeCoords.set_size( 3, 1 );

        aNodeCoords( 0 ) = ( 0.5 + tIJK[ 0 ]
                                   + mMySubDomain.mAuraIJK[ tLevel ][ 0 ][ 0 ] )
                                 * mElementLength[ tLevel ][ 0 ]
                         + mDomainOffset[ 0 ];

        aNodeCoords( 1 ) = ( 0.5 + tIJK[ 1 ]
                                   + mMySubDomain.mAuraIJK[ tLevel ][ 1 ][ 0 ] )
                                 * mElementLength[ tLevel ][ 1 ]
                         + mDomainOffset[ 1 ];

        aNodeCoords( 2 ) = ( 0.5 + tIJK[ 2 ]
                                   + mMySubDomain.mAuraIJK[ tLevel ][ 2 ][ 0 ] )
                                 * mElementLength[ tLevel ][ 2 ]
                         + mDomainOffset[ 2 ];
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::collect_coarsest_elements_on_side(
            uint                       aSideOrdinal,
            Cell< Background_Element_Base* >& aCoarsestElementsOnSide )
    {
        // clear output cell
        aCoarsestElementsOnSide.clear();

        // number of elements
        luint tNumberOfElementsI = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 0 ]
                                 - 2 * mParameters->get_padding_size();

        luint tNumberOfElementsJ = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 1 ]
                                 - 2 * mParameters->get_padding_size();

        luint tNumberOfElementsK = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ 2 ]
                                 - 2 * mParameters->get_padding_size();

        switch ( aSideOrdinal )
        {
            case ( 1 ):
            {
                // test if proc is on edge of domain
                if ( mCoarsestElements( 0 )->get_neighbor( 0 )->get_owner() == gNoProcID )
                {
                    // allocate cell
                    aCoarsestElementsOnSide.resize( tNumberOfElementsI * tNumberOfElementsK, nullptr );

                    luint tPivot = 0;
                    luint tCount = 0;

                    // loop over all coarsest elements
                    for ( luint k = 0; k < tNumberOfElementsK; ++k )
                    {
                        for ( luint i = 0; i < tNumberOfElementsI; ++i )
                        {
                            aCoarsestElementsOnSide( tCount++ ) = mCoarsestElements( tPivot++ );
                        }
                        tPivot += ( tNumberOfElementsJ - 1 ) * tNumberOfElementsI;
                    }
                }
                break;
            }
            case ( 2 ):
            {
                // test if proc is on edge of domain
                if ( mCoarsestElements( mCoarsestElements.size() - 1 )->get_neighbor( 1 )->get_owner() == gNoProcID )
                {

                    // allocate cell
                    aCoarsestElementsOnSide.resize( tNumberOfElementsJ * tNumberOfElementsK, nullptr );

                    luint tPivot = tNumberOfElementsI - 1;
                    luint tCount = 0;

                    // loop over all coarsest elements
                    for ( luint k = 0; k < tNumberOfElementsK; ++k )
                    {
                        for ( luint j = 0; j < tNumberOfElementsJ; ++j )
                        {
                            aCoarsestElementsOnSide( tCount++ ) = mCoarsestElements( tPivot );
                            tPivot += tNumberOfElementsI;
                        }
                    }
                }
                break;
            }
            case ( 3 ):
            {
                // test if proc is on edge of domain
                if ( mCoarsestElements( mCoarsestElements.size() - 1 )->get_neighbor( 2 )->get_owner() == gNoProcID )
                {

                    // allocate cell
                    aCoarsestElementsOnSide.resize( tNumberOfElementsI * tNumberOfElementsK, nullptr );
                    luint tJump  = tNumberOfElementsI * ( tNumberOfElementsJ - 1 );
                    luint tPivot = tJump;
                    luint tCount = 0;

                    // loop over all coarsest elements
                    for ( luint k = 0; k < tNumberOfElementsK; ++k )
                    {
                        for ( luint i = 0; i < tNumberOfElementsI; ++i )
                        {
                            aCoarsestElementsOnSide( tCount++ ) = mCoarsestElements( tPivot++ );
                        }
                        tPivot += tJump;
                    }
                }
                break;
            }
            case ( 4 ):
            {
                // test if proc is on edge of domain
                if ( mCoarsestElements( 0 )->get_neighbor( 3 )->get_owner() == gNoProcID )
                {

                    // allocate cell
                    aCoarsestElementsOnSide.resize( tNumberOfElementsJ * tNumberOfElementsK, nullptr );

                    luint tPivot = 0;
                    luint tCount = 0;

                    // loop over all coarsest elements
                    for ( luint k = 0; k < tNumberOfElementsK; ++k )
                    {
                        for ( luint j = 0; j < tNumberOfElementsJ; ++j )
                        {
                            aCoarsestElementsOnSide( tCount++ ) = mCoarsestElements( tPivot );
                            tPivot += tNumberOfElementsI;
                        }
                    }
                }
                break;
            }
            case ( 5 ):
            {
                // test if proc is on edge of domain
                if ( mCoarsestElements( 0 )->get_neighbor( 4 )->get_owner() == gNoProcID )
                {

                    // allocate cell
                    aCoarsestElementsOnSide.resize( tNumberOfElementsI * tNumberOfElementsJ, nullptr );

                    luint tCount = 0;

                    // loop over all coarsest elements
                    for ( luint j = 0; j < tNumberOfElementsJ; ++j )
                    {
                        for ( luint i = 0; i < tNumberOfElementsI; ++i )
                        {
                            aCoarsestElementsOnSide( tCount ) = mCoarsestElements( tCount );
                            tCount++;
                        }
                    }
                }
                break;
            }
            case ( 6 ):
            {
                // test if proc is on edge of domain
                if ( mCoarsestElements( mCoarsestElements.size() - 1 )->get_neighbor( 5 )->get_owner() == gNoProcID )
                {

                    // allocate cell
                    aCoarsestElementsOnSide.resize( tNumberOfElementsI * tNumberOfElementsJ, nullptr );

                    luint tCount = 0;
                    luint tPivot = tNumberOfElementsI * tNumberOfElementsJ * ( tNumberOfElementsK - 1 );
                    // loop over all coarsest elements
                    for ( luint j = 0; j < tNumberOfElementsJ; ++j )
                    {
                        for ( luint i = 0; i < tNumberOfElementsI; ++i )
                        {
                            aCoarsestElementsOnSide( tCount++ ) = mCoarsestElements( tPivot++ );
                        }
                    }
                }
                break;
            }
        }
    }

    //-------------------------------------------------------------------------------

    template<>
    inline void
    Background_Mesh< 3 >::collect_coarsest_elements_in_bounding_box(
            moris::Cell< Background_Element_Base* >& aBackgroundElements,
            luint                                    aBoundingBoxStartEndIJK[][ 2 ],
            uint                                     alevel )
    {
        aBackgroundElements.resize( ( aBoundingBoxStartEndIJK[ 0 ][ 1 ] - aBoundingBoxStartEndIJK[ 0 ][ 0 ] ) *            //
                                            ( aBoundingBoxStartEndIJK[ 1 ][ 1 ] - aBoundingBoxStartEndIJK[ 1 ][ 0 ] ) *    //
                                            ( aBoundingBoxStartEndIJK[ 2 ][ 1 ] - aBoundingBoxStartEndIJK[ 2 ][ 0 ] ),
                nullptr );

        luint tCounter = 0;

        for ( luint Ik = aBoundingBoxStartEndIJK[ 0 ][ 0 ]; Ik < aBoundingBoxStartEndIJK[ 0 ][ 1 ]; ++Ik )
        {
            for ( luint Ii = aBoundingBoxStartEndIJK[ 1 ][ 0 ]; Ii < aBoundingBoxStartEndIJK[ 1 ][ 1 ]; ++Ii )
            {
                for ( luint Ij = aBoundingBoxStartEndIJK[ 2 ][ 0 ]; Ij < aBoundingBoxStartEndIJK[ 2 ][ 1 ]; ++Ij )
                {
                    uint tId = calc_subdomain_id_of_element( alevel,
                            Ik,
                            Ii,
                            Ij );

                    aBackgroundElements( tCounter++ ) = mCoarsestElementsIncludingAura( tId );
                }
            }
        }
    }

    //-------------------------------------------------------------------------------

}    // namespace moris

#endif /* SRC_HMR_CL_HMR_BACKGROUND_MESH_3D_HPP_ */
