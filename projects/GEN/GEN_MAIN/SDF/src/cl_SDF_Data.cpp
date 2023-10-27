/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Data.cpp
 *
 */

#include "cl_SDF_Data.hpp"
namespace moris::sdf
{

    //-------------------------------------------------------------------------------
    Data::Data( Object &aObject )
            : mDimension( aObject.get_dimension() )
            , mFacets( aObject.get_facets() )
            , mVertices( aObject.get_vertices() )
            , mNumberOfFacets( mFacets.size() )
// #ifdef MORIS_USE_ARMA
//             , mFacetMinCoordsX( mNumberOfFacets, 1 )
//             , mFacetMinCoordsY( mNumberOfFacets, 1 )
//             , mFacetMinCoordsZ( mNumberOfFacets, 1 )
//             , mFacetMaxCoordsX( mNumberOfFacets, 1 )
//             , mFacetMaxCoordsY( mNumberOfFacets, 1 )
//             , mFacetMaxCoordsZ( mNumberOfFacets, 1 )
//             , mCandI( mNumberOfFacets )
//             , mCandJ( mNumberOfFacets )
//             , mCandK( mNumberOfFacets )

// #else
// #endif

    {
        // copy bounding box data
// #ifdef MORIS_USE_ARMA
//         for ( uint iFacetIndex = 0; iFacetIndex < mNumberOfFacets; iFacetIndex++ )
//         {
//             for ( uint iDimensionIndex = 0; iDimensionIndex < aObject.get_dimension(); iDimensionIndex++ )
//             {
//                 mFacetMinCoordsX( iFacetIndex ) = mFacets( iFacetIndex )->get_min_coord( 0 );
//                 mFacetMinCoordsY( iFacetIndex ) = mFacets( iFacetIndex )->get_min_coord( 1 );

//                 mFacetMaxCoordsX( iFacetIndex ) = mFacets( iFacetIndex )->get_min_coord( 0 );
//                 mFacetMaxCoordsY( iFacetIndex ) = mFacets( iFacetIndex )->get_min_coord( 1 );
//             }
//             if ( aObject.get_dimension() == 3 )
//             {
//                 for ( uint iDimensionIndex = 0; iDimensionIndex < aObject.get_dimension(); iDimensionIndex++ )
//                 {
//                     mFacetMinCoordsZ( iFacetIndex ) = mFacets( iFacetIndex )->get_min_coord( 2 );

//                     mFacetMaxCoordsZ( iFacetIndex ) = mFacets( iFacetIndex )->get_min_coord( 2 );
//                 }
//             }
//         }
// #else
// #endif
    }
} /* namespace moris::sdf */
