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
            : mFacets( aObject.get_facets() )
            , mVertices( aObject.get_vertices() )
            , mNumberOfFacets( mFacets.size() )
            , mFacetMinCoords( mNumberOfFacets, aObject.get_dimension() )
            , mFacetMaxCoords( mNumberOfFacets, aObject.get_dimension() )
            // #ifdef MORIS_USE_ARMA
            //             mCandI( mNumberOfFacets )
            //             , mCandJ( mNumberOfFacets )
            //             , mCandK( mNumberOfFacets )

            // #else
            , mCandJ( mNumberOfFacets, 1 )
            // #endif
            , mCandidateFacets( mNumberOfFacets, 1 )

    {
        // copy bounding box data
        for ( uint iFacetIndex = 0; iFacetIndex < mNumberOfFacets; iFacetIndex++ )
        {
            for ( uint iDimensionIndex = 0; iDimensionIndex < aObject.get_dimension(); iDimensionIndex++ )
            {
                mFacetMinCoords( iFacetIndex, iDimensionIndex ) = mFacets( iFacetIndex )->get_min_coord( iDimensionIndex );

                mFacetMaxCoords( iFacetIndex, iDimensionIndex ) = mFacets( iFacetIndex )->get_max_coord( iDimensionIndex );
            }
        }
    }
} /* namespace moris::sdf */
