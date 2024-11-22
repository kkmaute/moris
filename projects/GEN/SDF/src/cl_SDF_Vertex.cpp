/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Vertex.cpp
 *
 */

#include "cl_SDF_Vertex.hpp"
#include "cl_SDF_Cell.hpp"
#include "cl_SDF_Triangle.hpp"
#include "op_times.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------
    /**
     * constructor
     */
    Vertex::Vertex( const moris_index aIndex, const Matrix< DDRMat >& aNodeCoords )
            : mIndex( aIndex )
            , mNodeCoords( 1, 3 )
            , mOriginalNodeCoords( 3, 1 )
    {
        // convert dynamic array to fixed array
        for ( uint k = 0; k < 3; ++k )
        {
            mOriginalNodeCoords( k ) = aNodeCoords( k );
        }

        this->reset_coords();
    }

    //-------------------------------------------------------------------------------

    void
    Vertex::update_udf( Facet& aFacet )
    {
        // calculate distance to this point
        real tDistance = aFacet.get_distance_to_point( mNodeCoords );

        if ( tDistance < mSDF )
        {
            // remember value
            mSDF = tDistance;

            // set sdf flag
            mHasSDF = true;

            // remember triangle
            mClosestFacet = &aFacet;
        }
    }

    //-------------------------------------------------------------------------------

    void
    Vertex::insert_cell( Cell* aCell )
    {
        mCells( mCellCounter++ ) = aCell;
    }

    // -----------------------------------------------------------------------------

    uint
    Vertex::sweep()
    {
        bool tSwept = false;

        // loop over all neighbors
        for ( Vertex* tNeighbor : mNeighbors )
        {
            // get pointer to facet
            Facet* tFacet = tNeighbor->get_closest_facet();

            if ( tFacet != nullptr )
            {
                // get distance to triangle of neighbor
                real tDistance = tFacet->get_distance_to_point( mNodeCoords );

                if ( tDistance < mSDF )
                {
                    tSwept        = true;
                    mSDF          = tDistance;
                    mClosestFacet = tFacet;
                }
            }
        }

        if ( tSwept )
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }

    //-------------------------------------------------------------------------------

    void
    Vertex::rotate_coords( const Matrix< F33RMat >& aRotationMatrix )
    {
        mNodeCoords = aRotationMatrix * mOriginalNodeCoords;
    }

    // -----------------------------------------------------------------------------

    void
    Vertex::reset_coords()
    {
        mNodeCoords = mOriginalNodeCoords;
    }

    // -----------------------------------------------------------------------------

}    // namespace moris::sdf
