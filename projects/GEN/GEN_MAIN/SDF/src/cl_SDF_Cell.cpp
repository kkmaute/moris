/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Cell.cpp
 *
 */

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_SDF_Cell.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

    Cell::Cell(
            const moris_index aIndex,
            const moris_id    aID,
            const Matrix< IndexMat > & aIndices,
            moris::Vector< Vertex * >  & aAllVertices ) :
                    mIndex( aIndex ),
                    mID( aID )
    {
        // get number of vertices
        uint tNumberOfVertices = aIndices.length();

        // allocate member cell
        mVertices.resize( tNumberOfVertices, nullptr );

        // populate cell
        for( uint k=0; k<tNumberOfVertices; ++k )
        {
            mVertices( k ) = aAllVertices( aIndices( k ) );
        }
    }

//-------------------------------------------------------------------------------

        real
        Cell::get_buffer_diagonal()
        {
            // calculate min and max values of this cell
            Matrix< F31RMat > tMinValue  = mVertices( 0 )->get_coords();
            Matrix< F31RMat > tMaxValue = mVertices( 0 )->get_coords();

            // get number of nodes
            uint tNumberOfNodes = mVertices.size();

            // loop over all n
            for( uint k=1; k< tNumberOfNodes; ++k )
            {
                // get coordinates of node
                const Matrix< F31RMat > & tCoords =  mVertices( k )->get_coords();

                // get min and max coords
                for( uint i=0; i<3; ++i )
                {
                    tMinValue( i ) = std::min( tMinValue( i ), tCoords( i ) );
                    tMaxValue( i ) = std::max( tMaxValue( i ), tCoords( i ) );
                }
            }

            real aNorm = 0.0;
            for( uint i=0; i<3; ++i )
            {
                aNorm += std::pow( tMaxValue( i ) - tMinValue( i ), 2 );
            }
            return std::sqrt( aNorm );
        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

