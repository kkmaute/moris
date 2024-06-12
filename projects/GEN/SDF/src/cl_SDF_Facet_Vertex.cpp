/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Facet_Vertex.cpp
 *
 */

#include "cl_SDF_Facet_Vertex.hpp"
#include "op_times.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Facet_Vertex::Facet_Vertex(
                const moris_index       aIndex,
                const Matrix< DDRMat > &aNodeCoords )
                : mIndex( aIndex )
                , mNodeCoords( aNodeCoords )
                , mIterationNodeCoords( aNodeCoords )
        {
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::rotate_node_coords( const Matrix< DDRMat > &aRotationMatrix )
        {
            mNodeCoords = aRotationMatrix * mNodeCoords;
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::scale_node_coords( const moris::Vector< real > &aScaling )
        {
            for ( uint iAxis = 0; iAxis < mNodeCoords.numel(); iAxis++ )
            {
                mNodeCoords( iAxis ) *= aScaling( iAxis );
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::set_node_coords( const moris::Vector< real > &aCoordinates )
        {
            for ( uint iAxis = 0; iAxis < mNodeCoords.numel(); iAxis++ )
            {
                mNodeCoords( iAxis )          = aCoordinates( iAxis );
                mIterationNodeCoords( iAxis ) = aCoordinates( iAxis );
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::set_node_coords( const Matrix< DDRMat >& aCoordinates )
        {
            for ( uint iAxis = 0; iAxis < mNodeCoords.numel(); iAxis++ )
            {
                mNodeCoords( iAxis )          = aCoordinates( iAxis );
                mIterationNodeCoords( iAxis ) = aCoordinates( iAxis );
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::set_node_coord( const real aCoordinate, uint aDimension )
        {
            mNodeCoords( aDimension )          = aCoordinate;
            mIterationNodeCoords( aDimension ) = aCoordinate;
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::shift_node_coords_from_current( const moris::Vector< real > &aShift )
        {
            for ( uint iAxis = 0; iAxis < mNodeCoords.numel(); iAxis++ )
            {
                mNodeCoords( iAxis ) += aShift( iAxis );
            }
        }

        //-------------------------------------------------------------------------------

        void
        Facet_Vertex::reset_node_coords()
        {
            mNodeCoords = mIterationNodeCoords;
        }

        //-------------------------------------------------------------------------------

        Matrix< DDRMat >
        Facet_Vertex::get_coords() const
        {
            return mNodeCoords;
        }

    } /* namespace sdf */
} /* namespace moris */
