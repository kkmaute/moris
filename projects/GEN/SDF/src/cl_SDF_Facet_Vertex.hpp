/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Facet_Vertex.hpp
 *
 */

#pragma once

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        class Triangle;

        //-------------------------------------------------------------------------------
        class Facet_Vertex : public mtk::Vertex
        {

            const moris_index mIndex;

            //-------------------------------------------------------------------------------

            Matrix< DDRMat > mNodeCoords;

            const Matrix< DDRMat > mOriginalNodeCoords;
            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Facet_Vertex(
                    const moris_index       aIndex,
                    const Matrix< DDRMat >& aNodeCoords );

            //-------------------------------------------------------------------------------

            ~Facet_Vertex(){};

            //-------------------------------------------------------------------------------
            // Special SDF functions
            //-------------------------------------------------------------------------------

            void
            rotate_node_coords( const Matrix< F33RMat >& aRotationMatrix );

            //-------------------------------------------------------------------------------

            void
            scale_node_coords( const Matrix< DDRMat >& aScaling )
            {
                for ( uint iAxis = 0; iAxis < mNodeCoords.numel(); iAxis++ )
                {
                    mNodeCoords( iAxis ) *= aScaling( iAxis );
                }
            }

            //-------------------------------------------------------------------------------

            void
            shift_node_coords( const Matrix< DDRMat >& aShift )
            {
                mNodeCoords += aShift;
            }

            //-------------------------------------------------------------------------------

            void
            reset_node_coords();

            //-------------------------------------------------------------------------------
            // MTK API functions
            //-------------------------------------------------------------------------------

            real
            get_coord( uint aAxis ) const
            {
                MORIS_ASSERT( aAxis <= mNodeCoords.numel(), "SDF_Facet_Vertex::get_coord() - Provided axis of %u exceeds the dimension of the vertex.", aAxis );
                return mNodeCoords( aAxis );
            }

            //-------------------------------------------------------------------------------

            moris_id
            get_id() const
            {
                return mIndex + 1;
            }

            //-------------------------------------------------------------------------------

            moris_index
            get_index() const
            {
                return mIndex;
            }

            //-------------------------------------------------------------------------------

            moris_id
            get_owner() const
            {
                return 0;
            }

            //-------------------------------------------------------------------------------

            mtk::Vertex_Interpolation*
            get_interpolation( const uint aOrder )
            {
                MORIS_ERROR( false,
                        "get_interpolation() is not available for an SDF Vertex" );
                return nullptr;
            }

            //-------------------------------------------------------------------------------

            const mtk::Vertex_Interpolation*
            get_interpolation( const uint aOrder ) const
            {
                MORIS_ERROR( false,
                        "get_interpolation() is not available for an SDF Vertex" );
                return nullptr;
            }

            //-------------------------------------------------------------------------------
        };
        //-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */
