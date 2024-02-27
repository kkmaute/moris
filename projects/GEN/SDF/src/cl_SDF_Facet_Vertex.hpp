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
                , std::enable_shared_from_this< Facet_Vertex >
        {

            const moris_index mIndex;

            //-------------------------------------------------------------------------------

            Matrix< DDRMat > mNodeCoords;             // can be changed through coordinate transformations in raycasts and design iterations
            Matrix< DDRMat > mIterationNodeCoords;    // Can only be changed by design iterations

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
            scale_node_coords( const moris::Cell< real >& aScaling );


            //-------------------------------------------------------------------------------

            /**
             * Sets the coordinates of this node. Also changes mIterationNodeCoords. As such, this function should ONLY be called between design iterations.
             * 
             * @param aCoordinates new coordinates to be set to
             */
            void
            set_node_coords( const moris::Cell< real >& aCoordinates );
            
            //-------------------------------------------------------------------------------

            /**
             * Shifts only the node coordinates. Used for raycasts and coordinate transformations
             *
             * @param aShift Perturbation applied to the coordinates. size of mNodeCoords
             */
            void
            shift_node_coords_from_current( const moris::Cell< real >& aShift );

            //-------------------------------------------------------------------------------

            /**
             * Sets the coordinates of the node back to the coordinates at the start of the design iteration
             * 
             */
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

            Matrix< DDRMat > 
            get_coords() const override;

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

            uint
            get_dimension() const
            {
                return mNodeCoords.numel();
            }
        };
        //-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */
