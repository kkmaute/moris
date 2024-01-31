/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Object.hpp
 *
 */

#pragma once

#include <cl_SDF_Facet_Vertex.hpp>
#include <string>

#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_SDF_Triangle.hpp"
#include "cl_SDF_Line.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------
        class Object
        {
            const real mMeshHighPass = 1e-9;
            moris::Cell< std::shared_ptr< Facet > > mFacets;

            uint mDimension;
            uint mNumberOfFacets;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Object( const std::string&      aFilePath,
                    const moris::Cell< real >& aOffsets = { 0, 0, 0 },
                    const moris::Cell< real>& aScale = { 1.0, 1.0, 1.0 } );


            //-------------------------------------------------------------------------------

            /**
             * Performs a coordinate rotation of the object's facets and vertices
             * NOTE: This action itself cannot be undone without using reset_object_coordinates, which will also remove any applied scaling or translation.
             *
             * @param aRotationMatrix the direction cosine matrix defining the rotation
             */
            void
            rotate( const Matrix< DDRMat >& aRotationMatrix );

            //-------------------------------------------------------------------------------

            /**
             * Scales all the coordinates of the object.
             * NOTE: This action can be undone by calling scale_object( aScaling^-1 )
             *
             * @param aScaling factor to scale in each coordinate direction
             */
            void
            scale( const moris::Cell< real >& aScaling );

            //-------------------------------------------------------------------------------

            /**
             * Moves the object's spatial position.
             * NOTE: This action can be undone by calling translate_object( -aShift )
             *
             * @param aShift shift in each coordinate direction that is added to the objects coordinates.
             */
            void
            shift( const moris::Cell< real >& aShift );

            //-------------------------------------------------------------------------------

            /**
             * Resets the object back to its attitude when it was constructed,
             * removing any rotation, scaling, or translation that was applied
             *
             */
            void
            reset_coordinates();

            //-------------------------------------------------------------------------------

            Facet&
            get_facet( uint aFacetIndex )
            {
                MORIS_ASSERT( aFacetIndex >= 0, "SDF_Object: get_facet - aFacetIndex needs to be >= 0. Provided index: %u", aFacetIndex );
                return *mFacets( aFacetIndex );
            }

            //-------------------------------------------------------------------------------

            uint
            get_num_facets()
            {
                return mNumberOfFacets;
            }

            //-------------------------------------------------------------------------------

            uint
            get_dimension() const
            {
                return mDimension;
            }

            //-------------------------------------------------------------------------------

            real
            get_facet_min_coord(
                    uint aFacetIndex,
                    uint aAxis );


            //-------------------------------------------------------------------------------

            real
            get_facet_max_coord(
                    uint aFacetIndex,
                    uint aAxis );


            //-------------------------------------------------------------------------------
            // MTK
            //-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

          private:
            //-------------------------------------------------------------------------------

            /**
             * loads an ascii file and creates vertex and facet objects
             * Facets are either lines in 2D or triangles in 3D
             */
            void
            load_from_object_file( const std::string& aFilePath, const moris::Cell< real >& aOffsets, const moris::Cell< real >& aScale );

            //-------------------------------------------------------------------------------

            /**
             * loads an ascii file and creates vertex and facet objects
             * Facets are either lines in 2D or triangles in 3D
             */
            void
            load_from_stl_file( const std::string& aFilePath );

            //-------------------------------------------------------------------------------

            /**
             * loads an ASCII file into a buffer of strings.
             * Called through construction.
             */
            void
            load_ascii_to_buffer( const std::string& aFilePath,
                    moris::Cell< std::string >&      aBuffer );

            //-------------------------------------------------------------------------------
        };
        //-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */
