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

#include "typedefs.hpp"
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
            const real                   mMeshHighPass = 1e-9;
            moris::Cell< Facet_Vertex* > mVertices;
            moris::Cell< Facet* >        mFacets;
            Matrix< DDRMat >             mFacetMinCoords;
            Matrix< DDRMat >             mFacetMaxCoords;

            const Matrix< DDRMat >& mOffsets;
            uint                    mDimension;
            uint                    mNumberOfFacets;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Object( const std::string&      aFilePath,
                    const Matrix< DDRMat >& aOffsets = { { 0, 0, 0 } } );

            //-------------------------------------------------------------------------------

            ~Object();

            //-------------------------------------------------------------------------------

            /**
             * Performs a coordinate rotation of the object's facets and vertices
             *
             * @param aRotationMatrix the direction cosine matrix defining the rotation
             */
            void
            rotate_object( Matrix< DDRMat >& aRotationMatrix );

            //-------------------------------------------------------------------------------

            void
            undo_rotation();

            //-------------------------------------------------------------------------------

            moris::Cell< Facet* >&
            get_facets()
            {
                return mFacets;
            }

            //-------------------------------------------------------------------------------

            uint
            get_num_facets()
            {
                return mNumberOfFacets;
            }

            //-------------------------------------------------------------------------------

            moris::Cell< Facet_Vertex* >&
            get_vertices()
            {
                return mVertices;
            }

            //-------------------------------------------------------------------------------

            uint
            get_dimension()
            {
                return mDimension;
            }

            //-------------------------------------------------------------------------------

            Matrix< DDRMat >&
            get_facet_min_coords()
            {
                return mFacetMinCoords;
            }

            //-------------------------------------------------------------------------------

            Matrix< DDRMat >&
            get_facet_max_coords()
            {
                return mFacetMaxCoords;
            }

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
            load_from_object_file( const std::string& aFilePath );

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