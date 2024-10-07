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

#include "cl_MTK_Surface_Mesh.hpp"

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_SDF_Triangle.hpp"
#include "cl_SDF_Line.hpp"

namespace moris::sdf
{
        //-------------------------------------------------------------------------------
        class Object : public mtk::Surface_Mesh
        {
          protected:
            Vector< std::shared_ptr< Facet_Vertex > > mVertices  = {};    // vertices of all facets, can be modified by ADVs
            Vector< std::shared_ptr< Facet > >        mFacets = {};

          private:
            // real             mIntersectionTolerance = 1e-8;    // tolerance for interfaces when raycasting with this Object brendan moved
            const real mMeshHighPass = 1e-9;


        //-------------------------------------------------------------------------------

      public:
        //-------------------------------------------------------------------------------

        Object( const std::string&           aFilePath,
                real                         aIntersectionTolerance = 1e-8,
                const moris::Vector< real >& aOffsets               = { 0, 0, 0 },
                const moris::Vector< real >& aScale                 = { 1.0, 1.0, 1.0 } );

        //-------------------------------------------------------------------------------


        Facet&
        get_facet( uint aFacetIndex )
        {
            MORIS_ASSERT( aFacetIndex >= 0, "SDF_Object: get_facet - aFacetIndex needs to be >= 0. Provided index: %u", aFacetIndex );
            return *mFacets( aFacetIndex );
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

        protected:
        /**
         * Updates all member data for each facet of the object such as Hesse distance, normal, and center
         *
         */
        void
        update_all_facets();

        private:
        //-------------------------------------------------------------------------------

        /**
         * loads an ascii file and creates vertex and facet objects
         * Facets are either lines in 2D or triangles in 3D
         */
        void
        load_from_object_file( const std::string& aFilePath, const Vector< real >& aOffsets, const Vector< real >& aScale );

        //-------------------------------------------------------------------------------

        /**
         * loads an ascii file and creates vertex and facet objects
         * Facets are either lines in 2D or triangles in 3D
         */
        void
        load_from_stl_file( const std::string& aFilePath );

    };
    //-------------------------------------------------------------------------------
}    // namespace moris::sdf
