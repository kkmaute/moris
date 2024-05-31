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
#include "cl_Vector.hpp"
#include "cl_SDF_Triangle.hpp"
#include "cl_SDF_Line.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------
        class Object
        {
            const real                                mMeshHighPass = 1e-9;
            moris::Vector< std::shared_ptr< Facet > > mFacets;

            uint mNumberOfFacets;

            real mIntersectionTolerance = 1e-8;    // tolerance for interfaces when raycasting with this Object

          protected:
            uint                                      mDimension;
            Vector< std::shared_ptr< Facet_Vertex > > mVertices;    // vertices of all facets, can be modified by ADVs


            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Object( const std::string&           aFilePath,
                    real                         aIntersectionTolerance = 1e-8,
                    const moris::Vector< real >& aOffsets               = { 0, 0, 0 },
                    const moris::Vector< real >& aScale                 = { 1.0, 1.0, 1.0 } );


            //-------------------------------------------------------------------------------

            /**
             * Performs a coordinate rotation of the object's facets and vertices
             * NOTE: This action itself can be undone by calling reset_coordinates(), but will remove any scaling or shift that may have occurred
             *
             * @param aRotationMatrix the direction cosine matrix defining the rotation
             */
            void
            rotate( const Matrix< DDRMat >& aRotationMatrix );

            //-------------------------------------------------------------------------------

            /**
             * Scales all the coordinates of the object.
             * NOTE: This action can be undone by calling reset_coordinates(), but will remove any shift or rotation that may have occurred
             *
             * @param aScaling factor to scale in each coordinate direction
             */
            void
            scale( const Vector< real >& aScaling );

            //-------------------------------------------------------------------------------

            /**
             * Moves the object's spatial position.
             * NOTE: This action can be undone by calling reset_coordinates(), but will remove any scaling or rotation that may have occurred
             *
             * @param aShift shift in each coordinate direction that is added to the objects coordinates.
             */
            void
            shift( const Vector< real >& aShift );

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
            get_intersection_tolerance()
            {
                return mIntersectionTolerance;
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

            /**
             * Outputs the surface mesh to a .obj file
             * 
             * @param aFilePath the desired file name and location
             */
            void
            write_to_file( std::string aFilePath );


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

            //-------------------------------------------------------------------------------

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

            //-------------------------------------------------------------------------------

            /**
             * loads an ASCII file into a buffer of strings.
             * Called through construction.
             */
            void
            load_ascii_to_buffer( const std::string& aFilePath,
                    Vector< std::string >&           aBuffer );

            //-------------------------------------------------------------------------------
        };
        //-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */
