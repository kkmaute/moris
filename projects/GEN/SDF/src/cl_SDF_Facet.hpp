/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Facet.hpp
 *
 */

#pragma once

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "op_times.hpp"

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_SDF_Facet_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        /**
         * the facet class for the sdf generator
         */
        class Facet : public mtk::Cell, std::enable_shared_from_this< Facet >
        {
          protected:
            // index of this facet
            moris_index mIndex;

            // cells with vertex pointers
            moris::Cell< std::shared_ptr< Facet_Vertex > > mVertices;

            // container for center
            Matrix< DDRMat > mCenter;

            // container for normal
            Matrix< DDRMat > mNormal;

            moris::Cell< real > mMinCoord;
            moris::Cell< real > mMaxCoord;

            real mHesse;

            bool mFlag = false;

            real mIntersectionTolerance;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Facet(
                    moris_index                                     aIndex,
                    moris::Cell< std::shared_ptr< Facet_Vertex > >& aVertices,
                    uint                                            aDimension,
                    real aIntersectionTolerance = 1e-8 );

            //-------------------------------------------------------------------------------

            ~Facet(){};

            //-------------------------------------------------------------------------------

            virtual void
            update_data() = 0;

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief Returns the distance of a Point to the facet.
             * @param[in]  aPoint  point to be considered
             *
             */
            virtual real
            get_distance_to_point(
                    const Matrix< DDRMat >& aPoint ) = 0;

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
             * Resets all the transformed flags for each vertex of the facet
            */
            void
            reset_vertex_transformed_flags();

            //-------------------------------------------------------------------------------

            /**
             * Resets the object back to its attitude when it was constructed,
             * removing any rotation, scaling, or translation that was applied
             *
             */
            void
            reset_coordinates();

            //-------------------------------------------------------------------------------
            // SDF functions
            //-------------------------------------------------------------------------------

            /**
             * @brief returns the minimum coordinate of the triangle
             *
             * @param[in] aDimension   0: x-coordinate
             *                         1: y-coordinate
             *                         2: z-coordinate
             */
            real
            get_min_coord( uint aDimension ) const
            {
                return mMinCoord( aDimension );
            }

            //-------------------------------------------------------------------------------

            /**
             * @brief returns the maximum coordinate of the triangle
             *
             * @param[in] aDimension   0: x-coordinate
             *                         1: y-coordinate
             *                         2: z-coordinate
             *
             */
            real
            get_max_coord( uint aDimension ) const
            {
                return mMaxCoord( aDimension );
            }

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief returns the center of the facet
             *
             */
            const Matrix< DDRMat >&
            get_center() const
            {
                return mCenter;
            }

            /**
             * Get the normal vector to the facet
             *
             * @return Normal vector, mNormal set by child class
             */
            const Matrix< DDRMat >&
            get_normal() const
            {
                return mNormal;
            }

            /**
             * Get the Hesse normal distance. Dot product of the center vector and normal vector
             *
             * @return const real
             */
            real
            get_hesse() const
            {
                return mHesse;
            }

            //-------------------------------------------------------------------------------

            virtual bool
            check_edge(
                    const uint              aEdge,
                    const uint              aAxis,
                    const Matrix< DDRMat >& aPoint ) = 0;

            //-------------------------------------------------------------------------------


            /**
             * Intersects the Line in a coordinate axis direction from aPoint
             *
             *        g(i) = aPoint(i) + tParam*kronecker(i,aAxis)
             *
             * @param[in] aPoint the point in which to draw the coordinate axis line from
             * @param[in] aAxis the direction in which the intersection location will be computed\
             * @param[in] aCoordinate return value. The global location of the piercing point of the line
             * @param[in] aError return value. True if the line is close to parallel with the coordinate axis
             *
             */
            void
            intersect_with_coordinate_axis(
                    const Matrix< DDRMat >& aPoint,
                    const uint              aAxis,
                    real&                   aCoordinate,
                    bool&                   aError );

            //-------------------------------------------------------------------------------
            // MTK API functions
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

            uint
            get_number_of_vertices() const;

            //-------------------------------------------------------------------------------

            moris_id
            get_owner() const
            {
                return 0;
            }

            //-------------------------------------------------------------------------------

            moris::Cell< mtk::Vertex* >
            get_vertex_pointers() const;

            //-------------------------------------------------------------------------------

            void
            remove_vertex_pointer( moris_index aIndex );

            //-------------------------------------------------------------------------------

            Matrix< IdMat >
            get_vertex_ids() const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_vertex_inds() const;

            //-------------------------------------------------------------------------------

            /**
             * Gets the global coordinates of the vertices of the facet.
             *
             * @return Global coordinates. Each row is a vertex, each column is a dimension
             */
            Matrix< DDRMat >
            get_vertex_coords() const;

            //-------------------------------------------------------------------------------

            mtk::Geometry_Type
            get_geometry_type() const;

            //-------------------------------------------------------------------------------

            mtk::Interpolation_Order
            get_interpolation_order() const
            {
                return mtk::Interpolation_Order::LINEAR;
            }

            //-------------------------------------------------------------------------------
            // SDF Functions
            //-------------------------------------------------------------------------------

            /**
             * Sets flag marking this facet for further raycast consideration
             *
             */
            void
            flag()
            {
                mFlag = true;
            }

            //-------------------------------------------------------------------------------

            /**
             * sets flag to false when intersection computation is completed or this facet is no longer used for raycast
             *
             */
            void
            unflag()
            {
                mFlag = false;
            }

            //-------------------------------------------------------------------------------

            bool
            is_flagged() const
            {
                return mFlag;
            }

            //-------------------------------------------------------------------------------

            bool
            operator==( const Facet& aRHS ) const;

            //-------------------------------------------------------------------------------

          protected:
            /**
             * Computes the center of the facet and stores it in mCenter
             *
             */
            void
            compute_center();

            /**
             * computes the minimum and maximum coordinates for each coordinate direction
             * and stores them in mMinCoord and mMaxCoord respectively.
             *
             */
            void
            compute_min_and_max_coordinates();

            //-------------------------------------------------------------------------------

        };    // class Facet

        //-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */
