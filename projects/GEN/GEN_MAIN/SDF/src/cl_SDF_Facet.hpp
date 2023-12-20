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

#include "moris_typedefs.hpp"
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
        class Facet : public mtk::Cell
        {
          protected:
            // index of this triangle
            const moris_index mIndex;

            // cells with vertex pointers
            moris::Cell< Facet_Vertex* > mVertices;

            // container for node coordinates
            Matrix< DDRMat > mNodeCoords;

            // container for node indices
            Matrix< IndexMat > mNodeIndices;

            // container for center
            Matrix< DDRMat > mCenter;

            // container for normal
            Matrix< DDRMat > mNormal;

            moris::Cell< real > mMinCoord;
            moris::Cell< real > mMaxCoord;

            real mHesse;

            bool mFlag = false;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Facet(
                    moris_index                   aIndex,
                    moris::Cell< Facet_Vertex* >& aVertices,
                    uint                          aDimension );

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
            // SDF functions
            //-------------------------------------------------------------------------------

            /**
             * @brief returns the minimum coordinate of the triangle
             *
             * @param[in] aDimension   0: x-coordinate
             *                         1: y-coordinate
             *                         2: z-coordinate
             *
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
            get_normal()
            {
                return mNormal;
            }

            /**
             * Get the Hesse normal distance. Dot product of the center vector and normal vector
             *
             * @return const real
             */
            const real
            get_hesse()
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

            void
            flag()
            {
                mFlag = true;
            }

            //-------------------------------------------------------------------------------

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

          protected:
            /**
             * Sets mNodeCoords, mMinCoord, mMaxCoord, and mCenter from the coordinates of mVertices
             * Each column of mNodeCoords corresponds to a vertex, and each row a dimension
             *
             * @param aVertices Cell of Facet_Vertex's that each contain their global coordinates
             * @param aDimension Number of dimensions for the problem. 2 for line facets and 3 for triangle facets.
             */
            void
            copy_node_coords_and_inds( moris::Cell< Facet_Vertex* >& aVertices, uint aDimension );

            //-------------------------------------------------------------------------------

        };    // class Facet

        //-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */
