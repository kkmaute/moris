/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Triangle.hpp
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
        class Facet : public mtk::Cell
        {
          protected:
            // index of this triangle
            const moris_index mIndex;

            // cells with vertex pointers
            moris::Cell< Facet_Vertex* > mVertices;
            moris::Cell< Facet* >        mNeighbors;

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
             * @brief Returns the distance of a Point to the triangle.
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

            //-------------------------------------------------------------------------------

            virtual bool
            check_edge(
                    const uint              aEdge,
                    const uint              aAxis,
                    const Matrix< DDRMat >& aPoint ) = 0;

            //-------------------------------------------------------------------------------

            /**
             * @brief intersects the line
             *
             *        g(i) = aPoint(i) + tParam*kronecker(i,aAxis)
             *
             *        with the triangle and returns the coordinate of the axis
             *
             * @param[in] aPoint
             * @param[in] aAxis
             *
             */
            virtual void
            intersect_with_coordinate_axis(
                    const Matrix< DDRMat >& aPoint,
                    const uint              aAxis,
                    real&                   aCoordinate,
                    bool&                   aError ) = 0;

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
