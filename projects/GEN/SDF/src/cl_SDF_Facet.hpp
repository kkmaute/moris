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
#include "fn_eye.hpp"

#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_SDF_Facet_Vertex.hpp"
#include "SDF_Tools.hpp"

#include "GEN/src/PDV/cl_GEN_Intersection_Node_Surface_Mesh.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        /**
         * the facet class for the sdf generator
         */
        class Facet : public mtk::Cell
                , std::enable_shared_from_this< Facet >
        {
          protected:
            // index of this facet
            moris_index mIndex;

            // cells with vertex pointers
            Vector< std::shared_ptr< Facet_Vertex > > mVertices;

            // container for center
            Matrix< DDRMat > mCenter;

            // container for normal
            Matrix< DDRMat > mNormal;

            Vector< real > mMinCoord;
            Vector< real > mMaxCoord;

            real mHesse;

            bool mFlag = false;

            real mIntersectionTolerance;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            Facet(
                    moris_index                                aIndex,
                    Vector< std::shared_ptr< Facet_Vertex > >& aVertices,
                    uint                                       aDimension,
                    real                                       aIntersectionTolerance = 1e-8 );

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

            /**
             * Gets the specified axis coordinate of the specified vertex
             *
             * @param aVertexIndex Vertex to get the coord of
             * @param aAxis Component of the vertex to retrieve
             * @return global coordinate of the vertex
             */
            real
            get_vertex_coord( uint aVertexIndex, uint aAxis );

            //-------------------------------------------------------------------------------
            // GEN functions
            //-------------------------------------------------------------------------------

            /**
             * Computes the sensitivity of the local coordinate of a parent edge with respect to the facet vertices
             *
             * @param aXi Local coordinate of an intersection along a parent edge, NOT the local coordinate of the facet.
             * @param aIntersectionNode Intersection node to compute sensitivity for
             * @return Vector< real > Sensitivities of the local coordinate with respect to the facet vertices. Size <Object dimension x number of vertices>
             */
            virtual Matrix< DDRMat >
            compute_dxi_dvertices( moris::gen::Intersection_Node_Surface_Mesh& aIntersectionNode ) = 0;

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

            Vector< mtk::Vertex* >
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
