/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Triangle.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_TRIANGLE_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_TRIANGLE_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "op_times.hpp"

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_SDF_Triangle_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        /**
         * the triangle class for the sdf generator
         */
        class Triangle : public mtk::Cell
        {
            // index of this triangle
            const moris_index mIndex;

            // cells with vertex pointers
            moris::Cell< Triangle_Vertex* > mVertices;
            moris::Cell< Triangle* >        mNeighbors;

            struct BarycentricData
            {
                Matrix< F33RMat > mLocalEdgeDirectionVectors;
                Matrix< F31RMat > mLocalEdgeInverseMagnitudes;
                Matrix< F33RMat > mProjectionMatrix;
                Matrix< DDRMat >  mLocalNodeCoordsInPlane;
                real              mTwiceArea;
                real              mInvTwiceArea;
                BarycentricData()
                        : mLocalEdgeDirectionVectors( 3, 3 )
                        , mLocalEdgeInverseMagnitudes( 3, 1 )
                        , mProjectionMatrix( 3, 3 )
                        , mLocalNodeCoordsInPlane( 2, 3 ){};
                ~BarycentricData() = default;
            };

            BarycentricData mBarycentric;

            // container for node coordinates
            Matrix< F33RMat > mNodeCoords;

            // container for node indices
            Matrix< IndexMat > mNodeIndices;

            // container for center
            Matrix< F31RMat > mCenter;

            // container for normal
            Matrix< F31RMat > mNormal;

            real mHesse;

            Matrix< F33RMat > mPredictY;
            Matrix< F33RMat > mPredictYRA;
            Matrix< F33RMat > mPredictYRB;

            Matrix< F31RMat > mMinCoord;
            Matrix< F31RMat > mMaxCoord;

            bool mFlag = false;

            //-------------------------------------------------------------------------------

          public:
          
            //-------------------------------------------------------------------------------

            Triangle(
                    moris_index                      aIndex,
                    moris::Cell< Triangle_Vertex* >& aVertices );

            //-------------------------------------------------------------------------------

            ~Triangle(){};

            //-------------------------------------------------------------------------------

            void
            update_data();

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
             * @brief Returns the hesse distance of the plane describing the triangle
             *
             *
             */
            real
            get_hesse() const
            {
                return mHesse;
            }

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief returns the normal vector of the triangle
             *
             */
            const Matrix< F31RMat >&
            get_normal() const
            {
                return mNormal;
            }

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief returns the center of the triangle
             *
             */
            const Matrix< F31RMat >&
            get_center() const
            {
                return mCenter;
            }

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief returns the area of the triangle
             *
             */
            real
            get_area()
                    const
            {
                return 0.5 * mBarycentric.mTwiceArea;
            }

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
            void
            intersect_with_coordinate_axis(
                    const Matrix< F31RMat >& aPoint,
                    const uint               aAxis,
                    real&                    aCoordinate,
                    bool&                    aError );

            //-------------------------------------------------------------------------------

            bool
            check_edge(
                    const uint               aEdge,
                    const uint               aAxis,
                    const Matrix< F31RMat >& aPoint );

            //-------------------------------------------------------------------------------

            /**
             * @brief Projects a point on the local 2D coordinate system.
             *        The third entry contains a signed point-plane distance.
             * @param[in]  aPoint  point to project
             *
             */
            moris::Matrix< F31RMat >
            project_point_to_local_cartesian(
                    const moris::Matrix< F31RMat >& aPoint );

            //-------------------------------------------------------------------------------

            /**
             * @brief Returns the barycentric coordinates for a point.
             *       Point must be in local cartesian coordinates.
             *
             * @param[in]  aLocalPoint  point to project
             *
             */
            Matrix< F31RMat >
            get_barycentric_from_local_cartesian(
                    const Matrix< F31RMat >& aLocalPoint );

            //-------------------------------------------------------------------------------

            /**
             * @brief Returns the distance of a Point to an edge.
             *         Point must be in local cartesian coordinates.
             * @param[in]  aLocalPoint  point to be considered
             * @param[in]  aEdge        edge to be considered
             */
            real
            distance_point_to_edge_in_local_cartesian(
                    const Matrix< F31RMat >& aLocalPoint,
                    const uint               aEdge );

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief Returns the distance of a Point to the triangle.
             * @param[in]  aPoint  point to be considered
             *
             */
            real
            get_distance_to_point(
                    const Matrix< F31RMat >& aPoint );

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
            get_number_of_vertices() const
            {
                return 3;
            }

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
            get_geometry_type() const
            {
                return mtk::Geometry_Type::TRI;
            }

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

          private:

            //-------------------------------------------------------------------------------

            void
            copy_node_coords_and_inds( moris::Cell< Triangle_Vertex* >& aVertices );

            //-------------------------------------------------------------------------------

            void
            calculate_hesse_normal_form( Matrix< F31RMat >& aDirectionOfEdge );

            //-------------------------------------------------------------------------------

            void
            calculate_barycectric_data( const Matrix< F31RMat >& aDirectionOfEdge );

            //-------------------------------------------------------------------------------

            void
            calculate_prediction_helpers();

            //-------------------------------------------------------------------------------

        }; // class Triangle

        //-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_TRIANGLE_HPP_ */
