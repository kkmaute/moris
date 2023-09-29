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
#include "cl_SDF_Facet.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        /**
         * the triangle class for the sdf generator
         */
        class Triangle : public Facet
        {
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
                    moris_index                   aIndex,
                    moris::Cell< Facet_Vertex* >& aVertices );

            //-------------------------------------------------------------------------------

            ~Triangle(){};

            //-------------------------------------------------------------------------------

            void
            update_data() override;

            //-------------------------------------------------------------------------------
            // SDF functions
            //-------------------------------------------------------------------------------

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

            bool
            check_edge(
                    const uint              aEdge,
                    const uint              aAxis,
                    const Matrix< DDRMat >& aPoint ) override;

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
                    const Matrix< DDRMat >& aPoint ) override;

            //-------------------------------------------------------------------------------
            // SDF Functions
            //-------------------------------------------------------------------------------

          private:
            /**
             * Determines the hesse distance and the normal vector of the triangle
             * 
             * @param aDirectionOfEdge Helper vector for other functions. Direction vector pointing from the first vertex to the second vertex
             * Set in this function, not used in this function.
             */
            void
            calculate_hesse_normal_form( Matrix< DDRMat >& aDirectionOfEdge );

            //-------------------------------------------------------------------------------

            /**
             * Computes all of the values stored in the mBarycentric stuct
             * 
             * @param aDirectionOfEdge Direction vector pointing from the first vertex to the second vertex
             */
            void
            calculate_barycentric_data( const Matrix< DDRMat >& aDirectionOfEdge );

            //-------------------------------------------------------------------------------

            void
            calculate_prediction_helpers();

            //-------------------------------------------------------------------------------

        };    // class Triangle

        //-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */
