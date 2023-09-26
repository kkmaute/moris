/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Line.hpp
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

namespace moris::Sdf
{
    //-------------------------------------------------------------------------------

    /**
    * the line class for the sdf generator
    */
    class Line : public Facet
    {
        // index of this line
        const moris_index mIndex;

        // cells with vertex pointers
        moris::Cell< Facet_Vertex* > mVertices;
        moris::Cell< Line* >        mNeighbors;

        // container for node coordinates
        Matrix< DDRMat > mNodeCoords;

        // container for node indices
        Matrix< IndexMat > mNodeIndices;

        // container for center
        Matrix< DDRMat > mCenter;

        // container for normal
        Matrix< DDRMat > mNormal;

        Matrix< DDRMat > mPredictY;
        Matrix< DDRMat > mPredictYRA;
        Matrix< DDRMat > mPredictYRB;

        Matrix< DDRMat > mMinCoord;
        Matrix< DDRMat > mMaxCoord;

        bool mFlag = false;
            //-------------------------------------------------------------------------------

          public:
          
            //-------------------------------------------------------------------------------

            Line(
                    moris_index                      aIndex,
                    moris::Cell< Facet_Vertex* >& aVertices );

            //-------------------------------------------------------------------------------

            ~Line(){};

            //-------------------------------------------------------------------------------

            void
            update_data() override;

            //-------------------------------------------------------------------------------
            // SDF functions
            //-------------------------------------------------------------------------------

            /**
             *
             * @brief returns the normal vector of the line
             *
             */
            const Matrix< DDRMat >&
            get_normal() override
            {
                return mNormal;
            }

            //-------------------------------------------------------------------------------

            /**
             * @brief intersects the line
             *
             *        g(i) = aPoint(i) + tParam*kronecker(i,aAxis)
             *
             *        with the facet and returns the coordinate of the axis
             *
             * @param[in] aPoint
             * @param[in] aAxis
             *
             */
            void
            intersect_with_coordinate_axis(
                    const Matrix< DDRMat >& aPoint,
                    const uint               aAxis,
                    real&                    aCoordinate,
                    bool&                    aError );

            //-------------------------------------------------------------------------------

            bool
            check_edge(
                    const uint               aEdge,
                    const uint               aAxis,
                    const Matrix< DDRMat >& aPoint ) override;

            //-------------------------------------------------------------------------------

            /**
             * @brief Projects a point on the local 2D coordinate system.
             *        The third entry contains a signed point-plane distance.
             * @param[in]  aPoint  point to project
             *
             */
            moris::Matrix< DDRMat >
            project_point_to_local_cartesian(
                    const moris::Matrix< DDRMat >& aPoint );

            //-------------------------------------------------------------------------------

            /**
             * @brief Returns the distance of a Point to an edge.
             *         Point must be in local cartesian coordinates.
             * @param[in]  aLocalPoint  point to be considered
             * @param[in]  aEdge        edge to be considered
             */
            real
            distance_point_to_edge_in_local_cartesian(
                    const Matrix< DDRMat >& aLocalPoint,
                    const uint               aEdge );

            //-------------------------------------------------------------------------------

            /**
             *
             * @brief Returns the distance of a Point to the line.
             * @param[in]  aPoint  point to be considered
             *
             */
            real
            get_distance_to_point(
                    const Matrix< DDRMat >& aPoint ) override;

            //-------------------------------------------------------------------------------
            // SDF Functions
            //-------------------------------------------------------------------------------
            

        }; // class Line

        //-------------------------------------------------------------------------------

    } /* namespace moris::sdf */
