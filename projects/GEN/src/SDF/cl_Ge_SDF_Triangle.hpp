//
// Created by messe on 1/7/18.
//

#ifndef MORIS_CL_GETRIANGLE_HPP
#define MORIS_CL_GETRIANGLE_HPP

#include "assert.hpp"
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"

#include "typedefs.hpp" // COR/src
#include "GeUtilities.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"

namespace ge {

/*! The basic Triangle class for the SDF Generator. */
    class SDF_Triangle
    {
        struct BarycentricData
        {
            moris::Matrix< moris::DDRMat > mLocalEdgeDirectionVectors;
            moris::Matrix< moris::DDRMat > mLocalEdgeInverseMagnitudes;
            moris::Matrix< moris::DDRMat > mProjectionMatrix;
            moris::Matrix< moris::DDRMat > mLocalNodeCoordsInPlane;
            moris::real mTwiceArea;
            moris::real mInvTwiceArea;
            BarycentricData():mLocalEdgeDirectionVectors( 3, 3 ),
                              mLocalEdgeInverseMagnitudes( 3, 1 ),
                              mProjectionMatrix( 3, 3 ),
                              mLocalNodeCoordsInPlane( 2 , 3 )
            {
            };
            ~BarycentricData() = default;
        };
        moris::Matrix< moris::DDUMat > mNodeIDs;
        moris::Matrix< moris::DDRMat > mNodeCoords;

        BarycentricData mBarycentric;

        moris::Matrix< moris::DDRMat > mMinCoord;
        moris::Matrix< moris::DDRMat > mMaxCoord;

        moris::Matrix< moris::DDRMat > mCenter;
        moris::Matrix< moris::DDRMat > mNormal;

        moris::real mHesse;

        moris::Matrix< moris::DDRMat > mPredictY;
        moris::Matrix< moris::DDRMat > mPredictYRA;
        moris::Matrix< moris::DDRMat > mPredictYRB;

// =============================================================================
    public:
// =============================================================================
        /**
        * @brief default constructor
        *
        */
        SDF_Triangle(
                const moris::uint aNodeID0,
                const moris::uint aNodeID1,
                const moris::uint aNodeID2):
                mNodeIDs( 3 , 1 ),
                mNodeCoords( 3 , 3 ),
                mBarycentric(),
                mMinCoord( 3 , 1 ),
                mMaxCoord( 3, 1 ),
                mCenter( 3 , 1 ),
                mNormal( 3, 1 ),
                mPredictY( 3 , 3 ),
                mPredictYRA( 3, 3 ),
                mPredictYRB( 3, 3 )
        {

            mNodeIDs( 0 ) = aNodeID0;
            mNodeIDs( 1 ) = aNodeID1;
            mNodeIDs( 2 ) = aNodeID2;

        }

// -----------------------------------------------------------------------------

        ~SDF_Triangle() = default;

// -----------------------------------------------------------------------------
        /**
         * @brief Initializes the node coorditanes and various internal parameters.
         *
         *
         */
        void
        set_node_coords(const moris::Matrix< moris::DDRMat >& aAllNodeCoords);

// -----------------------------------------------------------------------------

        /**
         * @brief Returns the hesse distance of the plane describing the triangle
         *
         *
         */
        moris::real
        get_hesse() const
        {
            return mHesse;
        }

// -----------------------------------------------------------------------------

        /**
         * @brief returns the minimum coordinate of the triangle
         *
         * @param[in] aDimension   0: x-coordinate
         *                         1: y-coordinate
         *                         2: z-coordinate
         *
         */
        moris::real
        get_min_coord(moris::uint aDimension) const
        {
            return mMinCoord(aDimension);
        }

// -----------------------------------------------------------------------------

        /**
         *
         * @brief returns the center of the triangle
         *
         */
        moris::Matrix< moris::DDRMat >
        get_center() const
        {
            return mCenter;
        }

// -----------------------------------------------------------------------------

        /**
         *
         * @brief returns the normal vector of the triangle
         *
         */
        moris::Matrix< moris::DDRMat >
        get_normal() const
        {
            return mNormal;
        }

// -----------------------------------------------------------------------------

        /**
         *
         * @brief returns the area of the triangle
         *
         */
        moris::real get_area()
        const
        {
            return 0.5*mBarycentric.mTwiceArea;
        }

// -----------------------------------------------------------------------------

        /**
         * @brief returns the maximum coordinate of the triangle
         *
         * @param[in] aDimension   0: x-coordinate
         *                         1: y-coordinate
         *                         2: z-coordinate
         *
         */
        moris::real
        get_max_coord( moris::uint aDimension ) const
        {
            return mMaxCoord( aDimension );
        }

// -----------------------------------------------------------------------------

        /**
         * @brief intersects the line
         *
         *        g = aPoint + tParam*aDirection
         *
         *        with the triangle abd returns the crossing point.
         *
         * @param[in] aPoint
         * @param[in] aDirection
         *
         */
        moris::Matrix< moris::DDRMat >
        intersect_with_line(
                const moris::Matrix< moris::DDRMat >& aPoint,
                const moris::Matrix< moris::DDRMat >& aDirection);

// -----------------------------------------------------------------------------

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
                const moris::Matrix< moris::DDRMat >& aPoint,
                const moris::uint aAxis,
                moris::real & aCoordinate,
                bool        & aError )
        {

            if (ge::abs(mNormal(aAxis)) < MORIS_GE_EPSILON )
            {
                aCoordinate = 0;
                aError = true;
            }
            else
            {
                aCoordinate = aPoint( aAxis )+(mHesse - ge::dot( mNormal, aPoint )) / mNormal( aAxis );
                aError = false;
            }
            //MORIS_ASSERT(,
           //              "Division by zero in ge::Triangle::intersect_with_coordinate_axis");

            //return aPoint( aAxis )+(mHesse - ge::dot( mNormal, aPoint )) / mNormal( aAxis );
        }

// -----------------------------------------------------------------------------

        /**
         * @brief checks if a Ray originating from aPoint is between point and opposing
         *        edge of triangle
         *        opposing point
         *
         * @param[in] aEdge    edge to be considered
         * @param[in] aAxis    axis ray is parallel to:
         *                      0: x
         *                      1: y
         *                      2: z
         * @param[in] aPoint   point to be considered
         *
         */
        moris::bool_t
        check_edge(
                const moris::uint aEdge,
                const moris::uint aAxis,
                const moris::Matrix< moris::DDRMat >& aPoint );

// -----------------------------------------------------------------------------

        /**
         * @brief Projects a point on the local 2D coordinate system.
         *        The third entry contains a signed point-plane distance.
         * @param[in]  aPoint  point to project
         *
         */
        moris::Matrix< moris::DDRMat >
        project_point_to_local_cartesian(
                const moris::Matrix< moris::DDRMat >& aPoint)
        {
            return mBarycentric.mProjectionMatrix * ( aPoint - mCenter );
        }

// -----------------------------------------------------------------------------
        /**
        * @brief Returns the barycentric coordinates for a point.
        *       Point must be in local cartesian coordinates.
        *
        * @param[in]  aLocalPoint  point to project
        *
        */
        moris::Matrix< moris::DDRMat >
        get_barycentric_from_local_cartesian(
                const moris::Matrix< moris::DDRMat >& aLocalPoint);

// -----------------------------------------------------------------------------

        /**
         * @brief Returns the distance of a Point to an edge.
         *         Point must be in local cartesian coordinates.
         * @param[in]  aLocalPoint  point to be considered
         * @param[in]  aEdge        edge to be considered
         */
        moris::real distance_point_to_edge_in_local_cartesian(
                const moris::Matrix< moris::DDRMat >& aLocalPoint,
                const moris::uint aEdge);

// -----------------------------------------------------------------------------

        /**
         *
         * @brief Returns the distance of a Point to the triangle.
         * @param[in]  aPoint  point to be considered
         *
         */
        moris::real
        get_distance_to_point(
                const moris::Matrix< moris::DDRMat >& aPoint);

// -----------------------------------------------------------------------------

    };

}

#endif //MORIS_CL_GETRIANGLE_HPP
