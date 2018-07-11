//
// Created by messe on 1/7/18.
//

#ifndef MORIS_CL_GETRIANGLE_HPP
#define MORIS_CL_GETRIANGLE_HPP

#include "assert.hpp"
#include "cl_Mat.hpp" // LNA/src

#include "typedefs.hpp" // COR/src
#include "GeUtilities.hpp"

namespace ge {

/*! The basic Triangle class for the SDF Generator. */
    class SDF_Triangle
    {
        struct BarycentricData
        {
            moris::Mat< moris::real > mLocalEdgeDirectionVectors;
            moris::Mat< moris::real > mLocalEdgeInverseMagnitudes;
            moris::Mat< moris::real > mProjectionMatrix;
            moris::Mat< moris::real > mLocalNodeCoordsInPlane;
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
        moris::Mat< moris::uint> mNodeIDs;
        moris::Mat< moris::real > mNodeCoords;

        BarycentricData mBarycentric;

        moris::Mat< moris::real > mMinCoord;
        moris::Mat< moris::real > mMaxCoord;

        moris::Mat< moris::real > mCenter;
        moris::Mat< moris::real > mNormal;

        moris::real mHesse;

        moris::Mat< moris::real > mPredictY;
        moris::Mat< moris::real > mPredictYRA;
        moris::Mat< moris::real > mPredictYRB;

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
        set_node_coords(const moris::Mat< moris::real >& aAllNodeCoords);

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
        moris::Mat< moris::real >
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
        moris::Mat< moris::real >
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
        moris::Mat< moris::real >
        intersect_with_line(
                const moris::Mat< moris::real >& aPoint,
                const moris::Mat< moris::real >& aDirection);

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
                const moris::Mat< moris::real >& aPoint,
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
                const moris::Mat< moris::real >& aPoint );

// -----------------------------------------------------------------------------

        /**
         * @brief Projects a point on the local 2D coordinate system.
         *        The third entry contains a signed point-plane distance.
         * @param[in]  aPoint  point to project
         *
         */
        moris::Mat< moris::real >
        project_point_to_local_cartesian(
                const moris::Mat< moris::real >& aPoint)
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
        moris::Mat< moris::real >
        get_barycentric_from_local_cartesian(
                const moris::Mat< moris::real >& aLocalPoint);

// -----------------------------------------------------------------------------

        /**
         * @brief Returns the distance of a Point to an edge.
         *         Point must be in local cartesian coordinates.
         * @param[in]  aLocalPoint  point to be considered
         * @param[in]  aEdge        edge to be considered
         */
        moris::real distance_point_to_edge_in_local_cartesian(
                const moris::Mat< moris::real >& aLocalPoint,
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
                const moris::Mat< moris::real >& aPoint);

// -----------------------------------------------------------------------------

    };

}

#endif //MORIS_CL_GETRIANGLE_HPP
