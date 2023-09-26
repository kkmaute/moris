/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Line.cpp
 *
 */

#include "cl_SDF_Line.hpp"

#include "assert.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "SDF_Tools.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris::sdf
{
        //-------------------------------------------------------------------------------

        Line::Line(
                moris_index                      aIndex,
                moris::Cell< Facet_Vertex* >& aVertices )
                : Facet( aIndex, aVertices )
                , mIndex( aIndex )
                , mVertices( aVertices )
                , mNodeCoords( 2, 2 )
                , mNodeIndices( 2, 1 )
                , mCenter( 2, 1 )
                , mNormal( 2, 1 )
                , mPredictY( 2, 2 )
                , mPredictYRA( 2, 2 )
                , mPredictYRB( 2, 2 )
                , mMinCoord( 2, 1 )
                , mMaxCoord( 2, 1 )
        {
            this->update_data();
        }

        //-------------------------------------------------------------------------------

        void
        Line::update_data()
        {
            // step 1: copy node coordinates and determine center
            this->copy_node_coords_and_inds( mVertices, 2 );

            // help vector
            Matrix< DDRMat > tDirectionOfEdge( 3, 1 );

            // step 2: calculate hesse normal form of plane
            this->calculate_hesse_normal_form( tDirectionOfEdge );

            // step 3: calculate barycentric data
            this->calculate_barycentric_data( tDirectionOfEdge );

            // step 4: calculate helpers for cross prediction
            this->calculate_prediction_helpers();
        }

        //-------------------------------------------------------------------------------
        // SDF functions
        //-------------------------------------------------------------------------------

        void
        Line::intersect_with_coordinate_axis(
                const Matrix< DDRMat >& aPoint,
                const uint               aAxis,
                real&                    aCoordinate,
                bool&                    aError )
        {
            if ( std::abs( mNormal( aAxis ) ) < gSDFepsilon )
            {
                aCoordinate = 0;
                aError      = true;
            }
            else
            {
                aCoordinate = aPoint( aAxis ) + ( mHesse - dot( mNormal, aPoint ) ) / mNormal( aAxis );
                aError      = false;
            }
        }

        //-------------------------------------------------------------------------------

        bool
        Line::check_edge(
                const uint               aEdge,
                const uint               aAxis,
                const Matrix< DDRMat >& aPoint )
        {
            uint tI;
            uint tJ;
            uint tP;
            uint tQ;

            // permutation parameter for axis
            LinePermutation( aAxis, tI, tJ );

            // permutation parameter for edge
            LinePermutation( aEdge, tP, tQ );

            // R
            real tPredictYR = mPredictYRA( aEdge, aAxis ) * aPoint( tI ) + mPredictYRB( aEdge, aAxis );

            // check if point is within all three projected edges
            return ( ( mPredictY( aEdge, aAxis ) > mNodeCoords( tJ, aEdge ) )
                           && ( tPredictYR + gSDFepsilon > aPoint( tJ ) ) )
                || ( ( mPredictY( aEdge, aAxis ) < mNodeCoords( tJ, aEdge ) )
                        && ( tPredictYR - gSDFepsilon < aPoint( tJ ) ) )
                || ( std::abs( ( mNodeCoords( tJ, tP ) - mNodeCoords( tJ, tQ ) )
                               * ( mNodeCoords( tI, tP ) - aPoint( tI ) ) )
                        < gSDFepsilon );
        }

        //-------------------------------------------------------------------------------

        real
        Line::get_distance_to_point(
                const Matrix< DDRMat >& aPoint )
        {
            // step 1: Transform Point to in-plane coordinates
            Matrix< DDRMat > tLocalPointCoords = this->project_point_to_local_cartesian( aPoint );
            // step 2: calculate barycentric coordinates
            Matrix< DDRMat > tXi = this->get_barycentric_from_local_cartesian( tLocalPointCoords );

            // step 3: check if we are inside the triangle
            if ( ( tXi( 0 ) >= -gSDFepsilon )
                    && ( tXi( 1 ) >= -gSDFepsilon )
                    && ( tXi( 2 ) >= -gSDFepsilon ) )
            {
                // the absolute value of the local z-coordinate is the distance
                return std::abs( tLocalPointCoords( 2 ) );
            }
            else
            {
                if ( tXi( 0 ) > 0 )
                {
                    // this rules out edge 0
                    real tDist1 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords, 1 );
                    real tDist2 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords, 2 );
                    return std::min( tDist1, tDist2 );
                }
                else if ( tXi( 1 ) > 0 )
                {
                    // this rules out edge 1
                    real tDist0 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords, 0 );
                    real tDist2 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords, 2 );
                    return std::min( tDist0, tDist2 );
                }
                else
                {
                    // edge 2 must be the one to rule out
                    real tDist0 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords, 0 );
                    real tDist1 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords, 1 );
                    return std::min( tDist0, tDist1 );
                }
            }
        }

} /* namespace moris::sdf */
