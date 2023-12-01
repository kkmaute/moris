/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Triangle.cpp
 *
 */

#include "cl_SDF_Triangle.hpp"

#include "assert.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "SDF_Tools.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Triangle::Triangle(
                moris_index                      aIndex,
                moris::Cell< Facet_Vertex* >& aVertices )
                : Facet( aIndex, aVertices, 3 )
                , mPredictY( 3, 3 )
                , mPredictYRA( 3, 3 )
                , mPredictYRB( 3, 3 )
        {
            this->update_data();
        }

        //-------------------------------------------------------------------------------

        void
        Triangle::update_data()
        {
            // step 1: compute center of the triangle
            this->compute_center();
            
            // step 2: determine the triangle's bounding box
            this->compute_min_and_max_coordinates();

            // help vector
            Matrix< DDRMat > tDirectionOfEdge( 3, 1 );

            // step 3: calculate hesse normal form of plane
            this->calculate_hesse_normal_form( tDirectionOfEdge );

            // step 4: calculate barycentric data
            this->calculate_barycentric_data( tDirectionOfEdge );

            // step 5: calculate helpers for cross prediction
            this->calculate_prediction_helpers();
        }

        //-------------------------------------------------------------------------------

        void
        Triangle::calculate_hesse_normal_form( Matrix< DDRMat >& aDirectionOfEdge )
        {
            // step 2: calculate plane of triangle
            Matrix< DDRMat > tDirection02( 3, 1 );

            // help vectors: direction of sides 1 and 2
            for ( uint iAxis = 0; iAxis < 3; ++iAxis )
            {
                aDirectionOfEdge( iAxis ) = mVertices( 1 )->get_coord( iAxis ) - mVertices( 0 )->get_coord( iAxis );
                tDirection02( iAxis )     = mVertices( 2 )->get_coord( iAxis ) - mVertices( 0 )->get_coord( iAxis );
            }

            // norm of this triangle
            mNormal = cross( aDirectionOfEdge, tDirection02 );

            real tNorm = norm( mNormal );
            for ( uint i = 0; i < 3; ++i )
            {
                mNormal( i ) /= tNorm;
            }

            // Hesse parameter of this triangle
            mHesse = dot( mCenter, mNormal );
        }

        //-------------------------------------------------------------------------------

        void
        Triangle::calculate_barycentric_data( const Matrix< DDRMat >& aDirectionOfEdge )
        {

            // calculate direction orthogonal to plane and edge
            Matrix< DDRMat > tDirectionOrtho = cross( mNormal, aDirectionOfEdge );

            // normalize tDirection10
            real tNorm10 = norm( aDirectionOfEdge );

            // normalize tDirectionOrtho
            real tNormOrtho = norm( tDirectionOrtho );

            // Projection matrix
            for ( uint k = 0; k < 3; ++k )
            {
                mBarycentric.mProjectionMatrix( 0, k ) = aDirectionOfEdge( k ) / tNorm10;
                mBarycentric.mProjectionMatrix( 1, k ) = tDirectionOrtho( k ) / tNormOrtho;
                mBarycentric.mProjectionMatrix( 2, k ) = mNormal( k );
            }

            // node coordinates in triangle plane
            mBarycentric.mLocalNodeCoordsInPlane.fill( 0 );
            for ( uint k = 0; k < 3; ++k )
            {
                for ( uint j = 0; j < 3; ++j )
                {
                    for ( uint i = 0; i < 2; ++i )
                    {
                        mBarycentric.mLocalNodeCoordsInPlane( i, k ) +=
                                mBarycentric.mProjectionMatrix( i, j )
                                * ( mVertices( k )->get_coord( j ) - mCenter( j ) );
                    }
                }
            }

            // enlarge triangle
            mBarycentric.mLocalNodeCoordsInPlane( 0, 0 ) -= gSDFepsilon;
            mBarycentric.mLocalNodeCoordsInPlane( 1, 0 ) -= gSDFepsilon;
            mBarycentric.mLocalNodeCoordsInPlane( 0, 1 ) += gSDFepsilon;
            mBarycentric.mLocalNodeCoordsInPlane( 1, 1 ) -= gSDFepsilon;
            mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ) += gSDFepsilon;

            // twice the area
            mBarycentric.mTwiceArea = ( mBarycentric.mLocalNodeCoordsInPlane( 0, 0 )
                                              - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ) )
                                            * ( mBarycentric.mLocalNodeCoordsInPlane( 1, 1 )
                                                    - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ) )
                                    - ( mBarycentric.mLocalNodeCoordsInPlane( 1, 0 )
                                              - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ) )
                                              * ( mBarycentric.mLocalNodeCoordsInPlane( 0, 1 )
                                                      - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ) );

            
            // warn if the the triangle has a volume close to zero
            if( mBarycentric.mTwiceArea <= 2 * gSDFepsilon )
            {
                MORIS_LOG_WARNING( 
                        "TRI/TET with ID %i is potentially degenerate and has a volume of V = %e. ",
                        this->get_id(), 
                        mBarycentric.mTwiceArea );
                MORIS_LOG_WARNING( "Nodal coordinates = %s", ios::stringify_log( this->get_vertex_coords() ).c_str() );
            }

            // throw an error if an inverted triangle has been produced
            MORIS_ERROR( mBarycentric.mTwiceArea > -1.0 * MORIS_REAL_EPS,
                    "TRI/TET with ID %i is inverted and has a volume of V = %e. "
                    "Nodal coordinates of the element are: %s", 
                    this->get_id(), 
                    mBarycentric.mTwiceArea,
                    ios::stringify_log( this->get_vertex_coords() ).c_str() );

            // store the inverse of the jacobian
            mBarycentric.mInvTwiceArea = 1.0 / mBarycentric.mTwiceArea;

            // Edge directions
            for ( uint k = 0; k < 3; ++k )
            {
                uint i;
                uint j;
                TrianglePermutation( k, i, j );
                for ( uint l = 0; l < 2; ++l )
                {
                    mBarycentric.mLocalEdgeDirectionVectors( l, k ) = mBarycentric.mLocalNodeCoordsInPlane( l, j )
                                                                    - mBarycentric.mLocalNodeCoordsInPlane( l, i );
                }
                real tMagnitude = 0.0;
                for ( uint l = 0; l < 2; ++l )
                {
                    tMagnitude += mBarycentric.mLocalEdgeDirectionVectors( l, k )
                                * mBarycentric.mLocalEdgeDirectionVectors( l, k );
                }
                mBarycentric.mLocalEdgeInverseMagnitudes( k ) = 1.0 / tMagnitude;
            }
        }

        //-------------------------------------------------------------------------------

        void
        Triangle::calculate_prediction_helpers()
        {
            // values for cross prediction
            uint i;
            uint j;
            uint p;
            uint q;
            for ( moris::uint k = 0; k < 3; ++k )
            {
                TrianglePermutation( k, i, j );

                for ( uint r = 0; r < 3; ++r )
                {
                    TrianglePermutation( r, p, q );
                    real tDelta = mVertices( p )->get_coord( i ) - mVertices( q )->get_coord( i );
                    if ( std::abs( tDelta ) < gSDFepsilon )
                    {
                        if ( tDelta < 0 )
                            tDelta = -gSDFepsilon;
                        else
                            tDelta = gSDFepsilon;
                    }

                    mPredictYRA( r, k ) = ( mVertices( p )->get_coord( j ) - mVertices( q )->get_coord( j ) ) / tDelta;

                    mPredictY( r, k ) = mVertices( p )->get_coord( j )
                                      + mPredictYRA( r, k ) * ( mVertices( r )->get_coord( i ) - mVertices( p )->get_coord( i ) );

                    mPredictYRB( r, k ) = mVertices( p )->get_coord( j ) - mVertices( p )->get_coord( i ) * mPredictYRA( r, k );
                }
            }
        }

        //-------------------------------------------------------------------------------
        // SDF functions
        //-------------------------------------------------------------------------------

        bool
        Triangle::check_edge(
                const uint               aEdge,
                const uint               aAxis,
                const Matrix< DDRMat >& aPoint )
        {
            uint tI;
            uint tJ;
            uint tP;
            uint tQ;

            // permutation parameter for axis
            TrianglePermutation( aAxis, tI, tJ );

            // permutation parameter for edge
            TrianglePermutation( aEdge, tP, tQ );

            // R
            real tPredictYR = mPredictYRA( aEdge, aAxis ) * aPoint( tI ) + mPredictYRB( aEdge, aAxis );

            // check if point is within all three projected edges
            return ( ( mPredictY( aEdge, aAxis ) > mVertices( aEdge )->get_coord( tJ ) )
                           && ( tPredictYR + gSDFepsilon > aPoint( tJ ) ) )
                || ( ( mPredictY( aEdge, aAxis ) < mVertices( aEdge)->get_coord( tJ ) )
                        && ( tPredictYR - gSDFepsilon < aPoint( tJ ) ) )
                || ( std::abs( ( mVertices( tP )->get_coord( tJ ) - mVertices( tQ )->get_coord( tJ ) )
                               * ( mVertices( tP )->get_coord( tI ) - aPoint( tI ) ) )
                        < gSDFepsilon );
        }

        //-------------------------------------------------------------------------------

        Matrix< F31RMat >
        Triangle::get_barycentric_from_local_cartesian(
                const Matrix< F31RMat >& aLocalPoint )
        {
            Matrix< F31RMat > aXi( 3, 1 );

            // the first coordinate
            aXi( 0 ) = ( ( mBarycentric.mLocalNodeCoordsInPlane( 0, 1 )
                                 - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ) )
                                       * ( mBarycentric.mLocalNodeCoordsInPlane( 1, 2 )
                                               - aLocalPoint( 1 ) )
                               - ( mBarycentric.mLocalNodeCoordsInPlane( 1, 1 )
                                         - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ) )
                                         * ( mBarycentric.mLocalNodeCoordsInPlane( 0, 2 )
                                                 - aLocalPoint( 0 ) ) )
                     * mBarycentric.mInvTwiceArea;

            // the second coordinate
            aXi( 1 ) = ( ( mBarycentric.mLocalNodeCoordsInPlane( 1, 0 )
                                 - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ) )
                                       * ( mBarycentric.mLocalNodeCoordsInPlane( 0, 2 )
                                               - aLocalPoint( 0 ) )
                               - ( mBarycentric.mLocalNodeCoordsInPlane( 0, 0 )
                                         - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ) )
                                         * ( mBarycentric.mLocalNodeCoordsInPlane( 1, 2 )
                                                 - aLocalPoint( 1 ) ) )
                     * mBarycentric.mInvTwiceArea;

            // the third coordinate
            aXi( 2 ) = 1.0 - aXi( 0 ) - aXi( 1 );

            return aXi;
        }

        //-------------------------------------------------------------------------------

        real
        Triangle::distance_point_to_edge_in_local_cartesian(
                const Matrix< F31RMat >& aLocalPoint,
                const uint               aEdge )
        {
            real tParam = 0;
            uint i;
            uint j;

            // permutation parameter of current edge
            TrianglePermutation( aEdge, i, j );

            // calculate projection of point on edge

            // tParam = 0: orthogonal intersects with point i
            // tParam = 1: orthogonal intersects with point j

            for ( uint l = 0; l < 2; ++l )
            {
                tParam += ( aLocalPoint( l ) - mBarycentric.mLocalNodeCoordsInPlane( l, i ) )
                        * mBarycentric.mLocalEdgeDirectionVectors( l, aEdge );
            }
            tParam *= mBarycentric.mLocalEdgeInverseMagnitudes( aEdge );

            Matrix< F31RMat > aDirection( 3, 1 );

            if ( tParam < gSDFepsilon )
            {
                // snap to point i and set tParam = 0.0;
                aDirection( 0 ) = aLocalPoint( 0 ) - mBarycentric.mLocalNodeCoordsInPlane( 0, i );
                aDirection( 1 ) = aLocalPoint( 1 ) - mBarycentric.mLocalNodeCoordsInPlane( 1, i );
            }
            else if ( tParam > 1.0 - gSDFepsilon )
            {
                // snap to point j and set tParam = 1.0;
                aDirection( 0 ) = aLocalPoint( 0 ) - mBarycentric.mLocalNodeCoordsInPlane( 0, j );
                aDirection( 1 ) = aLocalPoint( 1 ) - mBarycentric.mLocalNodeCoordsInPlane( 1, j );
            }
            else
            {
                // find distance in plane
                aDirection( 0 ) = aLocalPoint( 0 ) - mBarycentric.mLocalNodeCoordsInPlane( 0, i )
                                - tParam * mBarycentric.mLocalEdgeDirectionVectors( 0, aEdge );
                aDirection( 1 ) = aLocalPoint( 1 ) - mBarycentric.mLocalNodeCoordsInPlane( 1, i )
                                - tParam * mBarycentric.mLocalEdgeDirectionVectors( 1, aEdge );
            }

            // add third dimension to distance
            aDirection( 2 ) = aLocalPoint( 2 );

            return norm( aDirection );
        }

        //-------------------------------------------------------------------------------

        Matrix< F31RMat >
        Triangle::project_point_to_local_cartesian(
                const moris::Matrix< F31RMat >& aPoint )
        {
            // fixme: times operator does not work with eigen
            // return mBarycentric.mProjectionMatrix * ( aPoint - mCenter ) ;
            Matrix< F31RMat > aOut( 3, 1 );
            aOut.fill( 0 );

            for ( uint k = 0; k < 3; ++k )
            {
                for ( uint i = 0; i < 3; ++i )
                {
                    aOut( k ) += mBarycentric.mProjectionMatrix( k, i )
                               * ( aPoint( i ) - mCenter( i ) );
                }
            }

            return aOut;
        }

        //-------------------------------------------------------------------------------

        real
        Triangle::get_distance_to_point(
                const Matrix< DDRMat >& aPoint )
        {
            // step 1: Transform Point to in-plane coordinates
            Matrix< F31RMat > tLocalPointCoords = this->project_point_to_local_cartesian( aPoint );
            // step 2: calculate barycentric coordinates
            Matrix< F31RMat > tXi = this->get_barycentric_from_local_cartesian( tLocalPointCoords );

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

        //-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */
