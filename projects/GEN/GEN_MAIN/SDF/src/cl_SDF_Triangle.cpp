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
//#include "fn_print.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Triangle::Triangle(
                moris_index aIndex,
                moris::Cell< Triangle_Vertex * > & aVertices  ) :
                        mIndex( aIndex ),
                        mVertices( aVertices ),
                        mNodeCoords( 3, 3 ),
                        mNodeIndices( 3, 1 ),
                        mCenter( 3, 1 ),
                        mNormal( 3, 1 ),
                        mPredictY( 3, 3 ),
                        mPredictYRA( 3, 3 ),
                        mPredictYRB( 3, 3 ),
                        mMinCoord( 3, 1 ),
                        mMaxCoord( 3, 1 )
        {
            this->update_data();
        }

//-------------------------------------------------------------------------------

        void
        Triangle::update_data()
        {
            // step 1: copy node coordinates and determine center
            this->copy_node_coords_and_inds( mVertices );

            // help vector
            Matrix< F31RMat > tDirectionOfEdge( 3, 1 );

            // step 2: calculate hesse normal form of plane
            this->calculate_hesse_normal_form( tDirectionOfEdge );

            // step 3: calculate barycentric data
            this->calculate_barycectric_data( tDirectionOfEdge );

            // step 4: calculate helpers for cross prediction
            this->calculate_prediction_helpers();
        }

//-------------------------------------------------------------------------------

        void
        Triangle::copy_node_coords_and_inds(
                moris::Cell< Triangle_Vertex * > & aVertices  )
        {
            // make sure that the length is correct
            MORIS_ASSERT( aVertices.size() >= 3,
                    "Triangle() needs at least three vertices as input");

            // step 1: copy node coordinates into member variables
            //         and calculate center

            // reset center
            mCenter.fill( 0 );

            // loop over all nodes
            for( uint k=0; k<3; ++k )
            {
                // get vertex coordinates
                auto tNodeCoords = aVertices( k )->get_coords();

                // copy coordinates into member matrix
                for( uint i=0; i<3; ++i )
                {
                    mNodeCoords( i, k ) = tNodeCoords( i );
                    mCenter( i ) += tNodeCoords( i );
                }

                // remember node indices
                mNodeIndices( k ) = aVertices( k )->get_index();
            }

            // identify minimum and maximum coordinate
            for (moris::uint i = 0; i < 3; ++i)
            {
                mMinCoord( i ) = min(mNodeCoords( i, 0 ),
                                     mNodeCoords( i, 1 ),
                                     mNodeCoords( i, 2 ) );

                mMaxCoord( i ) = max(mNodeCoords( i, 0 ),
                                     mNodeCoords( i, 1 ),
                                     mNodeCoords( i, 2 ) );
            }

            // divide center by three
            for( uint i=0; i<3; ++i )
            {
                mCenter( i ) /= 3.0;
            }
        }

//-------------------------------------------------------------------------------

        void
        Triangle::calculate_hesse_normal_form( Matrix< F31RMat > & aDirectionOfEdge )
        {
            // step 2: calculate plane of triangle
            Matrix< F31RMat > tDirection02( 3, 1 );

            // help vectors: direction of sides 1 and 2
            for ( uint i = 0; i < 3; ++i) {
                aDirectionOfEdge( i ) = mNodeCoords( i, 1 )-mNodeCoords( i, 0 );
                tDirection02( i ) = mNodeCoords( i, 2 )-mNodeCoords( i, 0 );
            }

            // norm of this triangle

            mNormal = cross( aDirectionOfEdge, tDirection02 );

            real tNorm = norm( mNormal );
            for ( uint i = 0; i < 3; ++i)
            {
                mNormal( i ) /= tNorm;
            }

            // Hesse parameter of this triangle
            mHesse = dot( mCenter, mNormal );

        }

//-------------------------------------------------------------------------------

        void
        Triangle::calculate_barycectric_data( const Matrix< F31RMat > & aDirectionOfEdge )
        {

            // calculate direction orthogonal to plane and edge
            Matrix< F31RMat > tDirectionOrtho = cross( mNormal, aDirectionOfEdge );

            // normalize tDirection10
            real tNorm10 = norm( aDirectionOfEdge );

            // normalize tDirectionOrtho
            real tNormOrtho = norm( tDirectionOrtho );

            // Projection matrix
            for ( uint k = 0; k < 3; ++k) {
                mBarycentric.mProjectionMatrix( 0, k ) = aDirectionOfEdge( k )/tNorm10;
                mBarycentric.mProjectionMatrix( 1, k ) = tDirectionOrtho( k )/tNormOrtho;
                mBarycentric.mProjectionMatrix( 2, k ) = mNormal( k );
            }

            // node coordinates in triangle plane
            mBarycentric.mLocalNodeCoordsInPlane.fill( 0 );
            for ( uint k = 0; k < 3; ++k ) {
                for ( uint j = 0; j < 3; ++j ) {
                    for ( uint i = 0; i < 2; ++i ) {
                        mBarycentric.mLocalNodeCoordsInPlane( i, k ) +=
                                mBarycentric.mProjectionMatrix( i, j )
                                *(mNodeCoords( j, k )-mCenter( j ));
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
            mBarycentric.mTwiceArea =  (mBarycentric.mLocalNodeCoordsInPlane( 0, 0 )
                    - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ))
                    *(mBarycentric.mLocalNodeCoordsInPlane( 1, 1 )
                    - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ))
                    -(mBarycentric.mLocalNodeCoordsInPlane( 1, 0 )
                    - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ))
                    *(mBarycentric.mLocalNodeCoordsInPlane( 0, 1 )
                    - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ));

            MORIS_ASSERT(mBarycentric.mTwiceArea > 2*gSDFepsilon,
                    "A degenerated triangle was found.");

            mBarycentric.mInvTwiceArea = 1.0/mBarycentric.mTwiceArea;

            // Edge directions
            for ( uint k = 0; k < 3; ++k)
            {
                uint i;
                uint j;
                TrianglePermutation(k,i,j);
                for ( uint l = 0; l < 2; ++l )
                {
                    mBarycentric.mLocalEdgeDirectionVectors( l, k ) = mBarycentric.mLocalNodeCoordsInPlane( l, j )
                            -mBarycentric.mLocalNodeCoordsInPlane( l, i );
                }
                real tMagnitude = 0.0;
                for ( uint l = 0; l < 2; ++l)
                {
                    tMagnitude +=  mBarycentric.mLocalEdgeDirectionVectors( l, k )
                                               *mBarycentric.mLocalEdgeDirectionVectors( l, k );
                }
                mBarycentric.mLocalEdgeInverseMagnitudes( k ) = 1.0/tMagnitude;
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
                TrianglePermutation(k, i, j);

                for ( uint r = 0; r < 3; ++r )
                {
                    TrianglePermutation(r, p, q);
                    real tDelta = mNodeCoords(i, p) - mNodeCoords(i, q);
                    if ( std::abs( tDelta ) < gSDFepsilon )
                    {
                        if (tDelta < 0)
                            tDelta = -gSDFepsilon;
                        else
                            tDelta = gSDFepsilon;
                    }

                    mPredictYRA( r, k ) =   ( mNodeCoords( j, p ) - mNodeCoords( j, q ) )/tDelta;

                    mPredictY( r, k )   =   mNodeCoords( j, p )
                                     + mPredictYRA( r, k ) * (mNodeCoords( i, r ) - mNodeCoords( i, p ));

                    mPredictYRB( r, k ) = mNodeCoords( j, p ) - mNodeCoords( i, p ) * mPredictYRA( r, k );
                }
            }
        }

//-------------------------------------------------------------------------------
// MTK Interface
//-------------------------------------------------------------------------------

        Cell< mtk::Vertex * >
        Triangle::get_vertex_pointers() const
        {
            moris::Cell< mtk::Vertex *  > aVertices( 3, nullptr );

            for( uint k=0; k<3; ++k )
            {
                aVertices( k ) = mVertices( k );

            }
            return aVertices;
        }

//-------------------------------------------------------------------------------

        Matrix< IdMat >
        Triangle::get_vertex_ids() const
        {

            Matrix< IdMat > aIDs( 3, 1 );
            for( uint k = 0; k<3; ++k )
            {
                aIDs( k ) = mVertices( k )->get_id();
            }

            return aIDs;
        }

//-------------------------------------------------------------------------------

        Matrix< IndexMat >
        Triangle::get_vertex_inds() const
        {

            Matrix< IndexMat > aINDs( 3, 1 );

            for( uint k = 0; k<3; ++k )
            {
                aINDs( k ) = mVertices( k )->get_index();
            }

            return aINDs;
        }

//-------------------------------------------------------------------------------

        Matrix< DDRMat >
        Triangle::get_vertex_coords() const
        {

            Matrix< DDRMat > aCoords( 3, 3 );

            for( uint k=0; k<3; ++k )
            {
                aCoords.set_row( k, mVertices( k )->get_coords() );
            }

            return aCoords;
        }

//-------------------------------------------------------------------------------
// SDF functions
//-------------------------------------------------------------------------------

        void
        Triangle::intersect_with_coordinate_axis(
                            const  Matrix< F31RMat > & aPoint,
                            const uint          aAxis,
                            real              & aCoordinate,
                            bool              & aError )
        {
            if (std::abs( mNormal( aAxis ) ) < gSDFepsilon )
            {
                aCoordinate = 0;
                aError = true;
            }
            else
            {
                aCoordinate = aPoint( aAxis )+( mHesse - dot( mNormal, aPoint )) / mNormal( aAxis );
                aError = false;
            }
        }

//-------------------------------------------------------------------------------

        bool
        Triangle::check_edge(
                const uint                aEdge,
                const uint                aAxis,
                const Matrix< F31RMat > & aPoint )
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

            // check of point is within all three projected edges
            return  (    (mPredictY( aEdge, aAxis ) > mNodeCoords( tJ, aEdge ) )
                      && (tPredictYR  + gSDFepsilon > aPoint( tJ ))) ||
                    (    (mPredictY( aEdge, aAxis ) < mNodeCoords( tJ, aEdge ) )
                      && (tPredictYR - gSDFepsilon  < aPoint( tJ ))) ||
                    ( std::abs((mNodeCoords( tJ, tP )-mNodeCoords( tJ, tQ ))
                            * (mNodeCoords( tI, tP )-aPoint( tI ))) < gSDFepsilon );
        }

//-------------------------------------------------------------------------------

        Matrix< F31RMat >
        Triangle::get_barycentric_from_local_cartesian(
                                   const  Matrix< F31RMat >& aLocalPoint )
        {
            Matrix< F31RMat > aXi( 3, 1 );

            // the first coordinate
            aXi( 0 ) =  ((  mBarycentric.mLocalNodeCoordsInPlane( 0, 1 )
                    -mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ))
                    *( mBarycentric.mLocalNodeCoordsInPlane( 1, 2 )
                            -aLocalPoint( 1 ))
                            -(  mBarycentric.mLocalNodeCoordsInPlane( 1, 1 )
                                    -mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ))
                                    *( mBarycentric.mLocalNodeCoordsInPlane( 0, 2 )
                                            -aLocalPoint( 0 )))
                                            * mBarycentric.mInvTwiceArea;

            // the second coordinate
            aXi( 1 ) =  ((  mBarycentric.mLocalNodeCoordsInPlane( 1, 0 )
                    -mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ))
                    *( mBarycentric.mLocalNodeCoordsInPlane( 0, 2 )
                            -aLocalPoint( 0 ))
                            -(  mBarycentric.mLocalNodeCoordsInPlane( 0, 0 )
                                    -mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ))
                                    *(mBarycentric.mLocalNodeCoordsInPlane( 1, 2 )
                                            -aLocalPoint( 1 )))
                                            *mBarycentric.mInvTwiceArea;

            // the third coordinate
            aXi( 2 ) =  1.0 - aXi( 0 ) - aXi( 1 );

            return aXi;
        }

//-------------------------------------------------------------------------------

        real
        Triangle::distance_point_to_edge_in_local_cartesian(
                            const Matrix< F31RMat >& aLocalPoint,
                            const uint aEdge)
        {
            real tParam = 0;
            uint i;
            uint j;

            // permutation parameter of current edge
            TrianglePermutation( aEdge, i, j );

            // calculate projection of point on edge

            // tParam = 0: orthogonal intersects with point i
            // tParam = 1: orthogonal intersects with point j

            for ( uint l=0; l<2; ++l )
            {
                tParam +=  ( aLocalPoint( l )-mBarycentric.mLocalNodeCoordsInPlane( l, i ) )
                                       *mBarycentric.mLocalEdgeDirectionVectors( l, aEdge );

            }
            tParam *= mBarycentric.mLocalEdgeInverseMagnitudes( aEdge );

            Matrix< F31RMat > aDirection( 3, 1 );

            if( tParam < gSDFepsilon )
            {
                // snap to point i and set tParam = 0.0;
                aDirection( 0 ) = aLocalPoint( 0 )-mBarycentric.mLocalNodeCoordsInPlane( 0, i );
                aDirection( 1 ) = aLocalPoint( 1 )-mBarycentric.mLocalNodeCoordsInPlane( 1, i );

            } else if( tParam > 1.0-gSDFepsilon ){
                // snap to point j and set tParam = 1.0;
                aDirection( 0 ) = aLocalPoint( 0 )-mBarycentric.mLocalNodeCoordsInPlane( 0, j );
                aDirection( 1 ) = aLocalPoint( 1 )-mBarycentric.mLocalNodeCoordsInPlane( 1, j );
            } else {
                // find distance in plane
                aDirection( 0 ) = aLocalPoint( 0 )-mBarycentric.mLocalNodeCoordsInPlane( 0, i )
                                   -tParam*mBarycentric.mLocalEdgeDirectionVectors( 0, aEdge );
                aDirection( 1 ) = aLocalPoint( 1 )-mBarycentric.mLocalNodeCoordsInPlane( 1, i )
                                   -tParam*mBarycentric.mLocalEdgeDirectionVectors( 1, aEdge );
            }

            // add third dimension to distance
            aDirection( 2 ) =  aLocalPoint( 2 );

            return norm( aDirection );
        }

//-------------------------------------------------------------------------------

        Matrix< F31RMat >
        Triangle::project_point_to_local_cartesian(
                const moris::Matrix< F31RMat >& aPoint)
        {
            // fixme: times operator does not work with eigen
            // return mBarycentric.mProjectionMatrix * ( aPoint - mCenter ) ;
        	Matrix< F31RMat > aOut( 3, 1 );
            aOut.fill( 0 );

            for( uint k=0; k<3; ++k )
            {
                for( uint i=0; i<3; ++i )
                {
                    aOut( k ) +=   mBarycentric.mProjectionMatrix( k, i )
                                 * ( aPoint( i ) - mCenter( i ) );
                }
            }

            return aOut;
        }

//-------------------------------------------------------------------------------

        real
        Triangle::get_distance_to_point(
                const Matrix< F31RMat > & aPoint )
        {
            // step 1: Transform Point to in-plane coordinates
            Matrix< F31RMat > tLocalPointCoords = this->project_point_to_local_cartesian( aPoint );
            // step 2: calculate barycentric coordinates
            Matrix< F31RMat > tXi = this->get_barycentric_from_local_cartesian( tLocalPointCoords  );

            // step 3: check if we are inside the triangle
            if (       (tXi( 0 ) >= -gSDFepsilon)
                    && (tXi( 1 ) >= -gSDFepsilon)
                    && (tXi( 2 ) >= -gSDFepsilon))
            {
                // the absolute value of the local z-coordinate is the distance
                return std::abs( tLocalPointCoords ( 2 ) );
            }
            else
            {
                if( tXi( 0 ) > 0 )
                {
                    // this rules out edge 0
                    real tDist1 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 1 );
                    real tDist2 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 2 );
                    return std::min( tDist1, tDist2 );
                }
                else if( tXi( 1 ) > 0 )
                {
                    // this rules out edge 1
                    real tDist0 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 0 );
                    real tDist2 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 2 );
                    return std::min( tDist0, tDist2 );
                }
                else
                {
                    // edge 2 must be the one to rule out
                    real tDist0 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 0 );
                    real tDist1 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 1 );
                    return std::min( tDist0, tDist1 );
                }
            }
        }
//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

