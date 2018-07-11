#include "cl_Ge_SDF_Triangle.hpp"

// =============================================================================
//  Public
// =============================================================================

void
ge::SDF_Triangle::set_node_coords (const moris::Mat< moris::real > &aAllNodeCoords)
{

    // copy node coordinates
    for (moris::uint k = 0; k < 3; ++k)
    {
        for (moris::uint i = 0; i < 3; ++i)
        {
            mNodeCoords( i, k ) = aAllNodeCoords( i, mNodeIDs(k) );
        }
    }

    // identify minimum and maximum coordinate
    for (moris::uint i = 0; i < 3; ++i)
    {
        mMinCoord( i ) = ge::min(mNodeCoords( i, 0 ),
                                 mNodeCoords( i, 1 ),
                                 mNodeCoords( i, 2 ));

        mMaxCoord( i ) = ge::max(mNodeCoords( i, 0 ),
                                 mNodeCoords( i, 1 ),
                                 mNodeCoords( i, 2 ));
    }

    // center of triangle
    for (moris::uint i = 0; i < 3; ++i)
    {
        mCenter( i ) = (  mNodeCoords( i, 0 )
                          + mNodeCoords( i, 1 )
                          + mNodeCoords( i, 2) ) / 3;
    }

    // help vectors: direction of sides 1 and 2
    moris::Mat< moris::real > tDirection01( 3 , 1 );
    moris::Mat< moris::real > tDirection02( 3 , 1 );

    for (moris::uint i = 0; i < 3; ++i) {
        tDirection01( i ) = mNodeCoords( i, 1 )-mNodeCoords( i, 0 );
        tDirection02( i ) = mNodeCoords( i, 2 )-mNodeCoords( i, 0 );
    }

    // Normal Vector
    mNormal = ge::cross( tDirection01, tDirection02 );
    moris::Mat< moris::real > tDirectionOrtho = ge::cross( mNormal, tDirection01 );

    // normalize vectors
    moris::real tNorm = ge::norm( mNormal );
    for (moris::uint i = 0; i < 3; ++i) {
        mNormal( i ) /= tNorm;
    }

    // normalize tDirection20
    moris::real tNorm10 = ge::norm( tDirection01 );
    // normalize tDirectionOrtho
    moris::real tNormOrtho = ge::norm( tDirectionOrtho );

    // Projection matrix
    for (moris::uint k = 0; k < 3; ++k) {
        mBarycentric.mProjectionMatrix( 0, k ) = tDirection01( k )/tNorm10;
        mBarycentric.mProjectionMatrix( 1, k ) = tDirectionOrtho( k )/tNormOrtho;
        mBarycentric.mProjectionMatrix( 2, k ) = mNormal( k );
    }

    // node coordinates in triangle plane
    mBarycentric.mLocalNodeCoordsInPlane.fill( 0 );
    for (moris::uint k = 0; k < 3; ++k) {
        for (moris::uint j = 0; j < 3; ++j) {
            for (moris::uint i = 0; i < 2; ++i) {
                mBarycentric.mLocalNodeCoordsInPlane( i, k ) +=
                        mBarycentric.mProjectionMatrix( i, j )
                        *(mNodeCoords( j, k )-mCenter( j ));
            }
        }
    }

    // enlarge triangle
    mBarycentric.mLocalNodeCoordsInPlane( 0, 0 ) -= MORIS_GE_EPSILON;
    mBarycentric.mLocalNodeCoordsInPlane( 1, 0 ) -= MORIS_GE_EPSILON;
    mBarycentric.mLocalNodeCoordsInPlane( 0, 1 ) += MORIS_GE_EPSILON;
    mBarycentric.mLocalNodeCoordsInPlane( 1, 1 ) -= MORIS_GE_EPSILON;
    mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ) += MORIS_GE_EPSILON;
    // twice the area
    mBarycentric.mTwiceArea =  (mBarycentric.mLocalNodeCoordsInPlane( 0, 0 )
                                - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ))
                               *(mBarycentric.mLocalNodeCoordsInPlane( 1, 1 )
                                 - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ))
                               -(mBarycentric.mLocalNodeCoordsInPlane( 1, 0 )
                                 - mBarycentric.mLocalNodeCoordsInPlane( 1, 2 ))
                                *(mBarycentric.mLocalNodeCoordsInPlane( 0, 1 )
                                  - mBarycentric.mLocalNodeCoordsInPlane( 0, 2 ));

    MORIS_ASSERT(mBarycentric.mTwiceArea > 2*MORIS_GE_EPSILON,
                 "A degenerated triangle was found.");

    mBarycentric.mInvTwiceArea = 1.0/mBarycentric.mTwiceArea;
    // Edge directions
    for (moris::uint k = 0; k < 3; ++k)
    {
        moris::uint i;
        moris::uint j;
        TrianglePermutation(k,i,j);
        for (moris::uint l = 0; l < 2; ++l)
        {
            mBarycentric.mLocalEdgeDirectionVectors( l, k ) = mBarycentric.mLocalNodeCoordsInPlane( l, j )
                                                              -mBarycentric.mLocalNodeCoordsInPlane( l, i );
        }
        moris::real tMagnitude = 0;
        for (moris::uint l = 0; l < 2; ++l)
        {
            tMagnitude +=  mBarycentric.mLocalEdgeDirectionVectors( l, k )
                           *mBarycentric.mLocalEdgeDirectionVectors( l, k );
        }
        mBarycentric.mLocalEdgeInverseMagnitudes( k ) = 1.0/tMagnitude;
    }

    // Hesse Parameter of triangle plane
    mHesse = ge::dot( mCenter, mNormal );
    // values for cross prediction
    moris::uint i;
    moris::uint j;
    moris::uint p;
    moris::uint q;
    for (moris::uint k = 0; k < 3; ++k)
    {
        TrianglePermutation(k, i, j);

        for (moris::uint r = 0; r < 3; ++r)
        {
            TrianglePermutation(r, p, q);
            moris::real tDelta = mNodeCoords(i, p) - mNodeCoords(i, q);
            if (abs(tDelta) < MORIS_GE_EPSILON)
            {
                if (tDelta < 0)
                    tDelta = -MORIS_GE_EPSILON;
                else
                    tDelta = MORIS_GE_EPSILON;
            }

            mPredictYRA( r, k ) =   (mNodeCoords( j, p ) - mNodeCoords( j, q ))/tDelta;

            mPredictY( r, k )   =   mNodeCoords( j, p )
                                    + mPredictYRA( r, k ) * (mNodeCoords( i, r ) - mNodeCoords( i, p ));

            mPredictYRB( r, k ) = mNodeCoords( j, p ) - mNodeCoords( i, p ) * mPredictYRA( r, k );

        }
    }
}

// -----------------------------------------------------------------------------

moris::Mat< moris::real >
ge::SDF_Triangle::intersect_with_line(
        const moris::Mat< moris::real >& aPoint,
        const moris::Mat< moris::real >& aDirection)
{

    // let the line be
    // g = aPoint + tParam*aDirection

    // get parameter of line
    moris::real tParam = (mHesse - ge::dot(mNormal, aPoint)) /
                         ge::dot( mNormal, aDirection );

    // calculate cross point
    moris::Mat< moris::real > aNormalPoint(3,1);
    for(moris::uint i=0; i<3; ++i)
    {
        aNormalPoint(i) = aPoint(i) + tParam*aDirection(i);
    }

    return aNormalPoint;
}

// -----------------------------------------------------------------------------

moris::bool_t ge::SDF_Triangle::check_edge(
        const moris::uint aEdge,
        const moris::uint aAxis,
        const moris::Mat< moris::real >& aPoint ){

    moris::uint tI;
    moris::uint tJ;
    moris::uint tP;
    moris::uint tQ;

    // permutation parameter for axis
    TrianglePermutation( aAxis, tI, tJ );

    // permutation parameter for edge
    TrianglePermutation( aEdge, tP, tQ );

    // R
    moris::real tPredictYR = mPredictYRA( aEdge, aAxis )*aPoint( tI ) + mPredictYRB( aEdge, aAxis );

    // check of point is within all three projected edges
    return  ( (mPredictY( aEdge, aAxis ) > mNodeCoords( tJ, aEdge ) )    && (tPredictYR  + MORIS_GE_EPSILON > aPoint( tJ ))) ||
               ( (mPredictY( aEdge, aAxis ) < mNodeCoords( tJ, aEdge ) ) && (tPredictYR - MORIS_GE_EPSILON  < aPoint( tJ ))) ||
               (ge::abs((mNodeCoords( tJ, tP )-mNodeCoords( tJ, tQ ))
                        * (mNodeCoords( tI, tP )-aPoint( tI ))) < MORIS_GE_EPSILON);
}

// -----------------------------------------------------------------------------

moris::Mat< moris::real >
ge::SDF_Triangle::get_barycentric_from_local_cartesian( const moris::Mat< moris::real >& aLocalPoint )
{
    moris::Mat< moris::real > aXi( 3, 1 );

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

// -----------------------------------------------------------------------------

moris::real
ge::SDF_Triangle::distance_point_to_edge_in_local_cartesian(
        const moris::Mat< moris::real >& aLocalPoint,
        const moris::uint aEdge)
{
    moris::real tParam = 0;
    moris::uint i;
    moris::uint j;

    // permutation parameter of current edge
    TrianglePermutation( aEdge, i, j );

    // calculate projection of point on edge

    // tParam = 0: orthogonal intersects with point i
    // tParam = 1: orthogonal intersects with point j
    for (moris::uint l=0; l<2; ++l)
    {
        tParam +=  (aLocalPoint( l )-mBarycentric.mLocalNodeCoordsInPlane( l, i ))
                   *mBarycentric.mLocalEdgeDirectionVectors( l, aEdge );

    }
    tParam *= mBarycentric.mLocalEdgeInverseMagnitudes( aEdge );

    moris::Mat< moris::real > aDirection( 3, 1 );

    if( tParam < MORIS_GE_EPSILON)
    {
        // snap to point i and set tParam = 0.0;
        for (moris::uint l=0; l<2; ++l)
        {
            aDirection( l ) = aLocalPoint( l )-mBarycentric.mLocalNodeCoordsInPlane( l, i );
        }
    } else if(tParam > 1.0-MORIS_GE_EPSILON){
        // snap to point j and set tParam = 1.0;
        for (moris::uint l=0; l<2; ++l)
        {
            aDirection( l ) = aLocalPoint( l )-mBarycentric.mLocalNodeCoordsInPlane( l, j );
        }
    } else {
        // find distance in plane
        for (moris::uint l=0; l<2; ++l){
            aDirection( l ) = aLocalPoint( l )-mBarycentric.mLocalNodeCoordsInPlane( l, i )
                              -tParam*mBarycentric.mLocalEdgeDirectionVectors( l, aEdge );
        }
    }

    // add third dimension to distance
    aDirection( 2 ) =  aLocalPoint( 2 );

    return ge::norm( aDirection );
}

// -----------------------------------------------------------------------------

moris::real
ge::SDF_Triangle::get_distance_to_point(
        const moris::Mat< moris::real >& aPoint)
{

    // step 1: Transform Point to in-plane coordinates
    moris::Mat< moris::real > tLocalPointCoords = project_point_to_local_cartesian( aPoint );

    // step 2: calculate barycentric coordinates
    moris::Mat< moris::real > tXi = get_barycentric_from_local_cartesian( tLocalPointCoords  );

    // step 3: check if we are inside the triangle
    if (    (tXi( 0 ) >= -MORIS_GE_EPSILON)
            && (tXi( 1 ) >= -MORIS_GE_EPSILON)
            && (tXi( 2 ) >= -MORIS_GE_EPSILON))
    {
        // the absolute value of the local z-coordinate is the distance
        return ge::abs(tLocalPointCoords ( 2 ));
    }
    else
    {
        if( tXi( 0 ) > 0 )
        {
            // this rules out edge 0
            moris::real tDist1 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 1 );
            moris::real tDist2 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 2 );
            return ge::min( tDist1, tDist2 );
        }
        else if( tXi( 1 ) > 0 )
        {
            // this rules out edge 1
            moris::real tDist0 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 0 );
            moris::real tDist2 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 2 );
            return ge::min( tDist0, tDist2 );
        }
        else
        {
            // edge 2 must be the one to rule out
            moris::real tDist0 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 0 );
            moris::real tDist1 = distance_point_to_edge_in_local_cartesian( tLocalPointCoords , 1 );
            return ge::min( tDist0, tDist1 );
        }
    }
}

// =============================================================================
//  Private
// =============================================================================
