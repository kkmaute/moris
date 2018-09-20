/*
 * cl_GE_SDF_Data.cpp
 *
 *  Created on: Mar 7, 2018
 *      Author: messe
 */

#include "cl_Ge_SDF_Data.hpp"

// -----------------------------------------------------------------------------

ge::SDF_Data::SDF_Data(moris::Matrix< moris::DDRMat >& aLocalSDF,
                       moris::BoostBitset& aSDFBitset,
                       const moris::Matrix< moris::DDUMat > &aTriangleTopology,
                       const moris::Matrix< moris::DDRMat > &aTriangleNodeCoords):
                       mLocalSDF(aLocalSDF),
                       mNumberOfTriangles(aTriangleTopology.n_cols()),
                       mTriangleMinCoordsX(mNumberOfTriangles, 1),
                       mTriangleMinCoordsY(mNumberOfTriangles, 1),
                       mTriangleMinCoordsZ(mNumberOfTriangles, 1),
                       mTriangleMaxCoordsX(mNumberOfTriangles, 1),
                       mTriangleMaxCoordsY(mNumberOfTriangles, 1),
                       mTriangleMaxCoordsZ(mNumberOfTriangles, 1),
#ifdef MORIS_USE_ARMA
                       mCandI(mNumberOfTriangles),
                       mCandJ(mNumberOfTriangles),
                       mCandK(mNumberOfTriangles),
#else
                       mCandJ(mNumberOfTriangles, 1),
#endif
                       mCandidateTriangles(mNumberOfTriangles, 1),
                       mIntersectedTriangles(mNumberOfTriangles, 1),
                       mCoordsK(mNumberOfTriangles, 1),
                       mLocalNodeHasSdfFlag(aSDFBitset)
{
    init_triangles(aTriangleTopology, aTriangleNodeCoords);
}

// -----------------------------------------------------------------------------

void
ge::SDF_Data::init_triangles
        (const moris::Matrix< moris::DDUMat > &aTriangleTopology,
         const moris::Matrix< moris::DDRMat > &aTriangleNodeCoords)
{
    // initialize triangle cell
    mTriangles.clear();
    mTriangles.reserve( mNumberOfTriangles );

    // create the triangle objects
    for( moris::uint k=0; k< mNumberOfTriangles; ++k ){
        ge::SDF_Triangle tTriangle( aTriangleTopology( 0, k ),
                                    aTriangleTopology( 1, k ),
                                    aTriangleTopology( 2, k ));

        // add triangle object to cell
        mTriangles.push_back( tTriangle );
    }

    // reset coordinates for triangles
    for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
    {
        mTriangles( k ).set_node_coords(aTriangleNodeCoords);
    }

    // copy triangle bounding box data
    for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
    {
        // minimum triangle coordinates for lower left point of bounding box
        mTriangleMinCoordsX( k )
                = mTriangles( k ).get_min_coord( 0 );

        mTriangleMinCoordsY( k )
                = mTriangles( k ).get_min_coord( 1 );

        mTriangleMinCoordsZ( k )
                = mTriangles( k ).get_min_coord( 2 );

        // maximum triangle coordinates for upper right point of bounding box
        mTriangleMaxCoordsX( k )
                = mTriangles( k ).get_max_coord( 0 );

        mTriangleMaxCoordsY( k )
                = mTriangles( k ).get_max_coord( 1 );

        mTriangleMaxCoordsZ ( k )
                = mTriangles( k ).get_max_coord( 2 );
    }
}

// -----------------------------------------------------------------------------
void
ge::SDF_Data::init_data_fields (const moris::uint aNumberOfNodes,
                                     const moris::uint aNumberOfElements)
{

    mLocalSDF.set_size(aNumberOfNodes, 1);

    // initialize list of unsure nodes
    mUnsureNodes.set_size(aNumberOfNodes, 1);
    for(moris::uint k=0; k<aNumberOfNodes; ++k)
    {
        mUnsureNodes(k) = k;
    }

    // reserve memory for new list
    mUnsureNodesNew.set_size(aNumberOfNodes, 1);

    // reserve memory for inside bitset and reset
    mLocalNodeInsideFlags.resize(aNumberOfNodes);
    mLocalNodeInsideFlags.reset();

    // reserve memory for inside bitset and reset
    mLocalNodeCandidateFlags.resize(aNumberOfNodes);
    mLocalNodeCandidateFlags.reset();

    // reserve memory of elements at surface
    mLocalElementsAtSurface.set_size(aNumberOfElements, 1);

    // reserve memory for elements in volume
    mLocalElementsInVolume.set_size(aNumberOfElements, 1);

    // reset buffer diagonal
    mBufferDiagonal = 0;

    // reserve memory for usf bitset
    mLocalNodeHasSdfFlag.resize(aNumberOfNodes);
    mLocalNodeHasSdfFlag.reset();
}
