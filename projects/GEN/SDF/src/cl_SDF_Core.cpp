#include "SDF_Tools.hpp"
#include "cl_SDF_Core.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Core::Core(
                const mtk::Mesh * aMesh,
                      Data      & aData ) :
                   mMesh( aMesh ),
                   mData( aData )
        {
            // fill unsure nodes list
            uint tNumberOfNodes = aMesh->get_num_nodes();

            // poputale unsure nodes list
            mData.mUnsureNodes.resize( tNumberOfNodes, 1 );
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                mData.mUnsureNodes( k ) = k;
            }
        }

//-------------------------------------------------------------------------------

        void
        Core::voxelize( const uint aAxis )
        {
            // reset unsure nodes counter
            mData.mUnsureNewNodesCount = 0;

            // get number of unsure nodes
            uint tNumberOfNodes = mData.mUnsureNodes.length();

            // loop over all nodes
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                // copy node coordinates into fixed vector
                auto tNodeCoords = mMesh->get_node_coordinate( k );

                // fixme: node coords are by default DDRMats
                // we must cast the node to a fixed size matrix
                Matrix< moris::F31RMat > tPoint;
                for( uint i=0; i<3; ++i )
                {
                    tPoint( i ) = tNodeCoords( i );
                }

                // preselect triangles for intersection test
                if(aAxis == 0)
                    this->preselect_triangles_x( tPoint );
                else if (aAxis == 1)
                    this->preselect_triangles_y( tPoint );
                else
                    this->preselect_triangles_z( tPoint );

                // from the candidate triangles, perform intersection
                if( mData.mCandidateTriangles.length() > 0 )
                {
                    this->intersect_triangles( aAxis, tPoint );

                    // intersect ray with triangles and check if node is inside
                    if( mData.mIntersectedTriangles.length() > 0 )
                    {

                    }
                }
            }
        }

//-------------------------------------------------------------------------------

        void
        Core::preselect_triangles_x( const Matrix< F31RMat >& aPoint )
        {
            // x: k = x, j = z, i = y
#ifdef MORIS_USE_ARMA

            // check bounding box in J-direction
            mData.mCandJ = arma::find(
                    (aPoint(2) - mData.mTriangleMinCoordsZ) %
                    (mData.mTriangleMaxCoordsZ - aPoint(2)) > -gSDFepsilon );

            // check bounding box in I-direction
            mData.mCandI = arma::find(
                    (aPoint(1) - mData.mTriangleMinCoordsY.elem(mData.mCandJ)) %
                    (mData.mTriangleMaxCoordsY.elem(mData.mCandJ) - aPoint(1)) > -gSDFepsilon );

            // help vector to be written in mData.mCandidateTriangles.data()
            mData.mCandK = mData.mCandJ.elem(mData.mCandI);
            // resize data object
            mData.mCandidateTriangles.resize(mData.mCandK.n_elem, 1);

            // link to current object
            arma::Mat<uint> &tCand = mData.mCandidateTriangles.matrix_data();

            // write data
            tCand = arma::conv_to<arma::Mat<uint> >::from(mData.mCandK);

#else
            // loop over all triangles in J-Direction
            uint tCountJ = 0;
            for (uint k = 0; k<mData.mNumberOfTriangles; ++k)
            {
                // check bounding box in J-direction
                if ( (aPoint(2) - mData.mTriangleMinCoordsZ(k)) *
                                     (mData.mTriangleMaxCoordsZ(k) - aPoint(2)) > -gSDFepsilon )
                {
                    // remember this triangle
                    mData.mCandJ(tCountJ) = k;

                    // increment counter
                    ++tCountJ;
                }
            }

            // counter for triangles
            uint tCount = 0;

            // reset candidate size
            mData.mCandidateTriangles.resize(mData.mNumberOfTriangles, 1);

            // loop over remaining triangles in I-direction
            for (uint k = 0; k<tCountJ; ++k)
            {
                // check bounding box in I-direction
                if((aPoint(1) - mData.mTriangleMinCoordsY(mData.mCandJ(k)))*
                (mData.mTriangleMaxCoordsY(mData.mCandJ(k)) - aPoint(1)) > -gSDFepsilon )
                {
                    mData.mCandidateTriangles(tCount) = mData.mCandJ(k);
                    ++tCount;
                }
            }

            mData.mCandidateTriangles.resize(tCount, 1);
#endif
        }

//-------------------------------------------------------------------------------

        void
        Core::preselect_triangles_y( const Matrix< F31RMat >& aPoint )
        {
            // y: k = y, j = x, i = z
#ifdef MORIS_USE_ARMA
            // check bounding box in J-direction
            mData.mCandJ = arma::find(
                    (aPoint(0) - mData.mTriangleMinCoordsX) %
                    (mData.mTriangleMaxCoordsX - aPoint(0)) > -gSDFepsilon );

            // check bounding box in I-direction
            mData.mCandI = arma::find(
                    (aPoint(2) - mData.mTriangleMinCoordsZ.elem(mData.mCandJ)) %
                    (mData.mTriangleMaxCoordsZ.elem(mData.mCandJ) - aPoint(2)) > -gSDFepsilon );

            // help vector to be written in mData.mCandidateTriangles.data()
            mData.mCandK = mData.mCandJ.elem(mData.mCandI);

            // resize data object
            mData.mCandidateTriangles.resize(mData.mCandK.n_elem, 1);

            // link to current object
            arma::Mat<uint> &tCand = mData.mCandidateTriangles.matrix_data();

            // write data
            tCand = arma::conv_to<arma::Mat<uint> >::from(mData.mCandK);

#else
            // loop over all triangles in J-Direction
            uint tCountJ = 0;
            for (uint k = 0; k<mData.mNumberOfTriangles; ++k)
            {
                // check bounding box in J-direction
                if ((aPoint(0) - mData.mTriangleMinCoordsX(k)) *
                        (mData.mTriangleMaxCoordsX(k) - aPoint(0)) > -gSDFepsilon )
                {
                    // remember this triangle
                    mData.mCandJ(tCountJ) = k;

                    // increment counter
                    ++tCountJ;
                }
            }

            // counter for triangles
            uint tCount = 0;

            // reset candidate size
            mData.mCandidateTriangles.resize(mData.mNumberOfTriangles, 1);

            // loop over remaining triangles in I-direction
            for (uint k = 0; k<tCountJ; ++k)
            {
                // check bounding box in I-direction
                if((aPoint(2) - mData.mTriangleMinCoordsZ(mData.mCandJ(k)))*
                        (mData.mTriangleMaxCoordsZ(mData.mCandJ(k)) - aPoint(2)) > -gSDFepsilon )
                {
                    mData.mCandidateTriangles(tCount) = mData.mCandJ(k);
                    ++tCount;
                }
            }

            mData.mCandidateTriangles.resize(tCount, 1);
#endif
        }

//-------------------------------------------------------------------------------

        void
        Core::preselect_triangles_z( const Matrix< F31RMat >& aPoint )
        {
            // z: k = z, j = y, i = x
#ifdef MORIS_USE_ARMA

            //bool_t tNothingFound = true;

            // check bounding box in J-direction
            mData.mCandJ = arma::find(
                    (aPoint(1) - mData.mTriangleMinCoordsY) %
                    (mData.mTriangleMaxCoordsY - aPoint(1)) > -gSDFepsilon );
            // check bounding box in I-direction
            mData.mCandI = arma::find(
                    (aPoint(0) - mData.mTriangleMinCoordsX.elem(mData.mCandJ)) %
                    (mData.mTriangleMaxCoordsX.elem(mData.mCandJ) - aPoint(0)) > -gSDFepsilon );

            // help vector to be written in mData.mCandidateTriangles.data()
            mData.mCandK = mData.mCandJ.elem(mData.mCandI);

            // resize data object
            mData.mCandidateTriangles.resize(mData.mCandK.n_elem, 1);

            // link to current object
            arma::Mat<uint> &tCand = mData.mCandidateTriangles.matrix_data();

            // write data
            tCand = arma::conv_to<arma::Mat<uint> >::from(mData.mCandK);
#else
            // loop over all triangles in J-Direction
            uint tCountJ = 0;
            for (uint k = 0; k<mData.mNumberOfTriangles; ++k)
            {
                // check bounding box in J-direction
                if( (aPoint(1) - mData.mTriangleMinCoordsY(k)) *
                        (mData.mTriangleMaxCoordsY(k) - aPoint(1)) > -gSDFepsilon )
                {
                    // remember this triangle
                    mData.mCandJ(tCountJ) = k;

                    // increment counter
                    ++tCountJ;
                }
            }

            // counter for triangles
            uint tCount = 0;

            // reset candidate size
            mData.mCandidateTriangles.resize(mData.mNumberOfTriangles, 1);

            // loop over remaining triangles in I-direction
            for (uint k = 0; k<tCountJ; ++k)
            {
                // check bounding box in I-direction
                if((aPoint(0) - mData.mTriangleMinCoordsX(mData.mCandJ(k)))*
                        (mData.mTriangleMaxCoordsX(mData.mCandJ(k)) - aPoint(0)) > -gSDFepsilon )
                {
                    mData.mCandidateTriangles(tCount) = mData.mCandJ(k);
                    ++tCount;
                }
            }

            mData.mCandidateTriangles.resize(tCount, 1);
#endif
        }
//-------------------------------------------------------------------------------

        void
        Core::intersect_triangles(
                const uint               aAxis,
                const Matrix< F31RMat >& aPoint )
        {
            // counter for intersected triangles
            uint tNumberOfTriangles = mData.mCandidateTriangles.length();

            // resize member container
            mData.mIntersectedTriangles.resize( tNumberOfTriangles, 1 );

            // reset counter
            uint tCount = 0;

            // loop over all candidate triangles, and check intersection criterion
            for ( uint k=0; k < tNumberOfTriangles; ++k)
            {
                // get index of current triangle
                uint tThisTriangle = mData.mCandidateTriangles( k );

                if ( mData.mTriangles( tThisTriangle )->check_edge( 0, aAxis, aPoint ) )
                {
                    if ( mData.mTriangles( tThisTriangle )->check_edge( 1, aAxis, aPoint ) )
                    {
                        if ( mData.mTriangles( tThisTriangle )->check_edge( 2, aAxis, aPoint ) )
                        {
                            mData.mIntersectedTriangles(tCount) = tThisTriangle;
                            ++tCount;
                        }
                    }
                }

            }

            // adapt size of list
            mData.mIntersectedTriangles.resize( tCount, 1 );
        }

//-------------------------------------------------------------------------------

        void
        Core::check_if_node_is_inside(
                const uint                aAxis,
                const uint                aLocalNodeInd,
                const Matrix< F31RMat > & aPoint )
        {
            uint tNumCoordsK = mData.mCoordsK.length();

            bool tNodeIsInside = false;

            // only even number of intersections is considered
            if( tNumCoordsK % 2 == 0 )
            {
                // loop over k-Coordinates.
                for ( uint k=0; k< tNumCoordsK / 2 ; ++k )
                {
                    tNodeIsInside = ( aPoint( aAxis ) > mData.mCoordsK( 2 * k ) ) &&
                                    ( aPoint( aAxis ) < mData.mCoordsK( 2 * k + 1 ) );

                    // break the loop if inside
                    if ( tNodeIsInside )
                    {
                        break;
                    }
                }

                std::cout << "set inside flag" << std::endl;
                // set the inside flag of this node to the corresponding value
                if( tNodeIsInside )
                {
                    // mData.mLocalNodeInsideFlags.set( aLocalNodeInd );
                }
                else
                {
                    // mData.mLocalNodeInsideFlags.reset( aLocalNodeInd );
                }
            }
            else
            {
                // add to unsure node list
                mData.mUnsureNodesNew( mData.mUnsureNewNodesCount ) = aLocalNodeInd;

                // increment unsure nodes counter
                ++mData.mUnsureNewNodesCount;
            }

        }

//-------------------------------------------------------------------------------

        void
        Core::calculate_candidate_points_and_buffer_diagonal()
        {

        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */
