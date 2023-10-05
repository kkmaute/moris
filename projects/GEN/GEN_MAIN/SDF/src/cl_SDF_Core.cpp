/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Core.cpp
 *
 */

#include <fstream>

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Tools.hpp"
#include "SDF_Tools.hpp"
#include "cl_SDF_Core.hpp"
#include "fn_sort.hpp"
#include "fn_print.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Core::Core( Mesh& aMesh, Data& aData, bool aVerbose )
                : mMesh( aMesh )
                , mData( aData )
                , mVerbose( aVerbose )
        {
            // fill unsure nodes list
            // uint tNumberOfNodes = aMesh.get_num_nodes();
        }
        //-------------------------------------------------------------------------------

        void
        Core::calculate_raycast(
                Matrix< IndexMat >& aElementsAtSurface,
                Matrix< IndexMat >& aElementsInVolume )
        {

            // call private routine
            this->calculate_raycast();

            // assign element containers
            aElementsAtSurface.set_size( mData.mSurfaceElements, 1 );

            aElementsInVolume.set_size( mData.mVolumeElements, 1 );

            // counters
            uint tSurfaceCount = 0;
            uint tVolumeCount  = 0;

            // get number of elements
            uint tNumberOfElements = mMesh.get_num_elems();

            // loop over all elements
            for ( uint k = 0; k < tNumberOfElements; ++k )
            {
                // get pointer to element
                Cell* tElement = mMesh.get_cell( k );

                if ( tElement->is_on_surface() )
                {
                    aElementsAtSurface( tSurfaceCount++ ) = k;
                }
                else if ( tElement->is_in_volume() )
                {
                    aElementsInVolume( tVolumeCount++ ) = k;
                }
            }

            // make sure that everything is OK
            MORIS_ASSERT( tSurfaceCount == mData.mSurfaceElements,
                    "Number of surface elements does not match. Surface Count after raycast: %d\t Struct Surface Elements: %d.",
                    tSurfaceCount,
                    mData.mSurfaceElements );

            MORIS_ASSERT( tVolumeCount == mData.mVolumeElements,
                    "Number of volume elements does not match. Volume Count after raycast: %d\t Struct Volume Elements: %d.",
                    tVolumeCount,
                    mData.mVolumeElements );
        }

        //-------------------------------------------------------------------------------

        void
        Core::calculate_raycast(
                Matrix< IndexMat >& aElementsAtSurface )
        {

            // call private routine
            this->calculate_raycast();

            // assign element containers
            aElementsAtSurface.set_size( mData.mSurfaceElements, 1 );

            // counters
            uint tSurfaceCount = 0;

            // get number of elements
            uint tNumberOfElements = mMesh.get_num_elems();

            // loop over all elements
            for ( uint k = 0; k < tNumberOfElements; ++k )
            {
                // get pointer to element
                Cell* tElement = mMesh.get_cell( k );

                if ( tElement->is_on_surface() )
                {
                    aElementsAtSurface( tSurfaceCount++ ) = k;
                }
            }

            // make sure that everything is OK
            MORIS_ASSERT( tSurfaceCount = mData.mSurfaceElements,
                    "Number of surface elements does not match" );
        }

        //-------------------------------------------------------------------------------

        void
        Core::calculate_raycast()
        {

            tic tTimer;

            // set unsure flag of all nodes to true
            uint tNumberOfNodes = mMesh.get_num_nodes();

            for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; ++iNodeIndex )
            {
                mMesh.get_vertex( iNodeIndex )->reset();
            }
            mData.mUnsureNodesCount = tNumberOfNodes;

            // flag that marks if rotation was called
            bool tRotation = false;

            while ( mData.mUnsureNodesCount > 0 )
            {
                // perform voxelizing algorithm in z-direction
                voxelize( 2 );
                if ( mData.mUnsureNodesCount > 0 )
                {
                    // perform voxelizing algorithm in y-direction
                    voxelize( 1 );
                    if ( mData.mUnsureNodesCount > 0 )
                    {
                        // perform voxelizing algorithm in x-direction
                        voxelize( 0 );
                    }
                }

                if ( mData.mUnsureNodesCount > 0 )
                {
                    tRotation = true;

                    this->random_rotation();
                }
            }

            if ( tRotation )
            {
                this->undo_rotation();
            }

            // remaining nodes are pushed outside
            this->force_unsure_nodes_outside();

            // identify elements in surface, volume and candidates
            this->calculate_candidate_points_and_buffer_diagonal();

            if ( mVerbose )
            {
                // stop the timer
                real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

                // print elapsed time
                if ( par_size() == 1 )
                {
                    std::fprintf( stdout, "Time for ray cast              : %5.3f [sec]\n", tElapsedTime / 1000.0 );
                }
                else
                {
                    std::fprintf( stdout, "Proc % i - Time for ray cast              : %5.3f [sec]\n", (int)par_rank(), tElapsedTime / 1000.0 );
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::calculate_raycast_and_sdf( Matrix< DDRMat >& aSDF )
        {
            this->calculate_raycast();

            moris::Cell< Vertex* > tCandidateList;          //========================================
            tCandidateList = this->set_candidate_list();    //===================================

            this->calculate_udf( tCandidateList );
            this->sweep();
            this->fill_sdf_with_values( aSDF );
        }

        //-------------------------------------------------------------------------------

        void
        Core::calculate_raycast_and_sdf(
                Matrix< DDRMat >&   aSDF,
                Matrix< IndexMat >& aElementsAtSurface,
                Matrix< IndexMat >& aElementsInVolume )
        {
            this->calculate_raycast( aElementsAtSurface, aElementsInVolume );

            moris::Cell< Vertex* > tCandidateList;          //========================================
            tCandidateList = this->set_candidate_list();    //===================================

            this->calculate_udf( tCandidateList );
            this->sweep();
            this->fill_sdf_with_values( aSDF );
        }

        //-------------------------------------------------------------------------------

        void
        Core::voxelize( const uint aAxis )
        {
            // reset unsure nodes counter
            mData.mUnsureNodesCount = 0;

            // get number of unsure nodes
            uint tNumberOfNodes = mMesh.get_num_nodes();

            // This loop is currently unnecessary, but may be needed if problems arise
            /*
            for( Facet* tFacet : mData.mFacets )
            {
                tFacet->unflag();
            }
            */

            // loop over all nodes
            for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; ++iNodeIndex )
            {
                if ( mMesh.get_vertex( iNodeIndex )->is_flagged() )
                {
                    // get node coordinate
                    const Matrix< DDRMat >& tPoint = mMesh.get_node_coordinate( iNodeIndex );

                    // preselect triangles for intersection test
                    if ( aAxis == 0 )
                        this->preselect_triangles_x( tPoint );
                    else if ( aAxis == 1 )
                        this->preselect_triangles_y( tPoint );
                    else
                        this->preselect_triangles_z( tPoint );

                    // from the candidate triangles, perform intersection
                    if ( mData.mCandidateFacets.length() > 0 )
                    {
                        this->intersect_triangles( aAxis, tPoint );

                        // intersect ray with triangles and check if node is inside
                        if ( mData.mIntersectedTriangles.size() > 0 )
                        {
                            this->intersect_ray_with_triangles( aAxis, tPoint, iNodeIndex );

                            this->check_if_node_is_inside( aAxis, iNodeIndex );
                        }
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------
        void
        Core::calculate_udf( moris::Cell< Vertex* >& aCandidateList )
        {
            tic tTimer;

            // get number of triangles
            uint tNumberOfFacets = mData.mFacets.size();
            std::cout << "number of triangles            : " << tNumberOfFacets << std::endl;    //======
            // loop over all triangles
            for ( uint k = 0; k < tNumberOfFacets; ++k )
            {
                // get pointer to triangle
                Facet* tFacet = mData.mFacets( k );

                // get nodes within triangle
                moris::Cell< Vertex* > tNodes;

                this->get_nodes_withing_bounding_box_of_triangle(
                        tFacet, tNodes, aCandidateList );

                // get number of nodes
                uint tNumberOfNodes = tNodes.size();

                // calculate distance of this point to the triangle
                // and update udf value if it is smaller
                for ( uint i = 0; i < tNumberOfNodes; ++i )
                {
                    // update UDF of this node
                    tNodes( i )->update_udf( tFacet );
                }

            }    // end loop over all triangles

            if ( mVerbose )
            {
                // stop the timer
                real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

                // print elapsed time
                if ( par_size() == 1 )
                {
                    std::fprintf( stdout, "Time for udf                   : %5.3f [sec]\n", tElapsedTime / 1000.0 );
                }
                else
                {
                    std::fprintf( stdout, "Proc % i - Time for udf                   : %5.3f [sec]\n", (int)par_rank(), tElapsedTime / 1000.0 );
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::preselect_triangles_x( const Matrix< F31RMat >& aPoint )
        {
            // x: k = x, j = z, i = y
            // #ifdef MORIS_USE_ARMA

            //             // check bounding box in J-direction
            //             mData.mCandJ = arma::find(
            //                     ( aPoint( 2 ) - mData.mTriangleMinCoordsZ ) % ( mData.mTriangleMaxCoordsZ - aPoint( 2 ) ) > -gSDFepsilon );

            //             // check bounding box in I-direction
            //             mData.mCandI = arma::find(
            //                     ( aPoint( 1 ) - mData.mFacetMinCoordsY.elem( mData.mCandJ ) ) % ( mData.mFacetMaxCoordsY.elem( mData.mCandJ ) - aPoint( 1 ) ) > -gSDFepsilon );

            //             // help vector to be written in mData.mCandidateFacets.data()
            //             mData.mCandK = mData.mCandJ.elem( mData.mCandI );
            //             // resize data object
            //             mData.mCandidateFacets.resize( mData.mCandK.n_elem, 1 );

            //             // link to current object
            //             arma::Mat< uint >& tCand = mData.mCandidateFacets.matrix_data();

            //             // write data
            //             tCand = arma::conv_to< arma::Mat< uint > >::from( mData.mCandK );

            // #else
            // loop over all triangles in K-Direction
            uint tCountJ = 0;
            for ( uint k = 0; k < mData.mNumberOfFacets; ++k )
            {
                // check bounding box in K-direction
                if ( ( aPoint( 2 ) - mData.mFacetMinCoords( k, 2 ) ) * ( mData.mFacetMaxCoords( k, 2 ) - aPoint( 2 ) ) > -gSDFepsilon )
                {
                    // remember this triangle
                    mData.mCandJ( tCountJ ) = k;

                    // increment counter
                    ++tCountJ;
                }
            }

            // counter for triangles
            uint tCount = 0;

            // reset candidate size
            mData.mCandidateFacets.resize( mData.mNumberOfFacets, 1 );

            // loop over remaining triangles in J-direction
            for ( uint k = 0; k < tCountJ; ++k )
            {
                // check bounding box in J-direction
                if ( ( aPoint( 1 ) - mData.mFacetMinCoords( mData.mCandJ( k ), 1 ) ) * ( mData.mFacetMaxCoords( mData.mCandJ( k ), 1 ) - aPoint( 1 ) ) > -gSDFepsilon )
                {
                    mData.mCandidateFacets( tCount ) = mData.mCandJ( k );
                    ++tCount;
                }
            }

            mData.mCandidateFacets.resize( tCount, 1 );
            // #endif
        }

        //-------------------------------------------------------------------------------

        void
        Core::preselect_triangles_y( const Matrix< F31RMat >& aPoint )
        {
            // y: k = y, j = x, i = z
            // #ifdef MORIS_USE_ARMA
            //             // check bounding box in J-direction
            //             mData.mCandJ = arma::find(
            //                     ( aPoint( 0 ) - mData.mFacetMinCoordsX ) % ( mData.mFacetMaxCoordsX - aPoint( 0 ) ) > -gSDFepsilon );

            //             // check bounding box in I-direction
            //             mData.mCandI = arma::find(
            //                     ( aPoint( 2 ) - mData.mTriangleMinCoordsZ.elem( mData.mCandJ ) ) % ( mData.mTriangleMaxCoordsZ.elem( mData.mCandJ ) - aPoint( 2 ) ) > -gSDFepsilon );

            //             // help vector to be written in mData.mCandidateFacets.data()
            //             mData.mCandK = mData.mCandJ.elem( mData.mCandI );

            //             // resize data object
            //             mData.mCandidateFacets.resize( mData.mCandK.n_elem, 1 );

            //             // link to current object
            //             arma::Mat< uint >& tCand = mData.mCandidateFacets.matrix_data();

            //             // write data
            //             tCand = arma::conv_to< arma::Mat< uint > >::from( mData.mCandK );

            // #else
            // loop over all triangles in I-Direction
            uint tCountJ = 0;
            for ( uint k = 0; k < mData.mNumberOfFacets; ++k )
            {
                // check bounding box in I-direction
                if ( ( aPoint( 0 ) - mData.mFacetMinCoords( k, 0 ) ) * ( mData.mFacetMaxCoords( k, 0 ) - aPoint( 0 ) ) > -gSDFepsilon )
                {
                    // remember this triangle
                    mData.mCandJ( tCountJ ) = k;

                    // increment counter
                    ++tCountJ;
                }
            }

            // counter for triangles
            uint tCount = 0;

            // reset candidate size
            mData.mCandidateFacets.resize( mData.mNumberOfFacets, 1 );

            // loop over remaining triangles in K-direction
            for ( uint k = 0; k < tCountJ; ++k )
            {
                // check bounding box in K-direction
                if ( ( aPoint( 2 ) - mData.mFacetMinCoords( mData.mCandJ( k ), 2 ) ) * ( mData.mFacetMaxCoords( mData.mCandJ( k ), 2 ) - aPoint( 2 ) ) > -gSDFepsilon )
                {
                    mData.mCandidateFacets( tCount ) = mData.mCandJ( k );
                    ++tCount;
                }
            }

            mData.mCandidateFacets.resize( tCount, 1 );
            // #endif
        }

        //-------------------------------------------------------------------------------

        void
        Core::preselect_triangles_z( const Matrix< F31RMat >& aPoint )
        {
            // z: k = z, j = y, i = x
            // #ifdef MORIS_USE_ARMA

            //             // bool_t tNothingFound = true;

            //             // check bounding box in J-direction
            //             mData.mCandJ = arma::find(
            //                     ( aPoint( 1 ) - mData.mFacetMinCoordsY ) % ( mData.mFacetMaxCoordsY - aPoint( 1 ) ) > -gSDFepsilon );
            //             // check bounding box in I-direction
            //             mData.mCandI = arma::find(
            //                     ( aPoint( 0 ) - mData.mFacetMinCoordsX.elem( mData.mCandJ ) ) % ( mData.mFacetMaxCoordsX.elem( mData.mCandJ ) - aPoint( 0 ) ) > -gSDFepsilon );

            //             // help vector to be written in mData.mCandidateFacets.data()
            //             mData.mCandK = mData.mCandJ.elem( mData.mCandI );

            //             // resize data object
            //             mData.mCandidateFacets.resize( mData.mCandK.n_elem, 1 );

            //             // link to current object
            //             arma::Mat< uint >& tCand = mData.mCandidateFacets.matrix_data();

            //             // write data
            //             tCand = arma::conv_to< arma::Mat< uint > >::from( mData.mCandK );
            // #else
            // loop over all triangles in J-Direction
            uint tCountJ = 0;
            for ( uint k = 0; k < mData.mNumberOfFacets; ++k )
            {
                // check bounding box in J-direction
                if ( ( aPoint( 1 ) - mData.mFacetMinCoords( k, 1 ) ) * ( mData.mFacetMaxCoords( k, 1 ) - aPoint( 1 ) ) > -gSDFepsilon )
                {
                    // remember this triangle
                    mData.mCandJ( tCountJ ) = k;

                    // increment counter
                    ++tCountJ;
                }
            }

            // counter for triangles
            uint tCount = 0;

            // reset candidate size
            mData.mCandidateFacets.resize( mData.mNumberOfFacets, 1 );

            // loop over remaining triangles in I-direction
            for ( uint k = 0; k < tCountJ; ++k )
            {
                // check bounding box in I-direction
                if ( ( aPoint( 0 ) - mData.mFacetMinCoords( mData.mCandJ( k ), 0 ) ) * ( mData.mFacetMaxCoords( mData.mCandJ( k ), 0 ) - aPoint( 0 ) ) > -gSDFepsilon )
                {
                    mData.mCandidateFacets( tCount ) = mData.mCandJ( k );
                    ++tCount;
                }
            }

            mData.mCandidateFacets.resize( tCount, 1 );
            // #endif
        }

        //-------------------------------------------------------------------------------

        void
        Core::preselect_lines(
                const uint              aDimension,
                const Matrix< DDRMat >& aPoint )
        {
            uint tCandidateCount = 0;
            // loop over all lines in the aAxis direction
            for ( uint iLineIndex = 0; iLineIndex < mData.mNumberOfFacets; iLineIndex++ )
            {
                // check bounding box of the line against the point (point is above min coord and below max coord)
                if ( ( mData.mFacetMaxCoords( iLineIndex, aDimension ) - aPoint( 0 ) ) * ( aPoint( 0 ) - mData.mFacetMinCoords( iLineIndex, aDimension ) )
                        < gSDFepsilon )
                {
                    // add line to candidate list
                    mData.mCandidateFacets( tCandidateCount ) = iLineIndex;
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::intersect_triangles( const uint aAxis, const Matrix< DDRMat >& aPoint )
        {
            // get number of candidate triangles
            uint tNumberOfFacets = mData.mCandidateFacets.length();

            // initialize counter for intersected triangles
            uint tCount = 0;

            // loop over all candidates
            for ( uint k = 0; k < tNumberOfFacets; ++k )
            {
                // get pointer to triangle
                Facet* tFacet = mData.mFacets( mData.mCandidateFacets( k ) );

                if ( tFacet->check_edge( 0, aAxis, aPoint ) )
                {
                    if ( tFacet->check_edge( 1, aAxis, aPoint ) )
                    {
                        if ( tFacet->check_edge( 2, aAxis, aPoint ) )
                        {
                            tFacet->flag();
                            ++tCount;
                        }
                    }
                }
            }

            // resize container with intersected triangles
            mData.mIntersectedTriangles.resize( tCount, nullptr );

            // reset counter
            tCount = 0;

            // loop over all candidates
            for ( uint k = 0; k < tNumberOfFacets; ++k )
            {
                // get pointer to triangle
                Facet* tFacet = mData.mFacets( mData.mCandidateFacets( k ) );

                if ( tFacet->is_flagged() )
                {
                    // add triangle to list
                    mData.mIntersectedTriangles( tCount++ ) = tFacet;

                    // unflag triangle
                    tFacet->unflag();
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::intersect_ray_with_triangles(
                const uint              aAxis,
                const Matrix< DDRMat >& aPoint,
                const uint              aNodeIndex )
        {
            // get number of triangles
            uint tNumberOfFacets = mData.mIntersectedTriangles.size();

            // initialize vector with coords in axis
            Matrix< DDRMat > tCoordsK( tNumberOfFacets, 1 );

            uint tCount = 0;

            bool tError;
            // loop over all intersected triangles and find intersection point
            for ( uint k = 0; k < tNumberOfFacets; ++k )
            {

                real tCoordK;

                // calculate intersection coordinate
                mData.mIntersectedTriangles( k )->intersect_with_coordinate_axis(
                        aPoint,
                        aAxis,
                        tCoordK,
                        tError );

                // tCoordK = std::max( std::min( tCoordK,  tMaxCoord ), tMinCoord );

                // error meant we would have divided by zero. This triangle is ignored
                // otherwise, the value is written into the result vector

                if ( !tError )
                {
                    tCoordsK( tCount++ ) = std::round( tCoordK / gSDFepsilon ) * gSDFepsilon;
                }
                else
                {
                    break;
                }
            }

            if ( tError )
            {
                // this way, the matrix is ignored
                mData.mCoordsK.set_size( 1, 1, 0.0 );
            }
            else
            {
                // resize coord array
                tCoordsK.resize( tCount, 1 );

                // sort array
                Matrix< DDRMat > tCoordsKSorted;
                sort( tCoordsK, tCoordsKSorted );

                // make result unique
                uint tCountUnique = 0;

                // set size of output array
                mData.mCoordsK.set_size( tCount, 1 );

                real tMinCoord = mMesh.get_min_coord( aAxis );
                real tMaxCoord = mMesh.get_max_coord( aAxis );

                // set first entry
                if ( tMinCoord < tCoordsKSorted( 0 ) )
                {
                    mData.mCoordsK( tCountUnique++ ) = tCoordsKSorted( 0 );
                }

                // find unique entries
                for ( uint k = 1; k < tCount; ++k )
                {
                    if ( tCoordsKSorted( k ) > tMinCoord && tCoordsKSorted( k ) < tMaxCoord )
                    {

                        if ( std::abs( tCoordsKSorted( k ) - tCoordsKSorted( k - 1 ) ) > 10 * gSDFepsilon )
                        {
                            mData.mCoordsK( tCountUnique++ ) = tCoordsKSorted( k );
                        }
                    }
                }

                // chop vector
                mData.mCoordsK.resize( tCountUnique, 1 );
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::check_if_node_is_inside(
                const uint aAxis,
                const uint aNodeIndex )
        {
            uint tNumCoordsK = mData.mCoordsK.length();

            bool tNodeIsInside = false;

            const Matrix< F31RMat >& aPoint = mMesh.get_node_coordinate( aNodeIndex );

            // only even number of intersections is considered
            if ( tNumCoordsK % 2 == 0 )
            {
                for ( uint k = 0; k < tNumCoordsK / 2; ++k )
                {
                    tNodeIsInside = ( aPoint( aAxis ) > mData.mCoordsK( 2 * k ) ) && ( aPoint( aAxis ) < mData.mCoordsK( 2 * k + 1 ) );

                    // break the loop if inside
                    if ( tNodeIsInside )
                    {
                        break;
                    }
                }

                // set the inside flag of this node to the corresponding value
                if ( tNodeIsInside )
                {
                    mMesh.get_vertex( aNodeIndex )->set_inside_flag();
                }
                else
                {
                    mMesh.get_vertex( aNodeIndex )->unset_inside_flag();
                }
                mMesh.get_vertex( aNodeIndex )->unflag();
            }
            else
            {
                // set unsure flag
                mMesh.get_vertex( aNodeIndex )->flag();

                // increment counter
                ++mData.mUnsureNodesCount;
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::calculate_candidate_points_and_buffer_diagonal()
        {
            // get number of elements
            uint tNumberOfElements = mMesh.get_num_elems();

            // counter for elements near surface
            mData.mSurfaceElements = 0;

            // counter for elements in volume
            mData.mVolumeElements = 0;

            // reset buffer diagonal
            mData.mBufferDiagonal = 0;

            // search all elements for sign change
            for ( uint e = 0; e < tNumberOfElements; ++e )
            {
                // unflag this element

                Cell* tElement = mMesh.get_cell( e );

                // reset flags of this element
                tElement->unflag();
                tElement->unset_surface_flag();
                tElement->unset_volume_flag();

                // get pointer to nodes
                const moris::Cell< Vertex* > tNodes = tElement->get_vertices();

                // get number of nodes
                uint tNumberOfNodes = tNodes.size();

                // get first sign
                bool tIsInside = tNodes( 0 )->is_inside();

                // assume element is not intersected
                bool tIsIntersected = false;

                // loop over all other nodes
                for ( uint k = 1; k < tNumberOfNodes; ++k )
                {
                    // check of sign is the same
                    if ( tNodes( k )->is_inside() != tIsInside )
                    {
                        // sign is not same
                        tIsIntersected = true;

                        // cancel loop
                        break;
                    }
                }

                // test if there is a sign change
                if ( tIsIntersected )
                {
                    // flag this element as surface element
                    tElement->set_surface_flag();
                    tElement->unset_volume_flag();

                    // increment counter
                    ++mData.mSurfaceElements;

                    // update buffer diagonal
                    mData.mBufferDiagonal = std::max(
                            mData.mBufferDiagonal,
                            tElement->get_buffer_diagonal() );

                    // flag to indicate that the buffer of this element
                    // has been calculated
                    tElement->flag();

                    // flag all nodes of this element as candidates
                    for ( uint k = 0; k < tNumberOfNodes; ++k )
                    {
                        tNodes( k )->set_candidate_flag();
                    }
                }
                else if ( tIsInside )
                {
                    // flag this element as volume element
                    tElement->unset_surface_flag();
                    tElement->set_volume_flag();

                    // increment counter
                    ++mData.mVolumeElements;
                }
                else
                {
                    // unflag element
                    tElement->unset_surface_flag();
                    tElement->unset_volume_flag();
                }
            }

            // add additional search depth
            for ( uint d = 1; d < mCandidateSearchDepth; ++d )
            {
                // loop over all elements
                for ( uint e = 0; e < tNumberOfElements; ++e )
                {
                    // get pointer to element
                    Cell* tElement = mMesh.get_cell( e );

                    // test if element is not flagged
                    if ( !tElement->is_flagged() )
                    {
                        // get pointer to nodes
                        const moris::Cell< Vertex* > tNodes = tElement->get_vertices();

                        // get number of nodes
                        uint tNumberOfNodes = tNodes.size();

                        bool tIsCandidate = false;

                        // test if at least one node of this element is flagged
                        // as candidate
                        for ( uint k = 0; k < tNumberOfNodes; ++k )
                        {
                            if ( tNodes( k )->is_candidate() )
                            {
                                tIsCandidate = true;
                                break;
                            }
                        }

                        // test if candidtae flag is set
                        if ( tIsCandidate )
                        {
                            // update buffer diagonal
                            mData.mBufferDiagonal = std::max(
                                    mData.mBufferDiagonal,
                                    tElement->get_buffer_diagonal() );

                            // flag this element
                            tElement->flag();

                            // flag all nodes of this element
                            for ( uint k = 0; k < tNumberOfNodes; ++k )
                            {
                                tNodes( k )->set_candidate_flag();
                            }
                        }
                    }    // end loop over all elements
                }
            }    // end candidate search depth loop
        }

        //-------------------------------------------------------------------------------

        moris::Cell< Vertex* >
        Core::set_candidate_list()
        {

            uint tNumberOfNodes = mMesh.get_num_nodes();
            //        	std::cout<<"number of nodes in mesh   : "<<tNumberOfNodes<<std::endl;
            moris::Cell< Vertex* > tCandidateVertices;

            for ( uint k = 0; k < tNumberOfNodes; k++ )
            {
                Vertex* tNode = mMesh.get_vertex( k );

                if ( tNode->is_candidate() )
                {
                    tCandidateVertices.push_back( tNode );
                }
                else
                {
                    continue;
                }
            }
            //        	std::cout<<"number of candidate nodes : "<<tCandidateVertices.size()<<std::endl;
            return tCandidateVertices;
        }

        //-------------------------------------------------------------------------------

        void
        Core::get_nodes_withing_bounding_box_of_triangle(
                Facet*                  aFacet,
                moris::Cell< Vertex* >& aNodes,
                moris::Cell< Vertex* >& aCandList )    //===========================================
        {
            // calculate minimum and maximum coordinate

            Matrix< F31RMat > tMinCoord( 3, 1 );
            Matrix< F31RMat > tMaxCoord( 3, 1 );

            for ( uint i = 0; i < 3; ++i )
            {
                tMinCoord( i ) = aFacet->get_min_coord( i ) - mData.mBufferDiagonal;
                tMaxCoord( i ) = aFacet->get_max_coord( i ) + mData.mBufferDiagonal;
            }

            // why is this necessary?

            for ( uint i = 0; i < 3; ++i )
            {
                tMinCoord( i ) = std::max( tMinCoord( i ), mMesh.get_min_coord( i ) );
                tMaxCoord( i ) = std::min( tMaxCoord( i ), mMesh.get_max_coord( i ) );
            }

            //            // get number of nodes on this mesh
            //            uint tNumberOfNodes = mMesh.get_num_nodes();

            // number of candidate nodes
            uint tNumberOfCandidates = aCandList.size();    //========================================

            // node counter
            uint tCount = 0;

            // loop over all nodes of this mesh
            //            for( uint k=0; k<tNumberOfNodes; ++k )
            //            {

            // loop over only the candidate nodes
            for ( uint k = 0; k < tNumberOfCandidates; k++ )    //========================================
            {
                // get pointer to node
                //                Vertex * tNode = mMesh.get_vertex( k );
                Vertex* tNode = aCandList( k );    //=================================================

                // unflag this node
                tNode->unflag();

                // test if node is a candidate
                //                if( tNode->is_candidate() )
                //                {
                // get coords of this node
                const Matrix< F31RMat >& tPoint = tNode->get_coords();

                // assume that node is in triangle
                bool tNodeIsWithinTriangle = true;

                for ( uint i = 0; i < 3; ++i )
                {
                    if ( tPoint( i ) < tMinCoord( i ) || tPoint( i ) > tMaxCoord( i ) )
                    {
                        tNodeIsWithinTriangle = false;
                        break;
                    }
                    /* tNodeIsWithinTriangle = tNodeIsWithinTriangle
                            && ( tPoint( i ) >= tMinCoord( i ) )
                            && ( tPoint( i ) <= tMaxCoord( i ) );
                    if( ! tNodeIsWithinTriangle )
                    {
                        break;
                    } */
                }

                // if node is in triangle
                if ( tNodeIsWithinTriangle )
                {
                    // flag this node
                    tNode->flag();

                    // increment counter
                    ++tCount;
                }
                //                } // end node is candidate
            }    // end loop over all nodes

            // reset output array
            aNodes.resize( tCount, nullptr );

            // reset counter
            tCount = 0;

            // loop over all nodes of this mesh
            //            for( uint k=0; k<tNumberOfNodes; ++k )
            //            {

            // loop over only the candidate nodes+
            for ( uint k = 0; k < tNumberOfCandidates; k++ )    //========================================
            {
                // get pointer to node
                //                Vertex * tNode = mMesh.get_vertex( k );
                Vertex* tNode = aCandList( k );    //==================================================

                // test if node is flagged
                if ( tNode->is_flagged() )
                {
                    aNodes( tCount++ ) = tNode;
                }
            }
        }

        //-------------------------------------------------------------------------------

        void
        Core::sweep()
        {

            uint tNumberOfVertices = mMesh.get_num_nodes();

            uint tSweepCount = 1;

            while ( tSweepCount != 0 )
            {
                tic tTimer;
                tSweepCount = 0;

                // loop over all vertices
                for ( uint k = 0; k < tNumberOfVertices; ++k )
                {
                    // get vertex
                    Vertex* tVertex = mMesh.get_vertex( k );

                    // test if node has sdf
                    if ( tVertex->has_sdf() )
                    {
                        // sweep this vertex
                        tSweepCount += tVertex->sweep();
                    }
                }

                if ( mVerbose )
                {
                    // stop the timer
                    real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

                    // print elapsed time
                    if ( par_size() == 1 )
                    {
                        std::fprintf( stdout, "Time for sweeping              : %5.3f [sec]\nSwept %i nodes\n", tElapsedTime / 1000.0, (int)tSweepCount );
                    }
                    else
                    {
                        std::fprintf( stdout, "Proc % i - Time for sweeping              : %5.3f [sec]\nSwept %i nodes\n", (int)par_rank(), tElapsedTime / 1000.0, (int)tSweepCount );
                    }
                }
            }
        }

        // -----------------------------------------------------------------------------

        void
        Core::fill_sdf_with_values( Matrix< DDRMat >& aSDF )
        {
            // get number of vertices
            uint tNumberOfVertices = mMesh.get_num_nodes();

            // min and max value
            real tMinSDF = -1e-12;    // std::numeric_limits<real>::max();
            real tMaxSDF = 1e-12;     // std::numeric_limits<real>::min();

            // allocate matrix
            aSDF.set_size( tNumberOfVertices, 1 );

            // loop over all nodes and write real values
            for ( uint k = 0; k < tNumberOfVertices; ++k )
            {
                // get pointer to vertex
                Vertex* tVertex = mMesh.get_vertex( k );

                // test if vertex has SDF
                if ( tVertex->has_sdf() )
                {
                    real tSDF = tVertex->get_sdf();
                    // write value
                    aSDF( tVertex->get_index() ) = tSDF;

                    if ( tVertex->is_inside() )
                    {
                        tMinSDF = std::min( tMinSDF, tSDF );
                    }
                    else
                    {
                        tMaxSDF = std::max( tMaxSDF, tSDF );
                    }
                }
            }

            // if parallel, synchronize min and max values for SDF
            if ( par_size() > 1 )
            {
                // container for min and max values
                Matrix< DDRMat > tValues;

                // communicate minimal value
                comm_gather_and_broadcast( tMinSDF, tValues );
                tMinSDF = tValues.min();

                // communicate maximal value
                comm_gather_and_broadcast( tMaxSDF, tValues );
                tMaxSDF = tValues.max();
            }

            // loop over all nodes and write fake values
            for ( uint k = 0; k < tNumberOfVertices; ++k )
            {
                // get pointer to vertex
                Vertex* tVertex = mMesh.get_vertex( k );

                // test if vertex does not have an SDF
                if ( !tVertex->has_sdf() )
                {
                    if ( tVertex->is_inside() )
                    {
                        aSDF( tVertex->get_index() ) = tMinSDF;
                    }
                    else
                    {
                        aSDF( tVertex->get_index() ) = tMaxSDF;
                    }
                }
            }
        }

        // -----------------------------------------------------------------------------

        void
        Core::save_to_vtk( const std::string& aFilePath )
        {
            // open the file
            std::ofstream tFile( aFilePath, std::ios::binary );

            // containers
            float tFChar = 0;
            int   tIChar = 0;

            Matrix< DDRMat > tSDF;

            this->fill_sdf_with_values( tSDF );

            tFile << "# vtk DataFile Version 3.0" << std::endl;
            tFile << "GO BUFFS!" << std::endl;
            tFile << "BINARY" << std::endl;
            // tFile << "ASCII" << std::endl;
            uint tNumberOfNodes = mMesh.get_num_nodes();

            // write node data
            tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            tFile << "POINTS " << tNumberOfNodes << " float" << std::endl;

            // loop over all nodes
            for ( luint k = 0; k < tNumberOfNodes; ++k )
            {
                // get coordinate from node
                const Matrix< F31RMat >& tCoords = mMesh.get_vertex( k )->get_coords();

                // write coordinates to mesh
                tFChar = swap_byte_endian( (float)tCoords( 0 ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
                tFChar = swap_byte_endian( (float)tCoords( 1 ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
                tFChar = swap_byte_endian( (float)tCoords( 2 ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
                // tFile << tCoords( 0 ) << " " << tCoords( 1 ) << " " << tCoords( 2 ) << std::endl;
            }
            tFile << std::endl;

            uint tNumberOfElements = mMesh.get_num_elems();

            // write header for cells
            tFile << "CELLS " << tNumberOfElements << " "
                  << 9 * tNumberOfElements << std::endl;

            // cell types
            Matrix< DDUMat > tCellTypes( tNumberOfElements, 1 );

            for ( luint k = 0; k < tNumberOfElements; ++k )
            {
                // get pointet to cell
                Cell* tCell = mMesh.get_cell( k );

                // get number of vertices
                uint tNumberOfCellVerts = tCell->get_number_of_vertices();

                // get vertex indices
                Matrix< DDUMat > tIndices( tNumberOfCellVerts, 1 );

                // VTK cell type
                uint tCellType = 0;

                switch ( tNumberOfCellVerts )
                {
                    case ( 4 ):    // TET 4
                    {
                        tCellType = 10;
                        break;
                    }
                    case ( 10 ):    // TET 10
                    {
                        tCellType = 24;
                        break;
                    }
                    case ( 8 ):    // HEX 8
                    {
                        tCellType = 12;
                        break;
                    }
                    case ( 27 ):    // HEX27
                    {
                        tCellType = 29;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "unknown cell type" );
                    }
                }

                // remember cell type
                tCellTypes( k ) = tCellType;

                if ( tCellType == 29 )
                {
                    // special case for HEX27
                    tIndices( 0 )  = tCell->get_vertex( 0 )->get_index();
                    tIndices( 1 )  = tCell->get_vertex( 1 )->get_index();
                    tIndices( 2 )  = tCell->get_vertex( 2 )->get_index();
                    tIndices( 3 )  = tCell->get_vertex( 3 )->get_index();
                    tIndices( 4 )  = tCell->get_vertex( 4 )->get_index();
                    tIndices( 5 )  = tCell->get_vertex( 5 )->get_index();
                    tIndices( 6 )  = tCell->get_vertex( 6 )->get_index();
                    tIndices( 7 )  = tCell->get_vertex( 7 )->get_index();
                    tIndices( 8 )  = tCell->get_vertex( 8 )->get_index();
                    tIndices( 9 )  = tCell->get_vertex( 9 )->get_index();
                    tIndices( 10 ) = tCell->get_vertex( 10 )->get_index();
                    tIndices( 11 ) = tCell->get_vertex( 11 )->get_index();
                    tIndices( 12 ) = tCell->get_vertex( 16 )->get_index();
                    tIndices( 13 ) = tCell->get_vertex( 17 )->get_index();
                    tIndices( 14 ) = tCell->get_vertex( 18 )->get_index();
                    tIndices( 15 ) = tCell->get_vertex( 19 )->get_index();
                    tIndices( 16 ) = tCell->get_vertex( 12 )->get_index();
                    tIndices( 17 ) = tCell->get_vertex( 13 )->get_index();
                    tIndices( 18 ) = tCell->get_vertex( 14 )->get_index();
                    tIndices( 19 ) = tCell->get_vertex( 15 )->get_index();
                    tIndices( 20 ) = tCell->get_vertex( 23 )->get_index();
                    tIndices( 21 ) = tCell->get_vertex( 24 )->get_index();
                    tIndices( 22 ) = tCell->get_vertex( 25 )->get_index();
                    tIndices( 23 ) = tCell->get_vertex( 26 )->get_index();
                    tIndices( 24 ) = tCell->get_vertex( 21 )->get_index();
                    tIndices( 25 ) = tCell->get_vertex( 22 )->get_index();
                    tIndices( 26 ) = tCell->get_vertex( 20 )->get_index();
                }
                else
                {
                    for ( uint i = 0; i < tNumberOfCellVerts; ++i )
                    {
                        tIndices( i ) = tCell->get_vertex( i )->get_index();
                    }
                }

                tIChar = swap_byte_endian( (int)tNumberOfCellVerts );
                tFile.write( (char*)&tIChar, sizeof( int ) );

                // write indices to file
                for ( uint i = 0; i < tNumberOfCellVerts; ++i )
                {
                    tIChar = swap_byte_endian( (int)tIndices( i ) );
                    tFile.write( (char*)&tIChar, sizeof( int ) );
                    // tFile << " " << tIndices( i );
                }
                // tFile << std::endl;
            }
            tFile << std::endl;
            // write cell types
            tFile << "CELL_TYPES " << tNumberOfElements << std::endl;

            for ( luint k = 0; k < tNumberOfElements; ++k )
            {
                tIChar = swap_byte_endian( tCellTypes( k ) );
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
            tFile << std::endl;

            /*
            // write element data
            tFile << "CELL_DATA " << tNumberOfElements << std::endl;
            tFile << "SCALARS ELEMENT_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k <  tNumberOfElements; ++k)
            {
                tIChar = swap_byte_endian( (int) mMesh.get_cell( k )->get_id() );
            }
            tFile << std::endl;

            tFile << "SCALARS ELEMENT_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k <  tNumberOfElements; ++k)
            {
                tIChar = swap_byte_endian( (int) mMesh.get_cell( k )->get_index() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            tFile << "SCALARS ELEMENT_IN_VOLUME int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k <  tNumberOfElements; ++k)
            {
                if( mMesh.get_cell( k )->is_in_volume() )
                {
                    tIChar = swap_byte_endian( (int) 1 );
                }
                else
                {
                    tIChar = swap_byte_endian( (int) 0 );
                }
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            tFile << "SCALARS ELEMENT_ON_SURFACE int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k <  tNumberOfElements; ++k)
            {
                if( mMesh.get_cell( k )->is_on_surface() )
                {
                    tIChar = swap_byte_endian( (int) 1 );
                }
                else
                {
                    tIChar = swap_byte_endian( (int) 0 );
                }
                tFile.write( (char*) &tIChar, sizeof(int));
            } */

            tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

            tFile << "SCALARS SDF float" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {

                tFChar = swap_byte_endian( (float)tSDF( k ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
            }

            tFile << "SCALARS HAS_SDF int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                if ( mMesh.get_vertex( k )->has_sdf() )
                {
                    tIChar = swap_byte_endian( (int)1 );
                }
                else
                {
                    tIChar = swap_byte_endian( (int)0 );
                }
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }

            tFile << "SCALARS IS_INSIDE int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                if ( mMesh.get_vertex( k )->is_inside() )
                {
                    tIChar = swap_byte_endian( (int)1 );
                }
                else
                {
                    tIChar = swap_byte_endian( (int)0 );
                }
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }

            tFile << "SCALARS VERTEX_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                tIChar = swap_byte_endian( (int)mMesh.get_vertex( k )->get_index() );
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }

            // close the output file
            tFile.close();
        }
        // -----------------------------------------------------------------------------

        void
        Core::save_unsure_to_vtk( const std::string& aFilePath )
        {
            // open the file
            std::ofstream tFile( aFilePath, std::ios::binary );

            // containers
            float tFChar = 0;
            int   tIChar = 0;

            Matrix< DDRMat > tSDF;

            this->fill_sdf_with_values( tSDF );

            tFile << "# vtk DataFile Version 3.0" << std::endl;
            tFile << "GO BUFFS!" << std::endl;
            tFile << "BINARY" << std::endl;
            // tFile << "ASCII" << std::endl;
            uint tNumberOfNodes = mMesh.get_num_nodes();

            // write node data
            tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            tFile << "POINTS " << tNumberOfNodes << " float" << std::endl;

            // loop over all nodes
            for ( luint k = 0; k < tNumberOfNodes; ++k )
            {
                // get coordinate from node
                const Matrix< F31RMat >& tCoords = mMesh.get_vertex( k )->get_coords();

                // write coordinates to mesh
                tFChar = swap_byte_endian( (float)tCoords( 0 ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
                tFChar = swap_byte_endian( (float)tCoords( 1 ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
                tFChar = swap_byte_endian( (float)tCoords( 2 ) );
                tFile.write( (char*)&tFChar, sizeof( float ) );
                // tFile << tCoords( 0 ) << " " << tCoords( 1 ) << " " << tCoords( 2 ) << std::endl;
            }
            tFile << std::endl;

            uint tNumberOfElements = mMesh.get_num_elems();

            // write header for cells
            tFile << "CELLS " << tNumberOfElements << " "
                  << 9 * tNumberOfElements << std::endl;

            // cell types
            Matrix< DDUMat > tCellTypes( tNumberOfElements, 1 );

            for ( luint k = 0; k < tNumberOfElements; ++k )
            {
                // get pointet to cell
                Cell* tCell = mMesh.get_cell( k );

                // get number of vertices
                uint tNumberOfCellVerts = tCell->get_number_of_vertices();

                // get vertex indices
                Matrix< DDUMat > tIndices( tNumberOfCellVerts, 1 );

                // VTK cell type
                uint tCellType = 0;

                switch ( tNumberOfCellVerts )
                {
                    case ( 4 ):    // TET 4
                    {
                        tCellType = 10;
                        break;
                    }
                    case ( 10 ):    // TET 10
                    {
                        tCellType = 24;
                        break;
                    }
                    case ( 8 ):    // HEX 8
                    {
                        tCellType = 12;
                        break;
                    }
                    case ( 27 ):    // HEX27
                    {
                        tCellType = 29;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "unknown cell type" );
                    }
                }

                // remember cell type
                tCellTypes( k ) = tCellType;

                if ( tCellType == 29 )
                {
                    // special case for HEX27
                    tIndices( 0 )  = tCell->get_vertex( 0 )->get_index();
                    tIndices( 1 )  = tCell->get_vertex( 1 )->get_index();
                    tIndices( 2 )  = tCell->get_vertex( 2 )->get_index();
                    tIndices( 3 )  = tCell->get_vertex( 3 )->get_index();
                    tIndices( 4 )  = tCell->get_vertex( 4 )->get_index();
                    tIndices( 5 )  = tCell->get_vertex( 5 )->get_index();
                    tIndices( 6 )  = tCell->get_vertex( 6 )->get_index();
                    tIndices( 7 )  = tCell->get_vertex( 7 )->get_index();
                    tIndices( 8 )  = tCell->get_vertex( 8 )->get_index();
                    tIndices( 9 )  = tCell->get_vertex( 9 )->get_index();
                    tIndices( 10 ) = tCell->get_vertex( 10 )->get_index();
                    tIndices( 11 ) = tCell->get_vertex( 11 )->get_index();
                    tIndices( 12 ) = tCell->get_vertex( 16 )->get_index();
                    tIndices( 13 ) = tCell->get_vertex( 17 )->get_index();
                    tIndices( 14 ) = tCell->get_vertex( 18 )->get_index();
                    tIndices( 15 ) = tCell->get_vertex( 19 )->get_index();
                    tIndices( 16 ) = tCell->get_vertex( 12 )->get_index();
                    tIndices( 17 ) = tCell->get_vertex( 13 )->get_index();
                    tIndices( 18 ) = tCell->get_vertex( 14 )->get_index();
                    tIndices( 19 ) = tCell->get_vertex( 15 )->get_index();
                    tIndices( 20 ) = tCell->get_vertex( 23 )->get_index();
                    tIndices( 21 ) = tCell->get_vertex( 24 )->get_index();
                    tIndices( 22 ) = tCell->get_vertex( 25 )->get_index();
                    tIndices( 23 ) = tCell->get_vertex( 26 )->get_index();
                    tIndices( 24 ) = tCell->get_vertex( 21 )->get_index();
                    tIndices( 25 ) = tCell->get_vertex( 22 )->get_index();
                    tIndices( 26 ) = tCell->get_vertex( 20 )->get_index();
                }
                else
                {
                    for ( uint i = 0; i < tNumberOfCellVerts; ++i )
                    {
                        tIndices( i ) = tCell->get_vertex( i )->get_index();
                    }
                }

                tIChar = swap_byte_endian( (int)tNumberOfCellVerts );
                tFile.write( (char*)&tIChar, sizeof( int ) );

                // write indices to file
                for ( uint i = 0; i < tNumberOfCellVerts; ++i )
                {
                    tIChar = swap_byte_endian( (int)tIndices( i ) );
                    tFile.write( (char*)&tIChar, sizeof( int ) );
                    // tFile << " " << tIndices( i );
                }
                // tFile << std::endl;
            }
            tFile << std::endl;
            // write cell types
            tFile << "CELL_TYPES " << tNumberOfElements << std::endl;

            for ( luint k = 0; k < tNumberOfElements; ++k )
            {
                tIChar = swap_byte_endian( tCellTypes( k ) );
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
            tFile << std::endl;

            tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

            tFile << "SCALARS RAYCAST int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                // test if vertex is determined
                if ( mMesh.get_vertex( k )->is_flagged() )
                {
                    tIChar = swap_byte_endian( (int)0 );
                }
                else if ( mMesh.get_vertex( k )->is_inside() )
                {
                    tIChar = swap_byte_endian( (int)-1 );
                }
                else
                {
                    tIChar = swap_byte_endian( (int)1 );
                }
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }

            tFile << "SCALARS VERTEX_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                tIChar = swap_byte_endian( (int)mMesh.get_vertex( k )->get_index() );
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }

            // close the output file
            tFile.close();
        }

        // -----------------------------------------------------------------------------

        void
        Core::force_unsure_nodes_outside()
        {
            // get number of nodes on mesh
            uint tNumberOfNodes = mMesh.get_num_nodes();

            // loop over all nodes
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                // get pointer to node
                mMesh.get_vertex( k )->unflag();
            }

            mData.mUnsureNodesCount = 0;
        }

        // -----------------------------------------------------------------------------

        void
        Core::random_rotation()
        {
            // determine the required dimensionality of the rotation (Number of entries in a vertex's coordinates)
            uint tNumDim = mData.mVertices( 0 )->get_coords().numel();

            // generate random angle
            real tAngle = random_angle();

            // generate random rotation matrix
            Matrix< DDRMat > tRotation;
            if ( tNumDim == 2 )
            {
                tRotation = rotation_matrix( tAngle );
            }
            else
            {
                // create random axis for cases larger than 2 dimensions
                Matrix< DDRMat > tAxis = random_axis( tNumDim );
                tRotation              = rotation_matrix( tAxis, tAngle );
            }

            // rotate all vertices of triangle mesh
            for ( Facet_Vertex* tVertex : mData.mVertices )
            {
                tVertex->rotate_node_coords( tRotation );
            }

            // update all triangles
            for ( Facet* tFacet : mData.mFacets )
            {
                tFacet->update_data();
            }

            // rotate unsure nodes
            uint tNumberOfNodes = mMesh.get_num_nodes();

            // loop over all nodes
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                // test if node is unsure
                if ( mMesh.get_vertex( k )->is_flagged() )
                {
                    mMesh.get_vertex( k )->rotate_coords( tRotation );
                }
            }
        }

        // -----------------------------------------------------------------------------

        void
        Core::undo_rotation()
        {
            // rotate all vertices of triangle mesh
            for ( Facet_Vertex* tVertex : mData.mVertices )
            {
                tVertex->reset_node_coords();
            }

            // update all triangles
            for ( Facet* tFacet : mData.mFacets )
            {
                tFacet->update_data();
            }

            // rotate unsure nodes
            uint tNumberOfNodes = mMesh.get_num_nodes();

            // loop over all nodes
            for ( uint k = 0; k < tNumberOfNodes; ++k )
            {
                mMesh.get_vertex( k )->reset_coords();
            }
        }

        // -----------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */
