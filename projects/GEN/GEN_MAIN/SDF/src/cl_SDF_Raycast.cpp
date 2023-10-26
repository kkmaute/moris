/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_SDF_Raycaster.cpp
 *
 */

#include "fn_sort.hpp"


#include "cl_SDF_Raycast.hpp"

namespace moris::sdf
{
    Raycast::Raycast(
            Object& aObject,
            Data&   aData )
            : mObject( aObject )
            , mData( aData )
            , mDimension( aObject.get_dimension() )
            , mPoint( mDimension, 1 )
            , mOriginalPoint( mDimension, 1 )
            , mFacetMinCoords( aObject.get_num_facets(), aObject.get_dimension() )
            , mFacetMaxCoords( aObject.get_num_facets(), aObject.get_dimension() )
            , mPointIsInside( 2 )
    {
        // Get facet information from the object
        Cell< Facet* > tFacets         = mObject.get_facets();
        uint           tNumberOfFacets = tFacets.size();

        // Determine and store the minimum and maximum coordinates of each facet
        for ( uint iFacetIndex = 0; iFacetIndex < tNumberOfFacets; iFacetIndex++ )
        {
            for ( uint iDimensionIndex = 0; iDimensionIndex < aObject.get_dimension(); iDimensionIndex++ )
            {
                mFacetMinCoords( iFacetIndex, iDimensionIndex ) = tFacets( iFacetIndex )->get_min_coord( iDimensionIndex );

                mFacetMaxCoords( iFacetIndex, iDimensionIndex ) = tFacets( iFacetIndex )->get_max_coord( iDimensionIndex );
            }
        }
    }

    void
    Raycast::calculate_raycast( const Matrix< DDRMat >& aPoint )
    {
        MORIS_ASSERT( aPoint.numel() == mDimension, "SDF_Raycast::calculate_raycast() - mDimension must match the number of coordinates of mPoint." );

        // copy point data
        // FIXME: mPoint and mOriginalPoint could be a reference to aPoint since it will exist elsewhere to avoid the copy
        mOriginalPoint = mPoint = aPoint;

        // flag that marks if rotation was called
        bool tRotation = false;

        switch ( mDimension )
        {
            case 2:
            {
                while ( mPointIsInside == 2 )
                {
                    // perform voxelizing algorithm in z-direction
                    voxelize_2D( 1 );
                    if ( mPointIsInside == 2 )
                    {
                        // perform voxelizing algorithm in y-direction
                        voxelize_2D( 0 );
                    }

                    if ( mPointIsInside == 2 )
                    {
                        tRotation = true;

                        this->random_rotation();
                    }
                }
                break;
            }
            case 3:
            {
                while ( mPointIsInside == 2 )
                {
                    // perform voxelizing algorithm in z-direction
                    voxelize_3D( 2 );
                    if ( mPointIsInside == 2 )
                    {
                        // perform voxelizing algorithm in y-direction
                        voxelize_3D( 1 );
                        if ( mPointIsInside == 2 )
                        {
                            // perform voxelizing algorithm in x-direction
                            voxelize_3D( 0 );
                        }
                    }

                    if ( mPointIsInside == 2 )
                    {
                        tRotation = true;

                        this->random_rotation();
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "SDF_Raycast::calculate_raycast() functionality not implemented for %dD problems.", mDimension );
            }
        }

        if ( tRotation )
        {
            this->undo_rotation();
        }

        // Force point outside if still unsure
        mPointIsInside = 0;
    }

    void
    Raycast::voxelize_2D( uint aAxis )
    {
        // preselect lines in the aAxis direction
        this->preselect_lines( aAxis );

        // compute intersection if the point is inside a line's bounding box
        if ( mData.mIntersectedFacets.size() > 0 )
        {
            this->intersect_ray_with_facets( aAxis );
        }

        // check if the node is inside the polygon
        this->check_if_node_is_inside_lines( aAxis );
    }


    void
    Raycast::voxelize_3D( uint aAxis )
    {
        // preselect triangles for intersection test
        if ( aAxis == 0 )
            this->preselect_triangles_x();
        else if ( aAxis == 1 )
            this->preselect_triangles_y();
        else
            this->preselect_triangles_z();

        // from the candidate triangles, perform intersection
        if ( mData.mCandidateFacets.size() > 0 )
        {
            this->intersect_triangles( aAxis );

            // intersect ray with triangles and check if node is inside
            if ( mData.mIntersectedFacets.size() > 0 )
            {
                this->intersect_ray_with_facets( aAxis );

                this->check_if_node_is_inside_triangles( aAxis );
            }
        }
    }

    void
    Raycast::preselect_triangles_x()
    {
        // x: k = x, j = z, i = y
        // #ifdef MORIS_USE_ARMA

        //             // check bounding box in J-direction
        //             mData.mCandJ = arma::find(
        //                     ( mPoint( 2 ) - mData.mFacetMinCoordsZ ) % ( mData.mFacetMaxCoordsZ - mPoint( 2 ) ) > -gSDFepsilon );

        //             // check bounding box in I-direction
        //             mData.mCandI = arma::find(
        //                     ( mPoint( 1 ) - mData.mFacetMinCoordsY.elem( mData.mCandJ ) ) % ( mData.mFacetMaxCoordsY.elem( mData.mCandJ ) - mPoint( 1 ) ) > -gSDFepsilon );

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
            if ( ( mPoint( 2 ) - mFacetMinCoords( k, 2 ) ) * ( mFacetMaxCoords( k, 2 ) - mPoint( 2 ) ) > -gSDFepsilon )
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
        mData.mCandidateFacets.resize( mData.mNumberOfFacets );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // check bounding box in I-direction
            if ( ( mPoint( 1 ) - mFacetMinCoords( mData.mCandJ( k ), 1 ) ) * ( mFacetMaxCoords( mData.mCandJ( k ), 1 ) - mPoint( 1 ) ) > -gSDFepsilon )
            {
                mData.mCandidateFacets( tCount ) = mData.mCandJ( k );
                ++tCount;
            }
        }

        mData.mCandidateFacets.resize( tCount );
        // #endif
    }

    //-------------------------------------------------------------------------------

    void
    Raycast::preselect_triangles_y()
    {
        // y: k = y, j = x, i = z
        // #ifdef MORIS_USE_ARMA
        //             // check bounding box in J-direction
        //             mData.mCandJ = arma::find(
        //                     ( mPoint( 0 ) - mData.mFacetMinCoordsX ) % ( mData.mFacetMaxCoordsX - mPoint( 0 ) ) > -gSDFepsilon );

        //             // check bounding box in I-direction
        //             mData.mCandI = arma::find(
        //                     ( mPoint( 2 ) - mData.mFacetMinCoordsZ.elem( mData.mCandJ ) ) % ( mData.mFacetMaxCoordsZ.elem( mData.mCandJ ) - mPoint( 2 ) ) > -gSDFepsilon );

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
            if ( ( mPoint( 0 ) - mFacetMinCoords( k, 0 ) ) * ( mFacetMaxCoords( k, 0 ) - mPoint( 0 ) ) > -gSDFepsilon )
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
        mData.mCandidateFacets.resize( mData.mNumberOfFacets );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // check bounding box in I-direction
            if ( ( mPoint( 2 ) - mFacetMinCoords( mData.mCandJ( k ), 2 ) ) * ( mFacetMaxCoords( mData.mCandJ( k ), 2 ) - mPoint( 2 ) ) > -gSDFepsilon )
            {
                mData.mCandidateFacets( tCount ) = mData.mCandJ( k );
                ++tCount;
            }
        }

        mData.mCandidateFacets.resize( tCount );
        // #endif
    }

    //-------------------------------------------------------------------------------

    void
    Raycast::preselect_triangles_z()
    {
        // z: k = z, j = y, i = x
        // #ifdef MORIS_USE_ARMA

        //             // bool_t tNothingFound = true;

        //             // check bounding box in J-direction
        //             mData.mCandJ = arma::find(
        //                     ( mPoint( 1 ) - mData.mFacetMinCoordsY ) % ( mData.mFacetMaxCoordsY - mPoint( 1 ) ) > -gSDFepsilon );
        //             // check bounding box in I-direction
        //             mData.mCandI = arma::find(
        //                     ( mPoint( 0 ) - mData.mFacetMinCoordsX.elem( mData.mCandJ ) ) % ( mData.mFacetMaxCoordsX.elem( mData.mCandJ ) - mPoint( 0 ) ) > -gSDFepsilon );

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
            if ( ( mPoint( 1 ) - mFacetMinCoords( k, 1 ) ) * ( mFacetMaxCoords( k, 1 ) - mPoint( 1 ) ) > -gSDFepsilon )
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
        mData.mCandidateFacets.resize( mData.mNumberOfFacets );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // check bounding box in I-direction
            if ( ( mPoint( 0 ) - mFacetMinCoords( mData.mCandJ( k ), 0 ) ) * ( mFacetMaxCoords( mData.mCandJ( k ), 0 ) - mPoint( 0 ) ) > -gSDFepsilon )
            {
                mData.mCandidateFacets( tCount ) = mData.mCandJ( k );
                ++tCount;
            }
        }

        mData.mCandidateFacets.resize( tCount );
        // #endif
    }

    //-------------------------------------------------------------------------------

    void
    Raycast::preselect_lines( uint aAxis )
    {
        moris::print( mPoint, "mPoint" );
        
        // Ensure the function is being called for the proper number of facets
        MORIS_ERROR( aAxis < mPoint.numel(),
                "SDF_ preselect_lines() mPoint is %lu while coordinate axis of %d specified.",
                mPoint.numel(),
                aAxis );
        MORIS_ASSERT( mPoint.numel() == 2,
                "SDF_ preselect_lines() should be called for 2D problems only. Query point dimension = %lu",
                mPoint.numel() );

        // reset candidate and intersected facet size
        mData.mCandidateFacets.resize( mData.mNumberOfFacets );
        mData.mIntersectedFacets.resize( mData.mNumberOfFacets );

        // get the other axis
        uint tOtherAxis = not aAxis;

        //  k = aAxis, j = other axis
        // #ifdef MORIS_USE_ARMA
        //             if ( aAxis == 0 )
        //             {
        //                 // check bounding box in J-direction
        //                 mData.mCandJ = arma::find(
        //                         ( mPoint( 0 ) - mData.mFacetMinCoordsX ) % ( mData.mFacetMaxCoordsX - mPoint( 0 ) ) > -gSDFepsilon );
        //             }
        //             else if ( aAxis == 1 )
        //             {
        //                 // check bounding box in J-direction
        //                 mData.mCandJ = arma::find(
        //                         ( mPoint( 1 ) - mData.mFacetMinCoordsY ) % ( mData.mFacetMaxCoordsY - mPoint( 1 ) ) > -gSDFepsilon );
        //             }

        //             // resize data object
        //             mData.mCandidateFacets.resize( mData.mCandJ.n_elem, 1 );

        //             // link to current object
        //             arma::Mat< uint >& tCand = mData.mCandidateFacets.matrix_data();

        //             // write data
        //             tCand = arma::conv_to< arma::Mat< uint > >::from( mData.mCandK );
        // #else

        uint tCandidateCount        = 0;
        uint tIntersectedFacetCount = 0;
        // loop over all lines in the aAxis direction
        for ( uint iLineIndex = 0; iLineIndex < mData.mNumberOfFacets; iLineIndex++ )
        {
            // check bounding box of the line against the point (point is above min coord and below max coord)
            if ( ( mFacetMaxCoords( iLineIndex, tOtherAxis ) - mPoint( tOtherAxis ) ) * ( mPoint( tOtherAxis ) - mFacetMinCoords( iLineIndex, tOtherAxis ) )
                    < gSDFepsilon )
            {
                // check if the point's !aAxis component is less the facet's minimum aAxis component. If so, the facet is intersected
                // NOTE: this makes the 2D raycast only cast in the positive axis direction
                if ( mFacetMinCoords( iLineIndex, aAxis ) - mPoint( aAxis ) > gSDFepsilon )
                {
                    mData.mCandidateFacets( tCandidateCount ) = iLineIndex;
                    tCandidateCount++;
                }
                // check the bounding box of the other axis to determine if the point is entirely in the bounding box
                else if ( ( mFacetMaxCoords( iLineIndex, aAxis ) - mPoint( aAxis ) )
                                  * ( mPoint( aAxis ) - mFacetMinCoords( iLineIndex, aAxis ) )
                          < gSDFepsilon )
                {
                    // if the point is totally in a line's bounding box, add line to candidate list
                    mData.mIntersectedFacets( tIntersectedFacetCount ) = mObject.get_facets()( iLineIndex );
                    tIntersectedFacetCount++;
                }
            }
        }

        // trim candidate and intersected matrix
        mData.mCandidateFacets.resize( tCandidateCount );
        mData.mIntersectedFacets.resize( tIntersectedFacetCount );
        // #endif
    }

    //-------------------------------------------------------------------------------

    void
    Raycast::intersect_triangles( uint aAxis )
    {
        // get number of candidate triangles
        uint tNumberOfCandidateFacets = mData.mCandidateFacets.size();

        // initialize counter for intersected triangles
        uint tCount = 0;

        // loop over all candidates
        for ( uint iCandidateFacetIndex = 0; iCandidateFacetIndex < tNumberOfCandidateFacets; ++iCandidateFacetIndex )
        {
            // get pointer to triangle
            Facet* tFacet = mObject.get_facets()( mData.mCandidateFacets( iCandidateFacetIndex ) );

            if ( tFacet->check_edge( 0, aAxis, mPoint ) )
            {
                if ( tFacet->check_edge( 1, aAxis, mPoint ) )
                {
                    if ( tFacet->check_edge( 2, aAxis, mPoint ) )
                    {
                        tFacet->flag();
                        ++tCount;
                    }
                }
            }
        }

        // resize container with intersected triangles
        mData.mIntersectedFacets.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // loop over all candidates
        for ( uint k = 0; k < tNumberOfCandidateFacets; ++k )
        {
            // get pointer to triangle
            Facet* tFacet = mObject.get_facets()( mData.mCandidateFacets( k ) );

            if ( tFacet->is_flagged() )
            {
                // add triangle to list
                mData.mIntersectedFacets( tCount++ ) = tFacet;

                // unflag triangle
                tFacet->unflag();
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Raycast::intersect_ray_with_facets( uint aAxis )
    {
        // get number of triangles
        uint tNumberOfFacets = mData.mIntersectedFacets.size();

        // initialize vector with coords in axis
        Matrix< DDRMat > tCoordsK( tNumberOfFacets, 1 );

        uint tCount = 0;

        bool tError;
        // loop over all intersected triangles and find intersection point
        for ( uint k = 0; k < tNumberOfFacets; ++k )
        {

            real tCoordK;

            // calculate intersection coordinate
            mData.mIntersectedFacets( k )->intersect_with_coordinate_axis(
                    mPoint,
                    aAxis,
                    tCoordK,
                    tError );

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
            // mData.mCoordsK.set_size( 1, 1, 0.0 );
            mData.mCoordsK.resize( 1, 0.0 );
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
            mData.mCoordsK.resize( tCount );


            // set first entry
            mData.mCoordsK( tCountUnique++ ) = tCoordsKSorted( 0 );

            // find unique entries
            for ( uint k = 1; k < tCount; ++k )
            {
                if ( std::abs( tCoordsKSorted( k ) - tCoordsKSorted( k - 1 ) ) > 10 * gSDFepsilon )
                {
                    mData.mCoordsK( tCountUnique++ ) = tCoordsKSorted( k );
                }
            }

            // chop vector
            mData.mCoordsK.resize( tCountUnique );
        }
    }

    //-------------------------------------------------------------------------------

    void
    Raycast::check_if_node_is_inside_triangles( uint aAxis )
    {
        MORIS_ASSERT( mDimension == 3,
                "SDF_Raycast::check_if_node_is_inside_triangles() only to be used for 3D raycasts. Current dimension = %d",
                mDimension );

        uint tNumCoordsK = mData.mCoordsK.size();

        bool tNodeIsInside = false;

        // only even number of intersections is considered
        if ( tNumCoordsK % 2 == 0 )
        {
            for ( uint k = 0; k < tNumCoordsK / 2; ++k )
            {
                tNodeIsInside = ( mPoint( aAxis ) > mData.mCoordsK( 2 * k ) ) && ( mPoint( aAxis ) < mData.mCoordsK( 2 * k + 1 ) );

                // break the loop if inside
                if ( tNodeIsInside )
                {
                    break;
                }
            }

            // set the inside flag of this node to the corresponding value
            mPointIsInside = tNodeIsInside;
        }
        else
        {
            mPointIsInside = 2;
        }
    }

    //-------------------------------------------------------------------------------


    void
    Raycast::check_if_node_is_inside_lines( uint aAxis )
    {
        MORIS_ASSERT( mDimension == 2,
                "SDF_Raycast::check_if_node_is_inside_lines() only to be used for 2D raycasts. Current dimension = %d",
                mDimension );

        // get length of the number of intersections computed
        uint tNumCoordsK = mData.mCoordsK.size();

        // check if the location of the intersection is greater than the location of the coordinate
        uint tIntersectionsRightOfPoint = 0;
        for ( uint iIntersectionIndex = 0; iIntersectionIndex < tNumCoordsK; iIntersectionIndex++ )
        {
            // BRENDAN, the axis indexing might not be correct, look here or other functions if problems arise
            if ( mData.mCoordsK( iIntersectionIndex ) - mPoint( aAxis ) > gSDFepsilon )
            {
                tIntersectionsRightOfPoint++;
            }
        }

        mPointIsInside = ( tIntersectionsRightOfPoint + mData.mIntersectedFacets.size() ) % 2;
    }

    void
    Raycast::random_rotation()
    {
        // generate random angle
        real tAngle = random_angle();

        // generate random rotation matrix
        Matrix< DDRMat > tRotation;
        if ( mDimension == 2 )
        {
            tRotation = rotation_matrix( tAngle );
        }
        else
        {
            // create random axis for cases larger than 2 dimensions
            Matrix< DDRMat > tAxis = random_axis( mDimension );
            tRotation              = rotation_matrix( tAxis, tAngle );
        }

        // rotate all vertices of triangle mesh
        for ( Facet_Vertex* tVertex : mObject.get_vertices() )
        {
            tVertex->rotate_node_coords( tRotation );
        }

        // update all facets
        for ( Facet* tFacet : mObject.get_facets() )
        {
            tFacet->update_data();
        }

        if ( mPointIsInside == 2 )
        {
            mPoint = tRotation * mPoint;
        }
    }

    // -----------------------------------------------------------------------------

    void
    Raycast::undo_rotation()
    {
        // rotate all vertices of triangle mesh
        for ( Facet_Vertex* tVertex : mObject.get_vertices() )
        {
            tVertex->reset_node_coords();
        }

        // update all facets
        for ( Facet* tFacet : mObject.get_facets() )
        {
            tFacet->update_data();
        }

        // Reset point coordinates
        mPoint = mOriginalPoint;
    }
}    // namespace moris::sdf