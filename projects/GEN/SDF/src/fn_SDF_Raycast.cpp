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

#include "fn_SDF_Raycast.hpp"

namespace moris::sdf
{
    Object_Region
    raycast_point(
            Object&           aObject,
            const Matrix< DDRMat >& aOriginalPoint )
    {
        // initialize return value to unsure
        Object_Region tPointIsInside = UNSURE;

        // Get information from the object
        uint tDimension = aObject.get_dimension();

        MORIS_ASSERT( aOriginalPoint.numel() == tDimension, "SDF_calculate_raycast() - Dimension mismatch. Point dimension = %lu\tObject Dimension = %d", aOriginalPoint.numel(), tDimension );

        // reset candidate and intersection points
        moris::Cell< uint >   tCandidateFacets;      // index of facets that lie within the bounding box, but might not actually be intersected
        moris::Cell< Facet* > tIntersectedFacets;    // intersection has been found, coordinate of intersection still to be computed

        // copy point data so it can be freely rotated
        Matrix< DDRMat > tPoint = aOriginalPoint;

        // flag that marks if rotation was called
        bool tRotation = false;

        switch ( tDimension )
        {
            case 2:
            {
                while ( tPointIsInside == UNSURE )
                {
                    // perform voxelizing algorithm in z-direction
                    tPointIsInside = voxelize_2D( aObject, tPoint, 1 );
                    if ( tPointIsInside == UNSURE )
                    {
                        // perform voxelizing algorithm in y-direction
                        tPointIsInside = voxelize_2D( aObject, tPoint, 0 );
                    }
                    if ( tPointIsInside == UNSURE )
                    {
                        tRotation = true;

                        random_rotation( aObject, tPoint );
                    }
                }
                break;
            }
            case 3:
            {
                while ( tPointIsInside == 2 )
                {
                    // perform voxelizing algorithm in z-direction
                    tPointIsInside = voxelize_3D( aObject, tPoint, 2 );
                    if ( tPointIsInside == UNSURE )
                    {
                        // perform voxelizing algorithm in y-direction
                        tPointIsInside = voxelize_3D( aObject, tPoint, 1 );

                        if ( tPointIsInside == UNSURE )
                        {
                            // perform voxelizing algorithm in x-direction
                            tPointIsInside = voxelize_3D( aObject, tPoint, 0 );
                        }
                    }

                    if ( tPointIsInside == UNSURE )
                    {
                        tRotation = true;

                        random_rotation( aObject, tPoint );
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "SDF_calculate_raycast() functionality not implemented for %dD problems.", tDimension );
            }
        }

        // reset the coordinates back to the orginal frame if they were rotated
        if ( tRotation )
        {
            aObject.undo_rotation();
        }

        return tPointIsInside;
    }

    Object_Region
    voxelize_2D(
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis )
    {
        // preselect lines in the aAxis direction
        moris::Cell< uint >   tCandidateFacets;
        moris::Cell< Facet* > tIntersectedFacets;
        preselect_lines( aObject, aPoint, aAxis, tCandidateFacets, tIntersectedFacets );

        // if there are no candidates and no intersected facets, the point is outside
        if ( tCandidateFacets.size() == 0 && tIntersectedFacets.size() == 0 )
        {
            return OUTSIDE;
        }

        // compute intersection if the point is inside a line's bounding box
        if ( tIntersectedFacets.size() > 0 )
        {
            moris::Cell< real > tIntersectionCoords = intersect_ray_with_facets( tIntersectedFacets, aPoint, aAxis );

            // check if the node is inside the polygon
            return check_if_node_is_inside_lines( tIntersectionCoords, tCandidateFacets, aPoint, aAxis );
        }

        // return unsure if there were no candidates
        return UNSURE;
    }


    Object_Region
    voxelize_3D(
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis )
    {
        // preselect triangles for intersection test
        moris::Cell< uint > tCandidateFacets = preselect_triangles( aObject, aPoint, aAxis );

        // from the candidate triangles, perform intersection
        if ( tCandidateFacets.size() > 0 )
        {
            moris::Cell< Facet* > tIntersectedFacets = intersect_triangles( tCandidateFacets, aObject, aPoint, aAxis );

            // intersect ray with triangles and check if node is inside
            if ( tIntersectedFacets.size() > 0 )
            {
                moris::Cell< real > tIntersectionCoords = intersect_ray_with_facets( tIntersectedFacets, aPoint, aAxis );

                return check_if_node_is_inside_triangles( tIntersectionCoords, aPoint, aAxis );
            }
        }
        // if there are no candidates, the point is outside
        return OUTSIDE;
    }

    moris::Cell< uint >
    preselect_triangles(
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis )
    {
        // select the axes that are not being cast in to preselect along (for a cast in the i-dir, preselect in j and k-dir)
        uint tFirstAxis;
        uint tSecondAxis;
        switch ( aAxis )
        {
            case 0:
            {
                tFirstAxis  = 2;
                tSecondAxis = 1;
                break;
            }
            case 1:
            {
                tFirstAxis  = 0;
                tSecondAxis = 2;
                break;
            }
            case 2:
            {
                tFirstAxis  = 1;
                tSecondAxis = 0;
                break;
            }
        }

        uint         tCountJ = 0;
        Cell< uint > tCandJ( aObject.get_num_facets() );
        for ( uint iFacetIndex = 0; iFacetIndex < aObject.get_num_facets(); ++iFacetIndex )
        {
            // check bounding box in J-direction
            if ( ( aPoint( tFirstAxis ) - aObject.get_facet_min_coords()( iFacetIndex, tFirstAxis ) ) * ( aObject.get_facet_max_coords()( iFacetIndex, tFirstAxis ) - aPoint( tFirstAxis ) ) > -gSDFepsilon )
            {
                // remember this triangle
                tCandJ( tCountJ ) = iFacetIndex;

                // increment counter
                ++tCountJ;
            }
        }

        // counter for triangles
        uint tCount = 0;

        // reset candidate size
        Cell< uint > tCandidateFacets( aObject.get_num_facets() );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // check bounding box in I-direction
            if ( ( aPoint( tSecondAxis ) - aObject.get_facet_min_coords()( tCandJ( k ), tSecondAxis ) ) * ( aObject.get_facet_max_coords()( tCandJ( k ), tSecondAxis ) - aPoint( tSecondAxis ) ) > -gSDFepsilon )
            {
                tCandidateFacets( tCount ) = tCandJ( k );
                ++tCount;
            }
        }

        tCandidateFacets.resize( tCount );

        return tCandidateFacets;
    }


    //-------------------------------------------------------------------------------

    void preselect_lines(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            moris::Cell< uint >&    aIntersectedFacets,
            moris::Cell< Facet* >&  aCandidateFacets )
    {
        // Ensure the function is being called for the proper number of facets
        MORIS_ERROR( aAxis < aPoint.numel(),
                "SDF_ preselect_lines() aPoint is %luD while coordinate axis %d specified.",
                aPoint.numel(),
                aAxis + 1 );
        MORIS_ASSERT( aPoint.numel() == 2,
                "SDF_ preselect_lines() should be called for 2D problems only. Query point dimension = %lu",
                aPoint.numel() );

        // reset candidate and intersected facet size
        aIntersectedFacets.resize( aObject.get_num_facets() );
        aCandidateFacets.resize( aObject.get_num_facets() );

        // get the other axis
        uint tOtherAxis = not aAxis;

        uint tCandidateCount        = 0;
        uint tIntersectedFacetCount = 0;
        // loop over all lines in the aAxis direction
        for ( uint iLineIndex = 0; iLineIndex < aObject.get_num_facets(); iLineIndex++ )
        {
            // check bounding box of the line against the point (point is above min coord and below max coord)
            if ( ( aObject.get_facet_max_coords()( iLineIndex, tOtherAxis ) - aPoint( tOtherAxis ) ) * ( aPoint( tOtherAxis ) - aObject.get_facet_min_coords()( iLineIndex, tOtherAxis ) )
                    > gSDFepsilon )
            {
                // check if the point's !aAxis component is less the facet's minimum aAxis component. If so, the facet is intersected
                // NOTE: this makes the 2D raycast only cast in the positive axis direction
                // FIXME: Candidate and intersected variable names are backward, need to be carefully considered to be switched back BRENDAN
                if ( aObject.get_facet_min_coords()( iLineIndex, aAxis ) - aPoint( aAxis ) > gSDFepsilon )
                {
                    aIntersectedFacets( tIntersectedFacetCount ) = iLineIndex;
                    tIntersectedFacetCount++;
                }
                // check the bounding box of the other axis to determine if the point is entirely in the bounding box
                else if ( ( aObject.get_facet_max_coords()( iLineIndex, aAxis ) - aPoint( aAxis ) )
                                  * ( aPoint( aAxis ) - aObject.get_facet_min_coords()( iLineIndex, aAxis ) )
                          > gSDFepsilon )
                {
                    // if the point is totally in a line's bounding box, add line to candidate list
                    aCandidateFacets( tCandidateCount ) = aObject.get_facets()( iLineIndex );
                    tCandidateCount++;
                }
            }
        }

        // trim candidate and intersected matrix
        aIntersectedFacets.resize( tIntersectedFacetCount );
        aCandidateFacets.resize( tCandidateCount );
    }

    //-------------------------------------------------------------------------------

    Cell< Facet* >
    intersect_triangles(
            Cell< uint >&     aCandidateFacets,
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis )
    {

        // get number of candidate triangles
        uint tNumberOfCandidateFacets = aCandidateFacets.size();

        // initialize counter for intersected triangles
        uint tCount = 0;

        // loop over all candidates
        for ( uint iCandidateFacetIndex = 0; iCandidateFacetIndex < tNumberOfCandidateFacets; ++iCandidateFacetIndex )
        {
            // get pointer to triangle
            Facet* tFacet = aObject.get_facets()( aCandidateFacets( iCandidateFacetIndex ) );

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
        moris::Cell< Facet* > tIntersectedFacets( tCount, nullptr );

        // reset counter
        tCount = 0;

        // loop over all candidates
        for ( uint iCandidateFacetIndex = 0; iCandidateFacetIndex < tNumberOfCandidateFacets; ++iCandidateFacetIndex )
        {
            // get pointer to triangle
            Facet* tFacet = aObject.get_facets()( aCandidateFacets( iCandidateFacetIndex ) );

            if ( tFacet->is_flagged() )
            {
                // add triangle to list
                tIntersectedFacets( tCount++ ) = tFacet;

                // unflag triangle
                tFacet->unflag();
            }
        }

        return tIntersectedFacets;
    }

    //-------------------------------------------------------------------------------

    moris::Cell< real >
    intersect_ray_with_facets(
            moris::Cell< Facet* >&  aIntersectedFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis )
    {
        // get number of triangles
        uint tNumberOfFacets = aIntersectedFacets.size();

        // initialize vector with coords in axis
        Matrix< DDRMat > tCoordsK( tNumberOfFacets, 1 );

        uint tCount = 0;

        bool tError;
        // loop over all intersected triangles and find intersection point
        for ( uint iFacetIndex = 0; iFacetIndex < tNumberOfFacets; ++iFacetIndex )
        {

            real tCoordK;

            // calculate intersection coordinate
            aIntersectedFacets( iFacetIndex )->intersect_with_coordinate_axis( aPoint, aAxis, tCoordK, tError );

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
            return { 0.0 };
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

            // initialize return Cell
            moris::Cell< real > tIntersectionCoords( tCount );

            // set first entry
            tIntersectionCoords( tCountUnique++ ) = tCoordsKSorted( 0 );

            // find unique entries
            for ( uint k = 1; k < tCount; ++k )
            {
                if ( std::abs( tCoordsKSorted( k ) - tCoordsKSorted( k - 1 ) ) > 10 * gSDFepsilon )
                {
                    tIntersectionCoords( tCountUnique++ ) = tCoordsKSorted( k );
                }
            }

            // chop vector
            tIntersectionCoords.resize( tCountUnique );

            return tIntersectionCoords;
        }
    }

    //-------------------------------------------------------------------------------

    Object_Region
    check_if_node_is_inside_triangles(
            moris::Cell< real >& aIntersectionCoords,
            Matrix< DDRMat >&    aPoint,
            uint                 aAxis )
    {
        uint tNumCoordsK = aIntersectionCoords.size();

        uint tNodeIsInside = 2;

        // only even number of intersections is considered
        if ( tNumCoordsK % 2 == 0 )
        {
            for ( uint iIntersectionIndex = 0; iIntersectionIndex < tNumCoordsK / 2; ++iIntersectionIndex )
            {
                tNodeIsInside = ( aPoint( aAxis ) > aIntersectionCoords( 2 * iIntersectionIndex ) ) && ( aPoint( aAxis ) < aIntersectionCoords( 2 * iIntersectionIndex + 1 ) );
            }
        }

        return static_cast< Object_Region >( tNodeIsInside );
    }

    //-------------------------------------------------------------------------------


    Object_Region
    check_if_node_is_inside_lines(
            const Cell< real >&     aIntersectionCoords,
            const Cell< uint >&     aCandidateFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis )
    {
        // get length of the number of intersections computed
        uint tNumCoordsK = aIntersectionCoords.size();

        // check if the location of the intersection is greater than the location of the coordinate
        uint tIntersectionsRightOfPoint = 0;
        for ( uint iIntersectionIndex = 0; iIntersectionIndex < tNumCoordsK; iIntersectionIndex++ )
        {
            if ( aIntersectionCoords( iIntersectionIndex ) - aPoint( aAxis ) > gSDFepsilon )
            {
                tIntersectionsRightOfPoint++;
            }
        }

        // if the ray passes through an even number of lines, the point is inside. otherwise, it is outside.
        return static_cast< Object_Region >( ( tIntersectionsRightOfPoint + aCandidateFacets.size() ) % 2 );
    }

    inline void
    random_rotation(
            Object&           aObject,
            Matrix< DDRMat >& aPoint )
    {
        // generate random angle
        real tAngle = random_angle();

        // generate random rotation matrix
        Matrix< DDRMat > tRotation;
        if ( aObject.get_dimension() == 2 )
        {
            tRotation = rotation_matrix( tAngle );
        }
        else
        {
            // create random axis for cases larger than 2 dimensions
            Matrix< DDRMat > tAxis = random_axis( aObject.get_dimension() );
            tRotation              = rotation_matrix( tAxis, tAngle );
        }

        aObject.rotate_object( tRotation );
        aPoint = tRotation * aPoint;
    }
    // -----------------------------------------------------------------------------

}    // namespace moris::sdf