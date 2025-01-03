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
#include "fn_trans.hpp"

#include "fn_SDF_Raycast.hpp"

namespace moris::sdf
{
    Object_Region
    raycast_point(
            Object&                 aObject,
            const Matrix< DDRMat >& aOriginalPoint )
    {
        // initialize return value to unsure
        Object_Region tPointIsInside = UNSURE;

        MORIS_ASSERT( aOriginalPoint.numel() == aObject.get_dimension(),
                "fn_SDF_Raycast: raycast_point() - Dimension mismatch. Point dimension = %lu\tObject Dimension = %d",
                aOriginalPoint.numel(),
                aObject.get_dimension() );

        // copy point data so it can be freely rotated
        Matrix< DDRMat > tPoint = aOriginalPoint;

        // flag that marks if rotation was called
        bool tRotation = false;
        while ( tPointIsInside == UNSURE )
        {
            // perform voxelizing algorithm in x-direction
            tPointIsInside = voxelize( aObject, tPoint, 0 );
            if ( tPointIsInside == UNSURE )
            {
                // perform voxelizing algorithm in y-direction
                tPointIsInside = voxelize( aObject, tPoint, 1 );

                // perform voxelizing algorithm in z-direction if 3D
                if ( tPointIsInside == UNSURE && aObject.get_dimension() == 3 )
                {
                    tPointIsInside = voxelize( aObject, tPoint, 2 );
                }
            }
            // if still unsure, rotate and cast again
            if ( tPointIsInside == UNSURE )
            {
                tRotation = true;

                random_rotation( aObject, tPoint );
            }
        }

        // reset the coordinates back to the orginal frame if they were rotated
        if ( tRotation )
        {
            aObject.reset_coordinates();
        }

        return tPointIsInside;
    }

    Vector< real >
    compute_distance_to_facets(
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis )
    {
        MORIS_ERROR( aObject.get_dimension() == aPoint.numel(),
                "SDF-Raycast::compute_distance_to_facets(): dimension mismatch. Object dimension: %d Point dimension: %lu",
                aObject.get_dimension(),
                aPoint.numel() );

        switch ( aObject.get_dimension() )
        {
            case 2:
            {
                // preselect lines in the aAxis direction
                Vector< uint >   tIntersectedFacets;
                Vector< Facet* > tCandidateFacets;
                bool                  tPreselectionSuccessful = preselect_lines( aObject, aPoint, aAxis, tIntersectedFacets, tCandidateFacets );

                if( tPreselectionSuccessful )
                {
                        
                }

                // get pointers to the facets from the candidates
                Vector< Facet* > tFacetsFromCandidates( tIntersectedFacets.size() );
                for ( uint iCandidate = 0; iCandidate < tIntersectedFacets.size(); iCandidate++ )
                {
                    tFacetsFromCandidates( iCandidate ) = &aObject.get_facet( tIntersectedFacets( iCandidate ) );
                }

                // append the candidates to the intersected facets to compute intersection
                tCandidateFacets.insert( tCandidateFacets.size(), tFacetsFromCandidates.begin(), tFacetsFromCandidates.end() );

                // compute intersection for all preselected facets
                return intersect_ray_with_facets( tCandidateFacets, aPoint, aAxis );
            }
            case 3:
            {
                // preselect triangles in positive and negative aAxis directions for intersection test
                Vector< uint > tCandidateFacets = preselect_triangles( aObject, aPoint, aAxis );

                // filter out facets that are definitely in the negative aAxis direction from the point
                // FIXME: a different preselection function could be written to avoid extra checks
                for ( uint iCandidate = 0; iCandidate < tCandidateFacets.size(); iCandidate++ )
                {
                    if ( aObject.get_facet_max_coord( iCandidate, aAxis ) < aPoint( aAxis ) )
                    {
                        tCandidateFacets.erase( iCandidate );
                    }
                }

                // from the candidate triangles, perform intersection
                Vector< Facet* > tIntersectedFacets = intersect_triangles( tCandidateFacets, aObject, aPoint, aAxis );

                // compute intersection locations
                Vector< real > tIntersectionCoordinates = intersect_ray_with_facets( tIntersectedFacets, aPoint, aAxis );

                // remove intersection locations that are behind the point
                for ( uint iIntersection = 0; iIntersection < tIntersectionCoordinates.size(); iIntersection++ )
                {
                    if ( tIntersectionCoordinates( iIntersection ) < aPoint( aAxis ) )
                    {
                        tIntersectionCoordinates.erase( iIntersection );
                    }
                }

                return tIntersectionCoordinates;
            }
            default:
            {
                MORIS_ERROR( false, "SDF-Raycast: Object dimension of %d provided. Not supported", aObject.get_dimension() );
                return {};
            }
        }
    }

    Object_Region
    voxelize(
            Object&           aObject,
            Matrix< DDRMat >& aPoint,
            uint              aAxis )
    {
        switch ( aObject.get_dimension() )
        {
            case 2:
            {
                // preselect lines in the aAxis direction
                Vector< uint >   tIntersectedFacets;
                Vector< Facet* > tCandidateFacets;
                bool                  tPreselectionSuccessful = preselect_lines( aObject, aPoint, aAxis, tIntersectedFacets, tCandidateFacets );

                // if the cast point was close to a vertex, cast again
                if ( !tPreselectionSuccessful )
                {
                    return UNSURE;
                }
                // if there are no candidates and no intersected facets, the point is outside
                else if ( tIntersectedFacets.size() == 0 && tCandidateFacets.size() == 0 )
                {
                    return OUTSIDE;
                }

                // compute intersection coordinates if the point is inside a line's bounding box
                Vector< real > tIntersectionCoords = intersect_ray_with_facets( tCandidateFacets, aPoint, aAxis );

                // check if the node is inside the polygon
                return check_if_node_is_inside_lines( tIntersectionCoords, tIntersectedFacets, aPoint, aAxis );
            }
            case 3:
            {
                // preselect triangles for intersection test
                Vector< uint > tCandidateFacets = preselect_triangles( aObject, aPoint, aAxis );

                // from the candidate triangles, perform intersection
                if ( tCandidateFacets.size() > 0 )
                {
                    Vector< Facet* > tIntersectedFacets = intersect_triangles( tCandidateFacets, aObject, aPoint, aAxis );

                    // intersect ray with triangles and check if node is inside
                    if ( tIntersectedFacets.size() > 0 )
                    {
                        Vector< real > tIntersectionCoords = intersect_ray_with_facets( tIntersectedFacets, aPoint, aAxis );

                        return check_if_node_is_inside_triangles( tIntersectionCoords, aPoint, aAxis );
                    }
                }
                // if there are no candidates, the point is outside
                return OUTSIDE;
            }
            default:
            {
                return UNSURE;
            }
        }
    }

    Vector< uint >
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
            default:
            {
                tFirstAxis  = 0;
                tSecondAxis = 0;
                MORIS_ERROR( false, "Raycasting axis must be given an index 0, 1, or 2. Given: %d", aAxis );
            }
        }

        uint         tCountJ = 0;
        Vector< uint > tCandJ( aObject.get_num_facets() );
        for ( uint iFacetIndex = 0; iFacetIndex < aObject.get_num_facets(); ++iFacetIndex )
        {
            // check bounding box in J-direction
            if ( ( aPoint( tFirstAxis ) - aObject.get_facet_min_coord( iFacetIndex, tFirstAxis ) )
                            * ( aObject.get_facet_max_coord( iFacetIndex, tFirstAxis ) - aPoint( tFirstAxis ) )
                    > -MORIS_REAL_EPS )
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
        Vector< uint > tCandidateFacets( aObject.get_num_facets() );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // check bounding box in I-direction
            if ( ( aPoint( tSecondAxis ) - aObject.get_facet_min_coord( tCandJ( k ), tSecondAxis ) )
                            * ( aObject.get_facet_max_coord( tCandJ( k ), tSecondAxis ) - aPoint( tSecondAxis ) )
                    > -MORIS_REAL_EPS )
            {
                tCandidateFacets( tCount ) = tCandJ( k );
                ++tCount;
            }
        }

        tCandidateFacets.resize( tCount );

        return tCandidateFacets;
    }

    //-------------------------------------------------------------------------------

    bool preselect_lines(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Vector< uint >&    aIntersectedFacets,
            Vector< Facet* >&  aCandidateFacets )
    {
        bool tPreselectionSuccessful = true;

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
            // get the difference of the cast point and the facet min and max coords in the !aAxis direction
            real tMaxCoordOffAxisDifference = aObject.get_facet_max_coord( iLineIndex, tOtherAxis ) - aPoint( tOtherAxis );
            real tMinCoordOffAxisDifference = aObject.get_facet_min_coord( iLineIndex, tOtherAxis ) - aPoint( tOtherAxis );

            // if the cast point is close to either vertex, return an error for preselection
            if ( std::abs( tMaxCoordOffAxisDifference ) < gSDFepsilon or std::abs( tMinCoordOffAxisDifference ) < gSDFepsilon )
            {
                aObject.get_facet( iLineIndex ).flag();
                tPreselectionSuccessful = false;
            }

            // check bounding box of the line against the point (point is above min coord and below max coord)
            if ( tMaxCoordOffAxisDifference * tMinCoordOffAxisDifference < MORIS_REAL_EPS )
            {
                // get the difference of the cast point and the facet min and max coords in the aAxis direction
                real tMaxCoordAxisDifference = aObject.get_facet_max_coord( iLineIndex, aAxis ) - aPoint( aAxis );
                real tMinCoordAxisDifference = aObject.get_facet_min_coord( iLineIndex, aAxis ) - aPoint( aAxis );

                // check if the point's !aAxis component is less the facet's minimum aAxis component. If so, the facet is intersected
                // NOTE: this makes the 2D raycast only cast in the positive axis direction
                if ( tMinCoordAxisDifference > MORIS_REAL_EPS )
                {
                    aIntersectedFacets( tIntersectedFacetCount ) = iLineIndex;
                    tIntersectedFacetCount++;
                }
                // check the bounding box of the other axis to determine if the point is entirely in the bounding box
                else if ( tMaxCoordAxisDifference * tMinCoordAxisDifference < MORIS_REAL_EPS )
                {
                    // if the point is totally in a line's bounding box, add line to candidate list
                    aCandidateFacets( tCandidateCount ) = &aObject.get_facet( iLineIndex );
                    tCandidateCount++;
                }
            }
        }

        // trim candidate and intersected matrix
        aIntersectedFacets.resize( tIntersectedFacetCount );
        aCandidateFacets.resize( tCandidateCount );

        return tPreselectionSuccessful;
    }

    //-------------------------------------------------------------------------------

    Vector< Facet* >
    intersect_triangles(
            Vector< uint >&     aCandidateFacets,
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
            // get reference to triangle
            Facet& tFacet = aObject.get_facet( aCandidateFacets( iCandidateFacetIndex ) );

            if ( tFacet.check_edge( 0, aAxis, aPoint ) )
            {
                if ( tFacet.check_edge( 1, aAxis, aPoint ) )
                {
                    if ( tFacet.check_edge( 2, aAxis, aPoint ) )
                    {
                        tFacet.flag();
                        ++tCount;
                    }
                }
            }
        }

        // resize container with intersected triangles
        Vector< Facet* > tIntersectedFacets( tCount );

        // reset counter
        tCount = 0;

        // loop over all candidates
        for ( uint iCandidateFacetIndex = 0; iCandidateFacetIndex < tNumberOfCandidateFacets; ++iCandidateFacetIndex )
        {
            // get reference to triangle
            Facet& tFacet = aObject.get_facet( aCandidateFacets( iCandidateFacetIndex ) );

            if ( tFacet.is_flagged() )
            {
                // add triangle to list
                tIntersectedFacets( tCount++ ) = &tFacet;

                // unflag triangle
                tFacet.unflag();
            }
        }

        return tIntersectedFacets;
    }

    //-------------------------------------------------------------------------------

    Vector< real >
    intersect_ray_with_facets(
            Vector< Facet* >&  aIntersectedFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis )
    {
        // return nothing if there are no intersected facets
        if ( aIntersectedFacets.size() == 0 )
        {
            return {};
        }

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
                tCoordsK( tCount++ ) = std::round( tCoordK / MORIS_REAL_EPS ) * MORIS_REAL_EPS;
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
            Vector< real > tIntersectionCoords( tCount );

            // set first entry
            tIntersectionCoords( tCountUnique++ ) = tCoordsKSorted( 0 );

            // find unique entries
            for ( uint k = 1; k < tCount; ++k )
            {
                if ( std::abs( tCoordsKSorted( k ) - tCoordsKSorted( k - 1 ) ) > 10 * MORIS_REAL_EPS )
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
            Vector< real >& aIntersectionCoords,
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
                tNodeIsInside = ( aPoint( aAxis ) > aIntersectionCoords( 2 * iIntersectionIndex ) )
                            and ( aPoint( aAxis ) < aIntersectionCoords( 2 * iIntersectionIndex + 1 ) );
            }
        }

        return static_cast< Object_Region >( tNodeIsInside );
    }

    //-------------------------------------------------------------------------------

    Object_Region
    check_if_node_is_inside_lines(
            const Vector< real >&     aIntersectionCoords,
            const Vector< uint >&     aCandidateFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis )
    {
        // get length of the number of intersections computed
        uint tNumCoordsK = aIntersectionCoords.size();

        // check if the location of the intersection is greater than the location of the coordinate
        uint tIntersectionsRightOfPoint = 0;
        for ( uint iIntersectionIndex = 0; iIntersectionIndex < tNumCoordsK; iIntersectionIndex++ )
        {
            if ( aIntersectionCoords( iIntersectionIndex ) - aPoint( aAxis ) > MORIS_REAL_EPS )
            {
                tIntersectionsRightOfPoint++;
            }
        }

        // if the ray passes through an even number of lines, the point is inside. otherwise, it is outside.
        return static_cast< Object_Region >( ( tIntersectionsRightOfPoint + aCandidateFacets.size() ) % 2 );
    }

    void
    random_rotation(
            Object&           aObject,
            Matrix< DDRMat >& aPoint )
    {
        MORIS_ASSERT( aPoint.n_cols() == 1 or aPoint.n_rows() == 1,
                "SDF-Raycast::random_rotation(): aPoint expected vector, but matrix received. Dimensions: %lux%lu",
                aPoint.n_rows(),
                aPoint.n_cols() );
        if ( aPoint.n_cols() != 1 )
        {
            aPoint = trans( aPoint );
        }

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

        aObject.rotate( tRotation );
        aPoint = tRotation * aPoint;
    }
    // -----------------------------------------------------------------------------

}    // namespace moris::sdf
