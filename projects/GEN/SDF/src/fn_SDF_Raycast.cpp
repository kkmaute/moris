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
            Object&          aObject,
            Matrix< DDRMat > aPoint )
    {
        Object_Region tPointIsInside = UNSURE;

        MORIS_ASSERT( aPoint.numel() == aObject.get_dimension(),
                "fn_SDF_Raycast: raycast_point() - Dimension mismatch. Point dimension = %lu\tObject Dimension = %d",
                aPoint.numel(),
                aObject.get_dimension() );

        // flag that marks if rotation was called
        bool tRotation = false;
        while ( tPointIsInside == UNSURE )
        {
            tPointIsInside = voxelize( aObject, aPoint, 0 );

            if ( tPointIsInside == UNSURE )
            {
                tPointIsInside = voxelize( aObject, aPoint, 1 );
                if ( aObject.get_dimension() == 3 )
                {
                    tPointIsInside = voxelize( aObject, aPoint, 2 );
                }

                if ( tPointIsInside == UNSURE )
                {
                    // if still unsure, rotate and cast again
                    tRotation = true;
                    random_rotation( aObject, aPoint );
                }
            }
        }

        // reset the coordinates back to the orginal frame if they were rotated
        if ( tRotation )
        {
            aObject.reset_coordinates();
        }

        return tPointIsInside;
    }

    // -------------------------------------------------------------------------------

    Vector< real >
    compute_distance_to_facets(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            const uint              aAxis,
            Vector< Facet* >&       aIntersectedFacets )
    {
        MORIS_ASSERT( aObject.get_dimension() == aPoint.numel(),
                "SDF-Raycast::compute_distance_to_facets(): Dimension mismatch. Object dimension: %d Point dimension: %lu",
                aObject.get_dimension(),
                aPoint.numel() );

        // Initialize preselection result and vector to store intersection locations
        Preselection_Result tPreselection;
        Vector< real >      tIntersectionCoordinates;

        switch ( aObject.get_dimension() )
        {
            case 2:
            {
                // preselect lines in the aAxis direction
                moris::Vector< uint > tIntersectedFacets;
                tPreselection = preselect_lines( aObject, aPoint, aAxis, tIntersectedFacets, aIntersectedFacets );

                // get pointers to all facets in tIntersectedFacets to compute intersection location
                moris::Vector< Facet* > tFacetsFromCandidates( tIntersectedFacets.size() );
                for ( uint iCandidate = 0; iCandidate < tIntersectedFacets.size(); iCandidate++ )
                {
                    tFacetsFromCandidates( iCandidate ) = &aObject.get_facet( tIntersectedFacets( iCandidate ) );
                }

                // append the candidates to the intersected facets to compute intersection
                aIntersectedFacets.insert( aIntersectedFacets.size(), tFacetsFromCandidates.begin(), tFacetsFromCandidates.end() );

                // compute intersection locations
                tIntersectionCoordinates = intersect_ray_with_facets( aIntersectedFacets, aPoint, aAxis, tPreselection, true );

                break;
            }
            case 3:
            {
                // preselect triangles in positive and negative aAxis directions for intersection test
                Vector< uint > tCandidateFacets;
                tPreselection = preselect_triangles( aObject, aPoint, aAxis, tCandidateFacets );

                tIntersectionCoordinates = intersect_triangles_moller_trumbore( tCandidateFacets, aObject, aPoint, aAxis, aIntersectedFacets );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "SDF-Raycast: Unsupported dimension of %d provided.", aObject.get_dimension() );
                return {};
            }
        }


        // --------------------------------------------------------------
        // remove intersections that are behind the cast point
        // --------------------------------------------------------------

        // Get the locations of the first entries of the coordinates and facets
        auto tItCoord = tIntersectionCoordinates.begin();
        auto tItFacet = aIntersectedFacets.begin();

        // Get the value of the point in the axis direction and the object's intersection tolerance
        real tPointAxisValue        = aPoint( aAxis );
        real tIntersectionTolerance = aObject.get_intersection_tolerance();

        while ( tItCoord != tIntersectionCoordinates.end() )
        {
            if ( *tItCoord + tIntersectionTolerance < tPointAxisValue )
            {
                tItCoord = tIntersectionCoordinates.erase( tItCoord );
                tItFacet = aIntersectedFacets.erase( tItFacet );
            }
            else
            {
                ++tItCoord;
                ++tItFacet;
            }
        }

        return tIntersectionCoordinates;
    }

    // -------------------------------------------------------------------------------

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
                Vector< uint >      tIntersectedFacets;
                Vector< Facet* >    tCandidateFacets;
                Preselection_Result tPreselection = preselect_lines( aObject, aPoint, aAxis, tIntersectedFacets, tCandidateFacets );

                switch ( tPreselection )
                {
                    case FAIL_CAST_TO_VERTEX:
                    {    // the ray will intersect a vertex, cast again in another direction
                        return UNSURE;
                    }
                    case FAIL_ON_VERTEX:    // the ray originates from a vertex, return interface
                    case FAIL_ON_EDGE:      // the ray originates from an edge, return interface
                        // the facet is along an axis, and the cast point is on it
                        return INTERFACE;
                    case SUCCESS:
                    {    // no pathological case was detected in preselection
                        // compute intersection coordinates if the point is inside a line's bounding box
                        moris::Vector< real > tIntersectionCoords = intersect_ray_with_facets( tCandidateFacets, aPoint, aAxis, tPreselection );

                        // check if the node is inside the polygon
                        return check_if_node_is_inside_lines( aObject, tIntersectionCoords, tIntersectedFacets, aPoint, aAxis );
                    }
                    default:
                    {
                        MORIS_ERROR( false, "Unexpected preselection result of %d returned from preselect_lines()", tPreselection );
                    }
                }
            }
            case 3:
            {
                // // preselect triangles for intersection test
                Vector< uint > tCandidateFacets;
                preselect_triangles_moller_trumbore( aObject, aPoint, aAxis, tCandidateFacets );

                // Initialize vector to store the facets the ray intersects
                Vector< Facet* > tIntersectedFacets;

                // Compute the intersection locations and their corresponding facets
                Vector< real > tIntersectionCoords = intersect_triangles_moller_trumbore( tCandidateFacets, aObject, aPoint, aAxis, tIntersectedFacets );

                return check_if_node_is_inside_triangles( tIntersectionCoords, aPoint, aAxis, aObject.get_intersection_tolerance() );

                // switch ( tPreselection )
                // {
                //     case FAIL_ON_VERTEX:    // the ray originates from a vertex, return interface
                //     case FAIL_ON_EDGE:      // the ray originates from an edge, return interface
                //     {
                //         return INTERFACE;
                //     }
                //     case SUCCESS:    // no pathological case was detected in preselection
                //     {
                //         // Initialize vector to store the facets the ray intersects
                //         Vector< Facet* > tIntersectedFacets;

                //         // Compute the intersection locations and their corresponding facets
                //         Vector< real > tIntersectionCoords = intersect_triangles_moller_trumbore( tCandidateFacets, aObject, aPoint, aAxis, tIntersectedFacets );
                //         // // From the candidate triangles, determine which will actually be intersected
                //         // Vector< Facet* > tIntersectedFacets = intersect_triangles( tCandidateFacets, aObject, aPoint, aAxis );

                //         // // intersect ray with triangles and check if node is inside
                //         // // FIXME: handle casting onto vertices by changing intersect_triangles to return a preselection result
                //         // Vector< real > tIntersectionCoords = intersect_ray_with_facets( tIntersectedFacets, aPoint, aAxis, tPreselection );

                //         return check_if_node_is_inside_triangles( tIntersectionCoords, aPoint, aAxis, aObject.get_intersection_tolerance() );
                //     }
                //     case FAIL_CAST_TO_VERTEX:    // the ray will intersect a vertex, cast again in another direction
                //     case FAIL_CAST_TO_EDGE:      // the ray will intersect an edge, cast again in another direction
                //     {
                //         return UNSURE;
                //     }
                //     default:
                //     {
                //         MORIS_ERROR( false, "Unexpected preselection result of %d returned from preselect_triangles()", tPreselection );
                //     }
                // }
            }
            default:
            {
                return UNSURE;
            }
        }
    }

    // -------------------------------------------------------------------------------

    Preselection_Result
    preselect_triangles(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Vector< uint >&         aCandidateFacets )
    {
        Preselection_Result tPreselection = SUCCESS;

        // select the axes that are not being cast in to preselect along (for a cast in the i-dir, preselect in j and k-dir)
        uint tFirstAxis;
        uint tSecondAxis;
        triangle_permutation( aAxis, tFirstAxis, tSecondAxis );

        // Get the intersection tolerance for this Object
        real tEpsilon        = aObject.get_intersection_tolerance();
        real tEpsilonScaling = 10.0;

        uint           tCountJ = 0;
        Vector< uint > tCandJ( aObject.get_num_facets() );
        for ( uint iFacetIndex = 0; iFacetIndex < aObject.get_num_facets(); ++iFacetIndex )
        {
            // Get the difference of the cast point and the facet min and max coords in the J direction
            real tJDirectionPreselection = ( aPoint( tFirstAxis ) - aObject.get_facet_min_coord( iFacetIndex, tFirstAxis ) )
                                         * ( aObject.get_facet_max_coord( iFacetIndex, tFirstAxis ) - aPoint( tFirstAxis ) );

            // check if the point is inside the bounding box of the triangle
            if ( tJDirectionPreselection > -tEpsilon )
            {
                // remember this triangle
                tCandJ( tCountJ ) = iFacetIndex;

                // increment counter
                ++tCountJ;

                // check that the point is not close to an edge or vertex
                if ( std::abs( tJDirectionPreselection ) < 2.0 * tEpsilon )
                {
                    const Vector< std::shared_ptr< sdf::Facet_Vertex > >& tVertices = aObject.get_facet( iFacetIndex ).get_facet_vertex_pointers();

                    // lambda function that determines if the point is on any of the vertices
                    auto tPointOnVertex = [ &aPoint, &tEpsilon, &tEpsilonScaling ]( std::shared_ptr< Facet_Vertex > aVertex ) {
                        const Matrix< DDRMat >& tVertexCoordinates = aVertex->get_coords_reference();
                        return std::equal( tVertexCoordinates.cbegin(), tVertexCoordinates.cend(), aPoint.cbegin(), [ & ]( real tCoord, real tPoint ) {
                            return std::abs( tCoord - tPoint ) < tEpsilon;
                        } ); };

                    // determine which is the case, edge or vertex
                    if ( std::any_of( tVertices.cbegin(), tVertices.cend(), tPointOnVertex ) )
                    {
                        aCandidateFacets.resize( 1 );
                        aCandidateFacets( 0 ) = iFacetIndex;
                        return FAIL_ON_VERTEX;
                    }
                    else if ( aObject.get_facet( iFacetIndex ).get_distance_to_point( aPoint ) < tEpsilon )
                    {
                        aCandidateFacets.resize( 1 );
                        aCandidateFacets( 0 ) = iFacetIndex;
                        return FAIL_ON_EDGE;
                    }
                    else
                    {
                        // mark that the cast hit an edge, but do not return since there are other facets to check
                        tPreselection = FAIL_CAST_TO_EDGE;
                    }
                }
            }
        }

        // counter for triangles
        uint tCount = 0;

        // reset candidate size
        aCandidateFacets.resize( aObject.get_num_facets() );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // Get the difference of the cast point and the facet min and max coords in the K direction
            real tKDirectionPreselection = ( aPoint( tSecondAxis ) - aObject.get_facet_min_coord( tCandJ( k ), tSecondAxis ) )
                                         * ( aObject.get_facet_max_coord( tCandJ( k ), tSecondAxis ) - aPoint( tSecondAxis ) );

            // see if the point is inside the bounding box of the triangle
            if ( tKDirectionPreselection > -tEpsilon )
            {
                aCandidateFacets( tCount ) = tCandJ( k );
                ++tCount;

                // check if the point is close to the edge
                if ( std::abs( tKDirectionPreselection ) < 2.0 * tEpsilon )
                {
                    const Vector< std::shared_ptr< sdf::Facet_Vertex > >& tVertices = aObject.get_facet( tCandJ( k ) ).get_facet_vertex_pointers();

                    // lambda function that determines if the point is on any of the vertices
                    auto tPointOnVertex = [ &aPoint, &tEpsilon, &tEpsilonScaling ]( std::shared_ptr< Facet_Vertex > aVertex ) {
                        const Matrix< DDRMat >& tVertexCoordinates = aVertex->get_coords_reference();
                        return std::equal( tVertexCoordinates.cbegin(), tVertexCoordinates.cend(), aPoint.cbegin(), [ & ]( real tCoord, real tPoint ) {
                            return std::abs( tCoord - tPoint ) < tEpsilon;
                        } ); };

                    // determine which is the case, edge or vertex
                    if ( std::any_of( tVertices.cbegin(), tVertices.cend(), tPointOnVertex ) )
                    {
                        aCandidateFacets.resize( 1 );
                        aCandidateFacets( 0 ) = tCandJ( k );
                        return FAIL_ON_VERTEX;
                    }
                    else if ( aObject.get_facet( tCandJ( k ) ).get_distance_to_point( aPoint ) < tEpsilon )
                    {
                        aCandidateFacets.resize( 1 );
                        aCandidateFacets( 0 ) = tCandJ( k );
                        return FAIL_ON_EDGE;
                    }
                    else
                    {
                        // mark that the cast hit an edge, but do not return since there are other facets to check
                        tPreselection = FAIL_CAST_TO_EDGE;
                    }
                }
            }
        }

        aCandidateFacets.resize( tCount );
        return tPreselection;
    }

    //-------------------------------------------------------------------------------

    void
    preselect_triangles_moller_trumbore(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Vector< uint >&         aCandidateFacets )
    {
        // select the axes that are not being cast in to preselect along (for a cast in the i-dir, preselect in j and k-dir)
        uint tFirstAxis;
        uint tSecondAxis;
        triangle_permutation( aAxis, tFirstAxis, tSecondAxis );

        // Get the intersection tolerance for this Object
        real tEpsilon = aObject.get_intersection_tolerance();

        uint           tCountJ = 0;
        Vector< uint > tCandJ( aObject.get_num_facets() );
        for ( uint iFacetIndex = 0; iFacetIndex < aObject.get_num_facets(); ++iFacetIndex )
        {
            // Get the difference of the cast point and the facet min and max coords in the J direction
            real tJDirectionPreselection = ( aPoint( tFirstAxis ) - aObject.get_facet_min_coord( iFacetIndex, tFirstAxis ) )
                                         * ( aObject.get_facet_max_coord( iFacetIndex, tFirstAxis ) - aPoint( tFirstAxis ) );

            // check if the point is inside the bounding box of the triangle
            if ( tJDirectionPreselection > -tEpsilon )
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
        aCandidateFacets.resize( aObject.get_num_facets() );

        // loop over remaining triangles in I-direction
        for ( uint k = 0; k < tCountJ; ++k )
        {
            // Get the difference of the cast point and the facet min and max coords in the K direction
            real tKDirectionPreselection = ( aPoint( tSecondAxis ) - aObject.get_facet_min_coord( tCandJ( k ), tSecondAxis ) )
                                         * ( aObject.get_facet_max_coord( tCandJ( k ), tSecondAxis ) - aPoint( tSecondAxis ) );

            // see if the point is inside the bounding box of the triangle
            if ( tKDirectionPreselection > -tEpsilon )
            {
                aCandidateFacets( tCount ) = tCandJ( k );
                ++tCount;
            }
        }

        aCandidateFacets.resize( tCount );
    }

    //-------------------------------------------------------------------------------

    Preselection_Result preselect_lines(
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Vector< uint >&         aIntersectedFacets,
            Vector< Facet* >&       aCandidateFacets )
    {
        Preselection_Result tPreselectionSuccessful = SUCCESS;

        // Ensure the function is being called for the proper number of facets
        MORIS_ASSERT( aPoint.numel() == 2,
                "SDF_preselect_lines() should be called for 2D problems only. Query point dimension = %lu",
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
            // Get a reference to this line and the Object's intersection tolerance
            real   tEpsilon = aObject.get_intersection_tolerance();
            Facet& tLine    = aObject.get_facet( iLineIndex );

            // get the difference of the cast point and the facet min and max coords in the !aAxis direction
            real tMaxCoordOffAxisDifference = tLine.get_max_coord( tOtherAxis ) - aPoint( tOtherAxis );
            real tMinCoordOffAxisDifference = tLine.get_min_coord( tOtherAxis ) - aPoint( tOtherAxis );

            // get the difference of the cast point and the facet min and max coords in the aAxis direction
            real tMaxCoordAxisDifference = tLine.get_max_coord( aAxis ) - aPoint( aAxis );
            real tMinCoordAxisDifference = tLine.get_min_coord( aAxis ) - aPoint( aAxis );

            Matrix< DDRMat > tLineCoords = tLine.get_vertex_coords();

            // the cast point is very close to a vertex
            if ( ( std::abs( tLineCoords( 0, aAxis ) - aPoint( aAxis ) ) < tEpsilon
                         and std::abs( tLineCoords( 0, tOtherAxis ) - aPoint( tOtherAxis ) ) < tEpsilon )
                    or ( std::abs( tLineCoords( 1, aAxis ) - aPoint( aAxis ) ) < tEpsilon
                            and std::abs( tLineCoords( 1, tOtherAxis ) - aPoint( tOtherAxis ) ) < tEpsilon ) )
            {
                // give only the facet whose vertex the cast point lies on
                aCandidateFacets( 0 ) = &aObject.get_facet( iLineIndex );
                aCandidateFacets.resize( 1 );
                aIntersectedFacets.resize( 0 );

                return FAIL_ON_VERTEX;
            }
            // the facet is along an axis, and the point is on it
            else if ( ( std::abs( tMaxCoordAxisDifference ) < tEpsilon
                              and std::abs( tMinCoordAxisDifference ) < tEpsilon and tMaxCoordOffAxisDifference * tMinCoordOffAxisDifference < MORIS_REAL_EPS )
                      or ( std::abs( tMinCoordOffAxisDifference ) < tEpsilon
                              and std::abs( tMaxCoordOffAxisDifference ) < tEpsilon and tMaxCoordAxisDifference * tMinCoordAxisDifference < MORIS_REAL_EPS ) )
            {
                // give only the facet whose vertex the cast point lies on
                aCandidateFacets( 0 ) = &aObject.get_facet( iLineIndex );
                aCandidateFacets.resize( 1 );
                aIntersectedFacets.resize( 0 );

                return FAIL_ON_EDGE;
            }
            // the ray will hit a vertex, but the cast point is not on a vertex
            else if ( std::abs( tMaxCoordOffAxisDifference ) < tEpsilon or std::abs( tMinCoordOffAxisDifference ) < tEpsilon )
            {
                tPreselectionSuccessful = FAIL_CAST_TO_VERTEX;
            }

            // check bounding box of the line against the point (point is above min coord and below max coord)
            if ( tMaxCoordOffAxisDifference * tMinCoordOffAxisDifference < MORIS_REAL_EPS )
            {
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
            Vector< uint >&         aCandidateFacets,
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis )
    {

        // get number of candidate triangles
        uint tNumberOfCandidateFacets = aCandidateFacets.size();

        // initialize counter for intersected triangles
        uint tCount = 0;

        // Vector to check each edge
        Vector< uint > tEdges = { 0, 1, 2 };

        // Determine how many facets are intersected and flag them accordingly
        for ( uint iCandidateFacetIndex : aCandidateFacets )
        {
            // get reference to triangle
            Facet& tFacet = aObject.get_facet( iCandidateFacetIndex );

            // check each edge of the triangle and flag it if the point is on the correct side for every edge
            if ( std::all_of( tEdges.cbegin(), tEdges.cend(), [ &tFacet, aAxis, &aPoint ]( uint iEdge ) { return tFacet.check_edge( iEdge, aAxis, aPoint ); } ) )
            {
                tFacet.flag();
                ++tCount;
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
    intersect_triangles_moller_trumbore(
            Vector< uint >&         aCandidateFacets,
            Object&                 aObject,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Vector< Facet* >&       aIntersectedFacets )
    {
        // get number of facets
        uint tNumberOfFacets = aObject.get_num_facets();

        // initialize vector with coords in axis
        Vector< real > tCoordsK( tNumberOfFacets );
        aIntersectedFacets.resize( tNumberOfFacets );

        // counter for number of valid intersections
        uint tCount = 0;

        // loop over all intersected triangles and find intersection point
        for ( uint iCandidateFacet : aCandidateFacets )
        {
            // return the intersection coordinate from moller trumbore algorithm
            real tIntersection = dynamic_cast< sdf::Triangle& >( aObject.get_facet( iCandidateFacet ) ).moller_trumbore( aAxis, aPoint );

            // if the ray intersects the triangle, add the intersection to the list
            if ( not std::isnan( tIntersection ) )
            {
                tCoordsK( tCount )             = tIntersection;
                aIntersectedFacets( tCount++ ) = &aObject.get_facet( iCandidateFacet );
            }
        }

        // check that there are valid intersection to process, and if not, return
        if ( tCount == 0 )
        {
            return {};
        }

        // resize coord array
        tCoordsK.resize( tCount );

        // Get the indices of the coordinate array
        Vector< uint > tCoordsIndices( tCoordsK.size() );
        std::iota( tCoordsIndices.begin(), tCoordsIndices.end(), 0 );

        // sort the indices of the array based on the intersection values
        std::sort( tCoordsIndices.begin(), tCoordsIndices.end(), [ &tCoordsK ]( uint i, uint j ) {
            return tCoordsK( i ) < tCoordsK( j );
        } );

        // rearrange both tCoordsK and the facets based on the sort
        Vector< real >   tCoordsKSorted( tCoordsK.size() );
        Vector< Facet* > tIntersectedFacetsSorted( aIntersectedFacets.size() );
        for ( uint iIntersectionIndex = 0; iIntersectionIndex < tCoordsK.size(); iIntersectionIndex++ )
        {
            tCoordsKSorted( iIntersectionIndex )           = tCoordsK( tCoordsIndices( iIntersectionIndex ) );
            tIntersectedFacetsSorted( iIntersectionIndex ) = aIntersectedFacets( tCoordsIndices( iIntersectionIndex ) );
        }

        // make result unique
        uint tCountUnique = 0;

        // initialize return vector
        Vector< real > tIntersectionCoords( tCount );

        // set first entry
        tIntersectionCoords( tCountUnique )  = tCoordsKSorted( 0 );
        aIntersectedFacets( tCountUnique++ ) = tIntersectedFacetsSorted( 0 );

        // find unique entries
        for ( uint k = 1; k < tCount; ++k )
        {
            if ( std::abs( tCoordsKSorted( k ) - tCoordsKSorted( k - 1 ) ) > 10 * MORIS_REAL_EPS )
            {
                tIntersectionCoords( tCountUnique ) = tCoordsKSorted( k );
                aIntersectedFacets( tCountUnique )  = tIntersectedFacetsSorted( k );

                tCountUnique++;
            }
        }

        // chop vector
        tIntersectionCoords.resize( tCountUnique );
        aIntersectedFacets.resize( tCountUnique );

        return tIntersectionCoords;
    }


    //-------------------------------------------------------------------------------

    Vector< real >
    intersect_ray_with_facets(
            Vector< Facet* >&       aIntersectedFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            Preselection_Result     aPreselection,
            bool                    aIgnoreErrors )
    {
        // return nothing if there are no intersected facets
        if ( aIntersectedFacets.size() == 0 )
        {
            return {};
        }

        // get number of facets
        uint tNumberOfFacets = aIntersectedFacets.size();

        // initialize vector with coords in axis
        Vector< real > tCoordsK( tNumberOfFacets );

        // counter for number of valid intersections
        uint tCount = 0;

        // return flag that is true if a facet is ever parallel to the ray, meaning our intersection coordinate is invalid
        bool tError;

        // check if a pathological case was determined
        switch ( aPreselection )
        {
            case Preselection_Result::FAIL_CAST_TO_VERTEX:
            case Preselection_Result::FAIL_CAST_TO_EDGE:
            case Preselection_Result::SUCCESS:
                // loop over all intersected triangles and find intersection point
                for ( uint iFacetIndex = 0; iFacetIndex < tNumberOfFacets; ++iFacetIndex )
                {
                    // return value for the following function call that is the intersection coordinate
                    real tCoordK;

                    // the ray will intersect the facet, or the facet is parallel to the ray. either case is handled here
                    aIntersectedFacets( iFacetIndex )->intersect_with_coordinate_axis( aPoint, aAxis, tCoordK, tError );

                    // error meant we would have divided by zero. This triangle is ignored
                    // otherwise, the value is written into the result vector
                    if ( not tError )
                    {
                        tCoordsK( tCount++ ) = std::round( tCoordK / MORIS_REAL_EPS ) * MORIS_REAL_EPS;
                    }
                    // if there is an error but we still want to compute intersections, make sure that this intersection is ignored
                    else if ( aIgnoreErrors )
                    {
                        // remove the facet from the list
                        aIntersectedFacets( iFacetIndex ) = nullptr;
                    }
                    // if there is an error and we want to know about it, return such that this ray will be ignored
                    else
                    {
                        // this way, the matrix is ignored
                        return { std::numeric_limits< real >::quiet_NaN() };
                    }
                }
                break;
            case FAIL_ON_VERTEX:
            case FAIL_ON_EDGE:
                // the cast point is on a vertex, return its coordinate
                return { aPoint( aAxis ) };
        }

        // move the nullptrs to the end of the vector and get the iterator for the index
        auto tItFacet = std::remove_if( aIntersectedFacets.begin(), aIntersectedFacets.end(), []( Facet* aFacet ) { return aFacet == nullptr; } );

        // erase the nullptrs
        aIntersectedFacets.erase_by_iterators( tItFacet, aIntersectedFacets.end() );

        // check that there are valid intersection to process, and if not, return
        if ( tCount == 0 )
        {
            return {};
        }

        // resize coord array
        tCoordsK.resize( tCount );

        // Get the indices of the coordinate array
        Vector< uint > tCoordsIndices( tCoordsK.size() );
        std::iota( tCoordsIndices.begin(), tCoordsIndices.end(), 0 );

        // sort the indices of the array based on the intersection values
        std::sort( tCoordsIndices.begin(), tCoordsIndices.end(), [ &tCoordsK ]( uint i, uint j ) {
            return tCoordsK( i ) < tCoordsK( j );
        } );

        // rearrange both tCoordsK and the facets based on the sort
        Vector< real >   tCoordsKSorted( tCoordsK.size() );
        Vector< Facet* > tIntersectedFacetsSorted( aIntersectedFacets.size() );
        for ( uint iIntersectionIndex = 0; iIntersectionIndex < tCoordsK.size(); iIntersectionIndex++ )
        {
            tCoordsKSorted( iIntersectionIndex )           = tCoordsK( tCoordsIndices( iIntersectionIndex ) );
            tIntersectedFacetsSorted( iIntersectionIndex ) = aIntersectedFacets( tCoordsIndices( iIntersectionIndex ) );
        }

        // make result unique
        uint tCountUnique = 0;

        // initialize return vector
        Vector< real > tIntersectionCoords( tCount );

        // set first entry
        tIntersectionCoords( tCountUnique )  = tCoordsKSorted( 0 );
        aIntersectedFacets( tCountUnique++ ) = tIntersectedFacetsSorted( 0 );

        // find unique entries
        for ( uint k = 1; k < tCount; ++k )
        {
            if ( std::abs( tCoordsKSorted( k ) - tCoordsKSorted( k - 1 ) ) > 10 * MORIS_REAL_EPS )
            {
                tIntersectionCoords( tCountUnique ) = tCoordsKSorted( k );
                aIntersectedFacets( tCountUnique )  = tIntersectedFacetsSorted( k );

                tCountUnique++;
            }
        }

        // chop vector
        tIntersectionCoords.resize( tCountUnique );
        aIntersectedFacets.resize( tCountUnique );

        return tIntersectionCoords;
    }

    //-------------------------------------------------------------------------------

    Object_Region
    check_if_node_is_inside_triangles(
            Vector< real >&         aIntersectionCoords,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis,
            const real&             aIntersectionTolerance )
    {
        uint tNumCoordsK = aIntersectionCoords.size();

        real tAxisCastPoint = aPoint( aAxis );

        // If the ray intersected no facets, the point is outside
        if ( tNumCoordsK == 0 )
        {
            return Object_Region::OUTSIDE;
        }
        // Check if any of the intersections are very close to the point
        else if ( std::any_of( aIntersectionCoords.begin(), aIntersectionCoords.end(), [ & ]( real aIntersection ) { return std::abs( aIntersection - tAxisCastPoint ) < 10.0 * aIntersectionTolerance; } ) )
        {
            return Object_Region::INTERFACE;
        }
        // only even number of intersections is considered
        else if ( tNumCoordsK % 2 == 0 )
        {
            uint tNodeIsInside;
            for ( uint iIntersectionIndex = 0; iIntersectionIndex < tNumCoordsK / 2; ++iIntersectionIndex )
            {
                // check if the point is in between two intersections
                tNodeIsInside = ( aPoint( aAxis ) > aIntersectionCoords( 2 * iIntersectionIndex ) )
                            and ( aPoint( aAxis ) < aIntersectionCoords( 2 * iIntersectionIndex + 1 ) );
            }

            return static_cast< Object_Region >( tNodeIsInside );
        }
        else
        {
            return Object_Region::UNSURE;
        }
    }

    //-------------------------------------------------------------------------------

    Object_Region
    check_if_node_is_inside_lines(
            Object&                 aObject,
            const Vector< real >&   aIntersectionCoords,
            const Vector< uint >&   aCandidateFacets,
            const Matrix< DDRMat >& aPoint,
            uint                    aAxis )
    {
        // get length of the number of intersections computed
        uint tNumCoordsK = aIntersectionCoords.size();

        // if the ray intersected no facets, the point is outside
        if ( tNumCoordsK == 0 and aCandidateFacets.size() == 0 )
        {
            return OUTSIDE;
        }

        // check if the location of the intersection is greater than the location of the coordinate
        uint tIntersectionsRightOfPoint = 0;
        for ( uint iIntersectionIndex = 0; iIntersectionIndex < tNumCoordsK; iIntersectionIndex++ )
        {
            if ( std::abs( aIntersectionCoords( iIntersectionIndex ) - aPoint( aAxis ) ) < 10.0 * aObject.get_intersection_tolerance() )
            {
                return INTERFACE;
            }
            else if ( aIntersectionCoords( iIntersectionIndex ) - aPoint( aAxis ) > MORIS_REAL_EPS )
            {
                tIntersectionsRightOfPoint++;
            }
        }

        // if the ray passes through an even number of lines, the point is inside. otherwise, it is outside.
        return static_cast< Object_Region >( ( tIntersectionsRightOfPoint + aCandidateFacets.size() ) % 2 );
    }

    void random_rotation(
            Object&           aObject,
            Matrix< DDRMat >& aPoint )
    {
        MORIS_ASSERT( aPoint.n_cols() == 1 or aPoint.n_rows() == 1,
                "SDF-Raycast::random_rotation(): aPoint expected vector, but matrix received. Dimensions: %lux%lu",
                aPoint.n_rows(),
                aPoint.n_cols() );

        // ensure aPoint is a column vector
        if ( aPoint.n_rows() == 1 )
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
