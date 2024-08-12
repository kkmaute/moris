/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper_Ray_ArborX.cpp
 *
 */

#include "cl_MTK_QuadraturePointMapper_Ray_ArborX.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "fn_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Ray_Line_Intersection.hpp"
#include "moris_typedefs.hpp"
#include "cl_Tracer.hpp"

namespace moris::mtk
{
    MappingResult QuadraturePointMapper_ArborX::map(
            moris_index             aSourceMeshIndex,
            Matrix< DDRMat > const &aParametricCoordinates,
            real                    aMaxNegativeRayLength,
            real                    aMaxPositiveRayLength ) const
    {
        Tracer                         tTracer( "Quadrature Point Mapper", "Map", "Map Quadrature Points" );
        Integration_Surface_Mesh const tSurfaceMesh = get_surface_meshes()( aSourceMeshIndex );
        Side_Set const *const          tSideSet     = get_side_sets()( aSourceMeshIndex );

        // skip, if the side set is empty
        if ( tSideSet->get_num_clusters_on_set() == 0 )
        {
            MORIS_LOG_WARNING( "Side set '%s' is empty. Skipping it in Contact Detection", tSideSet->get_set_name().c_str() );
            return { aSourceMeshIndex, tSideSet->get_spatial_dim(), 0 };
        }

        // initialize the mapping result with the correct size and the parametric coordinates and normals on each cell
        MappingResult tMappingResult = initialize_source_points( aSourceMeshIndex, aParametricCoordinates );
        auto const   &tBoxRayMap     = moris::mtk::arborx::map_rays_to_boxes( tMappingResult, get_target_surface_meshes( aSourceMeshIndex ) );

        // check the intersections of the rays with the target cells
        // since ArborX is only able to tell if a ray intersects with the bounding box of a cell, we need to check the intersection with the actual cell (e.g. line segment)
        check_cell_intersections( tMappingResult, aMaxNegativeRayLength, aMaxPositiveRayLength, tBoxRayMap );

        return tMappingResult;
    }
    void QuadraturePointMapper_ArborX::check_cell_intersections(
            MappingResult                  &tMappingResult,
            real                            aMaxNegativeRayLength,
            real                            aMaxPositiveRayLength,
            arborx::cell_locator_map const &tBoxRayMap ) const
    {
        Tracer tTracer( "Quadrature Point Mapper", "Map", "Check Cell Intersections" );
        for ( auto const &[ tTargetMeshIndex, tTargetCells ] : tBoxRayMap )
        {
            Integration_Surface_Mesh const &tTargetMesh = get_surface_meshes()( tTargetMeshIndex );
            for ( auto const &[ tTargetCellIndex, tRayIndices ] : tTargetCells )
            {
                // get the basic information from the cell like the vertex coordinates and calculate the origin and direction of the cell facet
                // i.e. the line segment between the first and second vertex which will be called "segment" in the following
                Matrix< DDRMat > tTargetCellCoordinates = tTargetMesh.get_all_vertex_coordinates_of_facet( tTargetCellIndex );

                /* Because the segments will always be oriented in opposing directions (e.g. the vertices of each triangle will be ordered counter-clockwise),
                 * the parametric coordinate will also be measured in opposing directions.
                 *              2\
                 *              │  \  Source
                 *              │    \
                 *  Source      │      \
                 *  Param.      0────────1
                 *  Direction-- ──────────► xi
                 *
                 *
                 *  Target ---- ◄────────── xi
                 *  Param.      1────────0 --Segment Origin
                 *  Direction   │      /
                 *              │    /
                 *              │  /  Target
                 *              2/
                 */
                Matrix< DDRMat > const tSegmentOrigin    = tTargetCellCoordinates.get_column( 0 );
                Matrix< DDRMat > const tSegmentDirection = tTargetCellCoordinates.get_column( 1 ) - tSegmentOrigin;
                Ray_Line_Intersection  tRayLineIntersection( tSegmentOrigin.n_rows() );
                tRayLineIntersection.set_target_origin( tSegmentOrigin );
                tRayLineIntersection.set_target_span( tSegmentDirection );
                for ( auto const &tRayIndex : tRayIndices )
                {
                    tRayLineIntersection.set_ray_origin( tMappingResult.mSourcePhysicalCoordinate.get_column( tRayIndex ) );
                    tRayLineIntersection.set_ray_direction( tMappingResult.mNormals.get_column( tRayIndex ) );
                    tRayLineIntersection.perform_raytracing();
                    if ( tRayLineIntersection.has_intersection()                                                               // check if the ray intersects the line segment
                            && tRayLineIntersection.get_signed_ray_length() > aMaxNegativeRayLength                            // check that the ray is not too long in the negative direction
                            && tRayLineIntersection.get_signed_ray_length() < aMaxPositiveRayLength                            //
                            && ( tRayLineIntersection.get_signed_ray_length() < tMappingResult.mSignedDistance( tRayIndex )    // check if the intersection is closer than the previous one
                                    || tMappingResult.mTargetCellIndices( tRayIndex ) == -1 ) )                                // or if the ray has not intersected anything before (initial distance is 0.0)
                    {
                        tMappingResult.mTargetParametricCoordinate.set_column( tRayIndex, tRayLineIntersection.get_intersection_parametric() );
                        tMappingResult.mTargetPhysicalCoordinate.set_column( tRayIndex, tRayLineIntersection.get_intersection_physical() );
                        tMappingResult.mSignedDistance( tRayIndex )       = tRayLineIntersection.get_signed_ray_length();
                        tMappingResult.mTargetSideSetIndices( tRayIndex ) = tTargetMeshIndex;
                        tMappingResult.mTargetCellIndices( tRayIndex )    = tTargetMesh.get_global_cell_index( tTargetCellIndex );
                        tMappingResult.mTargetClusterIndex( tRayIndex )   = tTargetMesh.get_cluster_of_cell( tTargetCellIndex );
                    }
                }
            }
        }
    }

    Vector< std::pair< moris_index, Surface_Mesh > > QuadraturePointMapper_ArborX::get_target_surface_meshes( moris_index aSourceMeshIndex ) const
    {
        Vector< std::pair< moris_index, Surface_Mesh > > tTargetSurfaceMeshes;
        for ( auto const &[ tSourceCandidateIndex, tTargetCandidateIndex ] : get_candidate_pairs() )
        {
            if ( tSourceCandidateIndex == aSourceMeshIndex )
            {
                tTargetSurfaceMeshes.push_back( { tTargetCandidateIndex, get_surface_meshes()( tTargetCandidateIndex ) } );
            }
        }
        return tTargetSurfaceMeshes;
    }

}    // namespace moris::mtk
