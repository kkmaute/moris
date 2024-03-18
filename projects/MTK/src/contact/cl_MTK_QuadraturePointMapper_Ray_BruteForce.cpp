//
// Created by frank on 11/28/23.
//

#include <deque>
#include <set>
#include "cl_MTK_QuadraturePointMapper_Ray_BruteForce.hpp"

#include "cl_MTK_QuadraturePointMapper.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_Vector.hpp"
#include "fn_assert.hpp"


#include <cl_MTK_Ray_Line_Intersection.hpp>

namespace moris::mtk
{


    MappingResult QuadraturePointMapper_Ray_BruteForce::map( moris_index const aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates ) const
    {
        Surface_Mesh const    tSurfaceMesh = get_surface_meshes()( aSourceSideSetIndex );
        Side_Set const *const tSideSet     = get_side_sets()( aSourceSideSetIndex );

        // skip, if the side set is empty
        if ( tSideSet->get_num_clusters_on_set() == 0 )
        {
            MORIS_LOG_WARNING( "Side set '%s' is empty. Skipping it in Contact Detection", tSideSet->get_set_name().c_str() );
            return { aSourceSideSetIndex, tSideSet->get_spatial_dim(), 0 };
        }

        // initialize the mapping result with the correct size and the parametric coordinates and normals on each cell
        MappingResult tMappingResult = initialize_source_points( aSourceSideSetIndex, aParametricCoordinates );

        Spatial_Indexing_Result const tSpatialIndexingResult = mSpatialIndexer.perform( aSourceSideSetIndex, 1.0 );

        // for each cell in the surface mesh (source mesh), interpolate the parametric coordinates and normals of the rays
        // and populate the results of the ray-casting in the MappingResult instance
        for ( moris_index iCell = 0; iCell < static_cast< moris_index >( tSurfaceMesh.get_number_of_cells() ); iCell++ )
        {
            this->raycast_cell( iCell, aParametricCoordinates.n_cols(), tMappingResult, tSpatialIndexingResult );
        }

        return tMappingResult;
    }

    void QuadraturePointMapper_Ray_BruteForce::raycast_cell(
            moris_index const              aSourceCellIndex,
            moris_index const              aNumberOfRays,
            MappingResult                 &aMappingResult,
            Spatial_Indexing_Result const &aSpatialIndexingResult ) const
    {
        // calculate the start index of the current cell i.e. the index at which the rays of the current cell start in the MappingResult
        uint const tResultOffset = aSourceCellIndex * aNumberOfRays;

        // use a deque to store unprocessed rays
        std::deque< moris_index > tUnprocessedRays( aNumberOfRays );
        std::iota( tUnprocessedRays.begin(), tUnprocessedRays.end(), 0 );

        // first, try to map the rays using the spatial indexing results.
        // In this case, the information about the closest vertices to each vertex of the current cell is used to find potential target cells.
        bool tBruteForce = false;
        process_rays(
                aSourceCellIndex,
                aSpatialIndexingResult,
                tResultOffset,
                aMappingResult,
                tUnprocessedRays,
                tBruteForce );

        if ( !tUnprocessedRays.empty() )
        {
            // if there are still unprocessed rays, try to map them using brute force
            // i.e. check all cells of all meshes. This is slow and inefficient and should be replaced by a better algorithm.
            tBruteForce = true;
            process_rays(
                    aSourceCellIndex,
                    aSpatialIndexingResult,
                    tResultOffset,
                    aMappingResult,
                    tUnprocessedRays,
                    tBruteForce );
        }

        if ( !tUnprocessedRays.empty() )
        {
            // MORIS_LOG_INFO( "Cell %d: %zu rays could not be mapped.", aSourceCellIndex, tUnprocessedRays.size() );
        }
    }

    void QuadraturePointMapper_Ray_BruteForce::process_rays(
            moris_index                    aSourceCellIndex,
            Spatial_Indexing_Result const &aSpatialIndexingResult,
            uint                           aResultOffset,
            MappingResult                 &aMappingResult,
            std::deque< moris_index >     &aUnprocessedRays,
            bool                           aBruteForce ) const
    {

        moris_index const &tSourceMeshIndex = aMappingResult.mSourceMeshIndex;

        std::set< moris_index > const tPotentialTargetMeshesIndices = get_potential_target_meshes( aBruteForce, tSourceMeshIndex, aSourceCellIndex, aSpatialIndexingResult );

        for ( auto const iTargetMeshIndex : tPotentialTargetMeshesIndices )
        {
            std::set< moris_index > tPotentialTargetCells = get_potential_target_cells( aBruteForce, tSourceMeshIndex, iTargetMeshIndex, aSourceCellIndex, aSpatialIndexingResult );
            std::set< moris_index > tCheckedCells;

            for ( moris_index iTargetCell : tPotentialTargetCells )
            {
                // if the cell has not been checked yet
                if ( tCheckedCells.find( iTargetCell ) == tCheckedCells.end() )
                {
                    tCheckedCells.insert( iTargetCell );
                    check_ray_cell_intersection(
                            aMappingResult,
                            aUnprocessedRays,
                            tSourceMeshIndex,
                            iTargetMeshIndex,
                            aSourceCellIndex,
                            iTargetCell,
                            aResultOffset );
                }
                // if all rays have been processed, stop the loop
                if ( aUnprocessedRays.empty() ) { break; }
            }    // end loop over neighboring cells
            // if all rays have been processed, stop the loop
            if ( aUnprocessedRays.empty() ) { break; }
        }    // end loop over source vertices
    }

    std::set< moris_index > QuadraturePointMapper_Ray_BruteForce::get_potential_target_cells(
            bool const                     aBruteForce,
            moris_index const              aSourceMeshIndex,
            moris_index const              tTargetMeshIndex,
            moris_index const              aSourceCellIndex,
            Spatial_Indexing_Result const &aSpatialIndexingResult ) const
    {
        Surface_Mesh const     &tTargetMesh = get_surface_meshes()( tTargetMeshIndex );
        std::set< moris_index > tPotentialTargetCells;
        if ( aBruteForce )
        {
            for ( moris_index iCell = 0; iCell < static_cast< moris_index >( tTargetMesh.get_number_of_cells() ); iCell++ ) { tPotentialTargetCells.insert( iCell ); }
        }
        else
        {
            // get the neighboring cells of the closest vertex to both vertices of the current cell on this target cell
            auto tSegmentVertices = get_surface_meshes()( aSourceMeshIndex ).get_vertices_of_cell( aSourceCellIndex );
            for ( auto const tVertex : tSegmentVertices )
            {
                // only consider the vertices that are actually closest to the target mesh
                if ( aSpatialIndexingResult[ tVertex ].mesh_index == tTargetMeshIndex )
                {
                    moris_index const     tClosestVertex    = aSpatialIndexingResult[ tVertex ].vertex;
                    Vector< moris_index > tNeighboringCells = tTargetMesh.get_cells_of_vertex( tClosestVertex );
                    tPotentialTargetCells.insert( tNeighboringCells.begin(), tNeighboringCells.end() );
                }
            }
        }
        return tPotentialTargetCells;
    }

    std::set< moris_index > QuadraturePointMapper_Ray_BruteForce::get_potential_target_meshes(
            bool const                     aBruteForce,
            moris_index const              aSourceMeshIndex,
            moris_index const              aSourceCellIndex,
            Spatial_Indexing_Result const &aSpatialIndexingResult ) const
    {
        std::set< moris_index > tPotentialTargetMeshesIndices;
        if ( aBruteForce )
        {
            for ( auto const &[ tSourceCandidate, tTargetCandidate ] : get_candidate_pairs() )
            {
                if ( tSourceCandidate == aSourceMeshIndex ) { tPotentialTargetMeshesIndices.insert( tTargetCandidate ); }
            }
        }
        else
        {
            // check the meshes that contain the closest vertices of all vertices of the current cell
            // e.g. a triangle with three vertices could have a different closest mesh for each of the three vertices
            auto tSegmentVertices = get_surface_meshes()( aSourceMeshIndex ).get_vertices_of_cell( aSourceCellIndex );
            for ( auto const tVertex : tSegmentVertices ) { tPotentialTargetMeshesIndices.insert( aSpatialIndexingResult[ tVertex ].mesh_index ); }
        }
        return tPotentialTargetMeshesIndices;
    }

    void QuadraturePointMapper_Ray_BruteForce::check_ray_cell_intersection(
            MappingResult             &aMappingResult,
            std::deque< moris_index > &aUnprocessedRays,
            moris_index const          aSourceMeshIndex,
            moris_index const          aTargetMeshIndex,
            moris_index const          aSourceCellIndex,
            moris_index const          aTargetCellIndex,
            uint const                 aResultOffset ) const
    {
        Surface_Mesh const &tTargetMesh = get_surface_meshes()( aTargetMeshIndex );

        // get the basic information from the cell like the vertex coordinates and calculate the origin and direction of the cell facet
        // i.e. the line segment between the first and second vertex which will be called "segment" in the following
        Matrix< DDRMat > tTargetCellCoordinates = tTargetMesh.get_vertex_coordinates_of_cell( aTargetCellIndex );

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

        // Keep the initial size of the deque because it might get smaller over time (because rays are removed from the deque if they got mapped)
        // If they are not mapped, they are pushed back into the deque. To prevent an infinite loop, the initial size is stored (i.e. each ray is processed once).
        uint const tNumUnprocessedRays = aUnprocessedRays.size();
        for ( uint tCounter = 0; tCounter < tNumUnprocessedRays; tCounter++ )
        {
            moris_index iRay = aUnprocessedRays.front();
            aUnprocessedRays.pop_front();

            tRayLineIntersection.set_ray_origin( aMappingResult.mSourcePhysicalCoordinate.get_column( aResultOffset + iRay ) );
            tRayLineIntersection.set_ray_direction( aMappingResult.mNormals.get_column( aResultOffset + iRay ) );
            tRayLineIntersection.perform_raytracing();

            if ( tRayLineIntersection.has_intersection() )
            {
                // it has an intersection: the ray information can be stored in the MappingResult
                // the results will be inserted relative to the start index of the current cell
                uint const tInsertIndex = aResultOffset + iRay;

                aMappingResult.mTargetParametricCoordinate.set_column( tInsertIndex, tRayLineIntersection.get_intersection_parametric() );
                aMappingResult.mTargetPhysicalCoordinate.set_column( tInsertIndex, tRayLineIntersection.get_intersection_physical() );

                aMappingResult.mDistances( tInsertIndex )            = tRayLineIntersection.get_ray_length();
                aMappingResult.mTargetSideSetIndices( tInsertIndex ) = aTargetMeshIndex;

                //                aMappingResult.mSourceCellIndex( tInsertIndex )   = tSourceMesh.get_global_cell_index( aSourceCellIndex ); // TOOD @ff remove if not needed
                aMappingResult.mTargetCellIndices( tInsertIndex ) = tTargetMesh.get_global_cell_index( aTargetCellIndex );

                //                aMappingResult.mSourceClusterIndex( tInsertIndex ) = tSourceMesh.get_cluster_of_cell( aSourceCellIndex ); // TOOD @ff remove if not needed
                aMappingResult.mTargetClusterIndex( tInsertIndex ) = tTargetMesh.get_cluster_of_cell( aTargetCellIndex );
            }
            else
            {
                // push back into the unprocessed deque because no intersection was found
                aUnprocessedRays.push_back( iRay );
            }
        }
    }
}    // namespace moris::mtk
