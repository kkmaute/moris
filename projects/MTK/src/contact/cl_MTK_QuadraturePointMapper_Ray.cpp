//
// Created by frank on 11/28/23.
//

#include <deque>
#include <set>
#include "cl_MTK_QuadraturePointMapper_Ray.hpp"

#include <math.h>
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_Json_Object.hpp"

namespace moris::mtk
{
    QuadraturePointMapper_Ray::QuadraturePointMapper_Ray(
            Integration_Mesh                                      *aIGMesh,
            Vector< Side_Set const * >                            &aSideSets,
            Vector< std::pair< moris_index, moris_index > > const &aCandidatePairs )
            : QuadraturePointMapper( aIGMesh, aSideSets, aCandidatePairs )
            , mSurfaceMeshes( initialize_surface_meshes( aIGMesh, aSideSets ) )
            , mSpatialIndexer( mSurfaceMeshes, mCandidatePairs )
    {
        Json  tSurfaceMeshes;
        auto &tMeshes = tSurfaceMeshes.put_child( "surface_meshes", Json() );
        for ( auto const &tSurfaceMesh : mSurfaceMeshes )
        {
            tMeshes.push_back( { "", tSurfaceMesh.to_json() } );
        }
        write_json( "surface_meshes.json", tSurfaceMeshes );
    }

    Vector< Surface_Mesh > QuadraturePointMapper_Ray::initialize_surface_meshes(
            Integration_Mesh const           *aIGMesh,
            Vector< Side_Set const * > const &aSideSets )
    {
        Vector< Surface_Mesh > tSurfaceMeshes;
        for ( auto const &tSideSet : aSideSets )
        {
            // initialize one surface mesh per side set
            Vector< mtk::Side_Set const * > tSideSetCast{ tSideSet };
            Surface_Mesh                    tSurfaceMesh( aIGMesh, tSideSetCast );
            tSurfaceMeshes.push_back( tSurfaceMesh );
        }
        return tSurfaceMeshes;
    }

    MappingResult QuadraturePointMapper_Ray::map( moris_index const aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates ) const
    {
        Side_Set const *const tSideSet = mSideSets( aSourceSideSetIndex );

        // skip, if the side set is empty
        if ( tSideSet->get_num_clusters_on_set() == 0 )
        {
            std::cout << "Side set " << tSideSet->get_set_name() << " is empty. Skipping it in Contact Detection" << std::endl;
            return MappingResult( aSourceSideSetIndex, tSideSet->get_spatial_dim(), 0 );
        }

        Interpolation_Rule const tInterpolationRule(
                tSideSet->get_integration_cell_geometry_type(),
                Interpolation_Type::LAGRANGE,
                Interpolation_Order::LINEAR,
                Interpolation_Type::UNDEFINED,
                Interpolation_Order::UNDEFINED );

        // using two Interpolator instances for faster processing of each point
        // (because the coefficients do not have to be changed between calculation of coordinates and normals each time
        Space_Interpolator            tCoordinateInterpolator( tInterpolationRule );
        Space_Interpolator            tNormalInterpolator( tInterpolationRule );
        Surface_Mesh const            tSurfaceMesh           = mSurfaceMeshes( aSourceSideSetIndex );
        Spatial_Indexing_Result const tSpatialIndexingResult = mSpatialIndexer.perform( aSourceSideSetIndex, 1.0 );

        // preallocate the MappingResult
        auto const tNumParametricCoordinates = static_cast< moris_index >( aParametricCoordinates.n_cols() );
        uint const tDim                      = tSideSet->get_spatial_dim();
        uint const tNumCells                 = tSurfaceMesh.get_number_of_cells();
        uint const tTotalNumPoints           = tNumCells * tNumParametricCoordinates;

        MappingResult tMappingResult( aSourceSideSetIndex, tDim, tTotalNumPoints );

        // for each cell in the surface mesh (source mesh), interpolate the parametric coordinates and normals of the rays
        // and populate the results of the ray-casting in the MappingResult instance
        for ( moris_index iCell = 0; iCell < static_cast< moris_index >( tNumCells ); iCell++ )
        {
            QuadraturePointMapper_Ray::interpolate_source_point(
                    iCell,
                    aParametricCoordinates,
                    tCoordinateInterpolator,
                    tNormalInterpolator,
                    tSurfaceMesh,
                    tMappingResult );


            this->raycast_cell(
                    iCell,
                    tNumParametricCoordinates,
                    tMappingResult,
                    tSpatialIndexingResult );
        }

        return tMappingResult;
    }

    void QuadraturePointMapper_Ray::interpolate_source_point(
            moris_index const       aCellIndex,
            Matrix< DDRMat > const &aParametricCoordinates,
            Space_Interpolator     &aCoordinateInterpolator,
            Space_Interpolator     &aNormalInterpolator,
            Surface_Mesh const     &aSurfaceMesh,
            MappingResult          &aMappingResult )
    {
        // initialize the coordinate-interpolator with the vertex coordinates
        uint const tNumRays    = aParametricCoordinates.n_cols();
        uint const tStartIndex = aCellIndex * tNumRays;

        Matrix< DDRMat > const tVertexCoordinates = aSurfaceMesh.get_vertex_coordinates_of_cell( aCellIndex );
        aCoordinateInterpolator.set_space_coeff( tVertexCoordinates );

        Matrix< DDRMat > const tVertexNormals = aSurfaceMesh.get_vertex_normals_of_cell( aCellIndex );
        aNormalInterpolator.set_space_coeff( tVertexNormals );

        for ( uint iPoint = 0; iPoint < tNumRays; iPoint++ )
        {
            // set the parametric coordinate of both interpolators
            aCoordinateInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
            aNormalInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );

            // map the parametric coordinate to the surface mesh
            aMappingResult.mSourcePhysicalCoordinate.set_column( tStartIndex + iPoint, tVertexCoordinates * trans( aCoordinateInterpolator.NXi() ) );
            aMappingResult.mNormals.set_column( tStartIndex + iPoint, tVertexNormals * trans( aNormalInterpolator.NXi() ) );
        }
    }

    void QuadraturePointMapper_Ray::raycast_cell(
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
            MORIS_LOG_INFO( "Cell %d: %zu rays could not be mapped.", aSourceCellIndex, tUnprocessedRays.size() );
        }
    }

    void QuadraturePointMapper_Ray::process_rays(
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

    std::set< moris_index > QuadraturePointMapper_Ray::get_potential_target_cells(
            bool const                     aBruteForce,
            moris_index const              aSourceMeshIndex,
            moris_index const              tTargetMeshIndex,
            moris_index const              aSourceCellIndex,
            Spatial_Indexing_Result const &aSpatialIndexingResult ) const
    {
        Surface_Mesh const     &tTargetMesh = mSurfaceMeshes( tTargetMeshIndex );
        std::set< moris_index > tPotentialTargetCells;
        if ( aBruteForce )
        {
            for ( moris_index iCell = 0; iCell < static_cast< moris_index >( tTargetMesh.get_number_of_cells() ); iCell++ ) { tPotentialTargetCells.insert( iCell ); }
        }
        else
        {
            // get the neighboring cells of the closest vertex to both vertices of the current cell on this target cell
            auto tSegmentVertices = mSurfaceMeshes( aSourceMeshIndex ).get_vertices_of_cell( aSourceCellIndex );
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

    std::set< moris_index > QuadraturePointMapper_Ray::get_potential_target_meshes(
            bool const                     aBruteForce,
            moris_index const              aSourceMeshIndex,
            moris_index const              aSourceCellIndex,
            Spatial_Indexing_Result const &aSpatialIndexingResult ) const
    {
        std::set< moris_index > tPotentialTargetMeshesIndices;
        if ( aBruteForce )
        {
            for ( auto const &[ tSourceCandidate, tTargetCandidate ] : mCandidatePairs )
            {
                if ( tSourceCandidate == aSourceMeshIndex ) { tPotentialTargetMeshesIndices.insert( tTargetCandidate ); }
            }
        }
        else
        {
            // check the meshes that contain the closest vertices of all vertices of the current cell
            // e.g. a triangle with three vertices could have a different closest mesh for each of the three vertices
            auto tSegmentVertices = mSurfaceMeshes( aSourceMeshIndex ).get_vertices_of_cell( aSourceCellIndex );
            for ( auto const tVertex : tSegmentVertices ) { tPotentialTargetMeshesIndices.insert( aSpatialIndexingResult[ tVertex ].mesh_index ); }
        }
        return tPotentialTargetMeshesIndices;
    }

    void QuadraturePointMapper_Ray::check_ray_cell_intersection(
            MappingResult             &aMappingResult,
            std::deque< moris_index > &aUnprocessedRays,
            moris_index const          aSourceMeshIndex,
            moris_index const          aTargetMeshIndex,
            moris_index const          aSourceCellIndex,
            moris_index const          aTargetCellIndex,
            uint const                 aResultOffset ) const
    {
        Surface_Mesh const &tSourceMesh = mSurfaceMeshes( aSourceMeshIndex );
        Surface_Mesh const &tTargetMesh = mSurfaceMeshes( aTargetMeshIndex );

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

        // Keep the initial size of the deque because it might get smaller over time (because rays are removed from the deque if they got mapped)
        // If they are not mapped, they are pushed back into the deque. To prevent an infinite loop, the initial size is stored (i.e. each ray is processed once).
        uint const tNumUnprocessedRays = aUnprocessedRays.size();
        for ( uint tCounter = 0; tCounter < tNumUnprocessedRays; tCounter++ )
        {
            moris_index iRay = aUnprocessedRays.front();
            aUnprocessedRays.pop_front();

            // the next call is the actual ray-casting. The resulting array contains the following information:
            // 0. bool: true if the ray intersects with the segment, false otherwise
            // 1. real: the distance between the ray origin and the intersection point
            // 2. Matrix< DDRMat >: the parametric coordinate of the intersection point
            // 3. Matrix< DDRMat >: the physical coordinate of the intersection point (mainly for debugging purposes)
            auto const &[ tHasIntersection, tDistance, tParamCoord, tPhysCoord ] =
                    calculate_ray_line_intersection(
                            aMappingResult.mSourcePhysicalCoordinate.get_column( aResultOffset + iRay ),
                            aMappingResult.mNormals.get_column( aResultOffset + iRay ),
                            tSegmentOrigin,
                            tSegmentDirection );

            if ( tHasIntersection )
            {
                // it has an intersection: the ray information can be stored in the MappingResult
                // the results will be inserted relative to the start index of the current cell
                uint const tInsertIndex = aResultOffset + iRay;

                aMappingResult.mTargetParametricCoordinate.set_column( tInsertIndex, tParamCoord );
                aMappingResult.mTargetPhysicalCoordinate.set_column( tInsertIndex, tPhysCoord );

                aMappingResult.mDistances( tInsertIndex )            = tDistance;
                aMappingResult.mTargetSideSetIndices( tInsertIndex ) = aTargetMeshIndex;

                aMappingResult.mSourceCellIndex( tInsertIndex )   = tSourceMesh.get_global_cell_index( aSourceCellIndex );
                aMappingResult.mTargetCellIndices( tInsertIndex ) = tTargetMesh.get_global_cell_index( aTargetCellIndex );

                aMappingResult.mSourceClusterIndex( tInsertIndex ) = tSourceMesh.get_cluster_of_cell( aSourceCellIndex );
                aMappingResult.mTargetClusterIndex( tInsertIndex ) = tTargetMesh.get_cluster_of_cell( aTargetCellIndex );
            }
            else
            {
                // push back into the unprocessed deque because no intersection was found
                aUnprocessedRays.push_back( iRay );
            }
        }
    }

    /**
     * @brief Implementation according to https://stackoverflow.com/a/2932601
     */
    std::tuple< bool, real, Matrix< DDRMat >, Matrix< DDRMat > > QuadraturePointMapper_Ray::calculate_ray_line_intersection(
            Matrix< DDRMat > const &aRayOrigin,
            Matrix< DDRMat > const &aRayDirection,
            Matrix< DDRMat > const &aSegmentOrigin,
            Matrix< DDRMat > const &aSegmentDirection )
    {
        // Nomenclature:
        //     aRayOrigin: r
        //     aRayDirection: dr
        //        -> p(u) = r + u * dr
        //     aSegmentOrigin: s
        //     aSegmentDirection: ds = e - s (between line segment vertices s and e)
        //        -> q(v) = s + v * ds

        // distance between origins r and s
        // dOrigins = s - r
        Matrix< DDRMat > dOrigins = aSegmentOrigin - aRayOrigin;

        // determinant of the matrix D = [dr, ds]
        real const detD = aSegmentDirection( 0 ) * aRayDirection( 1 ) - aSegmentDirection( 1 ) * aRayDirection( 0 );

        // For 2D faces (e.g. triangles), the parametric coordinate would be 2x1. We therefore store the result in a matrix and not as a scalar.
        // TODO: For memory efficiency, this could be changed...
        Matrix< DDRMat > tParametricCoordinate{ 1, 1 };
        Matrix< DDRMat > tPhsyicalCoordinate{ 2, 1 };
        bool             tHasIntersection = false;
        real             tGap             = 0.0;

        // if detD is zero, the lines are parallel and no calculation is necessary
        if ( std::abs( detD ) > 1e-16 )
        {
            // scaling factor u and v for the ray and the line segment, respectively
            real u = NAN;
            real v = NAN;

            // u = (dOrigins_y * ds_x - dOrigins_x * ds_y) / detD
            u = ( dOrigins( 1 ) * aSegmentDirection( 0 ) - dOrigins( 0 ) * aSegmentDirection( 1 ) ) / detD;

            // v = (dOrigins_y * dr_x - dOrigins_x * dr_y) / detD
            v = ( dOrigins( 1 ) * aRayDirection( 0 ) - dOrigins( 0 ) * aRayDirection( 1 ) ) / detD;

            // if v is between 0 and 1, the intersection point is on the line segment
            if ( v >= 0.0 && v <= 1.0 )
            {
                tHasIntersection = true;
                // q(v) = s + v * ds

                // calculate the gap between the ray and the line segment
                tGap = norm( u * aRayDirection );

                // calculate the physical coordinate of the intersection point
                tPhsyicalCoordinate = aRayOrigin + u * aRayDirection;

                // the parametric coordinate goes from -1 to 1 and has the center in the middle of the line segment
                tParametricCoordinate( 0 ) = 2.0 * ( v - 0.5 );
            }
        }
        return { tHasIntersection, tGap, tParametricCoordinate, tPhsyicalCoordinate };
    }
} // namespace moris::mtk
