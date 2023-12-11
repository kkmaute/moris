//
// Created by frank on 11/28/23.
//

#include <deque>
#include <set>
#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_Json_Object.hpp"

namespace moris::mtk
{
    QuadraturePointMapper_Ray::QuadraturePointMapper_Ray(
            mtk::Integration_Mesh                                      *aIGMesh,
            moris::Cell< mtk::Side_Set * >                             &aSideSets,
            moris::Cell< std::pair< moris_index, moris_index > > const &aCandidatePairs )
            : QuadraturePointMapper( aIGMesh, aSideSets, aCandidatePairs )
            , mSurfaceMeshes( initialize_surface_meshes( aIGMesh, aSideSets ) )
            , mSpatialIndexer( mSurfaceMeshes, mCandidatePairs )
    {
        Json  tSurfaceMeshes;
        auto &tMeshes = tSurfaceMeshes.put_child( "surface_meshes", Json() );
        for ( auto tSurfaceMesh : mSurfaceMeshes )
        {
            tMeshes.push_back( { "", tSurfaceMesh.to_json() } );
        }
        write_json( "surface_meshes.json", tSurfaceMeshes );
    }

    moris::Cell< Surface_Mesh > QuadraturePointMapper_Ray::initialize_surface_meshes( Integration_Mesh *aIGMesh, moris::Cell< mtk::Side_Set * > const &aSideSets )
    {
        moris::Cell< Surface_Mesh > tSurfaceMeshes;
        for ( auto &tSideSet : aSideSets )
        {
            // initialize one surface mesh per side set
            moris::Cell< mtk::Side_Set * > tSideSetCast{ tSideSet };
            Surface_Mesh                   tSurfaceMesh( aIGMesh, tSideSetCast );
            tSurfaceMeshes.push_back( tSurfaceMesh );
        }
        return tSurfaceMeshes;
    }

    MappingResult QuadraturePointMapper_Ray::map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates )
    {
        Interpolation_Rule tInterpolationRule(
                mSideSets( aSourceSideSetIndex )->get_integration_cell_geometry_type(),
                Interpolation_Type::LAGRANGE,
                Interpolation_Order::LINEAR,
                Interpolation_Type::UNDEFINED,
                Interpolation_Order::UNDEFINED );

        // using two Interpolator instances for faster processing of each point
        // (because the coefficients do not have to be changed between calculation of coordinates and normals each time
        Space_Interpolator      tCoordinateInterpolator( tInterpolationRule );
        Space_Interpolator      tNormalInterpolator( tInterpolationRule );
        Surface_Mesh            tSurfaceMesh           = mSurfaceMeshes( aSourceSideSetIndex );
        Spatial_Indexing_Result tSpatialIndexingResult = mSpatialIndexer.perform( aSourceSideSetIndex, 1.0 );

        // preallocate the MappingResult
        auto tNumParametricCoordinates = static_cast< moris_index >( aParametricCoordinates.n_cols() );
        uint tDim                      = mSideSets( aSourceSideSetIndex )->get_spatial_dim();
        uint tNumCells                 = tSurfaceMesh.get_number_of_cells();
        uint tTotalNumPoints           = tNumCells * tNumParametricCoordinates;

        MappingResult tMappingResult( tDim, tTotalNumPoints );

        // for each cell in the surface mesh
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
                    tSurfaceMesh.get_vertices_of_cell( iCell ),
                    tSpatialIndexingResult );
        }

        return tMappingResult;
    }

    void QuadraturePointMapper_Ray::interpolate_source_point(
            moris_index             aCellIndex,
            Matrix< DDRMat > const &aParametricCoordinates,
            Space_Interpolator     &aCoordinateInterpolator,
            Space_Interpolator     &aNormalInterpolator,
            Surface_Mesh const     &aSurfaceMesh,
            MappingResult          &aMappingResult )
    {    // initialize the coordinate-interpolator with the vertex coordinates
        uint tNumRays    = aParametricCoordinates.n_cols();
        uint tStartIndex = aCellIndex * tNumRays;

        Matrix< DDRMat > tVertexCoordinates = aSurfaceMesh.get_vertex_coordinates_of_cell( aCellIndex );
        aCoordinateInterpolator.set_space_coeff( tVertexCoordinates );

        Matrix< DDRMat > tVertexNormals = aSurfaceMesh.get_vertex_normals_of_cell( aCellIndex );
        aNormalInterpolator.set_space_coeff( tVertexNormals );

        for ( uint iPoint = 0; iPoint < tNumRays; iPoint++ )
        {
            // set the parametric coordinate of both interpolators
            aCoordinateInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
            aNormalInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );

            // map the parametric coordinate to the surface mesh
            aMappingResult.mSourcePhysicalCoordinate.set_column( tStartIndex + iPoint, tVertexCoordinates * trans( aCoordinateInterpolator.NXi() ) );
            aMappingResult.mNormal.set_column( tStartIndex + iPoint, tVertexNormals * trans( aNormalInterpolator.NXi() ) );
        }
    }

    void QuadraturePointMapper_Ray::raycast_cell(
            moris_index                       aCellIndex,
            moris_index                       aNumberOfRays,
            MappingResult                    &aMappingResult,
            moris::Cell< moris_index > const &aVertexIndices,
            Spatial_Indexing_Result const    &aSpatialIndexingResult )
    {
        // calculate the start index of the current cell i.e. the index at which the rays of the current cell start in the MappingResult
        uint tStartIndex = aCellIndex * aNumberOfRays;

        // use a deque to store unprocessed rays
        std::deque< moris_index > tUnprocessedRays( aNumberOfRays );
        std::iota( tUnprocessedRays.begin(), tUnprocessedRays.end(), 0 );

        std::set< moris_index > tCheckedCells;
        moris_index             tClosestVertex;
        moris_index             tClosestMeshIndex;

        for ( moris_index iSourceVertex : aVertexIndices )
        {
            tClosestVertex    = aSpatialIndexingResult[ iSourceVertex ].vertex;
            tClosestMeshIndex = aSpatialIndexingResult[ iSourceVertex ].mesh_index;

            // get the neighboring cells of the closest vertex
            Surface_Mesh const        &tSurfaceMesh   = mSurfaceMeshes( tClosestMeshIndex );
            moris::Cell< moris_index > tNeighborCells = tSurfaceMesh.get_cells_of_vertex( tClosestVertex );
            for ( moris_index iNeighborCell : tNeighborCells )
            {
                // if the cell has not been checked yet
                if ( tCheckedCells.find( iNeighborCell ) == tCheckedCells.end() )
                {
                    Matrix< DDRMat > tCellVertexCoordinates = tSurfaceMesh.get_vertex_coordinates_of_cell( iNeighborCell );
                    Matrix< DDRMat > tSegmentOrigin         = tCellVertexCoordinates.get_column( 0 );
                    Matrix< DDRMat > tSegmentDirection      = tCellVertexCoordinates.get_column( 1 ) - tSegmentOrigin;
                    tCheckedCells.insert( iNeighborCell );

                    // keep the initial size of the deque because it might get smaller over time
                    uint tNumUnprocessedRays = tUnprocessedRays.size();

                    // check all rays with the current cell
                    for ( uint tCounter = 0; tCounter < tNumUnprocessedRays; tCounter++ )
                    {
                        moris_index iRay = tUnprocessedRays.front();
                        tUnprocessedRays.pop_front();

                        std::tuple< bool, real, Matrix< DDRMat >, Matrix< DDRMat > > tIntersection =
                                calculate_ray_line_intersection(
                                        aMappingResult.mSourcePhysicalCoordinate.get_column( tStartIndex + iRay ),
                                        aMappingResult.mNormal.get_column( tStartIndex + iRay ),
                                        tSegmentOrigin,
                                        tSegmentDirection );

                        if ( std::get< 0 >( tIntersection ) )
                        {    // it has an intersection
                            aMappingResult.mDistances( tStartIndex + iRay ) = std::get< 1 >( tIntersection );
                            aMappingResult.mTargetParametricCoordinate.set_column( tStartIndex + iRay, std::get< 2 >( tIntersection ) );
                            aMappingResult.mTargetPhysicalCoordinate.set_column( tStartIndex + iRay, std::get< 3 >( tIntersection ) );
                            aMappingResult.mTargetCellIndices( tStartIndex + iRay )    = iNeighborCell;
                            aMappingResult.mTargetSideSetIndices( tStartIndex + iRay ) = tClosestMeshIndex;
                        }
                        else
                        {    // push back into the unprocessed deque
                            tUnprocessedRays.push_back( iRay );
                        }
                    }
                }
                // if all rays have been processed, stop the loop
                if ( tUnprocessedRays.empty() )
                {
                    break;
                }
            }    // end loop over neighboring cells
            // if all rays have been processed, stop the loop
            if ( tUnprocessedRays.empty() )
            {
                break;
            }
        }    // end loop over source vertices

        // get the neighboring cells of the closest vertex
        for ( uint iSurfaceMesh = 0; iSurfaceMesh < mSurfaceMeshes.size(); iSurfaceMesh++ )
        {
            Surface_Mesh const &tSurfaceMesh = mSurfaceMeshes( tClosestMeshIndex );
            for ( uint iCell = 0; iCell < tSurfaceMesh.get_number_of_cells(); iCell++ )
            {
                // if the cell has not been checked yet
                if ( tCheckedCells.find( iCell ) == tCheckedCells.end() )
                {
                    Matrix< DDRMat > tCellVertexCoordinates = tSurfaceMesh.get_vertex_coordinates_of_cell( iCell );
                    Matrix< DDRMat > tSegmentOrigin         = tCellVertexCoordinates.get_column( 0 );
                    Matrix< DDRMat > tSegmentDirection      = tCellVertexCoordinates.get_column( 1 ) - tSegmentOrigin;
                    tCheckedCells.insert( iCell );

                    // keep the initial size of the deque because it might get smaller over time
                    uint tNumUnprocessedRays = tUnprocessedRays.size();

                    // check all rays with the current cell
                    for ( uint tCounter = 0; tCounter < tNumUnprocessedRays; tCounter++ )
                    {
                        moris_index iRay = tUnprocessedRays.front();
                        tUnprocessedRays.pop_front();

                        std::tuple< bool, real, Matrix< DDRMat >, Matrix< DDRMat > > tIntersection =
                                calculate_ray_line_intersection(
                                        aMappingResult.mSourcePhysicalCoordinate.get_column( tStartIndex + iRay ),
                                        aMappingResult.mNormal.get_column( tStartIndex + iRay ),
                                        tSegmentOrigin,
                                        tSegmentDirection );

                        if ( std::get< 0 >( tIntersection ) )
                        {    // it has an intersection
                            aMappingResult.mDistances( tStartIndex + iRay ) = std::get< 1 >( tIntersection );
                            aMappingResult.mTargetParametricCoordinate.set_column( tStartIndex + iRay, std::get< 2 >( tIntersection ) );
                            aMappingResult.mTargetPhysicalCoordinate.set_column( tStartIndex + iRay, std::get< 3 >( tIntersection ) );
                            aMappingResult.mTargetCellIndices( tStartIndex + iRay )    = iCell;
                            aMappingResult.mTargetSideSetIndices( tStartIndex + iRay ) = iSurfaceMesh;
                        }
                        else
                        {    // push back into the unprocessed deque
                            tUnprocessedRays.push_back( iRay );
                        }
                    }
                }
                if ( tUnprocessedRays.empty() )
                {
                    break;
                }
            }
            if ( tUnprocessedRays.empty() )
            {
                break;
            }
        }

        if ( !tUnprocessedRays.empty() )
        {
            MORIS_LOG_INFO( "Cell %d: %zu rays could not be mapped.", aCellIndex, tUnprocessedRays.size() );
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
        // dO = s - r
        Matrix< DDRMat > dO = aSegmentOrigin - aRayOrigin;

        // determinant of the matrix D = [dr, ds]
        real detD = aSegmentDirection( 0 ) * aRayDirection( 1 ) - aSegmentDirection( 1 ) * aRayDirection( 0 );

        // For 2D faces (e.g. triangles), the parametric coordinate would be 2x1. We therefore store the result in a matrix and not as a scalar.
        // TODO: For memory efficiency, this could be changed...
        Matrix< DDRMat > tParametricCoordinate{ 1, 1 };
        Matrix< DDRMat > tPhsyicalCoordinate{ 2, 1 };
        bool             tHasIntersection = false;
        real             tGap             = 0.0;

        // if detD is zero, the lines are parallel and no calculation is necessary
        if ( std::abs( detD ) > 1e-12 )
        {
            // scaling factor u and v for the ray and the line segment, respectively
            real u;
            real v;

            // u = (dO_y * ds_x - dO_x * ds_y) / detD
            u = ( dO( 1 ) * aSegmentDirection( 0 ) - dO( 0 ) * aSegmentDirection( 1 ) ) / detD;

            // v = (dO_y * dr_x - dO_x * dr_y) / detD
            v = ( dO( 1 ) * aRayDirection( 0 ) - dO( 0 ) * aRayDirection( 1 ) ) / detD;

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
}    // namespace moris::mtk