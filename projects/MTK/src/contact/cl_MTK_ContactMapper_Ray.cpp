//
// Created by frank on 11/28/23.
//

#include <deque>
#include <set>
#include "cl_MTK_ContactMapper_Ray.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"

namespace moris::mtk
{
    ContactMapper_Ray::ContactMapper_Ray(
            mtk::Integration_Mesh                                      *aIGMesh,
            moris::Cell< mtk::Side_Set * >                             &aSideSets,
            moris::Cell< std::pair< moris_index, moris_index > > const &aCandidatePairs )
            : ContactMapper( aIGMesh, aSideSets, aCandidatePairs )
            , mSurfaceMeshes( initialize_surface_meshes( aIGMesh, aSideSets ) )
            , mSpatialIndexer( mSurfaceMeshes, mCandidatePairs )
    {
    }

    moris::Cell< Surface_Mesh > ContactMapper_Ray::initialize_surface_meshes( Integration_Mesh *aIGMesh, moris::Cell< mtk::Side_Set * > const &aSideSets )
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

    ContactMapper::MappingResult ContactMapper_Ray::map( Matrix< DDRMat > const &aParametricCoordinates )
    {
        Interpolation_Rule tInterpolationRule(
                mSideSets( 0 )->get_integration_cell_geometry_type(),
                Interpolation_Type::LAGRANGE,
                Interpolation_Order::LINEAR,
                Interpolation_Type::UNDEFINED,
                Interpolation_Order::UNDEFINED );

        // using two Interpolator instances for faster processing of each point
        // (because the coefficients do not have to be changed between calculation of coordinates and normals each time
        Space_Interpolator tCoordinateInterpolator( tInterpolationRule );
        Space_Interpolator tNormalInterpolator( tInterpolationRule );

        // preallocate the MappingResult
        uint tNumSurfaceMeshes         = mSurfaceMeshes.size();
        uint tNumParametricCoordinates = aParametricCoordinates.n_cols();
        uint tDim                      = mIGMesh->get_spatial_dim();
        uint tNumTotalCells            = std::accumulate(
                mSurfaceMeshes.begin(),
                mSurfaceMeshes.end(),
                0,
                []( uint const &aSum, Surface_Mesh const &aSurfaceMesh ) {
                    return aSum + aSurfaceMesh.get_number_of_cells();
                } );
        uint tTotalNumPoints = tNumTotalCells * tNumParametricCoordinates;

        moris::Cell< Spatial_Indexing_Result > tSpatialIndexingResults = mSpatialIndexer.perform( 1.0 );

        MappingResult tMappingResult{
            Matrix< DDRMat >( tDim - 1, tTotalNumPoints ),    // parametric coordinate is one dimension lower than physical coordinate
            moris::Cell< moris_index >( tNumSurfaceMeshes ),
            moris::Cell< moris_index >( tTotalNumPoints ),
            moris::Cell< moris_index >( tTotalNumPoints )
        };

        // for each source candidate side set
        moris_index tCurrentOffset = 0;
        for ( auto const &tCandidatePair : mCandidatePairs )
        {
            moris_index  iSurfaceMesh = tCandidatePair.first;
            Surface_Mesh tSurfaceMesh = mSurfaceMeshes( iSurfaceMesh );
            tMappingResult.mSideSetOffsets.push_back( tCurrentOffset );
            // for each cell in the surface mesh
            auto tNumCells = static_cast< moris_index >( tSurfaceMesh.get_number_of_cells() );
            for ( moris_index iCell = 0; iCell < tNumCells; iCell++ )
            {
                // initialize the coordinate-interpolator with the vertex coordinates
                Matrix< DDRMat > tVertexCoordinates = tSurfaceMesh.get_vertex_coordinates_of_cell( iCell );
                tCoordinateInterpolator.set_space_coeff( tVertexCoordinates );

                // initialize the normal-interpolator with the vertex normals
                Matrix< DDRMat > tVertexNormals = tSurfaceMesh.get_vertex_normals_of_cell( iCell );
                tNormalInterpolator.set_space_coeff( tVertexNormals );

                Matrix< DDRMat > tInterpolatedCoordinates{ tDim, tNumParametricCoordinates };
                Matrix< DDRMat > tInterpolatedNormals{ tDim, tNumParametricCoordinates };
                for ( size_t iPoint = 0; iPoint < aParametricCoordinates.n_cols(); iPoint++ )
                {
                    // set the parametric coordinate of both interpolators
                    tCoordinateInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
                    tNormalInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );

                    // map the parametric coordinate to the surface mesh
                    tInterpolatedCoordinates.set_column( iPoint, tVertexCoordinates * trans( tCoordinateInterpolator.NXi() ) );
                    tInterpolatedNormals.set_column( iPoint, tVertexNormals * trans( tNormalInterpolator.NXi() ) );
                }

                MappingResult tCellMappingResult = raycast_bundle(
                        iSurfaceMesh,
                        iCell,
                        tSurfaceMesh.get_vertices_of_cell( iCell ),
                        tSpatialIndexingResults( iSurfaceMesh ),
                        tInterpolatedCoordinates,
                        tInterpolatedNormals );

                // copy the results into the (global) MappingResult
                for ( uint i = 0; i < tNumParametricCoordinates; i++ )
                {
                    tMappingResult.mParametricCoordinates.set_column( i + tCurrentOffset, tCellMappingResult.mParametricCoordinates.get_column( i ) );
                    tMappingResult.mCellIndices( i + tCurrentOffset )    = tCellMappingResult.mCellIndices( i );
                    tMappingResult.mSideSetIndices( i + tCurrentOffset ) = tCellMappingResult.mSideSetIndices( i );
                }
                tCurrentOffset += tNumParametricCoordinates;
            }
        }

        return ContactMapper::MappingResult();
    }

    ContactMapper::MappingResult ContactMapper_Ray::raycast_bundle(
            moris_index                       aSourceMesh,
            moris_index                       aCellIndex,
            moris::Cell< moris_index > const &aVertexIndices,
            Spatial_Indexing_Result const    &aSpatialIndexingResult,
            Matrix< DDRMat > const           &aPhysicalCoordinates,
            Matrix< DDRMat > const           &aNormals )
    {
        // initialize the result variable
        size_t        tNumPoints = aPhysicalCoordinates.n_cols();
        size_t        tDim       = aPhysicalCoordinates.n_rows();
        MappingResult tMappingResult{
            Matrix< DDRMat >( tDim - 1, tNumPoints ),    // parametric coordinate is one dimension lower than physical coordinate
            moris::Cell< moris_index >( 0 ),
            moris::Cell< moris_index >( tNumPoints, -1 ),
            moris::Cell< moris_index >( tNumPoints, -1 )
        };

        // use a deque to store unprocessed rays
        std::deque< moris_index > tUnprocessedRays( aPhysicalCoordinates.n_cols() );
        std::iota( tUnprocessedRays.begin(), tUnprocessedRays.end(), 0 );

        std::set< moris_index > tCheckedCells;
        moris_index             tClosestVertex    = 0;
        moris_index             tClosestMeshIndex = 0;

        for ( moris_index iSourceVertex : aVertexIndices )
        {
            tClosestVertex    = aSpatialIndexingResult[ iSourceVertex ].vertex;
            tClosestMeshIndex = aSpatialIndexingResult[ iSourceVertex ].mesh_index;

            // get the neighboring cells of the closest vertex
            Surface_Mesh const        &tMesh          = mSurfaceMeshes( tClosestMeshIndex );
            moris::Cell< moris_index > tNeighborCells = tMesh.get_cells_of_vertex( tClosestVertex );
            for ( moris_index iNeighborCell : tNeighborCells )
            {
                // if the cell has not been checked yet
                if ( tCheckedCells.find( iNeighborCell ) == tCheckedCells.end() )
                {
                    Matrix< DDRMat > tCellVertexCoordinates = tMesh.get_vertex_coordinates_of_cell( iNeighborCell );
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

                        std::pair< bool, Matrix< DDRMat > > tIntersection = calculate_ray_line_intersection(
                                aPhysicalCoordinates.get_column( iRay ),
                                aNormals.get_column( iRay ),
                                tSegmentOrigin,
                                tSegmentDirection );

                        if ( tIntersection.first )
                        {    // it has an intersection
                            tMappingResult.mParametricCoordinates.set_column( iRay, tIntersection.second );
                            tMappingResult.mCellIndices( iRay )    = iNeighborCell;
                            tMappingResult.mSideSetIndices( iRay ) = tClosestMeshIndex;
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

        if ( tUnprocessedRays.size() > 0 )
        {
            std::cout << "Surface mesh: " << aSourceMesh
                      << ", cell: " << aCellIndex << ": "
                      << tUnprocessedRays.size() << " rays could not be mapped onto neighboring cells!" << std::endl;
        }
        return tMappingResult;
    }

    /**
     * @brief Implementation according to https://stackoverflow.com/a/2932601
     */
    std::pair< bool, Matrix< DDRMat > > ContactMapper_Ray::calculate_ray_line_intersection(
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

        // initialize results
        bool tHasIntersection = false;

        // For 2D faces (e.g. triangles), the parametric coordinate would be 2x1. We therefore store the result in a matrix and not as a scalar.
        // TODO: For memory efficiency, this could be changed...
        Matrix< DDRMat > tParametricCoordinate{ 1, 1 };

        // if detD is zero, the lines are parallel and no calculation is necessary
        if ( std::abs( detD ) > 1e-12 )
        {
            // scaling factor u and v for the ray and the line segment, respectively
            // real u           = 0.0;
            real v = 0.0;

            // u = (dO_y * ds_x - dO_x * ds_y) / detD
            // since only v is needed to find the intersection point and to check if the intersection point is on the line segment, u is not calculated:
            //            u = ( dO( 1 ) * aVertexDirection( 0 ) - dO( 0 ) * aVertexDirection( 1 ) ) / detD;

            // v = (dO_y * dr_x - dO_x * dr_y) / detD
            v = ( dO( 1 ) * aRayDirection( 0 ) - dO( 0 ) * aRayDirection( 1 ) ) / detD;

            // if v is between 0 and 1, the intersection point is on the line segment
            if ( v >= 0.0 && v <= 1.0 )
            {
                tHasIntersection = true;
                // q(v) = s + v * ds

                // the parametric coordinate goes from -1 to 1 and has the center in the middle of the line segment
                tParametricCoordinate( 0 ) = 2.0 * ( v - 0.5 );
            }
        }
        return { tHasIntersection, tParametricCoordinate };
    }

}    // namespace moris::mtk