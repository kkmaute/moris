//
// Created by frank on 11/28/23.
//

#include <deque>
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
            , mSpatialIndexer( mSurfaceMeshes, mCandidatePairs )
    {
        for ( auto &tSideSet : aSideSets )
        {
            // initialize one surface mesh per side set
            moris::Cell< mtk::Side_Set * > tSideSetCast{ tSideSet };
            Surface_Mesh                   tSurfaceMesh( aIGMesh, tSideSetCast );
            mSurfaceMeshes.push_back( tSurfaceMesh );
        }
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
        uint tNumSideSets              = mSideSets.size();
        uint tNumParametricCoordinates = aParametricCoordinates.n_cols();
        uint tDim                      = aParametricCoordinates.n_rows();
        uint tNumTotalCells            = std::accumulate(
                mSurfaceMeshes.begin(),
                mSurfaceMeshes.end(),
                0,
                []( uint const &aSum, Surface_Mesh const &aSurfaceMesh ) {
                    return aSum + aSurfaceMesh.get_number_of_cells();
                } );
        uint tTotalNumPoints = tNumTotalCells * tNumParametricCoordinates;

        moris::Cell< Spatial_Indexing_Result > tSpatialIndexingResult = mSpatialIndexer.perform( 0.0 );

        MappingResult tMappingResult{
            Matrix< DDRMat >( tDim, tTotalNumPoints ),
            moris::Cell< moris_index >( tNumSideSets ),
            moris::Cell< moris_index >( tTotalNumPoints ),
            moris::Cell< moris_index >( tTotalNumPoints )
        };


        //        Surface_Mesh tSurfaceMesh = mSurfaceMeshes( 0 );
        //        moris_index iCell = 1;
        //
        //        Matrix< DDRMat > tVertexCoordinates = tSurfaceMesh.get_vertex_coordinates_of_cell( iCell );
        //        tCoordinateInterpolator.set_space_coeff( tVertexCoordinates );
        //
        //        Matrix< DDRMat > tVertexNormals = tSurfaceMesh.get_vertex_normals_of_cell( iCell );
        //        tNormalInterpolator.set_space_coeff( tVertexNormals );
        //
        //        size_t iPoint = 0;
        //        // set the parametric coordinate of both interpolators
        //        tCoordinateInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
        //        tNormalInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
        //
        //        // map the parametric coordinate to the surface mesh
        //        Matrix< DDRMat > tInterpolatedCoordinate = tVertexCoordinates * trans( tCoordinateInterpolator.NXi() );
        //        Matrix< DDRMat > tInterpolatedNormal     = tVertexNormals * trans( tNormalInterpolator.NXi() );
        //
        //        std::string tMat;
        //        print_as_row_vector( tInterpolatedCoordinate, tMat );
        //        std::cout << tMat << std::endl;
        //
        //        print_as_row_vector( tInterpolatedNormal, tMat );
        //        std::cout << tMat << std::endl;
        //
        //        Surface_Mesh tTargetMesh = mSurfaceMeshes( 1 );
        //        moris_index  tTargetCell = 2;
        //
        //        Matrix< DDRMat > tTargetVertexCoordinates = tTargetMesh.get_vertex_coordinates_of_cell( tTargetCell );
        //
        //        Matrix< DDRMat > tSegmentOrigin    = tTargetVertexCoordinates.get_column( 0 );
        //        Matrix< DDRMat > tSegmentDirection = tTargetVertexCoordinates.get_column( 1 ) - tSegmentOrigin;
        //
        //        auto tResult = calculate_ray_line_intersection(
        //                tInterpolatedCoordinate,
        //                tInterpolatedNormal,
        //                tSegmentOrigin,
        //                tSegmentDirection );
        //        std::cout << "Has intersection: " << tResult.first << std::endl;


        // for each surface mesh
        moris_index tCurrentPointIndex = 0;
        for ( moris_index iSurfaceMesh = 0; iSurfaceMesh < mSurfaceMeshes.size(); iSurfaceMesh++ )
        {
            Surface_Mesh tSurfaceMesh = mSurfaceMeshes( iSurfaceMesh );
            tMappingResult.mSideSetOffsets.push_back( tCurrentPointIndex );
            // for each cell in the surface mesh
            auto tNumCells = static_cast< moris_index >( tSurfaceMesh.get_number_of_cells() );
            tCurrentPointIndex += tNumCells * static_cast< moris_index >( tNumParametricCoordinates );
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
                        tSpatialIndexingResult( iSurfaceMesh ),
                        tInterpolatedCoordinates,
                        tInterpolatedNormals );
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
            Matrix< DDRMat >( tDim, tNumPoints ),
            moris::Cell< moris_index >( 0 ),
            moris::Cell< moris_index >( tNumPoints ),
            moris::Cell< moris_index >( tNumPoints )
        };

        // use a deque to store unprocessed rays
        std::deque< moris_index > tUnprocessedRays;
        moris_index               tClosestVertex = 0;
        moris_index               tClosestMesh   = 0;

        for ( moris_index iSourceVertex : aVertexIndices )
        {
            tClosestVertex = aSpatialIndexingResult[ iSourceVertex ].vertex;
            tClosestMesh   = aSpatialIndexingResult[ iSourceVertex ].mesh_index;

            // get the neighboring cells of the closest vertex
            moris::Cell< moris_index > tNeighborCells = mSurfaceMeshes( tClosestMesh ).get_vertex_neighbors( tClosestVertex );
        }


        // map the points to the target mesh


        return tMappingResult;
    }

    /**
     * @brief Implementation according to https://stackoverflow.com/a/2932601
     */
    std::pair< bool, double > ContactMapper_Ray::calculate_ray_line_intersection(
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
        bool   tHasIntersection      = false;
        double tParametricCoordinate = 0.0;

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
                tParametricCoordinate = 2.0 * ( v - 0.5 );
            }
        }
        return { tHasIntersection, tParametricCoordinate };
    }
}    // namespace moris::mtk