/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper_Ray_ArborX_Details.cpp
 *
 */
#include "cl_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_Tracer.hpp"
#include <unordered_map>
#include <functional>

namespace moris::mtk::arborx
{
    template< typename MemorySpace, typename ExecutionSpace >
    QueryBoxes< MemorySpace > construct_query_boxes(
            ExecutionSpace const                                                      &aExecutionSpace,
            moris::Vector< std::pair< moris_index, moris::mtk::Surface_Mesh > > const &aTargetSurfaceMeshes )
    {
        uint const tNumCells = std::accumulate( aTargetSurfaceMeshes.begin(), aTargetSurfaceMeshes.end(), 0, []( auto a, const auto &b ) { return a + b.second.get_number_of_cells(); } );

        Kokkos::View< ArborX::Box *, MemorySpace > tBoxes( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:boxes" ), tNumCells );
        Kokkos::View< moris_index *, MemorySpace > tMeshIndices( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:mesh_indices" ), tNumCells );
        Kokkos::View< moris_index *, MemorySpace > tCellIndices( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:cell_indices" ), tNumCells );

        moris_index tBoxIndex = 0;
        for ( size_t iMeshIndex = 0; iMeshIndex < aTargetSurfaceMeshes.size(); ++iMeshIndex )
        {
            auto const &tSurfaceMesh = aTargetSurfaceMeshes( iMeshIndex ).second;
            for ( size_t iCellIndex = 0; iCellIndex < tSurfaceMesh.get_number_of_cells(); ++iCellIndex )
            {
                ArborX::Box                    tBox;
                moris::Matrix< moris::DDRMat > tVertices = tSurfaceMesh.get_vertex_coordinates_of_cell( iCellIndex );
                for ( size_t iVertexIndex = 0; iVertexIndex < tVertices.n_cols(); ++iVertexIndex )
                {
                    tBox += coordinate_to_arborx_point< ArborX::Point >( tVertices.get_column( iVertexIndex ) );
                }
                tBoxes( tBoxIndex )       = tBox;
                tMeshIndices( tBoxIndex ) = aTargetSurfaceMeshes( iMeshIndex ).first;
                tCellIndices( tBoxIndex ) = iCellIndex;
                ++tBoxIndex;
            }
        }

        return QueryBoxes< MemorySpace >{ tBoxes, tMeshIndices, tCellIndices };
    }

    template< typename MemorySpace, typename ExecutionSpace >
    QueryRays< MemorySpace > construct_query_rays( ExecutionSpace const &aExecutionSpace, moris::mtk::MappingResult const &aMappingResult )
    {
        uint const tNumPoints = aMappingResult.mSourcePhysicalCoordinate.n_cols();
        // since rays are directional, we need to double the number of rays to store rays pointing in both directions (positive and negative)
        uint const                                               tNumRays = 2 * tNumPoints;
        Kokkos::View< ArborX::Experimental::Ray *, MemorySpace > tRays( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:rays" ), tNumRays );
        Kokkos::View< moris_index *, MemorySpace >               tCellIndices( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:cell_indices" ), tNumRays );

        {    // Initialize the rays from the mapping result struct (preinitialized points and normals)
            Kokkos::parallel_for(
                    "initialize_rays",
                    Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumPoints ),
                    KOKKOS_LAMBDA( size_t const iRay ) {
                        size_t const tIndex = 2 * iRay;    // store the two rays of the first point pointing in positive, then negative direction. The second ray will be at index 2*iRay and so on...

                        auto const &tCoord  = aMappingResult.mSourcePhysicalCoordinate.get_column( iRay );
                        auto const &tNormal = aMappingResult.mNormals.get_column( iRay );

                        ArborX::Point const tOrigin = coordinate_to_arborx_point< ArborX::Point >( tCoord );

                        tRays( tIndex )     = ArborX::Experimental::Ray{ tOrigin, coordinate_to_arborx_point< ArborX::Experimental::Vector >( tNormal ) };
                        tRays( tIndex + 1 ) = ArborX::Experimental::Ray{ tOrigin, coordinate_to_arborx_point< ArborX::Experimental::Vector >( -1.0 * tNormal ) };

                        tCellIndices( tIndex )     = aMappingResult.mSourceCellIndex( iRay );    // ray in the positive direction
                        tCellIndices( tIndex + 1 ) = aMappingResult.mSourceCellIndex( iRay );    // ray in the negative direction
                    } );
        }

        return QueryRays< MemorySpace >{ tRays, tCellIndices };
    }

    cell_locator_map
    map_rays_to_boxes( moris::mtk::MappingResult const &aMappingResult, moris::Vector< std::pair< moris_index, moris::mtk::Surface_Mesh > > const &aTargetSurfaceMeshes )
    {
        Tracer tTracer( "Quadrature Point Mapper", "Map", "Perform Raytracing with ArborX" );
        using ExecutionSpace = Kokkos::DefaultExecutionSpace;
        using MemorySpace    = ExecutionSpace::memory_space;
        ExecutionSpace tExecutionSpace{};

        // Construct the query boxes from all cells in the target surface meshes
        QueryBoxes< MemorySpace > tQueryBoxes = construct_query_boxes< MemorySpace >( tExecutionSpace, aTargetSurfaceMeshes );

        // Construct the query rays from all the points that have been initialized in the mapping result (source side points)
        QueryRays< MemorySpace > tQueryRays = construct_query_rays< MemorySpace >( tExecutionSpace, aMappingResult );

        ArborX::BVH< MemorySpace > tBoundingVolumeHierarchy( tExecutionSpace, tQueryBoxes );
        //        ArborX::BruteForce< MemorySpace >                 tBoundingVolumeHierarchy( tExecutionSpace, tQueryBoxes ); // The brute force algorithm for comparison... much slower!
        Kokkos::View< QueryResult *, MemorySpace > tResults( "values", 0 );
        Kokkos::View< int *, MemorySpace >         tOffsets( "offsets", 0 );

        tBoundingVolumeHierarchy.query( tExecutionSpace, tQueryRays, IntersectionCallback< MemorySpace >{ tQueryBoxes, tQueryRays }, tResults, tOffsets );

        // return the results as unordered map
        // the key determines the cell < mesh index, cluster index, cell index >
        // the value is a vector of ray indices (columns in the mapping result) that intersected with the cell
        cell_locator_map tBoxRayMap;
        for ( size_t i = 0; i < tResults.extent( 0 ); ++i )
        {
            moris_index const tBoxIndex   = tResults( i ).mBoxIndex;
            moris_index const tPointIndex = tResults( i ).mPointIndex;    // index of the point from which a hitting ray originates (might be either the positive or negative ray)
            moris_index const tMeshIndex  = tQueryBoxes.mMeshIndices( tBoxIndex );
            moris_index const tCellIndex  = tQueryBoxes.mCellIndices( tBoxIndex );    // TODO @ff: This is trivial at the moment since box and cell indices are the same.
            tBoxRayMap[ tMeshIndex ][ tCellIndex ].push_back( tPointIndex );
        }

        return tBoxRayMap;
    }

    template< typename T >
    T coordinate_to_arborx_point( Matrix< moris::DDRMat > const &aMatrix )
    {
        MORIS_ASSERT( ( aMatrix.n_rows() == 3 || aMatrix.n_rows() == 2 ) && aMatrix.n_cols() == 1, "The input matrix must have 2 or 3 rows and 1 column." );
        float tZCoord = aMatrix.n_rows() == 3 ? aMatrix( 2, 0 ) : 0.0;
        return { static_cast< float >( aMatrix( 0, 0 ) ), static_cast< float >( aMatrix( 1, 0 ) ), tZCoord };
    }
}    // namespace moris::mtk::arborx
