/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper_Ray_ArborX.cpp
 *
 */

#include "cl_Logger.hpp"
#include "cl_Matrix.hpp"
#include "cl_Tracer.hpp"

#include "cl_MTK_QuadraturePointMapper_Ray_ArborX.hpp"
#include "fn_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_MTK_Ray_Line_Intersection.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Surface_Mesh.hpp"
#include "fn_MTK_Gather_Surface_Mesh_Arrays.hpp"
#include "fn_dot.hpp"


namespace moris::mtk
{
    namespace
    {
        // brendan delete all these?
        inline void assert_param_in_bounds_box( Matrix< DDRMat > const &aParam, const char *aContext )
        {
            if ( aParam.n_rows() >= 1 )
            {
                MORIS_ASSERT(
                        aParam( 0 ) >= -1.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                        "MTK parametric coordinate out of bounds (param[0]=%e, context=%s)",
                        aParam( 0 ),
                        aContext );
            }
            if ( aParam.n_rows() >= 2 )
            {
                MORIS_ASSERT(
                        aParam( 1 ) >= -1.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                        "MTK parametric coordinate out of bounds (param[1]=%e, context=%s)",
                        aParam( 1 ),
                        aContext );
            }
        }

        // inline void assert_param_in_bounds_simplex( Matrix< DDRMat > const &aParam, const char *aContext )
        // {
        //     if ( aParam.n_rows() >= 1 )
        //     {
        //         MORIS_ASSERT(
        //                 aParam( 0 ) >= 0.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
        //                 "MTK TRI parametric coordinate out of bounds (eta=%e, context=%s)",
        //                 aParam( 0 ),
        //                 aContext );
        //     }

        //     if ( aParam.n_rows() >= 2 )
        //     {
        //         MORIS_ASSERT(
        //                 aParam( 1 ) >= 0.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
        //                 "MTK TRI parametric coordinate out of bounds (zeta=%e, context=%s)",
        //                 aParam( 1 ),
        //                 aContext );

        //         MORIS_ASSERT(
        //                 aParam( 0 ) + aParam( 1 ) <= 1.0 + 1e-12,
        //                 "MTK TRI parametric coordinate out of bounds (eta+zeta=%e, context=%s)",
        //                 aParam( 0 ) + aParam( 1 ),
        //                 aContext );
        //     }
        // }

        // inline real dot3( Matrix< DDRMat > const &a, Matrix< DDRMat > const &b )
        // {
        //     return a( 0 ) * b( 0 ) + a( 1 ) * b( 1 ) + a( 2 ) * b( 2 );
        // }

        // inline Matrix< DDRMat > cross3( Matrix< DDRMat > const &a, Matrix< DDRMat > const &b )
        // {
        //     Matrix< DDRMat > c( 3, 1 );
        //     c( 0 ) = a( 1 ) * b( 2 ) - a( 2 ) * b( 1 );
        //     c( 1 ) = a( 2 ) * b( 0 ) - a( 0 ) * b( 2 );
        //     c( 2 ) = a( 0 ) * b( 1 ) - a( 1 ) * b( 0 );
        //     return c;
        // }

        // inline bool ray_triangle_intersect(
        //         Matrix< DDRMat > const &aOrigin,
        //         Matrix< DDRMat > const &aDirection,
        //         Matrix< DDRMat > const &aV0,
        //         Matrix< DDRMat > const &aV1,
        //         Matrix< DDRMat > const &aV2,
        //         real                   &aT,
        //         real                   &aU,
        //         real                   &aV )
        // {
        //     const real       tEps = 1e-12;
        //     Matrix< DDRMat > tE1  = aV1 - aV0;
        //     Matrix< DDRMat > tE2  = aV2 - aV0;
        //     Matrix< DDRMat > tP   = cross3( aDirection, tE2 );
        //     real             tDet = dot3( tE1, tP );
        //     if ( std::abs( tDet ) < tEps )
        //     {
        //         return false;
        //     }
        //     real             tInvDet = 1.0 / tDet;
        //     Matrix< DDRMat > tT      = aOrigin - aV0;
        //     aU                       = dot3( tT, tP ) * tInvDet;
        //     if ( aU < 0.0 || aU > 1.0 )
        //     {
        //         return false;
        //     }
        //     Matrix< DDRMat > tQ = cross3( tT, tE1 );
        //     aV                  = dot3( aDirection, tQ ) * tInvDet;
        //     if ( aV < 0.0 || ( aU + aV ) > 1.0 )
        //     {
        //         return false;
        //     }
        //     aT = dot3( tE2, tQ ) * tInvDet;
        //     return true;
        // }
    }    // namespace

    void QuadraturePointMapper_ArborX::update_ip_element_displacements( std::vector< std::tuple< moris_index, Matrix< DDRMat > > > const &aIPElementDisplacements )
    {
        QuadraturePointMapper_Ray::update_ip_element_displacements( aIPElementDisplacements );
        // mGatheredTargetMeshesDirty = true;
    }

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

        // Gather global arrays once per source mesh and refresh only when the surface-mesh state changes
        // (i.e. when displacement updates arrive). The gather uses the current deformed surface-mesh state,
        // so the resulting coordinates already include the IP-element displacements if they were stored.
        // auto tCachedIt = mCachedGatheredTargetMeshes.find( aSourceMeshIndex );
        // if ( mGatheredTargetMeshesDirty || tCachedIt == mCachedGatheredTargetMeshes.end() )
        // {
        //     Vector< arborx::GatheredSurfaceMesh > tGatheredTargetMeshes;
        //     for ( auto const &tPair : get_candidate_pairs() )
        //     {
        //         if ( tPair.first == aSourceMeshIndex )
        //         {
        //             arborx::GatheredSurfaceMesh tG;
        //             tG.mMeshIndex = tPair.second;

        //             mtk::gather_surface_mesh_arrays(
        //                     get_surface_meshes()( tPair.second ),
        //                     tG.mGlobalCells,
        //                     tG.mGlobalCellIds,
        //                     tG.mGlobalCellOwners,
        //                     tG.mGlobalVertexIds,
        //                     tG.mGlobalVertexCoords,
        //                     tG.mGlobalVertexDisplacements );

        //             tGatheredTargetMeshes.push_back( tG );
        //         }
        //     }
        //     mCachedGatheredTargetMeshes[ aSourceMeshIndex ] = tGatheredTargetMeshes;
        //     mGatheredTargetMeshesDirty                      = false;
        // }
        // const Vector< arborx::GatheredSurfaceMesh > &tGatheredTargetMeshes = mCachedGatheredTargetMeshes.at( aSourceMeshIndex );

        // auto const &tBoxRayMap = arborx::map_rays_to_boxes( tMappingResult, tGatheredTargetMeshes );

        const auto &tBoxRayMap = arborx::map_rays_to_boxes( tMappingResult, get_target_surface_meshes( aSourceMeshIndex ) );

        // // check the intersections of the rays with the target cells (use gathered arrays for geometry lookup)
        // std::unordered_map< moris_index, arborx::GatheredSurfaceMesh > tGatheredMap;
        // for ( auto const &tGM : tGatheredTargetMeshes )
        // {
        //     tGatheredMap[ tGM.mMeshIndex ] = tGM;
        // }
        check_cell_intersections( tMappingResult, aMaxNegativeRayLength, aMaxPositiveRayLength, tBoxRayMap );

        return tMappingResult;
    }

    void QuadraturePointMapper_ArborX::check_cell_intersections(
            MappingResult                  &tMappingResult,
            real                            aMaxNegativeRayLength,
            real                            aMaxPositiveRayLength,
            arborx::cell_locator_map const &aBoxRayMap ) const
    {
        // Refine ray/box hits with exact (possibly nonlinear) intersection checks.
        Tracer tTracer( "Quadrature Point Mapper", "Map", "Check Cell Intersections" );
        size_t tNumRaysChecked = 0;
        for ( auto const &[ tTargetMeshIndex, tTargetCells ] : aBoxRayMap )
        {
            Integration_Surface_Mesh const &tTargetMesh = get_surface_meshes()( tTargetMeshIndex );
            for ( auto const &[ tTargetCellIndex, tRayIndices ] : tTargetCells )
            {

                /*
                // ORIGINAL CODE (commented out):
                // moris_index tGlobalTargetCellIndex = tTargetMesh.get_global_cell_index( tTargetCellIndex );
                // // Use side set to get the cluster for this cell
                // const Side_Set     *tSideSet = get_side_sets()( tTargetMeshIndex );
                // const mtk::Cluster *tCluster = nullptr;
                // for ( uint iCluster = 0; iCluster < tSideSet->get_num_clusters_on_set(); ++iCluster )
                // {
                //     const mtk::Cluster *tCandidateCluster = tSideSet->get_clusters_by_index( iCluster );
                //     auto                tPrimaryCells     = tCandidateCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                //     for ( uint i = 0; i < tPrimaryCells.size(); ++i )
                //     {
                //         if ( tPrimaryCells( i )->get_index() == tGlobalTargetCellIndex )
                //         {
                //             tCluster = tCandidateCluster;
                //             break;
                //         }
                //     }
                //     if ( tCluster ) break;
                // }
                */

                moris_index tGlobalTargetCellIndex = -1;
                // if ( aGatheredMeshes != nullptr && aGatheredMeshes->find( tTargetMeshIndex ) != aGatheredMeshes->end() )
                // {
                //     auto const &tG         = aGatheredMeshes->at( tTargetMeshIndex );
                //     tGlobalTargetCellIndex = tG.mGlobalCellIds( tTargetCellIndex );
                // }
                // else
                // {
                tGlobalTargetCellIndex = tTargetMesh.get_global_cell_index( tTargetCellIndex );
                //}

                // Use side set to get the cluster for this cell
                const Side_Set     *tSideSet = get_side_sets()( tTargetMeshIndex );
                const mtk::Cluster *tCluster = nullptr;
                // int                 tClusterIndex = -1;
                for ( uint iCluster = 0; iCluster < tSideSet->get_num_clusters_on_set(); ++iCluster )
                {
                    const mtk::Cluster *tCandidateCluster = tSideSet->get_clusters_by_index( iCluster );
                    auto                tPrimaryCells     = tCandidateCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                    for ( uint iC = 0; iC < tPrimaryCells.size(); ++iC )
                    {
                        if ( tPrimaryCells( iC )->get_index() == tGlobalTargetCellIndex )
                        {
                            tCluster = tCandidateCluster;
                            // tClusterIndex = (int)iCluster;
                            break;
                        }
                    }
                    if ( tCluster ) break;
                }

                // Find the local index of the cell within the cluster
                moris_index aLeaderClusterLocalIndex = -1;
                if ( tCluster )
                {
                    auto tNumPrimaryCells = tCluster->get_num_primary_cells();
                    for ( uint iC = 0; iC < tNumPrimaryCells; ++iC )
                    {
                        Vector< moris::mtk::Cell const * > const &tPrimaryCellsInCluster = tCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                        if ( tPrimaryCellsInCluster( iC )->get_index() == tGlobalTargetCellIndex )
                        {
                            aLeaderClusterLocalIndex = iC;
                            break;
                        }
                    }
                }

                // Now use this index for local coordinates
                moris::Matrix< moris::DDRMat > tTargetLocalCoordinates;
                if ( tCluster )
                {
                    tTargetLocalCoordinates = tCluster->get_cell_local_coords_on_side_wrt_interp_cell( aLeaderClusterLocalIndex );
                }
                //  Get the interpolation cell from the cluster (skip if cluster not found)
                if ( tCluster == nullptr )
                {
                    continue;
                }
                const mtk::Cell &tIPElement      = tCluster->get_interpolation_cell();
                moris_index      tIPElementIndex = tIPElement.get_index();
                // Get the local coordinates of the cell vertices
                Matrix< DDRMat >      tIPElementVertices;
                const mtk::Cell_Info *tIPInfo = tIPElement.get_cell_info();
                tIPInfo->get_loc_coords_of_cell( tIPElementVertices );

                // Get the displacement at the IP element nodes
                Matrix< DDRMat >                                 tIPElementDisplacements;
                std::map< moris_index, Matrix< DDRMat > > const &tIPElementDisplacementsMap = tTargetMesh.get_ip_element_displacement();
                auto                                             it                         = tIPElementDisplacementsMap.find( tIPElementIndex );
                if ( it != tIPElementDisplacementsMap.end() )
                {
                    tIPElementDisplacements = it->second;

                    // sint tNiter = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
                    // std::cout << "MTK: Niter=" << tNiter << " IPelem=" << tIPElementIndex;
                    // PRINT( tIPElementDisplacements );

                    if ( tIPElementDisplacements.n_cols() > tIPElementVertices.n_cols() )
                    {
                        tIPElementDisplacements = tIPElementDisplacements( { 0, tIPElementDisplacements.n_rows() - 1 }, { 0, tIPElementVertices.n_cols() - 1 } );
                    }
                }
                else
                {
                    // If not found, fill with zeros
                    tIPElementDisplacements = Matrix< DDRMat >( tIPElementVertices.n_rows(), tIPElementVertices.n_cols(), 0.0 );
                }

                Matrix< DDRMat > tActualCellCoordinates = tIPElement.get_vertex_coords();

                // Geometry interpolator
                Interpolation_Rule tGeomInterpRule(
                        tIPElement.get_geometry_type(),
                        Interpolation_Type::LAGRANGE,
                        tIPElement.get_cell_info()->get_cell_interpolation_order(),
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );
                Space_Interpolator tGeomSpaceInterpolator( tGeomInterpRule );
                tGeomSpaceInterpolator.set_space_coeff( tActualCellCoordinates );      // Node physical coordinates
                tGeomSpaceInterpolator.set_space_param_coeff( tIPElementVertices );    // Node parametric coordinates

                // Field interpolator
                Interpolation_Rule tFieldInterpRule(
                        tIPElement.get_geometry_type(),
                        Interpolation_Type::LAGRANGE,
                        tIPElement.get_cell_info()->get_cell_interpolation_order(),
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );
                Space_Interpolator tFieldSpaceInterpolator( tFieldInterpRule );
                tFieldSpaceInterpolator.set_space_coeff( tIPElementDisplacements );     // Node physical coordinates
                tFieldSpaceInterpolator.set_space_param_coeff( tIPElementVertices );    // Node parametric coordinates

                Matrix< DDRMat > tTargetCellCoordinates;
                // if ( aGatheredMeshes != nullptr && aGatheredMeshes->find( tTargetMeshIndex ) != aGatheredMeshes->end() )
                // {
                //     auto const &tG        = aGatheredMeshes->at( tTargetMeshIndex );
                //     auto const &tCellConn = tG.mGlobalCells( tTargetCellIndex );
                //     uint const  tNumVerts = tCellConn.n_rows();
                //     tTargetCellCoordinates.set_size( tG.mGlobalVertexCoords.n_rows(), tNumVerts );
                //     for ( uint iv = 0; iv < tNumVerts; ++iv )
                //     {
                //         moris_index tCompIdx                    = tCellConn( iv );
                //         tTargetCellCoordinates.get_column( iv ) = tG.mGlobalVertexCoords.get_column( tCompIdx );
                //     }
                // }
                // else
                // {
                tTargetCellCoordinates = tTargetMesh.get_all_vertex_coordinates_of_facet( tTargetCellIndex );
                // }

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

                for ( size_t i = 0; i < tRayIndices.size(); ++i )
                {
                    tNumRaysChecked++;
                    auto             tRayIndex           = tRayIndices( i );
                    Matrix< DDRMat > tRayOrigin          = tMappingResult.mSourcePhysicalCoordinateNonlinear.get_column( tRayIndex );
                    Matrix< DDRMat > tRayDirectionSource = tMappingResult.mNormalsNonlinear.get_column( tRayIndex );

                    tRayLineIntersection.perform_nonlinear_raytracing(
                            tGeomSpaceInterpolator,
                            tFieldSpaceInterpolator,
                            tRayOrigin,
                            tRayDirectionSource,
                            tTargetLocalCoordinates );

                    if ( tRayLineIntersection.has_intersection()
                            && tRayLineIntersection.get_signed_ray_length() > aMaxNegativeRayLength
                            && tRayLineIntersection.get_signed_ray_length() < aMaxPositiveRayLength
                            && ( tRayLineIntersection.get_signed_ray_length() < tMappingResult.mSignedDistance( tRayIndex ) || tMappingResult.mTargetCellIndices( tRayIndex ) == -1 ) )
                    {
                        Matrix< DDRMat > tTargetParam = tRayLineIntersection.get_intersection_parametric();
                        assert_param_in_bounds_box( tTargetParam, "QuadraturePointMapper_ArborX::nonlinear_raytrace" );
                        tMappingResult.mTargetParametricCoordinate.set_column( tRayIndex, tTargetParam );

                        // Get the parametric coordinate of the Gauss point
                        Matrix< DDRMat > tTargetEtagp = tMappingResult.mTargetParametricCoordinate.get_column( tRayIndex );
                        assert_param_in_bounds_box( tTargetEtagp, "QuadraturePointMapper_ArborX::before_set_space" );

                        Interpolation_Rule tFacetInterpRule(
                                tSideSet->get_integration_cell_geometry_type(),
                                Interpolation_Type::LAGRANGE,
                                Interpolation_Order::LINEAR,
                                Interpolation_Type::UNDEFINED,
                                Interpolation_Order::UNDEFINED );

                        Space_Interpolator tFacetInterpolator( tFacetInterpRule );

                        // tTargetLocalCoordinates has rows = nodes, cols = parametric dims; transpose to (dim x nNodes)
                        Matrix< DDRMat > tTargetLocalCoordsTrans = trans( tTargetLocalCoordinates );
                        tFacetInterpolator.set_space_param_coeff( tTargetLocalCoordsTrans );
                        tFacetInterpolator.set_space( tTargetEtagp );

                        // Evaluate shape functions and interpolate the parametric coordinate along the side
                        Matrix< DDRMat > tTargetNrs  = trans( tFacetInterpolator.NXi() );
                        Matrix< DDRMat > tTargetRsgp = tTargetLocalCoordsTrans * tTargetNrs;
                        tGeomSpaceInterpolator.set_space( tTargetRsgp );
                        tFieldSpaceInterpolator.set_space( tTargetRsgp );
                        Matrix< DDRMat > tTargetYgp = trans( tGeomSpaceInterpolator.valx() + tFieldSpaceInterpolator.valx() );

                        tMappingResult.mTargetPhysicalCoordinate.set_column( tRayIndex, tTargetYgp );
                        tMappingResult.mSignedDistance( tRayIndex )       = tRayLineIntersection.get_signed_ray_length();
                        tMappingResult.mTargetSideSetIndices( tRayIndex ) = tTargetMeshIndex;
                        // if ( aGatheredMeshes != nullptr && aGatheredMeshes->find( tTargetMeshIndex ) != aGatheredMeshes->end() )
                        // {
                        //     tMappingResult.mTargetCellIndices( tRayIndex )  = tGlobalTargetCellIndex;
                        //     tMappingResult.mTargetClusterIndex( tRayIndex ) = tClusterIndex;
                        // }
                        // else
                        // {
                        tMappingResult.mTargetCellIndices( tRayIndex )  = tTargetMesh.get_global_cell_index( tTargetCellIndex );
                        tMappingResult.mTargetClusterIndex( tRayIndex ) = tTargetMesh.get_cluster_of_cell( tTargetCellIndex );
                        // }
                    }
                }
            }
        }
    }

    Vector< std::pair< moris_index, Integration_Surface_Mesh > > QuadraturePointMapper_ArborX::get_target_surface_meshes( moris_index aSourceMeshIndex ) const
    {
        Vector< std::pair< moris_index, Integration_Surface_Mesh > > tTargetSurfaceMeshes;
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
