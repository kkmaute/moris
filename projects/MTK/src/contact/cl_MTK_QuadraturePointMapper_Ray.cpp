/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper_Ray.cpp
 *
 */
#include <utility>
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Json_Object.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_QuadraturePointMapper.hpp"
#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_Tracer.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "fn_dot.hpp"
#include <iostream>

namespace
{
    inline void assert_param_in_bounds_box( moris::Matrix< moris::DDRMat > const &aParam, const char *aContext )
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

    inline void assert_param_in_bounds_simplex( moris::Matrix< moris::DDRMat > const &aParam, const char *aContext )
    {
        if ( aParam.n_rows() >= 1 )
        {
            MORIS_ASSERT(
                    aParam( 0 ) >= 0.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                    "MTK TRI parametric coordinate out of bounds (eta=%e, context=%s)",
                    aParam( 0 ),
                    aContext );
        }
        if ( aParam.n_rows() >= 2 )
        {
            MORIS_ASSERT(
                    aParam( 1 ) >= 0.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK TRI parametric coordinate out of bounds (zeta=%e, context=%s)",
                    aParam( 1 ),
                    aContext );
            MORIS_ASSERT(
                    aParam( 0 ) + aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK TRI parametric coordinate out of bounds (eta+zeta=%e, context=%s)",
                    aParam( 0 ) + aParam( 1 ),
                    aContext );
        }
    }
}    // namespace

namespace moris::mtk
{
    QuadraturePointMapper_Ray::QuadraturePointMapper_Ray(
            Integration_Mesh                                      *aIGMesh,
            Vector< Side_Set const * >                            &aSideSets,
            Vector< std::pair< moris_index, moris_index > > const &aCandidatePairs )
            : QuadraturePointMapper( aIGMesh, aSideSets, aCandidatePairs )
            , mSurfaceMeshes( initialize_surface_meshes( aIGMesh, aSideSets ) )
            , mReferenceSurfaceMeshes( mSurfaceMeshes )    // copy the surface meshes to the reference meshes
    {
    }

    void QuadraturePointMapper_Ray::write_surface_mesh_json() const
    {
        Json  tSurfaceMeshes;
        auto &tMeshes = tSurfaceMeshes.put_child( "surface_meshes", Json() );
        for ( auto const &tSurfaceMesh : mSurfaceMeshes )
        {
            tMeshes.push_back( { "", tSurfaceMesh.to_json() } );
        }
        uint const  tIteration = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve" );
        std::string tFileName  = "surface_meshes_" + std::to_string( tIteration ) + ".json";
        write_json( tFileName, tSurfaceMeshes );
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
            tSurfaceMesh.write_to_file( "surface_mesh_" + tSideSet->get_set_name() + std::to_string( gLogger.iter ) + ".obj" );
        }

        gLogger.iter++;
        return tSurfaceMeshes;
    }

    MappingResult QuadraturePointMapper_Ray::initialize_source_points( moris_index aSourceMeshIndex, Matrix< DDRMat > const &aParametricCoordinates ) const
    {
        Tracer tTracer( "Quadrature Point Mapper", "Map", "Initialize Source Points" );
        MORIS_ASSERT( aSourceMeshIndex < static_cast< moris_index >( get_surface_meshes().size() ), "QuadraturePointMapper_Ray::initialize_source_points: Source mesh index %d out of range.", aSourceMeshIndex );
        Surface_Mesh const &tSurfaceMesh          = get_surface_meshes()( aSourceMeshIndex );
        Surface_Mesh const &tReferenceSurfaceMesh = get_reference_surface_meshes()( aSourceMeshIndex );
        Side_Set const     *tSideSet              = get_side_sets()( aSourceMeshIndex );

        Interpolation_Rule const tInterpolationRule(
                tSideSet->get_integration_cell_geometry_type(),
                Interpolation_Type::LAGRANGE,
                Interpolation_Order::LINEAR,
                Interpolation_Type::UNDEFINED,
                Interpolation_Order::UNDEFINED );

        Space_Interpolator tInterpolator( tInterpolationRule );
        uint const         tDim            = tSideSet->get_spatial_dim();
        uint const         tNumCells       = tSurfaceMesh.get_number_of_cells();
        uint const         tNumRaysPerCell = aParametricCoordinates.n_cols();
        uint const         tTotalNumPoints = tNumCells * tNumRaysPerCell;

        MappingResult tMappingResult( aSourceMeshIndex, tDim, tTotalNumPoints );

        for ( moris_index iCell = 0; iCell < static_cast< moris_index >( tNumCells ); iCell++ )
        {    // Kurt should be const & and no vertex normals
            Matrix< DDRMat > const tVertexCoordinates      = tSurfaceMesh.get_vertex_coordinates_of_cell( iCell );
            Matrix< DDRMat > const tVertexNormals          = tSurfaceMesh.get_vertex_normals_of_cell( iCell );
            Matrix< DDRMat > const tReferenceVertexNormals = tReferenceSurfaceMesh.get_vertex_normals_of_cell( iCell );
            Matrix< DDRMat > const tNormals                = tSurfaceMesh.get_facet_normals();
            Matrix< DDRMat > const tReferenceNormals       = tReferenceSurfaceMesh.get_facet_normals();
            moris_index            tGlobalSourceCellIndex  = tSurfaceMesh.get_global_cell_index( iCell );
            const Side_Set        *tSideSet                = get_side_sets()( aSourceMeshIndex );
            const mtk::Cluster    *tCluster                = nullptr;

            // Find the cluster in which the examined cell (tGlobalSourceCellIndex) belongs
            for ( uint iCluster = 0; iCluster < tSideSet->get_num_clusters_on_set(); ++iCluster )
            {
                const mtk::Cluster *tCandidateCluster = tSideSet->get_clusters_by_index( iCluster );
                auto                tPrimaryCells     = tCandidateCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                for ( uint i = 0; i < tPrimaryCells.size(); ++i )
                {
                    if ( tPrimaryCells( i )->get_index() == tGlobalSourceCellIndex )
                    {
                        tCluster = tCandidateCluster;
                        break;
                    }
                }
                if ( tCluster ) break;
            }

            // Get the local index of the cell within the cluster
            moris_index tSourceCellLocalIndex = -1;
            auto        numPrimaryCells       = tCluster->get_num_primary_cells();
            for ( uint i = 0; i < numPrimaryCells; ++i )
            {
                Vector< moris::mtk::Cell const * > const &tPrimaryCellsInCluster = tCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                if ( tPrimaryCellsInCluster( i )->get_index() == tGlobalSourceCellIndex )
                {
                    tSourceCellLocalIndex = i;
                    break;
                }
            }

            Matrix< DDRMat > tSourceLocalCoordinates = tCluster->get_cell_local_coords_on_side_wrt_interp_cell( tSourceCellLocalIndex );

            // Get the physical coordinates of the cell vertices (x,y,z)
            const mtk::Cell      &tIPElement      = tCluster->get_interpolation_cell();
            moris_index           tIPElementIndex = tIPElement.get_index();
            Matrix< DDRMat >      tIPElementVertices;
            const mtk::Cell_Info *tIPInfo = tIPElement.get_cell_info();
            tIPInfo->get_loc_coords_of_cell( tIPElementVertices );
            Matrix< DDRMat > tActualCellCoordinates = tIPElement.get_vertex_coords();

            // Space Interpolator for the geometry of the source cell
            Interpolation_Rule tGeomInterpRule(
                    tIPElement.get_geometry_type(),
                    Interpolation_Type::LAGRANGE,
                    tIPElement.get_cell_info()->get_cell_interpolation_order(),
                    Interpolation_Type::UNDEFINED,
                    Interpolation_Order::UNDEFINED );
            Space_Interpolator tGeomSpaceInterpolator( tGeomInterpRule );
            tGeomSpaceInterpolator.set_space_coeff( tActualCellCoordinates );
            tGeomSpaceInterpolator.set_space_param_coeff( tIPElementVertices );

            Matrix< DDRMat > tIPElementDisplacements;
            auto             it = tSurfaceMesh.mIPElementDisplacements.find( tIPElementIndex );

            if ( it != tSurfaceMesh.mIPElementDisplacements.end() )
            {
                tIPElementDisplacements = it->second;
            }
            else
            {
                tIPElementDisplacements = Matrix< DDRMat >( tIPElementVertices.n_rows(), tIPElementVertices.n_cols(), 0.0 );
            }

            // Space Interpolator for the displacements of the source cell
            Interpolation_Rule tFieldInterpRule(
                    tIPElement.get_geometry_type(),
                    Interpolation_Type::LAGRANGE,
                    tIPElement.get_cell_info()->get_cell_interpolation_order(),
                    Interpolation_Type::UNDEFINED,
                    Interpolation_Order::UNDEFINED );
            Space_Interpolator tFieldSpaceInterpolator( tFieldInterpRule );
            tFieldSpaceInterpolator.set_space_coeff( tIPElementDisplacements );     // Node physical coordinates
            tFieldSpaceInterpolator.set_space_param_coeff( tIPElementVertices );    // Node parametric coordinates


            moris_index const tStartIndex = iCell * tNumRaysPerCell;
            for ( uint iPoint = 0; iPoint < tNumRaysPerCell; iPoint++ )
            {
                moris_index const tRayIndex = tStartIndex + iPoint;

                // Set the parametric coordinate of the interpolator and get the values of the shape functions
                tInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
                Matrix< DDRMat > const tNXi = trans( tInterpolator.NXi() );

                Matrix< DDRMat > tSourceEtagp = aParametricCoordinates.get_column( iPoint );

                // Build a Lagrange interpolator along the Source local coordinates on the side
                Interpolation_Rule tFacetInterpRule(
                        tSideSet->get_integration_cell_geometry_type(),
                        Interpolation_Type::LAGRANGE,
                        Interpolation_Order::LINEAR,
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );

                Space_Interpolator tFacetInterpolator( tFacetInterpRule );

                // tSourceLocalCoordinates has rows = nodes, cols = parametric dims; transpose to (dim x nNodes)
                Matrix< DDRMat > tSourceLocalCoordsTrans = trans( tSourceLocalCoordinates );
                tFacetInterpolator.set_space_param_coeff( tSourceLocalCoordsTrans );
                if ( tSideSet->get_integration_cell_geometry_type() == Geometry_Type::TRI )
                {
                    assert_param_in_bounds_simplex( tSourceEtagp, "QuadraturePointMapper_Ray::source_eta" );
                }
                else
                {
                    assert_param_in_bounds_box( tSourceEtagp, "QuadraturePointMapper_Ray::source_eta" );
                }
                tFacetInterpolator.set_space( tSourceEtagp );

                // Evaluate shape functions and interpolate the parametric coordinate along the side
                Matrix< DDRMat > tSourceNRsgp = trans( tFacetInterpolator.NXi() );
                Matrix< DDRMat > tSourceRsgp  = tSourceLocalCoordsTrans * tSourceNRsgp;
                // tSourceRsgp is in interpolation-cell parametric coordinates (r,s[,t]),
                // so box bounds are appropriate even when side geometry is TRI.
                assert_param_in_bounds_box( tSourceRsgp, "QuadraturePointMapper_Ray::source_rs_interp" );
                tGeomSpaceInterpolator.set_space_time( tSourceRsgp );
                tFieldSpaceInterpolator.set_space( tSourceRsgp );
                Matrix< DDRMat > tSourceYgp = tGeomSpaceInterpolator.valx() + tFieldSpaceInterpolator.valx();

                // Compute chain rule factors in matrix form
                Matrix< DDRMat > tSourcedNrsdXi = trans( tFacetInterpolator.dNdXi() );
                Matrix< DDRMat > tSourcedrsdXi  = tSourceLocalCoordsTrans * tSourcedNrsdXi;

                // Derivatives with respect to r and s
                Matrix< DDRMat > tSourcedNXgpdrs = tGeomSpaceInterpolator.dNdXi();
                Matrix< DDRMat > tSourceIGNodes  = tGeomSpaceInterpolator.get_space_coeff();

                // dX/dr and dX/ds as column vectors (nDim x 1)
                Matrix< DDRMat > tSourcedXgpdrs = tSourcedNXgpdrs * tSourceIGNodes;
                Matrix< DDRMat > tSourcedXgpdXi = tSourcedXgpdrs * tSourcedrsdXi;

                // // Displacement derivatives
                Matrix< DDRMat > tSourcedNUgpdr = tFieldSpaceInterpolator.dNdXi();
                Matrix< DDRMat > tSourceUHat    = tFieldSpaceInterpolator.get_space_coeff();
                Matrix< DDRMat > tSourcedUgpdrs = tSourcedNUgpdr * tSourceUHat;
                Matrix< DDRMat > tSourcedUgpdXi = trans( trans( tSourcedrsdXi ) * tSourcedUgpdrs );

                // compute derivatives along IG element side
                // const Matrix< DDRMat > tSourcedXgpdXi = trans( tSourcedNXgpdXi * tSourceLocalCoordinates );
                // const Matrix< DDRMat > tSourcedRgpdXi = tSourcedNXgpdXi * tSourceIGNodesParam;
                // const Matrix< DDRMat > tSourcedUgpdr  = tSourcedNUgpdr * tSourceUHat;
                // const Matrix< DDRMat > tSourcedUgpdXi = trans( tSourcedRgpdXi * tSourcedUgpdr );

                // Compute normal based on dimension
                Matrix< DDRMat > tSourceNormalTilde;
                Matrix< DDRMat > aSourceRefNormal;
                uint tSpaceDim = tSourcedXgpdXi.n_rows();

                if ( tSpaceDim == 2 )
                {
                    // For 2D, the normal is the rotated tangent
                    const Matrix< DDRMat > RotMat = { { 0, 1 }, { -1, 0 } };
                    tSourceNormalTilde            = RotMat * ( tSourcedXgpdXi + tSourcedUgpdXi );
                    aSourceRefNormal              = RotMat * tSourcedXgpdXi;
                    aSourceRefNormal              = 1.0 / norm( aSourceRefNormal ) * aSourceRefNormal;
                }
                else    // 3D
                {
                    // For 3D, compute normal as cross product
                    Matrix< DDRMat > tDxDxi  = tSourcedXgpdXi.get_column( 0 );
                    Matrix< DDRMat > tDxDeta = ( tSourcedXgpdXi.n_cols() > 1 ) ? tSourcedXgpdXi.get_column( 1 ) : Matrix< DDRMat >( tSpaceDim, 1, 0.0 );
                    aSourceRefNormal         = { { tDxDxi( 1 ) * tDxDeta( 2 ) - tDxDxi( 2 ) * tDxDeta( 1 ) },
                        { tDxDxi( 2 ) * tDxDeta( 0 ) - tDxDxi( 0 ) * tDxDeta( 2 ) },
                        { tDxDxi( 0 ) * tDxDeta( 1 ) - tDxDxi( 1 ) * tDxDeta( 0 ) } };
                    aSourceRefNormal         = 1.0 / norm( aSourceRefNormal ) * aSourceRefNormal;

                    Matrix< DDRMat > tDxDxiDef  = tDxDxi + tSourcedUgpdXi.get_column( 0 );
                    Matrix< DDRMat > tDxDetaDef = ( tSourcedUgpdXi.n_cols() > 1 ) ? tDxDeta + tSourcedUgpdXi.get_column( 1 ) : tDxDeta;
                    tSourceNormalTilde          = { { tDxDxiDef( 1 ) * tDxDetaDef( 2 ) - tDxDxiDef( 2 ) * tDxDetaDef( 1 ) },
                        { tDxDxiDef( 2 ) * tDxDetaDef( 0 ) - tDxDxiDef( 0 ) * tDxDetaDef( 2 ) },
                        { tDxDxiDef( 0 ) * tDxDetaDef( 1 ) - tDxDxiDef( 1 ) * tDxDetaDef( 0 ) } };
                }

                const real       tNtildeSquared   = dot( tSourceNormalTilde, tSourceNormalTilde );
                const real       tNtildeOm12      = 1.0 / std::sqrt( tNtildeSquared );
                Matrix< DDRMat > tNormalNonlinear = tNtildeOm12 * tSourceNormalTilde;

                // DEBUG OUTPUT
                // sint              tNiter    = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
                // if ( x0( 0 ) == 4.117391304347825e-01 && tNiter == 36 )
                // {
                //     PRINT( tSourceUHat );
                //     PRINT( tSourcedXgpdXi );
                //     PRINT( tSourcedNUgpdr );
                //     PRINT( tSourcedUgpdXi );
                //     PRINT( tSourceNormalTilde );
                //     PRINT( aSourceRefNormal );
                //     std::cerr << tNtildeSquared << std::endl;
                //     std::cerr << tNtildeOm12 << std::endl;
                //     PRINT( ray_direction );
                //     PRINT( tFieldSpaceInterpolator.dNdXi() );
                // }

                tMappingResult.mSourceCellIndex( tRayIndex )    = tSurfaceMesh.get_global_cell_index( iCell );
                tMappingResult.mSourceClusterIndex( tRayIndex ) = tSurfaceMesh.get_cluster_of_cell( iCell );

                // get the interpolated coordinates and normals of the parametric point
                tMappingResult.mSourcePhysicalCoordinate.set_column( tRayIndex, tVertexCoordinates * tNXi );
                tMappingResult.mNormals.set_column( tStartIndex + iPoint, tNormals.get_column( iCell ) );    // use tRayIndex
                tMappingResult.mReferenceNormals.set_column( tStartIndex + iPoint, tReferenceNormals.get_column( iCell ) );
                tMappingResult.mNormalsNonlinear.set_column( tRayIndex, tNormalNonlinear );
                tMappingResult.mSourcePhysicalCoordinateNonlinear.set_column( tRayIndex, trans( tSourceYgp ) );
            }
        }
        return tMappingResult;
    }

    void QuadraturePointMapper_Ray::update_ip_element_displacements( std::vector< std::tuple< moris_index, Matrix< DDRMat > > > const &aIPElementDisplacements )
    {
        // Cache IP element displacement matrices on each surface mesh.
        Tracer tTracer( "Quadrature Point Mapper", "Update Facet Displacements", "Update" );

        // std::cout << "=== QuadraturePointMapper_Ray::update_facet_displacements ===" << std::endl;
        // std::cout << "Received " << aIPElementDisplacements.size() << " facets with param coords and displacements" << std::endl;
        // std::cout << "Number of surface meshes: " << mSurfaceMeshes.size() << std::endl;

        for ( auto &tSurfaceMesh : mSurfaceMeshes )
        {
            for ( const auto &[ iCellIndex, iDispMat ] : aIPElementDisplacements )
            {
                // std::cout << "  Calling set_facet_displacement for cell=" << iCellIndex
                //           << ", side=" << iSideOrdinal << std::endl;
                // PRINT(iDispMat);
                // Store the displacement matrix (we ignore param coords as they're regenerated in get_all_facet_data_with_param_coords)
                tSurfaceMesh.set_ip_element_displacement( iCellIndex, iDispMat );
            }

            tSurfaceMesh.refresh_derived_quantities();
        }
    }
}    // namespace moris::mtk
