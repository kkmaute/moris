//
// Created by frank on 3/17/24.
//

#include <utility>
#include <set>
#include "cl_MTK_Interpolation_Rule.hpp"
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
        //        write_surface_mesh_json(); // TODO @ff: Remove! Only for Debug!
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
        }
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
        {
            Matrix< DDRMat > const tVertexCoordinates      = tSurfaceMesh.get_vertex_coordinates_of_cell( iCell );
            Matrix< DDRMat > const tVertexNormals          = tSurfaceMesh.get_vertex_normals_of_cell( iCell );
            Matrix< DDRMat > const tReferenceVertexNormals = tReferenceSurfaceMesh.get_vertex_normals_of_cell( iCell );
            Matrix< DDRMat > const tNormals                = tSurfaceMesh.get_facet_normals();
            Matrix< DDRMat > const tReferenceNormals       = tReferenceSurfaceMesh.get_facet_normals();

            moris_index const tStartIndex = iCell * tNumRaysPerCell;
            for ( uint iPoint = 0; iPoint < tNumRaysPerCell; iPoint++ )
            {
                moris_index const tRayIndex = tStartIndex + iPoint;

                // set the parametric coordinate of the interpolator and get the values of the shape functions
                tInterpolator.set_space( aParametricCoordinates.get_column( iPoint ) );
                Matrix< DDRMat > const tNxi = trans( tInterpolator.NXi() );

                tMappingResult.mSourceCellIndex( tRayIndex )    = tSurfaceMesh.get_global_cell_index( iCell );
                tMappingResult.mSourceClusterIndex( tRayIndex ) = tSurfaceMesh.get_cluster_of_cell( iCell );

                // get the interpolated coordinates and normals of the parametric point
                tMappingResult.mSourcePhysicalCoordinate.set_column( tRayIndex, tVertexCoordinates * tNxi );
                tMappingResult.mNormals.set_column( tStartIndex + iPoint, tNormals.get_column( iCell ) );
                tMappingResult.mReferenceNormals.set_column( tStartIndex + iPoint, tReferenceNormals.get_column( iCell ) );
            }
        }
        return tMappingResult;
    }

    void QuadraturePointMapper_Ray::update_displacements( std::unordered_map< moris_index, Vector< real > > const &aSetDisplacements )
    {
        Tracer tTracer( "Quadrature Point Mapper", "Update Displacements", "Update" );
        for ( auto &tSurfaceMesh : mSurfaceMeshes )
        {
            uint const tNumSurfaceVertices = tSurfaceMesh.get_number_of_vertices();
            uint const tNumDimensions      = tSurfaceMesh.get_spatial_dimension();
            // rows are dimensions, columns are vertices
            Matrix< DDRMat > tDisplacements( tNumDimensions, tNumSurfaceVertices, 0.0 );
            for ( size_t iLocalVertexIndex = 0; iLocalVertexIndex < tNumSurfaceVertices; ++iLocalVertexIndex )
            {
                moris_index const tGlobalVertexIndex = tSurfaceMesh.get_global_vertex_index( iLocalVertexIndex );
                auto const       &tDisplacement      = aSetDisplacements.find( tGlobalVertexIndex );
                MORIS_ASSERT( tDisplacement != aSetDisplacements.end(), "QuadraturePointMapper_Ray::update_displacements: Vertex %d not found in displacement map.", tGlobalVertexIndex );
                MORIS_ASSERT( tDisplacement->second.size() == tNumDimensions, "QuadraturePointMapper_Ray::update_displacements: Vertex %d: Displacement vector has wrong size.", tGlobalVertexIndex );
                for ( size_t iDimension = 0; iDimension < tNumDimensions; ++iDimension )
                {
                    tDisplacements( iDimension, iLocalVertexIndex ) = tDisplacement->second( iDimension );
                }
            }
            tSurfaceMesh.set_displacement( tDisplacements );
        }
        //        write_surface_mesh_json(); // TODO @ff Remove! Only for Debug!
    }
}    // namespace moris::mtk