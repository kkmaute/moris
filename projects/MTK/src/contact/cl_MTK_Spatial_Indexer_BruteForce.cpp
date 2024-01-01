//
// Created by frank on 11/13/23.
//

#include "cl_MTK_Spatial_Indexer_BruteForce.h"

namespace moris::mtk
{
    Spatial_Indexing_Result Spatial_Indexer_BruteForce::perform( moris_index const aSourceMeshIndex, real epsilon ) const
    {
        Spatial_Indexing_Result tResult;
        for ( auto& tPair : mCandidatePairs )
        {
            if ( tPair.first == aSourceMeshIndex )
            {
                Spatial_Indexing_Result tNewResult = perform_on_mesh_pair( tPair.first, tPair.second, epsilon );
                tResult.merge( tNewResult );
            }
        }
        return tResult;
    }

    Spatial_Indexing_Result Spatial_Indexer_BruteForce::perform_on_mesh_pair( moris_index aSourceMeshIndex, moris_index aTargetMeshIndex, real aEpsilon ) const
    {
        Spatial_Indexing_Result tResult;
        mtk::Surface_Mesh       tSourceMesh               = mSurfaceMeshes( aSourceMeshIndex );
        mtk::Surface_Mesh       tTargetMesh               = mSurfaceMeshes( aTargetMeshIndex );
        bool                    tIsSelfIntersectionSearch = aSourceMeshIndex == aTargetMeshIndex;
        Matrix< DDRMat >        tSourceCoordinates        = tSourceMesh.get_vertex_coordinates();
        Matrix< DDRMat >        tTargetCoordinates        = tTargetMesh.get_vertex_coordinates();
        Matrix< DDRMat >        tSourceNormals            = tSourceMesh.get_vertex_normals();
        Matrix< DDRMat >        tTargetNormals            = tTargetMesh.get_vertex_normals();

        // if the meshes are the same, we can use the neighbor information to exclude neighboring vertices from the search
        Vector< Vector< moris_index > > mNeighbors;
        if ( tIsSelfIntersectionSearch ) mNeighbors = tSourceMesh.get_vertex_neighbors();

        moris_index tNumSourceVertices = tSourceCoordinates.n_cols();
        moris_index tNumTargetVertices = tTargetCoordinates.n_cols();

        // outer loop over all source coordinates
        for ( moris_index source_vertex = 0; source_vertex < tNumSourceVertices; source_vertex++ )
        {
            moris_index tClosestVertex = MORIS_INDEX_MAX;
            real        tMinDistance   = MORIS_REAL_MAX;

            // inner loop to compare all other coordinates to the first one
            for ( moris_index target_vertex = 0; target_vertex < tNumTargetVertices; target_vertex++ )
            {

                if ( tIsSelfIntersectionSearch )
                {
                    // skip the same vertex if both meshes are the same (i.e. self-intersection)
                    if ( source_vertex == target_vertex ) continue;

                    // check if the second vertex is a neighbor of the first one
                    bool tIsNeighbor = false;
                    for ( auto& tNeighbor : mNeighbors( source_vertex ) )
                    {
                        if ( tNeighbor == target_vertex )
                        {
                            tIsNeighbor = true;
                            break;
                        }
                    }

                    // if the second vertex is a neighbor, skip it
                    if ( tIsNeighbor ) continue;
                }

                // check if the distance between the two vertices is smaller than the current distance
                auto tDistance = norm( tSourceCoordinates.get_column( source_vertex ) - tTargetCoordinates.get_column( target_vertex ) );
                if ( tDistance < tMinDistance )
                {
                    // check that vertex normals are pointing in different directions
                    auto tDot = dot( tSourceNormals.get_column( source_vertex ), tTargetNormals.get_column( target_vertex ) );
                    if ( tDot > 0.0 ) continue;
                    tMinDistance   = tDistance;
                    tClosestVertex = target_vertex;
                }
            }    // end inner loop over all coordinates

            if ( tMinDistance < aEpsilon )
            {
                tResult[ source_vertex ] = { tClosestVertex, aTargetMeshIndex, tMinDistance };
            }
        }    // end outer loop over all coordinates
        return tResult;
    }
}    // namespace moris::mtk
