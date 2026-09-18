/*
 * fn_MTK_Gather_Surface_Mesh_Arrays.cpp
 */

#include "fn_MTK_Gather_Surface_Mesh_Arrays.hpp"
#include "cl_Communication_Tools.hpp"

#include <unordered_map>

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    void
    gather_surface_mesh_arrays(
            Integration_Surface_Mesh const & aLocalSurfaceMesh,
            Vector< Matrix< IndexMat > >&    aGlobalCells,
            Matrix< IndexMat >&              aGlobalCellIds,
            Matrix< IndexMat >&              aGlobalCellOwners,
            Matrix< IndexMat >&              aGlobalVertexIds,
            Matrix< DDRMat >&                aGlobalVertexCoords,
            Matrix< DDRMat >&                aGlobalVertexDisplacements,
            moris_index                      aRootProc )
    {

        uint tNumLocalVerts = aLocalSurfaceMesh.get_number_of_vertices();
        uint tNumLocalCells = aLocalSurfaceMesh.get_number_of_facets();
        uint tDim           = aLocalSurfaceMesh.get_spatial_dimension();

        // local vertex coordinates: (d x nLocalVerts), already in this format
        Matrix< DDRMat > tLocalVertexCoords        = aLocalSurfaceMesh.get_all_vertex_coordinates();
        Matrix< DDRMat > tLocalVertexDisplacements = aLocalSurfaceMesh.get_vertex_displacements();
        if ( tLocalVertexDisplacements.n_rows() != tLocalVertexCoords.n_rows() || tLocalVertexDisplacements.n_cols() != tLocalVertexCoords.n_cols() )
        {
            tLocalVertexDisplacements.set_size( tLocalVertexCoords.n_rows(), tLocalVertexCoords.n_cols() );
            tLocalVertexDisplacements.fill( 0.0 );
        }

        // local -> mesh-wide global vertex id, one row per local vertex
        Matrix< DDSMat > tLocalVertexGlobalIds( tNumLocalVerts, 1 );
        for ( uint iV = 0; iV < tNumLocalVerts; ++iV )
        {
            tLocalVertexGlobalIds( iV ) = aLocalSurfaceMesh.get_global_vertex_index( (moris_index)iV );
        }

        // flattened cell -> global-vertex-id connectivity, plus a parallel
        // array recording how many vertices each cell
        Matrix< DDSMat > tLocalCellVertexCounts( tNumLocalCells, 1 );
        Matrix< DDSMat > tLocalCellGlobalIds( tNumLocalCells, 1 );

        Vector< moris_index > tFlatConnectivity;
        tFlatConnectivity.reserve( tNumLocalCells * 4 );

        for ( moris_index iC = 0; iC < (moris_index)tNumLocalCells; ++iC )
        {
            Vector< moris_index > tCellVerts = aLocalSurfaceMesh.get_facets_vertex_indices( iC );

            tLocalCellVertexCounts( iC ) = (sint)tCellVerts.size();
            tLocalCellGlobalIds( iC )    = aLocalSurfaceMesh.get_global_cell_index( iC );

            for ( moris_index tLocalVertIndex : tCellVerts )
            {
                tFlatConnectivity.push_back(
                        aLocalSurfaceMesh.get_global_vertex_index( tLocalVertIndex ) );
            }
        }

        // Convert to matrix for gatherv_mats
        Matrix< DDSMat > tLocalCellVertexGlobalIds( tFlatConnectivity.size(), 1 );
        std::transform( tFlatConnectivity.begin(), tFlatConnectivity.end(), tLocalCellVertexGlobalIds.begin(), []( moris_index a ) { return (sint)a; } );

        // ---------------------------------------------------------------
        // 2) gather the four local arrays onto the root processor
        // ---------------------------------------------------------------

        Vector< Matrix< DDRMat > > tGatheredVertexCoords;
        Vector< Matrix< DDSMat > > tGatheredVertexGlobalIds;
        Vector< Matrix< DDSMat > > tGatheredCellVertexCounts;
        Vector< Matrix< DDSMat > > tGatheredCellVertexGlobalIds;
        Vector< Matrix< DDSMat > > tGatheredCellGlobalIds;

        moris::gatherv_mats( tLocalVertexCoords, tGatheredVertexCoords, aRootProc );
        Vector< Matrix< DDRMat > > tGatheredVertexDisplacements;
        moris::gatherv_mats( tLocalVertexDisplacements, tGatheredVertexDisplacements, aRootProc );
        moris::gatherv_mats( tLocalVertexGlobalIds, tGatheredVertexGlobalIds, aRootProc );
        moris::gatherv_mats( tLocalCellVertexCounts, tGatheredCellVertexCounts, aRootProc );
        moris::gatherv_mats( tLocalCellVertexGlobalIds, tGatheredCellVertexGlobalIds, aRootProc );
        moris::gatherv_mats( tLocalCellGlobalIds, tGatheredCellGlobalIds, aRootProc );


        // flattened form of aGlobalCells (built on root,
        // reconstructed on every rank after the broadcast in step 4)
        Matrix< DDSMat > tGlobalCellVertexCounts;
        Matrix< DDSMat > tGlobalCellVertexIndices;

        if ( moris::par_rank() == aRootProc )
        {
            std::unordered_map< moris_index, moris_index > tGlobalIdToCompactIndex;

            // first pass: figure out the unique set of vertices and their coordinates
            Vector< moris_index >      tUniqueGlobalIds;
            Vector< Matrix< DDRMat > > tUniqueCoordColumns;
            Vector< Matrix< DDRMat > > tUniqueDispColumns;

            uint tNumProcs = tGatheredVertexGlobalIds.size();
            for ( uint p = 0; p < tNumProcs; ++p )
            {
                Matrix< DDSMat > const & tProcVertIds = tGatheredVertexGlobalIds( p );
                Matrix< DDRMat > const & tProcCoords  = tGatheredVertexCoords( p );

                for ( uint iV = 0; iV < tProcVertIds.numel(); ++iV )
                {
                    moris_index tGlobalId = tProcVertIds( iV );

                    if ( tGlobalIdToCompactIndex.find( tGlobalId ) == tGlobalIdToCompactIndex.end() )
                    {
                        moris_index tCompactIndex            = (moris_index)tUniqueGlobalIds.size();
                        tGlobalIdToCompactIndex[ tGlobalId ] = tCompactIndex;

                        tUniqueGlobalIds.push_back( tGlobalId );
                        tUniqueCoordColumns.push_back( Matrix< DDRMat >( tProcCoords.get_column( iV ) ) );
                        // Also need disp for this vertex - find it in gathered disps
                        tUniqueDispColumns.push_back( Matrix< DDRMat >( tDim, 1, 0.0 ) );
                    }
                }
            }

            // Fill disp columns by looking up each unique global id in gathered disps
            {
                std::unordered_map< moris_index, Matrix< DDRMat > > tIdToDisp;
                for ( uint p = 0; p < tNumProcs; ++p )
                {
                    Matrix< DDSMat > const & tProcIdsDisp = tGatheredVertexGlobalIds( p );
                    Matrix< DDRMat > const & tProcDisps   = tGatheredVertexDisplacements( p );
                    for ( uint iV = 0; iV < tProcIdsDisp.numel(); ++iV )
                    {
                        moris_index tGId = tProcIdsDisp( iV );
                        if ( tIdToDisp.find( tGId ) == tIdToDisp.end() )
                        {
                            Matrix< DDRMat > tCol( tDim, 1 );
                            for ( uint d = 0; d < tDim; ++d )
                            {
                                tCol( d, 0 ) = tProcDisps( d, iV );
                            }
                            tIdToDisp[ tGId ] = tCol;
                        }
                    }
                }
                for ( uint i = 0; i < tUniqueGlobalIds.size(); ++i )
                {
                    auto it = tIdToDisp.find( tUniqueGlobalIds( i ) );
                    if ( it != tIdToDisp.end() )
                    {
                        tUniqueDispColumns( i ) = it->second;
                    }
                }
            }

            uint tNumUniqueVerts = tUniqueGlobalIds.size();

            aGlobalVertexIds.set_size( tNumUniqueVerts, 1 );
            aGlobalVertexCoords.set_size( tDim, tNumUniqueVerts );
            aGlobalVertexDisplacements.set_size( tDim, tNumUniqueVerts );

            for ( uint iV = 0; iV < tNumUniqueVerts; ++iV )
            {
                aGlobalVertexIds( iV ) = tUniqueGlobalIds( iV );
                aGlobalVertexCoords.set_column( iV, tUniqueCoordColumns( iV ) );
                aGlobalVertexDisplacements.set_column( iV, tUniqueDispColumns( iV ) );
            }

            Vector< moris_index > tFlatCounts;
            Vector< moris_index > tFlatIndices;
            Vector< moris_index > tFlatCellIds;
            Vector< moris_index > tFlatCellOwners;

            for ( uint p = 0; p < tNumProcs; ++p )
            {
                Matrix< DDSMat > const & tProcCounts  = tGatheredCellVertexCounts( p );
                Matrix< DDSMat > const & tProcConn    = tGatheredCellVertexGlobalIds( p );
                Matrix< DDSMat > const & tProcCellIds = tGatheredCellGlobalIds( p );

                uint tOffset = 0;
                for ( uint iC = 0; iC < tProcCounts.numel(); ++iC )
                {
                    uint tNumVertsInCell = (uint)tProcCounts( iC );
                    tFlatCounts.push_back( (moris_index)tNumVertsInCell );
                    tFlatCellIds.push_back( tProcCellIds( iC ) );
                    tFlatCellOwners.push_back( (moris_index)p );

                    for ( uint k = 0; k < tNumVertsInCell; ++k )
                    {
                        moris_index tGlobalId = tProcConn( tOffset + k );
                        tFlatIndices.push_back( tGlobalIdToCompactIndex.at( tGlobalId ) );
                    }
                    tOffset += tNumVertsInCell;
                }
            }

            tGlobalCellVertexCounts.set_size( tFlatCounts.size(), 1 );
            for ( uint i = 0; i < tFlatCounts.size(); ++i )
            {
                tGlobalCellVertexCounts( i ) = tFlatCounts( i );
            }

            tGlobalCellVertexIndices.set_size( tFlatIndices.size(), 1 );
            for ( uint i = 0; i < tFlatIndices.size(); ++i )
            {
                tGlobalCellVertexIndices( i ) = tFlatIndices( i );
            }

            aGlobalCellIds.set_size( tFlatCellIds.size(), 1 );
            for ( uint i = 0; i < tFlatCellIds.size(); ++i )
            {
                aGlobalCellIds( i ) = tFlatCellIds( i );
            }

            aGlobalCellOwners.set_size( tFlatCellOwners.size(), 1 );
            for ( uint i = 0; i < tFlatCellOwners.size(); ++i )
            {
                aGlobalCellOwners( i ) = tFlatCellOwners( i );
            }
        }


        moris::broadcast_mat( aGlobalVertexIds, aRootProc );
        moris::broadcast_mat( aGlobalVertexCoords, aRootProc );
        moris::broadcast_mat( aGlobalVertexDisplacements, aRootProc );
        moris::broadcast_mat( aGlobalCellIds, aRootProc );
        moris::broadcast_mat( aGlobalCellOwners, aRootProc );
        moris::broadcast_mat( tGlobalCellVertexCounts, aRootProc );
        moris::broadcast_mat( tGlobalCellVertexIndices, aRootProc );

        // every rank (including root) rebuilds aGlobalCells from the
        // now-identical flattened connectivity
        aGlobalCells.clear();

        uint tOffset = 0;
        for ( uint iC = 0; iC < tGlobalCellVertexCounts.numel(); ++iC )
        {
            uint tNumVertsInCell = (uint)tGlobalCellVertexCounts( iC );

            Matrix< IndexMat > tCellConn( tNumVertsInCell, 1 );
            for ( uint k = 0; k < tNumVertsInCell; ++k )
            {
                tCellConn( k ) = tGlobalCellVertexIndices( tOffset + k );
            }

            aGlobalCells.push_back( tCellConn );
            tOffset += tNumVertsInCell;
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::mtk