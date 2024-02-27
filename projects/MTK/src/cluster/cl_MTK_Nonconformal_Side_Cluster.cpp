//
// Created by frank on 11/30/23.
//

#include "cl_MTK_Nonconformal_Side_Cluster.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_IntegrationPointPairs.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"

#include <fn_sum.hpp>
#include <map>

namespace moris::mtk
{
    Vector< mtk::Cell const * >
    Nonconformal_Side_Cluster::get_nonconforming_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        Vector< moris_index >       tCellMask     = get_cell_local_indices( aIsLeader );
        Vector< mtk::Cell const * > tPrimaryCells = get_primary_cells_in_cluster( aIsLeader );
        Vector< mtk::Cell const * > tNonconformingPrimaryCells;
        tNonconformingPrimaryCells.reserve( tCellMask.size() );
        for ( auto const &tCellIndex : tCellMask )
        {
            tNonconformingPrimaryCells.push_back( tPrimaryCells( tCellIndex ) );
        }
        return tNonconformingPrimaryCells;
    }

    moris::Matrix< moris::IndexMat > Nonconformal_Side_Cluster::get_nonconforming_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const
    {
        Vector< moris_index >    tCellMask     = get_cell_local_indices( aIsLeader );
        Matrix< IndexMat > const tSideOrdinals = get_cell_side_ordinals( aIsLeader );

        moris::Matrix< moris::IndexMat > tCellOrdinals( 1, tCellMask.size() );
        for ( size_t iIndex = 0; iIndex < tCellMask.size(); ++iIndex )
        {
            // the tCellLocalIndex is the index at which the cell is stored in the vector of primary cells in this cluster, it is not the index of the cell itself!
            moris_index const tCellLocalIndex = tCellMask( iIndex );
            tCellOrdinals( 0, iIndex )        = tSideOrdinals( tCellLocalIndex );
        }
        return tCellOrdinals;
    }

    Vector< MappingPointPairs > const &Nonconformal_Side_Cluster::get_integration_point_pairs() const
    {
        return mIntegrationPointPairs;
    }

    Vector< moris_index > Nonconformal_Side_Cluster::get_cell_local_indices( const mtk::Leader_Follower aIsLeader ) const
    {
        // this map stores the index of the cell in the cluster to the index of the cell in the vector of primary cells
        std::map< moris_index, moris_index > tCellIndexToStoringIndex;
        Vector< mtk::Cell const * >          tPrimaryCells = get_primary_cells_in_cluster( aIsLeader );
        for ( size_t iStoringIndex = 0; iStoringIndex < tPrimaryCells.size(); ++iStoringIndex )
        {
            tCellIndexToStoringIndex[ tPrimaryCells( iStoringIndex )->get_index() ] = iStoringIndex;
        }

        Vector< moris_index > tCellMask;
        tCellMask.reserve( mIntegrationPointPairs.size() );
        for ( auto const &tIntegrationPointPair : mIntegrationPointPairs )
        {
            moris_index const tCellIndex =
                    ( aIsLeader == mtk::Leader_Follower::LEADER )
                            ? tIntegrationPointPair.get_leader_cell_index()
                            : tIntegrationPointPair.get_follower_cell_index();
            tCellMask.push_back( tCellIndexToStoringIndex[ tCellIndex ] );
        }

        return tCellMask;
    }

}    // namespace moris::mtk
