//
// Created by frank on 11/30/23.
//

#include "cl_MTK_Nonconformal_Side_Cluster.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_PointPairs.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"

#include <fn_sum.hpp>
#include <map>

namespace moris::mtk
{
    std::pair< moris_index, moris_index > get_leader_follower_index_pair(
            std::map< moris_index, moris_index > aLeaderCellIndexToStoringIndex,
            std::map< moris_index, moris_index > aFollowerCellIndexToStoringIndex,
            MappingPointPairs                    aMappingPointPairs )
    {
        moris_index const tLeaderCellIndex    = aMappingPointPairs.get_leader_cell_index();
        moris_index const tFollowerCellIndex  = aMappingPointPairs.get_follower_cell_index();
        moris_index const tLocalLeaderIndex   = aLeaderCellIndexToStoringIndex[ tLeaderCellIndex ];
        moris_index const tLocalFollowerIndex = aFollowerCellIndexToStoringIndex[ tFollowerCellIndex ];
        return std::make_pair( tLocalLeaderIndex, tLocalFollowerIndex );
    }

    std::map< moris_index, moris_index > get_cell_index_to_storing_index( Vector< mtk::Cell const * > aListOfCells )
    {
        std::map< moris_index, moris_index > tCellIndexToStoringIndex;
        for ( size_t iStoringIndex = 0; iStoringIndex < aListOfCells.size(); ++iStoringIndex )
        {
            tCellIndexToStoringIndex[ aListOfCells( iStoringIndex )->get_index() ] = iStoringIndex;
        }
        return tCellIndexToStoringIndex;
    }

    Nonconformal_Side_Cluster::NonconformalCellPairing Nonconformal_Side_Cluster::get_nonconformal_cell_pairing() const
    {
        std::map< moris_index, moris_index > tLeaderCellIndexToStoringIndex   = get_cell_index_to_storing_index( get_primary_cells_in_cluster( Leader_Follower::LEADER ) );
        std::map< moris_index, moris_index > tFollowerCellIndexToStoringIndex = get_cell_index_to_storing_index( get_primary_cells_in_cluster( Leader_Follower::FOLLOWER ) );

        NonconformalCellPairing tNonconformalCellPairing;
        for ( size_t iIPPIndex = 0; iIPPIndex < mIntegrationPointPairs.size(); iIPPIndex++ )
        {
            auto const &tIntegrationPointPair = mIntegrationPointPairs( iIPPIndex );

            std::pair< moris_index, moris_index > tLocalLeaderFollowerIndexPair = get_leader_follower_index_pair(
                    tLeaderCellIndexToStoringIndex, tFollowerCellIndexToStoringIndex, tIntegrationPointPair );

            tNonconformalCellPairing[ tLocalLeaderFollowerIndexPair ] = std::make_pair( iIPPIndex, -1 );
        }

        for ( size_t iNPPIndex = 0; iNPPIndex < mNodalPointPairs.size(); iNPPIndex++ )
        {
            auto const &tNodalPointPair = mNodalPointPairs( iNPPIndex );

            std::pair< moris_index, moris_index > const tLocalLeaderFollowerIndexPair = get_leader_follower_index_pair(
                    tLeaderCellIndexToStoringIndex, tFollowerCellIndexToStoringIndex, tNodalPointPair );

            if ( tNonconformalCellPairing.find( tLocalLeaderFollowerIndexPair ) != tNonconformalCellPairing.end() )
            {
                tNonconformalCellPairing[ tLocalLeaderFollowerIndexPair ].second = iNPPIndex;
            }
            else
            {
                tNonconformalCellPairing[ tLocalLeaderFollowerIndexPair ] = std::make_pair( -1, iNPPIndex );
            }
        }
        return tNonconformalCellPairing;
    }

}    // namespace moris::mtk
