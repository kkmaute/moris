//
// Created by frank on 1/3/24.
//

#include "cl_FEM_Element_Nonconformal_Sideset.hpp"

namespace moris::fem
{
    Matrix< DDRMat > Element_Nonconformal_Sideset::get_leader_integration_point( uint const aGPIndex ) const
    {
        return mLeaderIntegrationPoints.get_column( aGPIndex );
    }

    Matrix< DDRMat > Element_Nonconformal_Sideset::get_follower_integration_point( uint const aGPIndex ) const
    {
        return mFollowerIntegrationPoints.get_column( aGPIndex );
    }

    real Element_Nonconformal_Sideset::get_integration_weight( uint aGPIndex ) const
    {
        return mIntegrationPointWeights( aGPIndex );
    }

    uint Element_Nonconformal_Sideset::get_number_of_integration_points() const
    {
        return mIntegrationPointWeights.size();
    }

    void Element_Nonconformal_Sideset::compute_jacobian_and_residual()
    {
        std::cout << "NCS: L: "
                  << mCluster->get_mesh_cluster()->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER )( get_leader_local_cell_index() )->get_index()
                  << " F: "
                  << mCluster->get_mesh_cluster()->get_primary_cells_in_cluster( mtk::Leader_Follower::FOLLOWER )( get_follower_local_cell_index() )->get_index();

        for ( auto const &tCoordinate : mLeaderIntegrationPoints.get_row( 0 ) )
        {
            std::cout << " " << std::setw( 4 ) << tCoordinate;
        }

        std::cout << "\n";
        Element_Double_Sideset::compute_jacobian_and_residual();
    }

    moris_index Element_Nonconformal_Sideset::get_follower_local_cell_index() const
    {
        return mFollowerCellIndexInCluster;
    }
}    // namespace moris::fem