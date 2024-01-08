//
// Created by frank on 1/3/24.
//

#include "cl_FEM_Element_Nonconformal_Sideset.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Element_Double_Sideset.hpp"
#include <iostream>
#include "cl_MTK_Enums.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"

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

    moris_index Element_Nonconformal_Sideset::get_follower_local_cell_index() const
    {
        return mFollowerCellIndexInCluster;
    }

    void Element_Nonconformal_Sideset::initialize_leader_follower_ig_interpolator(
            mtk::Cell const           *aCell,
            moris_index const          aSideOrdinal,
            moris_index const          aLocalCellIndex,
            mtk::Leader_Follower const aLeaderFollowerType ) const
    {
        Element_Double_Sideset::initialize_leader_follower_ig_interpolator( aCell, aSideOrdinal, aLocalCellIndex, aLeaderFollowerType );

        // the nonconformal sideset element can provide precomputed integration point distances
        Geometry_Interpolator *tIGInterpolator = mSet->get_field_interpolator_manager( aLeaderFollowerType )->get_IG_geometry_interpolator();
        tIGInterpolator->set_integration_point_distances( mIntegrationPointDistances );    // TODO: find different way to access ray length in the IWGs/IQIs
    }
}    // namespace moris::fem