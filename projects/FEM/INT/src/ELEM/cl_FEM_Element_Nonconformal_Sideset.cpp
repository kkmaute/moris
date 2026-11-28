/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Nonconformal_Sideset.cpp
 *
 */

#include "cl_FEM_Element_Nonconformal_Sideset.hpp"
#include "cl_FEM_Element_Double_Sideset.hpp"
#include "cl_MTK_Enums.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"

#include <cl_FEM_Model.hpp>

namespace moris::fem
{
    Matrix< DDRMat > Element_Nonconformal_Sideset::get_leader_integration_point( uint const aGPIndex ) const
    {
        return mIntegrationPointPairs.get_leader_coordinates().get_column( aGPIndex );
    }

    //----------------------------------------------------------------

    Matrix< DDRMat > Element_Nonconformal_Sideset::get_follower_integration_point( uint const aGPIndex ) const
    {
        return mIntegrationPointPairs.get_follower_coordinates().get_column( aGPIndex );
    }

    //----------------------------------------------------------------

    real Element_Nonconformal_Sideset::get_integration_weight( uint aGPIndex ) const
    {
        return mIntegrationPointPairs.get_integration_weights()( aGPIndex );
    }

    //----------------------------------------------------------------

    uint Element_Nonconformal_Sideset::get_number_of_integration_points() const
    {
        return mIntegrationPointPairs.get_integration_weights().size();
    }

    //----------------------------------------------------------------

    void Element_Nonconformal_Sideset::compute_jacobian_and_residual()
    {
        Element_Double_Sideset::compute_jacobian_and_residual();

        // undo the increment from the double sided sideset
        mSet->mFemModel->mDoubleSidedSideSetsGaussPoints -= get_number_of_integration_points();

        mSet->mFemModel->mNonconformalSideSetsGaussPoints += get_number_of_integration_points();
    }

    //----------------------------------------------------------------

    moris_index Element_Nonconformal_Sideset::get_follower_local_cell_index() const
    {
        return mFollowerCellIndexInCluster;
    }

    //----------------------------------------------------------------

    void Element_Nonconformal_Sideset::initialize_quadrature_point( uint iGP )
    {
        // get local integration point for the leader and follower integration cell
        const Matrix< DDRMat >& tLeaderLocalIntegPoint   = get_leader_integration_point( iGP );
        const Matrix< DDRMat >& tFollowerLocalIntegPoint = get_follower_integration_point( iGP );

        // set evaluation point for leader and follower interpolators
        mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )
                ->set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

        mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                ->set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );

        if ( mSet->get_time_continuity() )
        {
            mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER )
                    ->set_space_time_from_local_IG_point( tLeaderLocalIntegPoint );

            mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::FOLLOWER )
                    ->set_space_time_from_local_IG_point( tFollowerLocalIntegPoint );
        }
    }

}    // namespace moris::fem
