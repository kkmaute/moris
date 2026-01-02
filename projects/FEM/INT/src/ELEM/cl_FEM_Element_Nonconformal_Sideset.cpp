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

    void Element_Nonconformal_Sideset::initialize_leader_follower_ig_interpolator( mtk::Leader_Follower const aLeaderFollowerType ) const
    {
        // call the parent class function to
        Element_Double_Sideset::initialize_leader_follower_ig_interpolator( aLeaderFollowerType );

        if ( mSet->get_time_continuity() )
        {
          initialize_leader_follower_ig_interpolator_previous_time( aLeaderFollowerType );
        }
    }

    //------------------------------------------------------------------------------

    void Element_Nonconformal_Sideset::initialize_leader_follower_ig_interpolator_previous_time( mtk::Leader_Follower const aLeaderFollowerType ) const
    {
        // get treated side ordinal on the leader and on the follower
        moris_index      tSideOrd;
        moris_index      tLocalCellIndex;
        const mtk::Cell* tCell;
        if ( aLeaderFollowerType == mtk::Leader_Follower::LEADER )
        {
            tCell           = mLeaderCell;
            tLocalCellIndex = get_leader_local_cell_index();
            tSideOrd        = mCluster->mLeaderListOfSideOrdinals( tLocalCellIndex );
        }
        else
        {
            tCell           = mFollowerCell;
            tLocalCellIndex = get_follower_local_cell_index();
            tSideOrd        = mCluster->mFollowerListOfSideOrdinals( tLocalCellIndex );
        }

        // get previous time IG geometry interpolator
        Geometry_Interpolator* tPreviousIGGI = mSet->get_field_interpolator_manager_previous_time( aLeaderFollowerType )->get_IG_geometry_interpolator();

        // physical coefficients
        tPreviousIGGI->set_space_coeff( tCell->get_cell_physical_coords_on_side_ordinal( tSideOrd ) );
        tPreviousIGGI->set_time_coeff( mCluster->mInterpolationElement->get_previous_time() );

        // parametric coefficients
        tPreviousIGGI->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( tLocalCellIndex, tSideOrd, aLeaderFollowerType ) );
        tPreviousIGGI->set_time_param_coeff( { { -1.0 }, { 1.0 } } );    // FIXME: not true if time is not linear
    }

    //------------------------------------------------------------------------------

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
