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


#include <cl_FEM_Model.hpp>

namespace moris::fem
{
    Matrix< DDRMat > Element_Nonconformal_Sideset::get_leader_integration_point( uint const aGPIndex ) const
    {
        return mIntegrationPointPairs.get_leader_coordinates().get_column( aGPIndex );
    }

    Matrix< DDRMat > Element_Nonconformal_Sideset::get_follower_integration_point( uint const aGPIndex ) const
    {
        return mIntegrationPointPairs.get_follower_coordinates().get_column( aGPIndex );
    }

    real Element_Nonconformal_Sideset::get_integration_weight( uint aGPIndex ) const
    {
        return mIntegrationPointPairs.get_integration_weights()( aGPIndex );
    }

    uint Element_Nonconformal_Sideset::get_number_of_integration_points() const
    {
        return mIntegrationPointPairs.get_integration_weights().size();
    }

    void Element_Nonconformal_Sideset::compute_jacobian_and_residual()
    {
        Element_Double_Sideset::compute_jacobian_and_residual();
        mSet->mFemModel->mDoubleSidedSideSetsGaussPoints -= get_number_of_integration_points();    // undo the increment from the double sided sideset
        mSet->mFemModel->mNonconformalSideSetsGaussPoints += get_number_of_integration_points();
    }

    moris_index Element_Nonconformal_Sideset::get_follower_local_cell_index() const
    {
        return mFollowerCellIndexInCluster;
    }
}    // namespace moris::fem