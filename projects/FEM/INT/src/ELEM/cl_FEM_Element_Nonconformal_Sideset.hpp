/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Nonconformal_Sideset.hpp
 *
 */

#ifndef CL_FEM_ELEMENT_NONCONFORMAL_SIDESET_HPP
#define CL_FEM_ELEMENT_NONCONFORMAL_SIDESET_HPP
#include <cl_FEM_Element_Double_Sideset.hpp>

#include <cl_MTK_PointPairs.hpp>

namespace moris::fem
{
    class Element_Nonconformal_Sideset : public Element_Double_Sideset
    {
      public:
        Element_Nonconformal_Sideset(
                mtk::Cell const                  *aLeftIGCell,
                mtk::Cell const                  *aRightIGCell,
                Set                              *aSet,
                Cluster                          *aCluster,
                moris_index                       aLeaderCellIndexInCluster,
                moris_index                       aFollowerCellIndexInCluster,
                mtk::IntegrationPointPairs const &aIntegrationPointPairs,
                mtk::NodalPointPairs const       &aNodalPointPairs )
                : Element_Double_Sideset( aLeftIGCell, aRightIGCell, aSet, aCluster, aLeaderCellIndexInCluster )
                , mFollowerCellIndexInCluster( aFollowerCellIndexInCluster )
                , mIntegrationPointPairs( aIntegrationPointPairs )
                , mNodalPointPairs( aNodalPointPairs )
        {
        }

        ~Element_Nonconformal_Sideset() override = default;

        Matrix< DDRMat > get_leader_integration_point( uint const aGPIndex ) const override;
        Matrix< DDRMat > get_follower_integration_point( uint const aGPIndex ) const override;
        real             get_integration_weight( uint const aGPIndex ) const override;
        uint             get_number_of_integration_points() const override;
        void             compute_jacobian_and_residual() override;

      protected:
        moris_index get_follower_local_cell_index() const override;

      private:
        moris_index                mFollowerCellIndexInCluster;
        mtk::IntegrationPointPairs mIntegrationPointPairs;
        mtk::NodalPointPairs       mNodalPointPairs;
    };
}    // namespace moris::fem

#endif    // CL_FEM_ELEMENT_NONCONFORMAL_SIDESET_HPP
