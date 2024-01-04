//
// Created by frank on 1/3/24.
//

#ifndef CL_FEM_ELEMENT_NONCONFORMAL_SIDESET_HPP
#define CL_FEM_ELEMENT_NONCONFORMAL_SIDESET_HPP
#include <cl_FEM_Element_Double_Sideset.hpp>

#include <cl_MTK_IntegrationPointPairs.hpp>

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
                moris::moris_index                aLeaderCellIndexInCluster,
                moris::moris_index                aFollowerCellIndexInCluster,
                mtk::IntegrationPointPairs const &aIntegrationPointPairs )
                : Element_Double_Sideset( aLeftIGCell, aRightIGCell, aSet, aCluster, aLeaderCellIndexInCluster )
                , mFollowerCellIndexInCluster( aFollowerCellIndexInCluster )
                , mLeaderIntegrationPoints( aIntegrationPointPairs.get_leader_coordinates() )
                , mFollowerIntegrationPoints( aIntegrationPointPairs.get_follower_coordinates() )
                , mIntegrationPointWeights( aIntegrationPointPairs.get_integration_weights() )
        {
        }

        ~Element_Nonconformal_Sideset() override = default;

        Matrix< DDRMat > get_leader_integration_point( uint const aGPIndex ) const override;
        Matrix< DDRMat > get_follower_integration_point( uint const aGPIndex ) const override;
        real             get_integration_weight( uint const aGPIndex ) const override;
        uint             get_number_of_integration_points() const override;

        void compute_jacobian_and_residual() override;

      protected:
        moris_index get_follower_local_cell_index() const override;

      private:
        moris_index      mFollowerCellIndexInCluster;
        Matrix< DDRMat > mLeaderIntegrationPoints;
        Matrix< DDRMat > mFollowerIntegrationPoints;
        Vector< real >   mIntegrationPointWeights;
    };
}    // namespace moris::fem

#endif    // CL_FEM_ELEMENT_NONCONFORMAL_SIDESET_HPP
