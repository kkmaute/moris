//
// Created by frank on 11/30/23.
//

#ifndef MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP
#define MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP

#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Integration_Rule.hpp"

namespace moris::mtk
{
    class Nonconformal_Side_Cluster : public Double_Side_Cluster
    {
      private:
        std::shared_ptr< Integration_Rule > mIntegrationRule;

      public:
        Nonconformal_Side_Cluster(
                moris::mtk::Cluster const                       *aLeaderSideCluster,
                moris::mtk::Cluster const                       *aFollowerSideCluster,
                moris::Cell< moris::mtk::Vertex const * > const &aLeaderToFollowerVertexPair,
                std::shared_ptr< Integration_Rule >              aIntegrationRule )
                : Double_Side_Cluster( aLeaderSideCluster, aFollowerSideCluster, aLeaderToFollowerVertexPair )
                , mIntegrationRule( aIntegrationRule ){};
    };


}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP
