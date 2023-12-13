//
// Created by frank on 11/30/23.
//

#ifndef MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP
#define MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP

#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_IntegrationPointPairs.hpp"

namespace moris::mtk
{
    class Nonconformal_Side_Cluster : public Double_Side_Cluster
    {
      public:
        Nonconformal_Side_Cluster(
                Cluster const                              *aFollowerSideCluster,
                Cluster const                              *aLeaderSideCluster,
                moris::Cell< IntegrationPointPairs > const &aIntegrationPointPairs )
                : Double_Side_Cluster( aLeaderSideCluster, aFollowerSideCluster, {} )
                , mIntegrationPointPairs( aIntegrationPointPairs ){};

      private:
        moris::Cell< IntegrationPointPairs > mIntegrationPointPairs;
    };


}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_NONCONFORMAL_SIDE_CLUSTER_HPP
