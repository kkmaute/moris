//
// Created by frank on 11/13/23.
//

#ifndef MORIS_CL_MTK_SPATIAL_INDEXER_BRUTEFORCE_H
#define MORIS_CL_MTK_SPATIAL_INDEXER_BRUTEFORCE_H


#include "cl_MTK_Spatial_Indexer.h"

namespace moris::mtk
{
    class Spatial_Indexer_BruteForce : public Spatial_Indexer
    {
      public:
        Spatial_Indexer_BruteForce( moris::Cell< mtk::Surface_Mesh > aSurfaceMeshes, moris::Cell< std::pair< moris_index, moris_index > > aSurfacePairs )
                : Spatial_Indexer( aSurfaceMeshes, aSurfacePairs ){};

        moris::Cell< Spatial_Indexing_Result > perform( real epsilon ) override;

      private:
        Spatial_Indexing_Result perform_on_mesh_pair( moris_index aSourceMeshIndex, moris_index aTargetMeshIndex, real aEpsilon );
    };
}    // namespace moris::mtk

#endif    // MORIS_CL_MTK_SPATIAL_INDEXER_BRUTEFORCE_H
