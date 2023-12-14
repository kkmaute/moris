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
      public:    // constructors
        Spatial_Indexer_BruteForce( moris::Cell< mtk::Surface_Mesh > const &aSurfaceMeshes, const moris::Cell< std::pair< moris_index, moris_index > > &aSurfacePairs )
                : Spatial_Indexer( aSurfaceMeshes, aSurfacePairs ){};

        // methods
        Spatial_Indexing_Result perform( moris_index aSourceMeshIndex, real epsilon ) const override;

      private:    // methods
        Spatial_Indexing_Result perform_on_mesh_pair( moris_index aSourceMeshIndex, moris_index aTargetMeshIndex, real aEpsilon ) const;
    };
}    // namespace moris::mtk

#endif    // MORIS_CL_MTK_SPATIAL_INDEXER_BRUTEFORCE_H
