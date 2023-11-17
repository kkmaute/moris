//
// Created by frank on 11/13/23.
//

#ifndef MORIS_CL_MIG_SPATIAL_INDEXER_BRUTEFORCE_H
#define MORIS_CL_MIG_SPATIAL_INDEXER_BRUTEFORCE_H


#include "cl_MIG_Spatial_Indexer_Base.h"

namespace moris::mig
{
    class Spatial_Indexer_BruteForce : public Spatial_Indexer_Base
    {


      public:
        Spatial_Indexer_BruteForce(
                const Matrix< DDRMat >                          &aCoordinates,
                const moris::Cell< moris::Cell< moris_index > > &aNeighbors,
                const Matrix< DDRMat >                          &aDisplacements,
                const Matrix< DDRMat >                          &aVertexNormals );

        Spatial_Indexer_BruteForce( mtk::Surface_Mesh const &aSurfaceMesh, Matrix< DDRMat > const &aDisplacements );


        moris::map< moris_index, moris_index > perform( real epsilon ) override;
    };
}    // namespace moris::mig

#endif    // MORIS_CL_MIG_SPATIAL_INDEXER_BRUTEFORCE_H
