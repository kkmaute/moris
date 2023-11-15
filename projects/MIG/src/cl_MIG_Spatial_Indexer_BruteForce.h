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
                const moris::Matrix< DDRMat >                   &aCoordinates,
                const moris::Cell< moris::Cell< moris_index > > &aNeighbors,
                const moris::Matrix< moris::DDRMat >            &aDisplacements )
                : Spatial_Indexer_Base( aCoordinates, aNeighbors, aDisplacements )
        {
        }

        moris::map< moris_index, moris_index > perform( real epsilon ) override;
    };
}    // namespace moris::mig

#endif    // MORIS_CL_MIG_SPATIAL_INDEXER_BRUTEFORCE_H
