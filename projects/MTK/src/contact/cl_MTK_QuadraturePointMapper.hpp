//
// Created by frank on 11/27/23.
//

#ifndef MORIS_CL_MTK_QUADRATUREPOINTMAPPER_HPP
#define MORIS_CL_MTK_QUADRATUREPOINTMAPPER_HPP

#include <ostream>
#include "cl_MTK_Side_Set.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_MappingResult.hpp"

namespace moris::mtk
{
    class Integration_Mesh;
    class QuadraturePointMapper
    {
      public:
        QuadraturePointMapper(
                Integration_Mesh                                           *aIGMesh,
                moris::Cell< Side_Set * >                                  &aSideSets,
                moris::Cell< std::pair< moris_index, moris_index > > const &aCandidatePairs );

        virtual ~QuadraturePointMapper() = default;

        virtual auto map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinate ) -> MappingResult = 0;

      protected:
        Integration_Mesh                                    *mIGMesh;
        moris::Cell< mtk::Side_Set * >                       mSideSets;
        moris::Cell< std::pair< moris_index, moris_index > > mCandidatePairs;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_QUADRATUREPOINTMAPPER_HPP
