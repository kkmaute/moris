/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper.hpp
 *
 */
#ifndef MORIS_CL_MTK_QUADRATUREPOINTMAPPER_HPP
#define MORIS_CL_MTK_QUADRATUREPOINTMAPPER_HPP

#include <ostream>
#include "cl_MTK_Side_Set.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_MappingResult.hpp"

namespace moris::mtk
{
    class Integration_Mesh;
    class QuadraturePointMapper
    {
      public:
        QuadraturePointMapper(
                Integration_Mesh                                      *aIGMesh,
                Vector< Side_Set const * >                            &aSideSets,
                Vector< std::pair< moris_index, moris_index > > const &aCandidatePairs );

        virtual ~QuadraturePointMapper() = default;

        virtual MappingResult map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinate, real aMaxNegativeRayLength, real aMaxPositiveRayLength ) const = 0;

        virtual void update_displacements( std::unordered_map< moris_index, Vector< real > > const &aSetDisplacements ) = 0;

      protected:
        Vector< Side_Set const * > const                      &get_side_sets() const { return mSideSets; }
        Vector< std::pair< moris_index, moris_index > > const &get_candidate_pairs() const { return mCandidatePairs; }

      private:
        Integration_Mesh                               *mIGMesh;
        Vector< Side_Set const * >                      mSideSets;
        Vector< std::pair< moris_index, moris_index > > mCandidatePairs;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_QUADRATUREPOINTMAPPER_HPP
