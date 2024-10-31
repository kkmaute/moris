/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper_Ray_Dummy.hpp
 *
 */
#pragma once

#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_MTK_Integration_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris::mtk
{
    class QuadraturePointMapper_Ray_Dummy : public QuadraturePointMapper_Ray
    {
      public:
        ~QuadraturePointMapper_Ray_Dummy() override = default;

        QuadraturePointMapper_Ray_Dummy(
                mtk::Integration_Mesh                                 *aIGMesh,
                Vector< Side_Set const * >                            &aSideSets,
                const Vector< std::pair< moris_index, moris_index > > &aCandidatePairs )
                : QuadraturePointMapper_Ray( aIGMesh, aSideSets, aCandidatePairs ){};

        MappingResult map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates, real aMaxNegativeRayLength, real aMaxPositiveRayLength ) const override
        {
            MORIS_ERROR( false, "To use nonconformal side sets with ray-based mapping, please install ArborX and activate the CMAKE Flag MORIS_HAVE_ARBORX!" );
            return { 0, 0, 0 };
        };
    };
}    // namespace moris::mtk
