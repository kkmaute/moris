//
// Created by frank on 3/17/24.
//

#pragma once
#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp"

namespace moris::mtk
{
    class QuadraturePointMapper_ArborX : public QuadraturePointMapper_Ray
    {
      public:
        QuadraturePointMapper_ArborX(
                mtk::Integration_Mesh                                 *aIGMesh,
                Vector< Side_Set const * >                            &aSideSets,
                const Vector< std::pair< moris_index, moris_index > > &aCandidatePairs )
                : QuadraturePointMapper_Ray( aIGMesh, aSideSets, aCandidatePairs ){};

        MappingResult map( moris_index aSourceMeshIndex, Matrix< DDRMat > const &aParametricCoordinates, real aMaxNegativeRayLength, real aMaxPositiveRayLength ) const override;

      private:
        Vector< std::pair< moris_index, Surface_Mesh > > get_target_surface_meshes( moris_index aSourceMeshIndex ) const;
        void                                             check_cell_intersections( MappingResult &tMappingResult, real aMaxNegativeRayLength, real aMaxPositiveRayLength, arborx::cell_locator_map const &tBoxRayMap ) const;
    };
}    // namespace moris::mtk
