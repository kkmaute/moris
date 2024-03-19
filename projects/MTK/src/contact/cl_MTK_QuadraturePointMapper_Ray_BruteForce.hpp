//
// Created by frank on 11/28/23.
//

#pragma once

#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Spatial_Indexer_BruteForce.h"
#include "cl_Json_Object.hpp"
#include <deque>
#include <set>


namespace moris::mtk
{
    class QuadraturePointMapper_Ray_BruteForce : public QuadraturePointMapper_Ray
    {
      public:
        ~QuadraturePointMapper_Ray_BruteForce() override = default;

        QuadraturePointMapper_Ray_BruteForce(
                mtk::Integration_Mesh                                 *aIGMesh,
                Vector< Side_Set const * >                            &aSideSets,
                const Vector< std::pair< moris_index, moris_index > > &aCandidatePairs )
                : QuadraturePointMapper_Ray( aIGMesh, aSideSets, aCandidatePairs )
                , mSpatialIndexer( get_surface_meshes(), aCandidatePairs ){};

        MappingResult map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates, real aMaxNegativeRayLength ) const override;

      private:
        // methods
        void                    raycast_cell( moris_index aSourceCellIndex, moris_index aNumberOfRays, MappingResult &aMappingResult, Spatial_Indexing_Result const &aSpatialIndexingResult ) const;
        void                    check_ray_cell_intersection( MappingResult &aMappingResult, std::deque< moris_index > &aUnprocessedRays, moris_index aSourceMeshIndex, moris_index aTargetMeshIndex, moris_index aSourceCellIndex, moris_index aTargetCellIndex, uint aResultOffset ) const;
        void                    process_rays( moris_index aSourceCellIndex, Spatial_Indexing_Result const &aSpatialIndexingResult, uint aResultOffset, MappingResult &aMappingResult, std::deque< moris_index > &aUnprocessedRays, bool aBruteForce ) const;
        std::set< moris_index > get_potential_target_meshes( bool aBruteForce, moris_index aSourceMeshIndex, moris_index aSourceCellIndex, Spatial_Indexing_Result const &aSpatialIndexingResult ) const;
        std::set< moris_index > get_potential_target_cells( bool aBruteForce, moris_index aSourceMeshIndex, moris_index tTargetMeshIndex, moris_index aSourceCellIndex, Spatial_Indexing_Result const &aSpatialIndexingResult ) const;

        // members
        Spatial_Indexer_BruteForce mSpatialIndexer;
    };
}    // namespace moris::mtk