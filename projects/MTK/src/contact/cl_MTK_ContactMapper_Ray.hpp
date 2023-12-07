//
// Created by frank on 11/28/23.
//

#ifndef MORIS_CL_MTK_CONTACTMAPPER_RAY_HPP
#define MORIS_CL_MTK_CONTACTMAPPER_RAY_HPP

#include "cl_MTK_ContactMapper.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Spatial_Indexer_BruteForce.h"

namespace moris::mtk
{
    class ContactMapper_Ray : public ContactMapper
    {
      public:
        ContactMapper_Ray(
                mtk::Integration_Mesh                                      *aIGMesh,
                moris::Cell< mtk::Side_Set * >                             &aSideSets,
                const moris::Cell< std::pair< moris_index, moris_index > > &aCandidatePairs );

        ContactMapper::MappingResult map( Matrix< DDRMat > const &aParametricCoordinates ) override;

      private:    // methods
        ContactMapper::MappingResult raycast_bundle(
                moris_index                       aSourceMesh,
                moris_index                       aCellIndex,
                const moris::Cell< moris_index > &aVertexIndices,
                Spatial_Indexing_Result const    &aSpatialIndexingResult,
                Matrix< DDRMat > const           &aPhysicalCoordinates,
                Matrix< DDRMat > const           &aNormals );

        static std::pair< bool, double > calculate_ray_line_intersection( Matrix< DDRMat > const &aRayOrigin, Matrix< DDRMat > const &aRayDirection, Matrix< DDRMat > const &aSegmentOrigin, Matrix< DDRMat > const &aSegmentDirection );

      private:    // data
        moris::Cell< Surface_Mesh > mSurfaceMeshes{};
        Spatial_Indexer_BruteForce  mSpatialIndexer;
    };
}    // namespace moris::mtk

#endif    // MORIS_CL_MTK_CONTACTMAPPER_RAY_HPP
