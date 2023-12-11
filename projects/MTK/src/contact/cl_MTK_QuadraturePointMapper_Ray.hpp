//
// Created by frank on 11/28/23.
//

#ifndef MORIS_CL_MTK_QUADRATUREPOINTMAPPER_RAY_HPP
#define MORIS_CL_MTK_QUADRATUREPOINTMAPPER_RAY_HPP

#include "cl_MTK_QuadraturePointMapper.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Spatial_Indexer_BruteForce.h"
#include "cl_Json_Object.hpp"


namespace moris::mtk
{
    class QuadraturePointMapper_Ray : public QuadraturePointMapper
    {
      public:
        QuadraturePointMapper_Ray(
                mtk::Integration_Mesh                                      *aIGMesh,
                moris::Cell< mtk::Side_Set * >                             &aSideSets,
                const moris::Cell< std::pair< moris_index, moris_index > > &aCandidatePairs );

        MappingResult map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates ) override;

      private:    // methods
        void raycast_cell(
                moris_index                       aCellIndex,
                moris_index                       aNumberOfRays,
                MappingResult                    &aMappingResult,
                moris::Cell< moris_index > const &aVertexIndices,
                Spatial_Indexing_Result const    &aSpatialIndexingResult );

        static std::tuple< bool, real, Matrix< DDRMat >, Matrix< DDRMat > > calculate_ray_line_intersection(
                Matrix< DDRMat > const &aRayOrigin,
                Matrix< DDRMat > const &aRayDirection,
                Matrix< DDRMat > const &aSegmentOrigin,
                Matrix< DDRMat > const &aSegmentDirection );

      private:    // methods
        /**
         * @brief
         * @param aIGMesh
         * @param aSideSets
         * @return
         */
        static moris::Cell< Surface_Mesh > initialize_surface_meshes( Integration_Mesh *aIGMesh, moris::Cell< mtk::Side_Set * > const &aSideSets );

        /**
         * @brief
         * @param aCellIndex
         * @param aParametricCoordinates
         * @param aCoordinateInterpolator
         * @param aNormalInterpolator
         * @param aSurfaceMesh
         * @param aMappingResult
         */
        static void interpolate_source_point( moris_index aCellIndex, Matrix< DDRMat > const &aParametricCoordinates, Space_Interpolator &aCoordinateInterpolator, Space_Interpolator &aNormalInterpolator, Surface_Mesh const &aSurfaceMesh, MappingResult &aMappingResult ) ;

      private:    // data
        moris::Cell< Surface_Mesh > mSurfaceMeshes{};
        Spatial_Indexer_BruteForce  mSpatialIndexer;
    };
}    // namespace moris::mtk

#endif    // MORIS_CL_MTK_QUADRATUREPOINTMAPPER_RAY_HPP
