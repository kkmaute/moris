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
#include <deque>
#include <set>


namespace moris::mtk
{
    class QuadraturePointMapper_Ray : public QuadraturePointMapper
    {
      public:
        ~QuadraturePointMapper_Ray() override = default;

        QuadraturePointMapper_Ray(
                mtk::Integration_Mesh                                      *aIGMesh,
                Vector< Side_Set const * >                            &aSideSets,
                const Vector< std::pair< moris_index, moris_index > > &aCandidatePairs );

        MappingResult map( moris_index aSourceSideSetIndex, Matrix< DDRMat > const &aParametricCoordinates ) const override;

      private:
        // methods
        void raycast_cell( moris_index aSourceCellIndex, moris_index aNumberOfRays, MappingResult &aMappingResult, Spatial_Indexing_Result const &aSpatialIndexingResult ) const;

        /**
         * @brief
         * @param aIGMesh
         * @param aSideSets
         * @return
         */
        static auto initialize_surface_meshes(
                Integration_Mesh const                     *aIGMesh,
                Vector< mtk::Side_Set const * > const &aSideSets ) -> Vector< Surface_Mesh >;

        /**
         * @brief
         * @param aCellIndex
         * @param aParametricCoordinates
         * @param aCoordinateInterpolator
         * @param aNormalInterpolator
         * @param aSurfaceMesh
         * @param aMappingResult
         */
        static void interpolate_source_point( moris_index aCellIndex, Matrix< DDRMat > const &aParametricCoordinates, Space_Interpolator &aCoordinateInterpolator, Space_Interpolator &aNormalInterpolator, Surface_Mesh const &aSurfaceMesh, MappingResult &aMappingResult );

        // data
        Vector< Surface_Mesh > mSurfaceMeshes;
        Spatial_Indexer_BruteForce  mSpatialIndexer;
        void                        check_ray_cell_intersection( MappingResult &aMappingResult, std::deque< moris_index > &aUnprocessedRays, moris_index aSourceMeshIndex, moris_index aTargetMeshIndex, moris_index aSourceCellIndex, moris_index aTargetCellIndex, uint aResultOffset ) const;
        void                        process_rays( moris_index aSourceCellIndex, Spatial_Indexing_Result const &aSpatialIndexingResult, uint aResultOffset, MappingResult &aMappingResult, std::deque< moris_index > &aUnprocessedRays, bool aBruteForce ) const;
        std::set< moris_index >     get_potential_target_meshes( bool aBruteForce, moris_index aSourceMeshIndex, moris_index aSourceCellIndex, Spatial_Indexing_Result const &aSpatialIndexingResult ) const;
        std::set< moris_index >     get_potential_target_cells( bool aBruteForce, moris_index aSourceMeshIndex, moris_index tTargetMeshIndex, moris_index aSourceCellIndex, Spatial_Indexing_Result const &aSpatialIndexingResult ) const;
    };
}    // namespace moris::mtk

#endif    // MORIS_CL_MTK_QUADRATUREPOINTMAPPER_RAY_HPP
