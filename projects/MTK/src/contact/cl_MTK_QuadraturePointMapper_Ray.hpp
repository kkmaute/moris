/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_QuadraturePointMapper_Ray.hpp
 *
 */
#pragma once

#include <deque>
#include <set>
#include "cl_MTK_QuadraturePointMapper.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_Json_Object.hpp"

namespace moris::mtk
{
    class QuadraturePointMapper_Ray : public QuadraturePointMapper
    {

      public:
        ~QuadraturePointMapper_Ray() override = default;

        QuadraturePointMapper_Ray(
                mtk::Integration_Mesh                                 *aIGMesh,
                Vector< Side_Set const * >                            &aSideSets,
                const Vector< std::pair< moris_index, moris_index > > &aCandidatePairs );

        void update_displacements( std::unordered_map< moris_index, Vector< real > > const &aSetDisplacements ) override;

      protected:
        Vector< Surface_Mesh > const &get_surface_meshes() const { return mSurfaceMeshes; }
        Vector< Surface_Mesh > const &get_reference_surface_meshes() const { return mReferenceSurfaceMeshes; }
        MappingResult                 initialize_source_points( moris_index aSourceMeshIndex, Matrix< DDRMat > const &aParametricCoordinates ) const;

      private:
        static auto initialize_surface_meshes(
                Integration_Mesh const                *aIGMesh,
                Vector< mtk::Side_Set const * > const &aSideSets ) -> Vector< Surface_Mesh >;

        void write_surface_mesh_json() const;

        // data
        Vector< Surface_Mesh > mSurfaceMeshes;             // stores the surface meshes in their current (possibly deformed) state
        Vector< Surface_Mesh > mReferenceSurfaceMeshes;    // stores the surface meshes in their original state (not deformed)
    };
}    // namespace moris::mtk
