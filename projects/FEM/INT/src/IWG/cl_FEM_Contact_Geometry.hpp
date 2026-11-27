/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#ifndef MORIS_CL_FEM_CONTACT_GEOMETRY_HPP
#define MORIS_CL_FEM_CONTACT_GEOMETRY_HPP

#include <memory>
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris::fem
{
    class Field_Interpolator_Manager;
    struct GapData;

    /** Stateless leader-to-follower ray mapping for contact assembly.
     * Managers and gap data are borrowed. The caller owns and invalidates its point
     * cache and supplies previous-time managers when tangential motion is needed.
     * Mapping updates follower interpolation points; no quadrature weight is
     * applied here. The existing invalid-intersection sentinel is preserved.
     */
    class Contact_Geometry
    {
      public:
        static Matrix< DDRMat > remap_nonconformal_rays(
                const bool                     aUseDeformedGeometryForGap,
                const Vector< MSI::Dof_Type >& aDisplDofTypes,
                Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
                Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
                Field_Interpolator_Manager*    aLeaderPreviousFieldInterpolatorManager,
                Field_Interpolator_Manager*    aFollowerPreviousFieldInterpolatorManager,
                std::unique_ptr< GapData >&    aGapData );

        static Matrix< DDRMat > remap_nonconformal_rays_deformed_geometry(
                const Vector< MSI::Dof_Type >& aDisplDofTypes,
                Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
                Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
                Field_Interpolator_Manager*    aLeaderPreviousFieldInterpolatorManager,
                Field_Interpolator_Manager*    aFollowerPreviousFieldInterpolatorManager,
                std::unique_ptr< GapData >&    aGapData );

        static Matrix< DDRMat > remap_nonconformal_rays_undeformed_geometry(
                Field_Interpolator_Manager* aLeaderFieldInterpolatorManager,
                Field_Interpolator_Manager* aFollowerFieldInterpolatorManager );
    };
}    // namespace moris::fem

#endif
