/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Crosswind_Stabilization_Tools.hpp
 *
 */

#ifndef SRC_FEM_FN_FEM_IWG_CROSSWIND_STABILIZATION_TOOLS_HPP_
#define SRC_FEM_FN_FEM_IWG_CROSSWIND_STABILIZATION_TOOLS_HPP_

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Property.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        /**
         * compute crosswind projection of gradient of target dof
         * C gradx(w),
         * C = crosswind projection operator
         * w = target dof type for stabilization
         * @param[ in ] aVelocityDofGroup  list of velocity dof types
         * @param[ in ] aTargetDofGroup    list of target dof types
         * @param[ in ] aLeaderFIManager   field interpolator manager
         * @param[ out ] cgradxw
         */
        void compute_cgradxw(
                const Vector< MSI::Dof_Type > & aVelocityDofGroup,
                const Vector< MSI::Dof_Type > & aTargetDofGroup,
                Field_Interpolator_Manager         * aLeaderFIManager,
                const real                         & aSpaceDim,
                const real                         & aEpsilon,
                bool                               & aIsCrosswind,
                Matrix< DDRMat >                   & acgradxw );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of crosswind projection of gradient of target dof
         * wrt to a dof type
         * @param[ in ] aVelocityDofGroup  list of velocity dof types
         * @param[ in ] aTargetDofGroup    list of target dof types
         * @param[ in ] aLeaderFIManager   field interpolator manager
         * @param[ in ] aDofTypes          a list of dof type wrt which
         *                                 the derivative is requested
         * @param[ in ] adcgradxwdu        a matrix to fill with derivative
         */
        void compute_dcgradxwdu(
                const Vector< MSI::Dof_Type > & aVelocityDofGroup,
                const Vector< MSI::Dof_Type > & aTargetDofGroup,
                Field_Interpolator_Manager         * aLeaderFIManager,
                const Vector< MSI::Dof_Type > & aDofTypes,
                const real                         & aSpaceDim,
                const real                         & aEpsilon,
                bool                               & aIsCrosswind,
                Matrix< DDRMat >                   & adcgradxwdu );

    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_IWG_CROSSWIND_STABILIZATION_TOOLS_HPP_ */

