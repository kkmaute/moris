/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_simple_mesh.hpp
 *
 */

#pragma once

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris::ge
{
    /**
     * Creates a simple HMR interpolation mesh for GEN testing.
     *
     * @param aNumXElements Number of elements in the X direction
     * @param aNumYElements Number of elements in the Y direction
     * @param aLagrangeOrder Order of the lagrange mesh
     * @param aBSplineOrder Order of the B-spline mesh
     * @return HMR mesh pointer
     */
    mtk::Interpolation_Mesh* create_simple_mesh(
            uint aNumXElements,
            uint aNumYElements,
            uint aLagrangeOrder = 1,
            uint aBSplineOrder = 1,
            uint aRefinement = 0);
}
