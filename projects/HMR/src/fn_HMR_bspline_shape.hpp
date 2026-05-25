/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_bspline_shape.hpp
 *
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"

namespace moris::hmr
{
    /**
     * Calculates a B-spline shape function at a given point
     *
     * @param aOrder Polynomial order
     * @param aBasisNumber Basis number
     * @param aXi Parametric coordinate
     * @return Evaluated shape function
     */
    real bspline_shape(
            uint aOrder,
            uint aBasisNumber,
            real aXi );

    //-------------------------------------------------------------------------------------
    
    /**
     * Calculates an extended B-spline shape function at a given point
     *
     * @param aOrder Polynomial order
     * @param aBasisNumber Basis number
     * @param aXi Parametric coordinate
     * @return Evaluated shape function
     */
    real bspline_shape_extended(
            uint aOrder,
            uint aBasisNumber,
            real aXi );

    //-------------------------------------------------------------------------------------

    real
    eval_spline(
            const uint aNumDims, // d
            const uint aOrder,   // p
            Vector< moris_index > const & aRelativeIJK,
            Matrix< DDRMat > const & aXi );

    real
    eval( 
            const uint aNumDims,
            const uint aOrder,
            const uint aBfLevel,
            const luint* aBfIJK,
            const uint aElementLevel,
            const luint* aElementIJK,
            const Matrix< DDRMat > & aXi );

    //-------------------------------------------------------------------------------------
}
