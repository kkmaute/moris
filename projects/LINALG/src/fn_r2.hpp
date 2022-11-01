/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_r2.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_R2_HPP_
#define PROJECTS_LINALG_SRC_FN_R2_HPP_

// MORIS library header files.
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "op_minus.hpp"
#include "fn_sum.hpp"
#include "fn_norm.hpp"
#include "op_div.hpp"

namespace moris
{
    /**
     * calculates the R2 coefficient of determination from values
     * of an evaluated function with respect to given samples.
     *
     * See also
     *
     * https://en.wikipedia.org/wiki/Coefficient_of_determination
     *
     * @param[in] aFunctionValues   values of evaluated function
     * @param[in] aSamples          data samples or exact solution
     */
    inline real
    r2(
            const Matrix< DDRMat >& aFunctionValues,
            const Matrix< DDRMat >& aSamples )
    {

        // calculate average of samples

        double tValue = moris::sum( aSamples );
        // real tAverage = 1.0;

        double tAverage = tValue / (double)aSamples.length();

        // sum of square residuals
        real tRootOfSSres = norm( aSamples - aFunctionValues );

        // total sum of squares
        real tRootOfSStot = norm( aSamples - tAverage );

        // return R2
        return 1.0 - std::pow( tRootOfSSres / tRootOfSStot, 2 );
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_R2_HPP_ */
