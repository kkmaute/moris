/*
 * fn_R2.hpp
 *
 *  Created on: Aug 16, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LNA_SRC_FN_R2_HPP_
#define PROJECTS_LNA_SRC_FN_R2_HPP_

#include "typedefs.hpp"
#include "cl_Mat.hpp"
#include "fn_sum.hpp"
#include "op_minus.hpp"

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
    real
    r2( const Mat< real > & aFunctionValues, const Mat< real > & aSamples )
    {
        // get number of samples
        uint tNumberOfSamples = aSamples.length();

        // calculate average of samples
        real tAverage = sum( aSamples ) / tNumberOfSamples;

        // sum of square residuals
        Mat< real > tDelta = aSamples - aFunctionValues;
        real tRootOfSSres = tDelta.norm();

        // total sum of squares
        tDelta = aSamples - tAverage;
        real tRootOfSStot = tDelta.norm();

        // return R2
        return 1.0 - std::pow( tRootOfSSres / tRootOfSStot , 2 );
    }
}



#endif /* PROJECTS_LNA_SRC_FN_R2_HPP_ */
