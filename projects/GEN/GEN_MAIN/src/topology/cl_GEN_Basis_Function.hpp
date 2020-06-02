/*
 * cl_XTK_Basis_Function.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_RIPPED_ADDITIONAL_CL_GE_ANALYTIC_HPP_
#define PROJECTS_GEN_SRC_RIPPED_ADDITIONAL_CL_GE_ANALYTIC_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
namespace ge
{
class Basis_Function
{
public:
    ~Basis_Function()
    {

    }

    /**
     * Evaluates the values of the basis functions at a give local or parametric coordinate [-1,1]
     */
    virtual void evaluate_basis_function(moris::Matrix< moris::DDRMat > const & aLocalCoordinate,
                                         moris::Matrix< moris::DDRMat > & aBasisFunctionValues) const = 0;
};

}
}

#endif /* PROJECTS_GEN_SRC_RIPPED_ADDITIONAL_CL_GE_ANALYTIC_HPP_ */
