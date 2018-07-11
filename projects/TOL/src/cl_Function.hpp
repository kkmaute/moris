/*
 * cl_Function.hpp
 *
 *  Created on: Oct 28, 2016
 *      Author: Keenan Doble
 */

#ifndef SRC_TOOLS_CL_FUNCTION_HPP_
#define SRC_TOOLS_CL_FUNCTION_HPP_

#include "ios.hpp"
#include "core.hpp"
#include "assert.hpp"
#include "cl_Mat.hpp" // LNA/src

namespace moris
{
    namespace tools
    {
        class Function
        {
        // Pure virtual function class
        public:

            Function();

            virtual ~Function() = default;

            moris::real
            evaluate_function(moris::Mat<moris::real> & aCoords);

            moris::Mat<moris::real>
            evaluate_dxdp(moris::Mat<moris::real> & aCoords);

        private:
            virtual moris::real compute_func_val(const moris::Mat<moris::real>    & aIndepVars)  = 0;
            virtual moris::Mat<moris::real> compute_func_dxdp(const moris::Mat<moris::real>    & aIndepVars) = 0;
        };
    }
}



#endif /* SRC_TOOLS_CL_FUNCTION_HPP_ */
