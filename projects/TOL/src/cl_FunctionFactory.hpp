/*
 * cl_tools_Function.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: doble
 *
 *
 */

#ifndef SRC_TOOLS_CL_FUNCTIONFACTORY_HPP_
#define SRC_TOOLS_CL_FUNCTIONFACTORY_HPP_


#include "ios.hpp"
#include "core.hpp"

#include "cl_Enums.hpp"
#include "cl_Function_Levelset.hpp"


namespace moris
{

namespace tools
{
    class FunctionFactory
    {
    public:
        FunctionFactory(){};

        ~FunctionFactory(){};

        moris::tools::Function*
        create_explicit_function(enum FunctionType                       aFunctionType,
                                 moris::Mat<moris::real>                 aFunctionParameters);



    };

}
}


#endif /* SRC_TOOLS_CL_FUNCTIONFACTORY_HPP_ */
