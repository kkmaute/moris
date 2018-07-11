/*
 * cl_tools_Function.cpp
 *
 *  Created on: Oct 20, 2016
 *      Author: doble
 */
#include "cl_FunctionFactory.hpp"


/* ------------------------------------------------------------------------- */

moris::tools::Function*
moris::tools::FunctionFactory::create_explicit_function(enum FunctionType          aFunctionType,
                                                        moris::Mat<moris::real>    aFunctionParameters)
{

    moris::tools::Function* pFunction = NULL;
    switch (aFunctionType)
    {
        case (FunctionType::LEVELSET_SPHERE):
            {

                pFunction = new FunctionLevelset(aFunctionParameters, aFunctionType);
                break;
            }
        case (FunctionType::LEVELSET_PLANE):
            {

                pFunction = new FunctionLevelset(aFunctionParameters, aFunctionType);
                break;
            }
        case (FunctionType::LEVELSET_RANDOM):
            {
                pFunction = new FunctionLevelset(aFunctionParameters, aFunctionType);
                break;
            }
//        case (FunctionType::BezierSpline):
//        {
//            pFunction = new BezierSpline(aFunctionParameters,poly)
//            break;
//        }
    }

    return pFunction;
}

/* ------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------- */


