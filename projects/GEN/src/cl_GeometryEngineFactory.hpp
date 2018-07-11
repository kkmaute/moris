/*
 * cl_GeometryEngineFactory.hpp
 *
 *  Created on: Nov 8, 2016
 *      Author: doble
 */

#ifndef SRC_GEOMENG_CL_GEOMETRYENGINEFACTORY_HPP_
#define SRC_GEOMENG_CL_GEOMETRYENGINEFACTORY_HPP_

#include "ios.hpp"
#include "core.hpp"

#include "cl_FunctionFactory.hpp" // TOL/src

#include "cl_GeometryEngine_Levelset.hpp"
#include "cl_GeometryEngine.hpp"

namespace moris
{
class GeometryEngineFactory
{
public:
    GeometryEngineFactory(){};

    ~GeometryEngineFactory(){};

    moris::GeometryEngine*
    create_geometry_engine(enum FunctionType            aFunction,
                           moris::Mat<moris::real>      aFunctionParameters);

};
}
#endif /* SRC_GEOMENG_CL_GEOMETRYENGINEFACTORY_HPP_ */
