/*
 * cl_GeometryEngineFactory.cpp
 *
 *  Created on: Nov 8, 2016
 *      Author: doble
 */
#include "cl_GeometryEngineFactory.hpp"

moris::GeometryEngine*
moris::GeometryEngineFactory::create_geometry_engine(enum FunctionType         aFunctionType,
                                                     moris::Mat<moris::real>   aFunctionParameters)
{
    moris::tools::FunctionFactory tFcnFactory;
    moris::tools::Function* pFunction = tFcnFactory.create_explicit_function(aFunctionType,aFunctionParameters);

    moris::GeometryEngine* pGeometryEngine = NULL;
    switch (aFunctionType)
    {
    case(FunctionType::LEVELSET_SPHERE):
        {
        pGeometryEngine = new LevelsetEngine(pFunction);
        break;
        }
    case(FunctionType::LEVELSET_PLANE):
        {
        pGeometryEngine = new LevelsetEngine(pFunction);
        break;
        }
    case(FunctionType::LEVELSET_RANDOM):
        {
        pGeometryEngine = new LevelsetEngine(pFunction);
        break;
        }
    }

    return pGeometryEngine;
}

