/*
 * cl_tools_Function_Levelset.cpp
 *
 *  Created on: Oct 20, 2016
 *      Author: doble
 */

#include "cl_Function_Levelset.hpp"

moris::tools::FunctionLevelset::~FunctionLevelset()
{

}

moris::tools::FunctionLevelset::FunctionLevelset( moris::Mat<moris::real>         aFunctionParameters,
                                                  enum FunctionType               aLevelsetType )
{
    mFuncParameters = aFunctionParameters;
    mFuncType       = aLevelsetType;
}

//-------------------------------------------------------------------------------

moris::real
moris::tools::FunctionLevelset::compute_func_val(const moris::Mat<moris::real>        & aCoords)
{
    moris::real tFuncVal = 0.0;

    switch (mFuncType)
    {
        case(FunctionType::LEVELSET_SPHERE):
            {

            tFuncVal = (aCoords(0,0)-mFuncParameters(1,0))*(aCoords(0,0)-mFuncParameters(1,0))
                               +(aCoords(0,1)-mFuncParameters(2,0))*(aCoords(0,1)-mFuncParameters(2,0))
                               +(aCoords(0,2)-mFuncParameters(3,0))*(aCoords(0,2)-mFuncParameters(3,0))
                               -(mFuncParameters(0,0)*mFuncParameters(0,0));
            break;
            }
        case(FunctionType::LEVELSET_PLANE):
            {
            // x - center = mFuncParameters(0,0)
            // y - center = mFuncParameters(0,1)
            // z - center = mFuncParameters(0,2)
            // x - normal = mFuncParameters(0,3)
            // y - normal = mFuncParameters(0,4)
            // z - normal = mFuncParameters(0,5)

            tFuncVal = mFuncParameters(0,3)*(aCoords(0,0)-mFuncParameters(0,0))+
                    mFuncParameters(0,4)*(aCoords(0,1)-mFuncParameters(0,1))+
                    mFuncParameters(0,5)*(aCoords(0,2)-mFuncParameters(0,2));
            //
            break;
            }

        case(FunctionType::LEVELSET_RANDOM):
            {
            tFuncVal = std::rand() % RAND_MAX;
            tFuncVal = tFuncVal-RAND_MAX/2;
            break;
            }
        default:
        {
            MORIS_ASSERT(0 , "Function type not implemented." );
        }
    }

    return tFuncVal;
}


moris::Mat<moris::real>
moris::tools::FunctionLevelset::compute_func_dxdp(const moris::Mat<moris::real>        & aCoords)
{
    moris::Mat<moris::real> tFuncVal(mFuncParameters.n_rows(),3);

    switch (mFuncType)
    {
        case(FunctionType::LEVELSET_SPHERE):
            {
            moris::real sign = INT_MAX;


            // dx/dr
            moris::real tSqrt =    mFuncParameters(0,0)*mFuncParameters(0,0)
                                - (aCoords(0,1)-mFuncParameters(2,0))*(aCoords(0,1)-mFuncParameters(2,0))
                                - (aCoords(0,2)-mFuncParameters(3,0))*(aCoords(0,2)-mFuncParameters(3,0));
            if( tSqrt <0)
            {
                sign  = -1;
            }
            else if(tSqrt>0)
                {
                    sign = 1;
                }
            else
                MORIS_LOG_ERROR<<"zero denominator detected";

            tFuncVal(0,0)  = sign * mFuncParameters(0,0)/std::sqrt(abs(tSqrt));

            //dy/dr
            tSqrt =    mFuncParameters(0,0)*mFuncParameters(0,0)
                    - (aCoords(0,0)-mFuncParameters(1,0))*(aCoords(0,0)-mFuncParameters(1,0))
                    - (aCoords(0,2)-mFuncParameters(3,0))*(aCoords(0,2)-mFuncParameters(3,0));
            if( tSqrt <0)
            {
                sign  = -1;
            }
            else if(tSqrt>0)
                {
                    sign = 1;
                }
            else
                MORIS_LOG_ERROR<<"zero denominator detected";

            tFuncVal(0,1)  = mFuncParameters(0,0)/std::sqrt(abs(tSqrt));

            //dz/dr
            tSqrt =    mFuncParameters(0,0)*mFuncParameters(0,0)
                    - (aCoords(0,0)-mFuncParameters(1,0))*(aCoords(0,0)-mFuncParameters(1,0))
                    - (aCoords(0,1)-mFuncParameters(2,0))*(aCoords(0,1)-mFuncParameters(2,0));
            if( tSqrt <0)
            {
                sign  = -1;
            }
            else if(tSqrt>0)
                {
                    sign = 1;
                }
            else
                MORIS_LOG_ERROR<<"zero denominator detected";

            tFuncVal(0,2)  = sign*mFuncParameters(0,0)/std::sqrt(abs(tSqrt));

            // dx/dxc
            tFuncVal(1,0) = 1;
            // dy/dxc
            tFuncVal(1,1) = 0;
            // dz/dxc
            tFuncVal(1,2) = 0;

            // dx/dyc
            tFuncVal(2,0) = 0;
            // dy/dyc
            tFuncVal(2,1) = 1;
            // dz/dyc
            tFuncVal(2,2) = 0;

            // dx/dzc
            tFuncVal(3,0) = 0;
            // dy/dzc
            tFuncVal(3,1) = 0;
            // dz/dzc
            tFuncVal(3,2) = 1;

            break;
            }
        case(FunctionType::LEVELSET_PLANE):
            {
            // x - center = mFuncParameters(0,0)
            // y - center = mFuncParameters(0,1)
            // z - center = mFuncParameters(0,2)
            // x - normal = mFuncParameters(0,3)
            // y - normal = mFuncParameters(0,4)
            // z - normal = mFuncParameters(0,5)

            MORIS_LOG_ERROR<<"Sensitivity Not implemented.";
            //
            break;
            }

        case(FunctionType::LEVELSET_RANDOM):
            {
            tFuncVal(0,0)  = UINT_MAX;
            MORIS_LOG_ERROR<<"Sensitivity not implemented and does not make sense ";
            break;
            }
        default:
        {
            tFuncVal(0,0)  = UINT_MAX;
            MORIS_ASSERT(0 , "Function type not implemented." );
        }
    }

    return tFuncVal;
}
