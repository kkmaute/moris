/*
 * cl_Interpolation.cpp
 *
 *  Created on: Mar 20, 2017
 *      Author: doble
 */

#include "cl_Interpolation.hpp"

moris::Mat<moris::real>
moris::Interpolation::linear_interpolation_location(const moris::Mat<moris::real>  & aInterpVars,
                                                    const moris::Mat<moris::real>  & aLocation)
{
    moris::real xi        = aLocation(0,0);
    moris::uint tNumInt = aInterpVars.n_cols();
    moris::Mat<moris::real> tInterpResult(1,tNumInt,UINT_MAX);
    moris::Mat<moris::real> tTmpVar(1,1,UINT_MAX);
    for (moris::uint  i = 0; i<tNumInt;i++)
    {
        tTmpVar = aInterpVars.col(i);

        tInterpResult(0,i) = (tTmpVar(0,0)*(1-xi)+
                              tTmpVar(1,0)*(1+xi))/2;
    }
    return tInterpResult;
}

moris::Mat<moris::real>
moris::Interpolation::linear_interpolation_value(const moris::Mat<moris::real>  & aInterpVars,
                                                 const moris::Mat<moris::real>  & aValue)
{
    moris::real val        = aValue(0,0);
    moris::uint tNumInt = aInterpVars.n_cols();
    moris::Mat<moris::real> tInterpResult(1,tNumInt,UINT_MAX);
    moris::Mat<moris::real> tTmpLoc(1,1,UINT_MAX);
    for (moris::uint  i = 0; i<tNumInt;i++)
    {
        tTmpLoc = aInterpVars.col(i);

        tInterpResult(0,i) = (2*val-tTmpLoc(1,0)-tTmpLoc(0,0))/(tTmpLoc(1,0)-tTmpLoc(0,0));
    }
    return tInterpResult;
}

moris::Mat<moris::real>
moris::Interpolation::bilinear_interpolation(const moris::Mat<moris::real>  & aInterpVars,
                                             const moris::Mat<moris::real>    & aValue)
{
    moris::real xi        = aValue(0,0);
    moris::real eta       = aValue(0,1);
    moris::uint tNumInt = aInterpVars.n_cols();
    moris::Mat<moris::real> tInterpResult(1,tNumInt,UINT_MAX);
    moris::Mat<moris::real> tTmpVar(1,1,UINT_MAX);
    for (moris::uint  i = 0; i<tNumInt;i++)
    {
        tTmpVar = aInterpVars.col(i);

        tInterpResult(0,i) = (tTmpVar(0,0)*(1-xi)*(1-eta)+
                              tTmpVar(1,0)*(1+xi)*(1-eta)+
                              tTmpVar(2,0)*(1+xi)*(1+eta)+
                              tTmpVar(3,0)*(1-xi)*(1+eta))/4;
    }
    return tInterpResult;
}

moris::Mat<moris::real>
moris::Interpolation::trilinear_interpolation(const moris::Mat<moris::real>  & aInterpVars,
                                              const moris::Mat<moris::real>  & aValue)
{
    moris::real xi          = aValue(0,0);
    moris::real eta         = aValue(0,1);
    moris::real mu          = aValue(0,2);
    moris::uint tNumInt     = aInterpVars.n_cols();
    moris::Mat<moris::real>   tInterpResult(1,tNumInt,UINT_MAX);
    moris::Mat<moris::real>   tTmpVar(1,1,UINT_MAX);

    for (moris::uint  i = 0; i<tNumInt;i++)
    {
        tTmpVar = aInterpVars.col(i);

        tInterpResult(0,i) = (tTmpVar(0,0)*(1-xi)*(1-eta)*(1-mu)+
                              tTmpVar(1,0)*(1+xi)*(1-eta)*(1-mu)+
                              tTmpVar(2,0)*(1+xi)*(1+eta)*(1+mu)+
                              tTmpVar(3,0)*(1-xi)*(1+eta)*(1-mu)+
                              tTmpVar(4,0)*(1-xi)*(1-eta)*(1+mu)+
                              tTmpVar(5,0)*(1+xi)*(1-eta)*(1+mu)+
                              tTmpVar(6,0)*(1+xi)*(1+eta)*(1+mu)+
                              tTmpVar(7,0)*(1-xi)*(1+eta)*(1+mu))/8;
    }
    return tInterpResult;
}
