/*
 * cl_GE_Geometry_Library.hpp
 *
 *  Created on: Apr 5, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_GEOMETRY_LIBRARY_HPP_
#define PROJECTS_GEN_SRC_CL_GE_GEOMETRY_LIBRARY_HPP_

#include "HDF5_Tools.hpp"

#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris{
namespace ge{
//------------------------------------------------------------------------------
/*
 * *****************************************************************************
 * list of standard functions
 * *****************************************************************************
*/
//------------------------------------------------------------------------------
real
circle_function( const Matrix< DDRMat > & aPoint,
                       Cell< real >       aInputs )
{
    // aInputs(0) = radius
    // aInputs(1) = x location of center
    // aInputs(2) = y location of center

    Matrix< DDRMat > tCenter(1,2);
    tCenter(0) = aInputs(1); tCenter(1) = aInputs(2);

    return norm(aPoint-tCenter) - aInputs(0);   // using this form linearizes the circle equation
//    return aInputs(0) - norm(aPoint-tCenter);   // using this form linearizes the circle equation
//    return std::pow((aPoint(0,0) - aInputs(0)),2) + std::pow((aPoint(0,1) - aInputs(1)),2) - std::pow(aInputs(2),1);  // using this form results in a non-linear circle equaiton
}
//------------------------------------------------------------------------------
Matrix< DDRMat >
circle_function_dphi_dp( const Matrix< DDRMat > & aPoint,
                               Cell< real >       aInputs )
{
    real tX = aPoint(0);
    real tY = aPoint(1);

    real tR  = aInputs(0);
    real tXc = aInputs(1);
    real tYc = aInputs(2);

    moris::Matrix< DDRMat > tSensitivity( 3, 2, 0.0 );
    moris::real tSign = 0.0;

    // dx/dr = r/sqrt(r^2 - (y-yc)^2)
    moris::real tSqrt = tR*tR - std::pow((tY-tYc),2);
    if(tSqrt < 0.0)
    {
        tSign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        tSign = 1.0;
    }
    else
    {
        std::cout << "zero denominator detected";
    }
    tSensitivity(0,0) = tSign * tR / std::sqrt(std::abs(tSqrt));

    // dy/dr = r/sqrt(r^2 - (x-xc)^2)
    tSqrt = tR*tR - std::pow((tX-tXc),2);
    if(tSqrt < 0.0)
    {
        tSign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        tSign = 1.0;
    }
    else
    {
        std::cout << "zero denominator detected";
    }
    tSensitivity(0,1) = tSign * tR / std::sqrt(std::abs(tSqrt));

    // fill remaining values in tSensitivity
    tSensitivity(1,0) = 1.0; // dx/dxc
    tSensitivity(1,1) = 0.0; // dy/dxc
    tSensitivity(2,0) = 0.0; // dx/dyc
    tSensitivity(2,1) = 1.0; // dy/dyc

    return tSensitivity;
}
//------------------------------------------------------------------------------
real
sphere_function( const Matrix< DDRMat > & aCoordinate,
                       Cell< real >       aInputs )
{   /* aCoordinate = point vector to determine value at (x,y,z)
     * aInputs(0)  = radius of sphere;
     * aInputs(1)  = x location of center;
     * aInputs(2)  = y location of center;
     * aInputs(3)  = z location of center */
    Matrix< DDRMat > tCenterVec(1,3);   // (x,y,z)
    tCenterVec(0,0) = aInputs(1);
    tCenterVec(0,1) = aInputs(2);
    tCenterVec(0,2) = aInputs(3);
    return norm( aCoordinate - tCenterVec ) - aInputs(0);     //linearized form of the sphere equation

//    real tFuncVal = (aCoordinate(0,0) - aInputs(0))*(aCoordinate(0,0) - aInputs(0)) +
//                    (aCoordinate(0,1) - aInputs(1))*(aCoordinate(0,1) - aInputs(1)) +
//                    (aCoordinate(0,2) - aInputs(2))*(aCoordinate(0,2) - aInputs(2)) -
//                    (aInputs(3)*aInputs(3));
//    return tFuncVal;
}
//------------------------------------------------------------------------------
Matrix< DDRMat >
sphere_function_dphi_dp( const Matrix< DDRMat > & aCoordinate,
                               Cell< real>        aInputs)
{
    real tX = aCoordinate(0);
    real tY = aCoordinate(1);
    real tZ = aCoordinate(2);

    real tR  = aInputs(0);
    real tXc = aInputs(1);
    real tYc = aInputs(2);
    real tZc = aInputs(3);
    moris::Matrix< moris::DDRMat > tSensitivityDxDp(4, 3, 0.0);

    moris::real sign = 0.0;

    // dx/dr
    moris::real tSqrt = tR * tR - (tY - tYc) * (tY - tYc)
                                - (tZ - tZc) * (tZ - tZc);

    if(tSqrt < 0.0)
    {
        sign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        sign = 1.0;
    }
    else
    {
        std::cout << "zero denominator detected";
    }

    tSensitivityDxDp(0, 0) = sign * tR/ std::sqrt(std::abs(tSqrt));

    //dy/dr
    tSqrt = tR * tR - (tX - tXc) * (tX - tXc)
                    - (tZ - tZc) * (tZ - tZc);

    if(tSqrt < 0.0)
    {
        sign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        sign = 1.0;
    }
    else
    {
        std::cout << "zero denominator detected";
    }

    tSensitivityDxDp(0, 1) = tR / std::sqrt(std::abs(tSqrt));

    //dz/dr
    tSqrt = tR*tR - (tX - tXc) * (tX - tXc)
                  - (tY - tYc) * (tY - tYc);
    if(tSqrt < 0.0)
    {
        sign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        sign = 1.0;
    }
    else
    {
        std::cout << "zero denominator detected";
    }
    tSensitivityDxDp(0, 2) = sign*tR / std::sqrt(std::abs(tSqrt));

    tSensitivityDxDp(1, 0) = 1.0; // dx/dxc
    tSensitivityDxDp(1, 1) = 0.0; // dy/dxc
    tSensitivityDxDp(1, 2) = 0.0; // dz/dxc
    tSensitivityDxDp(2, 0) = 0.0; // dx/dyc
    tSensitivityDxDp(2, 1) = 1.0; // dy/dyc
    tSensitivityDxDp(2, 2) = 0.0; // dz/dyc
    tSensitivityDxDp(3, 0) = 0.0; // dx/dzc
    tSensitivityDxDp(3, 1) = 0.0; // dy/dzc
    tSensitivityDxDp(3, 2) = 1.0; // dz/dzc

    return tSensitivityDxDp;
}
//------------------------------------------------------------------------------
real
plane_function( const Matrix< DDRMat > & aCoordinates,
                      Cell< real>        aInputs)
{
    real tXc = aInputs(0);
    real tYc = aInputs(1);
    real tZc = aInputs(2);
    real tXn = aInputs(3);
    real tYn = aInputs(4);
    real tZn = aInputs(5);

    real tDist = tXn*(aCoordinates(0)-tXc) + tYn*(aCoordinates(1)-tYc) + tZn*(aCoordinates(2)-tZc);
    // moris::real tDist = mXn*(aCoordinates(aRowIndex,0)-mXc) + mYn*(aCoordinates(aRowIndex,1)-mYc) + mZn*(aCoordinates(aRowIndex,2)-mZc);
    return tDist;
}
//------------------------------------------------------------------------------
real
spiral_function( const Matrix< DDRMat > & aCoordinates,
                       Cell< real>        aInputs )
{
    uint const & tNumSpheres  = (uint)aInputs(0);
    real const & tFiberRadius = aInputs(1);
    real const & tCoilRadius  = aInputs(2);
    real const & tNumCoils    = aInputs(3);
    real const & tStemLength  = aInputs(4);
    real const & tPitchLength = aInputs(5);
    real const & tFiberExpo   = aInputs(6);
    real const & tFiberXctr   = aInputs(7);
    real const & tFiberYctr   = aInputs(8);
    real const & tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real coilLength  = 2.0*MATH_PI*tCoilRadius*tNumCoils;

    real totalLength = coilLength+tStemLength+tCoilRadius;

    real numSpheresCoil = std::ceil(coilLength/totalLength*tNumSpheres)+1;
    real numSpheresStem = std::ceil(tStemLength/totalLength*tNumSpheres)+1;
    real numSpheresBar  = std::ceil(tCoilRadius/totalLength*tNumSpheres)+1;

    real ofrx = std::pow(1.0/tFiberRadius,tFiberExpo);

    real xci,yci,zci;

    // coils

    for (uint i=0;i<numSpheresCoil;i++)
    {
        real tc = 2*MATH_PI*tNumCoils*(i-1)/(numSpheresCoil-1);
        real tz = (i-1)/(numSpheresCoil-1);

        xci = tFiberXctr + tCoilRadius*std::cos(tc);
        yci = tFiberYctr + tCoilRadius*std::sin(tc);
        zci = tFiberZctr + tPitchLength*tz + tStemLength;

        LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
        + ofrx*std::pow(ycoord-yci, tFiberExpo)
        + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
    }

    // straight section
    for (uint i=0;i<numSpheresStem;i++)
    {
        real tz = (i-1)/(numSpheresStem-1);

        xci = tFiberXctr;
        yci = tFiberYctr;
        zci = tFiberZctr + tStemLength*tz;

        LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                + ofrx*std::pow(ycoord-yci, tFiberExpo)
                + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
    }


    // connector
    for (uint i=0;i<numSpheresBar;i++)
    {
        real tz = (i-1)/(numSpheresBar-1);

        xci = tFiberXctr + tCoilRadius*tz;
        yci = tFiberYctr;
        zci = tFiberZctr + tStemLength;

        LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                + ofrx*std::pow(ycoord-yci, tFiberExpo)
                + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
    }

    return LSval;
}
//------------------------------------------------------------------------------
real
gyroid_function(const Matrix< DDRMat > & aCoordinates,
                      Cell< real>        aInputs  =  0.0)
{
    real tFuncVal = std::sin(aCoordinates(0))*std::cos(aCoordinates(1))+
                    std::sin(aCoordinates(1))*std::cos(aCoordinates(2))+
                    std::sin(aCoordinates(2))*std::cos(aCoordinates(0));
    return tFuncVal;
}
//------------------------------------------------------------------------------
real
composite_fiber_function( const Matrix< DDRMat > & aCoordinates,
                                Cell< real>        aInputs )
{
    real tFiberRadius         = aInputs(0);
    moris::size_t tNumSpheres = (moris::size_t)aInputs(1);
    real tFiberFrq            = aInputs(2);
    real tFiberExpo           = aInputs(3);
    real tFiberDelX           = aInputs(4);
    real tFiberDelY           = aInputs(5);
    real tFiberDelZ           = aInputs(6);
    real tFiberAmp            = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax           = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr           = aInputs(7);
    real tFiberYctr           = aInputs(8);
    real tFiberZctr           = aInputs(9);

    moris::real MATH_PI = 3.14159265359;
    moris::real xcoord  = aCoordinates(0);
    moris::real ycoord  = aCoordinates(1);
    moris::real zcoord  = aCoordinates(2);

    moris::real LSval = -1e99;

    moris::real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    moris::real xci,yci,zci;

    for (moris::size_t i=0;i<tNumSpheres;i++)
    {
        // wavy fiber
        for (moris::size_t k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+tFiberAmp*std::cos(tFiberFrq*MATH_PI*i/tNumSpheres);
            yci   = tFiberYctr+k*6.0*tFiberRadius;
            zci   = tFiberZctr+tFiberDelZ*i/tNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

        // straight fibers
        for (moris::size_t k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+2.0/MATH_PI*tFiberAmp-tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*(2.0*k-0.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
            xci   = tFiberXctr-2.0/MATH_PI*tFiberAmp+tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*((2.0*k-1.0)-0.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

        // wavy fiber
        for (moris::size_t k=0;k<tFiberKmax;k++)//------------------------------------------------------------------------------
        {
            xci   = tFiberXctr+tFiberAmp*std::cos(tFiberFrq*MATH_PI*i/tNumSpheres+2.0/3.0*MATH_PI);
            yci   = tFiberYctr+2.0*tFiberRadius+k*6.0*tFiberRadius;
            zci   = tFiberZctr+tFiberDelZ*i/tNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }


        // straight fibers
        for (moris::size_t k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+2.0/MATH_PI*tFiberAmp-tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*(2.0*k-2.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));

            xci   = tFiberXctr-2.0/MATH_PI*tFiberAmp+tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*((2.0*k-1.0)-2.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

        // wavy fiber
        for (moris::size_t k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+tFiberAmp*std::cos(tFiberFrq*MATH_PI*i/tNumSpheres+4.0/3.0*MATH_PI);
            yci   = tFiberYctr+4.0*tFiberRadius+k*6.0*tFiberRadius;
            zci   = tFiberZctr+tFiberDelZ*i/tNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

        // straight fibers
        for (moris::size_t k=0;k<tFiberKmax;k++)
        {
            xci=tFiberXctr+2.0/MATH_PI*tFiberAmp-tFiberRadius;
            yci=tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci=tFiberZctr+tFiberDelZ*(2.0*k-4.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));

            xci=tFiberXctr-2.0/MATH_PI*tFiberAmp+tFiberRadius;
            yci=tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci=tFiberZctr+tFiberDelZ*((2.0*k-1)-4.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }
    }

    return LSval;
}
//------------------------------------------------------------------------------
real
composite_fiber_wave_1_function( const Matrix< DDRMat > & aCoordinates,
                                       Cell< real>        aInputs )
{
    real tFiberRadius = aInputs(0);
    uint tNumSpheres  = (uint)aInputs(1);
    real tFiberFrq    = aInputs(2);
    real tFiberExpo   = aInputs(3);
    real tFiberDelX   = aInputs(4);
    real tFiberDelY   = aInputs(5);
    real tFiberDelZ   = aInputs(6);
    real tFiberAmp    = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax   = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr   = aInputs(7);
    real tFiberYctr   = aInputs(8);
    real tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<tNumSpheres;i++)
    {
        // wavy fiber
        for (uint k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+tFiberAmp*std::cos(tFiberFrq*MATH_PI*i/tNumSpheres);
            yci   = tFiberYctr+k*6.0*tFiberRadius;
            zci   = tFiberZctr+tFiberDelZ*i/tNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }
    }


    return LSval;
}
//------------------------------------------------------------------------------
real
composite_fiber_wave_2_function( const Matrix< DDRMat > & aCoordinates,
                                       Cell< real>        aInputs )
{
    real tFiberRadius = aInputs(0);
    uint tNumSpheres  = (uint)aInputs(1);
    real tFiberFrq    = aInputs(2);
    real tFiberExpo   = aInputs(3);
    real tFiberDelX   = aInputs(4);
    real tFiberDelY   = aInputs(5);
    real tFiberDelZ   = aInputs(6);
    real tFiberAmp    = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax   = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr   = aInputs(7);
    real tFiberYctr   = aInputs(8);
    real tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<tNumSpheres;i++)
    {
        // wavy fiber
        for (uint k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+tFiberAmp*std::cos(tFiberFrq*MATH_PI*i/tNumSpheres+2.0/3.0*MATH_PI);
            yci   = tFiberYctr+2.0*tFiberRadius+k*6.0*tFiberRadius;
            zci   = tFiberZctr+tFiberDelZ*i/tNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

    }

    return LSval;
}
//------------------------------------------------------------------------------
real
composite_fiber_wave_3_function( const Matrix< DDRMat > & aCoordinates,
                                       Cell< real>        aInputs )
{
    real tFiberRadius = aInputs(0);
    uint tNumSpheres  = (uint)aInputs(1);
    real tFiberFrq    = aInputs(2);
    real tFiberExpo   = aInputs(3);
    real tFiberDelX   = aInputs(4);
    real tFiberDelY   = aInputs(5);
    real tFiberDelZ   = aInputs(6);
    real tFiberAmp    = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax   = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr   = aInputs(7);
    real tFiberYctr   = aInputs(8);
    real tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<tNumSpheres;i++)
    {
        // wavy fiber
        for (uint k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+tFiberAmp*std::cos(tFiberFrq*MATH_PI*i/tNumSpheres+4.0/3.0*MATH_PI);
            yci   = tFiberYctr+4.0*tFiberRadius+k*6.0*tFiberRadius;
            zci   = tFiberZctr+tFiberDelZ*i/tNumSpheres;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

    }

    return LSval;
}
//------------------------------------------------------------------------------
real
composite_fiber_straight_1_function( const Matrix< DDRMat > & aCoordinates,
                                           Cell< real>        aInputs )
{
    real tFiberRadius = aInputs(0);
    uint tNumSpheres  = (uint)aInputs(1);
    real tFiberFrq    = aInputs(2);
    real tFiberExpo   = aInputs(3);
    real tFiberDelX   = aInputs(4);
    real tFiberDelY   = aInputs(5);
    real tFiberDelZ   = aInputs(6);
    real tFiberAmp    = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax   = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr   = aInputs(7);
    real tFiberYctr   = aInputs(8);
    real tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<tNumSpheres;i++)
    {
        // straight fibers
        for (uint k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+2.0/MATH_PI*tFiberAmp-tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*(2.0*k-0.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
            xci   = tFiberXctr-2.0/MATH_PI*tFiberAmp+tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*((2.0*k-1.0)-0.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }
    }

    return LSval;
}
//------------------------------------------------------------------------------
real
composite_fiber_straight_2_function( const Matrix< DDRMat > & aCoordinates,
                                           Cell< real>        aInputs )
{
    real tFiberRadius = aInputs(0);
    uint tNumSpheres  = (uint)aInputs(1);
    real tFiberFrq    = aInputs(2);
    real tFiberExpo   = aInputs(3);
    real tFiberDelX   = aInputs(4);
    real tFiberDelY   = aInputs(5);
    real tFiberDelZ   = aInputs(6);
    real tFiberAmp    = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax   = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr   = aInputs(7);
    real tFiberYctr   = aInputs(8);
    real tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<tNumSpheres;i++)
    {

        // straight fibers
        for (uint k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+2.0/MATH_PI*tFiberAmp-tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*(2.0*k-2.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));

            xci   = tFiberXctr-2.0/MATH_PI*tFiberAmp+tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*((2.0*k-1.0)-2.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }

    }


    return LSval;

}
//------------------------------------------------------------------------------
real
composite_fiber_straight_3_function( const Matrix< DDRMat > & aCoordinates,
                                           Cell< real>        aInputs )
{
    real tFiberRadius = aInputs(0);
    uint tNumSpheres  = (uint)aInputs(1);
    real tFiberFrq    = aInputs(2);
    real tFiberExpo   = aInputs(3);
    real tFiberDelX   = aInputs(4);
    real tFiberDelY   = aInputs(5);
    real tFiberDelZ   = aInputs(6);
    real tFiberAmp    = tFiberDelX/2.0 - 1.1*tFiberRadius;
    real tFiberKmax   = std::floor(tFiberDelY/(6.0*tFiberRadius))+1;
    real tFiberXctr   = aInputs(7);
    real tFiberYctr   = aInputs(8);
    real tFiberZctr   = aInputs(9);

    real MATH_PI = 3.14159265359;
    real xcoord  = aCoordinates(0);
    real ycoord  = aCoordinates(1);
    real zcoord  = aCoordinates(2);

    real LSval = -1e99;

    real ofrx  = std::pow(1.0/tFiberRadius, tFiberExpo);

    real xci,yci,zci;

    for (uint i=0;i<tNumSpheres;i++)
    {
       // straight fibers
        for (uint k=0;k<tFiberKmax;k++)
        {
            xci   = tFiberXctr+2.0/MATH_PI*tFiberAmp-tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*(2.0*k-4.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));

            xci   = tFiberXctr-2.0/MATH_PI*tFiberAmp+tFiberRadius;
            yci   = tFiberYctr+tFiberDelY*i/tNumSpheres;
            zci   = tFiberZctr+tFiberDelZ*((2.0*k-1)-4.0/3.0)/tFiberFrq;
            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, tFiberExpo)
                    + ofrx*std::pow(ycoord-yci, tFiberExpo)
                    + ofrx*std::pow(zcoord-zci, tFiberExpo), 1.0/tFiberExpo));
        }
    }

    return LSval;
}
//------------------------------------------------------------------------------
// multi_cylinder_function() uses this (below)
real
getSingleCylLSVal(Cell<moris::real> const              & aCenter,
                  Cell<moris::real> const              & aAxis,
                  moris::real const                    & aRad,
                  moris::real const                    & aLength,
                  moris::Matrix< moris::DDRMat > const & aPointPosition)
   {
       MORIS_ASSERT(aCenter.size() == 3,"Centers need to have length 3");
       MORIS_ASSERT(aAxis.size() == 3, "axis need to have length 3");
       MORIS_ASSERT(aPointPosition.n_cols() == 3, "pointPosition need to have length 3");

       Cell<moris::real> relativePosition = {(aPointPosition(0) - aCenter(0)),(aPointPosition(1) - aCenter(1)),(aPointPosition(2) - aCenter(2))};
       moris::real lsFromLeft             = (relativePosition(0)*(-aAxis(0)) + relativePosition(1)*(-aAxis(1))+ relativePosition(2)*(-aAxis(2))) - aLength/2.0;
       moris::real lsFromRight            = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2))) - aLength/2.0;

       moris::real axialCrd     = (relativePosition(0)*(aAxis(0)) + relativePosition(1)*(aAxis(1))+ relativePosition(2)*(aAxis(2)));
       Cell<moris::real> radDir = {(relativePosition(0) - aAxis(0)*axialCrd), (relativePosition(1) - aAxis(1)*axialCrd),(relativePosition(2) - aAxis(2)*axialCrd)};
       moris::real radDist      = std::pow(radDir(0)*radDir(0)+radDir(1)*radDir(1)+radDir(2)*radDir(2), 0.5);
       moris::real lsFromRad    = radDist - aRad;

       return std::max(std::max(lsFromLeft, lsFromRight), lsFromRad);
   }
//------------------------------------------------------------------------------
real
multi_cylinder_function( const Matrix< DDRMat >        & aCoordinates,
                               Cell<Cell<moris::real>> & aCenter,
                               Cell<moris::real>       & aRadius,
                               Cell<moris::real>       & aLength,
                               Cell<Cell<moris::real>> & aAxis  )
{
    Cell< Cell< real > > tCenters = aCenter;
    Cell< real > tRadii           = aRadius;
    Cell< real > tLengths         = aLength;
    Cell< Cell< real > > tAxes    = aAxis;

    real lsVal = getSingleCylLSVal(tCenters(0),tAxes(0),tRadii(0),tLengths(0),aCoordinates);

    for (size_t cInd = 1; cInd < tCenters.size(); ++cInd)
    {
        real thisLsVal  =  getSingleCylLSVal(tCenters(cInd),tAxes(cInd),tRadii(cInd),tLengths(cInd),aCoordinates);

        lsVal = std::min(thisLsVal, lsVal);
    }

    return lsVal;
}
//------------------------------------------------------------------------------
//----------- functions for the generation of a toroid structure ---------------
//------------------------------------------------------------------------------
real CSGUnion(real tA, real tB)
{
  return std::min(tA, tB);
}
//------------------------------------------------------------------------------
real CSGUnionRound(real tA, real tB, real r)
{
  auto vc0 = r - tA;
  auto vc1 = r - tB;
  auto u0 = std::max(vc0, 0.0);
  auto u1 = std::max(vc1, 0.0);
  auto len = std::sqrt(std::pow((u0),2) + std::pow((u1),2));

  return std::max(r, std::min(tA, tB)) - len;
}
//------------------------------------------------------------------------------
real CSGIntersect(real tA, real tB)
{
  return std::max(tA, tB);
}
//------------------------------------------------------------------------------
real CSGSubtract(real tA, real tB)
{
  return CSGIntersect(tA, -tB);
}
//------------------------------------------------------------------------------
real CSGIntersectionRound(real tA, real tB, real r)
{
  auto vc0 = r + tA;
  auto vc1 = r + tB;
  auto u0 = std::max(vc0, 0.0);
  auto u1 = std::max(vc1, 0.0);
  auto lenU = std::sqrt(pow((u0),2) + pow((u1),2));

  return std::min(-r, std::max(tA, tB)) + lenU;
}
//------------------------------------------------------------------------------
real offset_tor(real t, real r)
{
  return t - r;
}
//------------------------------------------------------------------------------
real clearance_tor(real tA, real tB, real r)
{
  return CSGSubtract(tA, offset_tor(tB, r));
}
//------------------------------------------------------------------------------
real shell_tor(real t, real r)
{
  return clearance_tor(t, t, r);
}
//------------------------------------------------------------------------------
real box_tor(const moris::Matrix< moris::DDRMat > & pt, const std::vector<real> & lower, const std::vector<real> & upper)
{
  return std::max(std::max(
          std::max(lower[0] - pt(0),
                pt(0) - upper[0]),
          std::max(lower[1] - pt(0),
                pt(1) - upper[1])),
          std::max(lower[2] - pt(2),
                pt(2) - upper[2])   );
}
//------------------------------------------------------------------------------
real sphere_tor(const moris::Matrix< moris::DDRMat > & aPoint, real r, const moris::Matrix< moris::DDRMat > & center)
{
  return std::sqrt(std::pow((aPoint(0) - center(0)),2) +

          std::pow((aPoint(1) - center(1)),2) +

          std::pow((aPoint(2) - center(2)),2)) - r;
}
//------------------------------------------------------------------------------
real getDistanceToGyroids_SphereBox(moris::Matrix< moris::DDRMat > pt)
{
real scale =.2; // scale of the gyroid period
real thickness = .1; // thickness of the gyroids
auto gyroidSrf =
    sin(pt(0) / scale) * cos(pt(1) / scale) +
    sin(pt(1) / scale) * cos(pt(2) / scale) +
    sin(pt(2) / scale) * cos(pt(0) / scale);
  auto gyroid = shell_tor(gyroidSrf, -thickness);

  moris::Matrix<DDRMat> tTemp(1,3);
  tTemp(0) = 2.0; tTemp(1) = 2.0; tTemp(2) = 1.0;

  auto sphere1 = sphere_tor(pt, 3.0, tTemp);

  auto boxB = box_tor(pt, { -4.0,-4.0,-1.5 }, { 2.0,2.0,1.5 });

  auto boxC = box_tor(pt, { -4.125,-4.125,-1.75 }, { -3.875,3.,1.75 });

  auto boxD = box_tor(pt, { -4.125,-4.125,-1.75 }, { 3.0,-3.875,1.75 });

  auto negBox = CSGUnion(boxC, boxD);

  auto blnBoxSphere = CSGUnionRound(boxB, sphere1, .5);

  auto shelledBoxSphere = shell_tor(blnBoxSphere, .12);

  auto shelledBox = CSGSubtract(shelledBoxSphere, negBox);

  auto sphereGyroid = CSGIntersectionRound(blnBoxSphere, gyroid, .1);

  sphereGyroid = CSGSubtract(sphereGyroid, negBox);

  sphereGyroid = CSGSubtract(sphereGyroid, shelledBox);

  auto gyroidWithShell = CSGUnionRound(shelledBox, sphereGyroid, .1);

  gyroidWithShell = CSGIntersect(gyroidWithShell, blnBoxSphere);

  gyroidWithShell = CSGSubtract(gyroidWithShell, negBox);

  return CSGUnion(gyroidWithShell, shelledBox);
}
//------------------------------------------------------------------------------
// --- main function for the gyroid ---
//real getDistanceToGyroidsMassive( const moris::Matrix< moris::DDRMat > & aPoint, const moris::Cell< moris::real > aConst )
real getDistanceToGyroidsMassive( const moris::Matrix< moris::DDRMat > & aPoint )
{
//real scale = 1.0/2.0/M_PI; // scale of the gyroid period
    real scale = 1;
real thickness = 0.75; // thickness of the gyroids
auto gyroidSrf =
    sin(aPoint(0) / scale) * cos(aPoint(1) / scale) +
    sin(aPoint(1) / scale) * cos(aPoint(2) / scale) +
    sin(aPoint(2) / scale) * cos(aPoint(0) / scale);


auto gyroid = shell_tor(gyroidSrf, -thickness);

  moris::Matrix<DDRMat> ttTemp(1,3, 0.0);

  auto sphere1 = sphere_tor(aPoint, 5, ttTemp);

  return CSGIntersectionRound(gyroid, sphere1, 0.75);

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
} // end ge namespace
} // end moris namespace


#endif /* PROJECTS_GEN_SRC_CL_GE_GEOMETRY_LIBRARY_HPP_ */
