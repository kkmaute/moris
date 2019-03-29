/*
 * cl_GE_Geometry.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_MTK_GE_SRC_CL_GE_BASE_HPP_
#define PROJECTS_MTK_GE_SRC_CL_GE_BASE_HPP_

// GE includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Enums.hpp"
//------------------------------------------------------------------------------
// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"
//------------------------------------------------------------------------------
// other includes
#include <cmath>
//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
	class Geometry{

	public:
		Geometry(){};

		~Geometry(){};
        //------------------------------------------------------------------------------
		virtual
		void set_analytical_function( real ( *mFuncAnalytic )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant ) )
		{
			std::cout<<"set_analytical_function(): please specify your own analytic function"<<std::endl;
		};
		//------------------------------------------------------------------------------
		virtual
        void set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat >        & aCoordinates,
                                                                   Cell<Cell<moris::real>> & aCenter,
                                                                   Cell<moris::real>       & aRadius,
                                                                   Cell<moris::real>       & aLength,
                                                                   Cell<Cell<moris::real>> & aAxis ) )
        {
		    std::cout<<"set_analytical_function(): please choose a valid function"<<std::endl;
        };
        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function( type aGeomType )
        {
            std::cout<<"set_analytical_function(): please choose a valid function"<<std::endl;
        };

        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function_dphi_dx( Matrix< DDRMat > ( *mFuncAnalyticDphiDx )( const Matrix< DDRMat > & aPoint, Cell< real > aConst ) )
        {
            std::cout<<"set_analytical_function_dphi_dx(): please specify your own analytic function dphi/dx"<<std::endl;
        };

        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function_dphi_dx( type aGeomType )
        {
            std::cout<<"set_analytical_function_dphi_dx(): please choose a valid dphi/dx function"<<std::endl;
        };

        //------------------------------------------------------------------------------

        virtual
        real get_field_val_at_coordinate( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConst )
        {
            std::cout<<"get_field_val_at_coordinate(): function not implemented"<<std::endl;
            return 0.0;
        };

        //------------------------------------------------------------------------------
        /* used when the function does not require any constants as inputs ( e.g. gyroid_function() ) */
        virtual
        real get_field_val_at_coordinate( const Matrix< DDRMat > & aPoint )
        {
            std::cout<<"get_field_val_at_coordinate(): function not implemented"<<std::endl;
            return 0.0;
        };

        //------------------------------------------------------------------------------

        virtual
        Matrix< DDRMat >
        get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint,
                                                     moris::Cell< real > aConst)
        {
            Matrix< DDRMat > tSensitivityDxDp(4, 3, 0.0);
            return tSensitivityDxDp;
        };

        //------------------------------------------------------------------------------


	//------------------------------------------------------------------------------
    private:

	//------------------------------------------------------------------------------
    protected:
		uint mNumEle = 0;

	};

//------------------------------------------------------------------------------
/*
 * *****************************************************************************
 * list of standard functions
 * *****************************************************************************
*/
//------------------------------------------------------------------------------
real
sphere_function( const Matrix< DDRMat > & aCoordinate,
                       Cell< real >       aInputs )
{   /* aCoordinate = point vector to determine value at (x,y,z)
     * aInputs(0)  = x location of center;      aInputs(1)  = y location of center
     * aInputs(2)  = z location of center;      aInputs(3)  = radius of sphere */
    real tFuncVal = (aCoordinate(0,0) - aInputs(0))*(aCoordinate(0,0) - aInputs(0)) +
                    (aCoordinate(0,1) - aInputs(1))*(aCoordinate(0,1) - aInputs(1)) +
                    (aCoordinate(0,2) - aInputs(2))*(aCoordinate(0,2) - aInputs(2)) -
                    (aInputs(3)*aInputs(3));
    return tFuncVal;
}
//------------------------------------------------------------------------------
Matrix< DDRMat >
sphere_function_dphi_dx( const Matrix< DDRMat > & aCoordinate,
                               Cell< real>        aInputs)
{
    moris::Matrix< moris::DDRMat > tSensitivityDxDp(4, 3, 0.0);

    moris::real sign = 0.0;

    // dx/dr
    moris::real tSqrt = aInputs(3) * aInputs(3) - (aCoordinate(1) - aInputs(1)) * (aCoordinate(1) - aInputs(1))
                                  - (aCoordinate(2) - aInputs(2)) * (aCoordinate(2) - aInputs(2));

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

    (tSensitivityDxDp)(0, 0) = sign * aInputs(3) / std::sqrt(std::abs(tSqrt));

    //dy/dr
    tSqrt = aInputs(3) * aInputs(3) - (aCoordinate(0) - aInputs(0)) * (aCoordinate(0) - aInputs(0))
                             - (aCoordinate(2) - aInputs(2)) * (aCoordinate(2) - aInputs(2));

    if(tSqrt < 0.0)
    {
        sign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        sign = 1.0;
    }
    else
        std::cout << "zero denominator detected";

    tSensitivityDxDp(0, 1) = aInputs(3) / std::sqrt(std::abs(tSqrt));

    //dz/dr
    tSqrt = aInputs(3) * aInputs(3) - (aCoordinate(0) - aInputs(0)) * (aCoordinate(0) - aInputs(0))
                             - (aCoordinate(1) - aInputs(1)) * (aCoordinate(1) - aInputs(1));
    if(tSqrt < 0.0)
    {
        sign = -1.0;
    }
    else if(tSqrt > 0.0)
    {
        sign = 1.0;
    }
    else
        std::cout << "zero denominator detected";

    tSensitivityDxDp(0, 2) = sign * aInputs(3) / std::sqrt(std::abs(tSqrt));

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
{   /* aInputs(0) = aXc;     aInputs(1) = aYc;
     * aInputs(2) = aZc;     aInputs(3) = aXn;
     * aInputs(4) = aYn;     aInputs(5) = aZn;     */

    real tDist = aInputs(3)*(aCoordinates(0)-aInputs(0)) + aInputs(4)*(aCoordinates(1)-aInputs(1)) + aInputs(5)*(aCoordinates(2)-aInputs(2));
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
        for (moris::size_t k=0;k<tFiberKmax;k++)
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
    Cell< real > tRadiuses        = aRadius;
    Cell< real > tLengths         = aLength;
    Cell< Cell< real > > tAxes    = aAxis;

    real lsVal = getSingleCylLSVal(tCenters(0),tAxes(0),tRadiuses(0),tLengths(0),aCoordinates);

    for (size_t cInd = 1; cInd < tCenters.size(); ++cInd)
    {
        real thisLsVal  =  getSingleCylLSVal(tCenters(cInd),tAxes(cInd),tRadiuses(cInd),tLengths(cInd),aCoordinates);

        lsVal = std::min(thisLsVal, lsVal);
    }

    return lsVal;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_GEN_SRC_CL_GE_BASE_HPP_ */
