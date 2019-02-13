/*
 * cl_Composite_Fiber_Wave_3.hpp
 *
 *  Created on: Jul 17, 2018
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_COMPOSITE_FIBER_WAVE_3_HPP_
#define SRC_GEOMETRY_CL_COMPOSITE_FIBER_WAVE_3_HPP_

#include <cmath>

#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"

// XTKL: Matrix Include


namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Composite_Fiber_Wave_3: public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Composite_Fiber_Wave_3(Real    const & aFiberRadius,
                             Integer const & aNumSpheres,
                             Real    const & aFiberFrq,
                             Real    const & aFiberExpo,
                             Real    const & aFiberDelX,
                             Real    const & aFiberDelY,
                             Real    const & aFiberDelZ,
                             Real    const & aFiberXctr,
                             Real    const & aFiberYctr,
                             Real    const & aFiberZctr):
                                 mFiberRadius(aFiberRadius),
                                 mNumSpheres(aNumSpheres),
                                 mFiberFrq(aFiberFrq),
                                 mFiberExpo(aFiberExpo),
                                 mFiberDelX(aFiberDelX),
                                 mFiberDelY(aFiberDelY),
                                 mFiberDelZ(aFiberDelZ),
                                 mFiberAmp(aFiberDelX/2.0-1.1*aFiberRadius),
                                 mFiberKmax(std::floor(mFiberDelY/(6.0*mFiberRadius))+1),
                                 mFiberXctr(aFiberXctr),
                                 mFiberYctr(aFiberYctr),
                                 mFiberZctr(aFiberZctr)

    {

    }

    bool is_analytic() const
    {
        return true;
    }

    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 1;
    }

    Real evaluate_field_value_with_coordinate(Integer const & aRowIndex,
                                              moris::Matrix< Real_Matrix > const & aCoordinates) const
    {

        Real MATH_PI = 3.14159265359;
        Real xcoord  = aCoordinates(aRowIndex,0);
        Real ycoord  = aCoordinates(aRowIndex,1);
        Real zcoord  = aCoordinates(aRowIndex,2);

        Real LSval = -1e99;

        Real ofrx  = std::pow(1.0/mFiberRadius, mFiberExpo);

        Real xci,yci,zci;

        for (Integer i=0;i<mNumSpheres;i++)
        {
           // wavy fiber
            for (Integer k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres+4.0/3.0*MATH_PI);
                yci=mFiberYctr+4.0*mFiberRadius+k*6.0*mFiberRadius;
                zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }

        }


        return LSval;
    }


private:
    // design parameters
    const Real mFiberRadius;    // radius of fibers
    const Integer mNumSpheres;     // number of spheres
    const Real mFiberFrq;     // frequency of weave
    const Real mFiberExpo;     // cuboid exponent

    const Real mFiberDelX;
    const Real mFiberDelY;
    const Real mFiberDelZ;

    const Real mFiberAmp;
    const Real mFiberKmax;

    const Real mFiberXctr;
    const Real mFiberYctr;
    const Real mFiberZctr;
};
}



#endif /* SRC_GEOMETRY_CL_COMPOSITE_FIBER_WAVE_3_HPP_ */
