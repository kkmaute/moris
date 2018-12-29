/*
 * cl_Analytic_Fibers.hpp
 *
 *  Created on: Jul 27, 2017
 *      Author: ktdoble
 */

#include <cmath>

#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"

// XTKL: Matrix Include


namespace xtk
{
class Composite_Fiber: public Geometry
{
public:
    Composite_Fiber(moris::real   const & aFiberRadius,
                    moris::size_t const & aNumSpheres,
                    moris::real   const & aFiberFrq,
                    moris::real   const & aFiberExpo,
                    moris::real   const & aFiberDelX,
                    moris::real   const & aFiberDelY,
                    moris::real   const & aFiberDelZ,
                    moris::real   const & aFiberXctr,
                    moris::real   const & aFiberYctr,
                    moris::real   const & aFiberZctr):
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

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 1;
    }

    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                              moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {

        moris::real MATH_PI = 3.14159265359;
        moris::real xcoord  = aCoordinates(aRowIndex,0);
        moris::real ycoord  = aCoordinates(aRowIndex,1);
        moris::real zcoord  = aCoordinates(aRowIndex,2);

        moris::real LSval = -1e99;

        moris::real ofrx  = std::pow(1.0/mFiberRadius, mFiberExpo);

        moris::real xci,yci,zci;

        for (moris::size_t i=0;i<mNumSpheres;i++)
        {
            // wavy fiber
            for (moris::size_t k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres);
                yci=mFiberYctr+k*6.0*mFiberRadius;
                zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                    + ofrx*std::pow(ycoord-yci, mFiberExpo)
                    + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }

            // straight fibers
            for (moris::size_t k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+2.0/MATH_PI*mFiberAmp-mFiberRadius;
                yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
                zci=mFiberZctr+mFiberDelZ*(2.0*k-0.0)/mFiberFrq;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                    + ofrx*std::pow(ycoord-yci, mFiberExpo)
                    + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
                xci=mFiberXctr-2.0/MATH_PI*mFiberAmp+mFiberRadius;
                yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
                zci=mFiberZctr+mFiberDelZ*((2.0*k-1.0)-0.0)/mFiberFrq;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }

            // wavy fiber
            for (moris::size_t k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres+2.0/3.0*MATH_PI);
                yci=mFiberYctr+2.0*mFiberRadius+k*6.0*mFiberRadius;
                zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }


            // straight fibers
            for (moris::size_t k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+2.0/MATH_PI*mFiberAmp-mFiberRadius;
                yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
                zci=mFiberZctr+mFiberDelZ*(2.0*k-2.0/3.0)/mFiberFrq;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));

                xci=mFiberXctr-2.0/MATH_PI*mFiberAmp+mFiberRadius;
                yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
                zci=mFiberZctr+mFiberDelZ*((2.0*k-1.0)-2.0/3.0)/mFiberFrq;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }

            // wavy fiber
            for (moris::size_t k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+mFiberAmp*std::cos(mFiberFrq*MATH_PI*i/mNumSpheres+4.0/3.0*MATH_PI);
                yci=mFiberYctr+4.0*mFiberRadius+k*6.0*mFiberRadius;
                zci=mFiberZctr+mFiberDelZ*i/mNumSpheres;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }

            // straight fibers
            for (moris::size_t k=0;k<mFiberKmax;k++)
            {
                xci=mFiberXctr+2.0/MATH_PI*mFiberAmp-mFiberRadius;
                yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
                zci=mFiberZctr+mFiberDelZ*(2.0*k-4.0/3.0)/mFiberFrq;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));

                xci=mFiberXctr-2.0/MATH_PI*mFiberAmp+mFiberRadius;
                yci=mFiberYctr+mFiberDelY*i/mNumSpheres;
                zci=mFiberZctr+mFiberDelZ*((2.0*k-1)-4.0/3.0)/mFiberFrq;
                LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
                + ofrx*std::pow(ycoord-yci, mFiberExpo)
                + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
            }
        }


        return LSval;
    }

private:
    // design parameters
    const moris::real mFiberRadius;    // radius of fibers
    const moris::size_t mNumSpheres;     // number of spheres
    const moris::real mFiberFrq;     // frequency of weave
    const moris::real mFiberExpo;     // cuboid exponent

    const moris::real mFiberDelX;
    const moris::real mFiberDelY;
    const moris::real mFiberDelZ;

    const moris::real mFiberAmp;
    const moris::real mFiberKmax;

    const moris::real mFiberXctr;
    const moris::real mFiberYctr;
    const moris::real mFiberZctr;
};
}
