/*
 * cl_Analytic_Level_Set_Spiral.hpp
 *
 *  Created on: Aug 18, 2017
 *      Author: ktdoble
 */

#ifndef SRC_GEOMETRY_CL_SPIRAL_HPP_
#define SRC_GEOMETRY_CL_SPIRAL_HPP_


#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "geometry/cl_Geometry.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Spiral : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Spiral(Integer const & aNumSphere  ,
           Real    const & aFiberRadius,
           Real    const & aCoilRadius ,
           Real    const & aNumCoils   ,
           Real    const & aStemLength ,
           Real    const & aPitchLength,
           Real    const & aFiberExpo  ,
           Real    const & aFiberXctr  ,
           Real    const & aFiberYctr  ,
           Real    const & aFiberZctr  ) :
               mNumSpheres  ( aNumSphere ),
               mFiberRadius( aFiberRadius ),
               mCoilRadius ( aCoilRadius ),
               mNumCoils   ( aNumCoils ),
               mStemLength ( aStemLength ),
               mPitchLength( aPitchLength ),
               mFiberExpo  ( aFiberExpo ),
               mFiberXctr  ( aFiberXctr ),
               mFiberYctr  ( aFiberYctr ),
               mFiberZctr  ( aFiberZctr )




    {
    }

    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
    }

    bool is_analytic() const
    {
        return true;
    }

    void get_dx_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 0;
        aNumCols = 0;
    }

    Real evaluate_field_value_with_coordinate(Integer const & aRowIndex,
                                              moris::Mat_New<Real,Real_Matrix> const & aCoordinates) const
    {

        Real MATH_PI = 3.14159265359;
        Real xcoord = aCoordinates(aRowIndex,0);
        Real ycoord = aCoordinates(aRowIndex,1);
        Real zcoord = aCoordinates(aRowIndex,2);

        Real LSval = -1e99;

        Real coilLength  = 2.0*MATH_PI*mCoilRadius*mNumCoils;

        Real totalLength = coilLength+mStemLength+mCoilRadius;

        Real numSpheresCoil = std::ceil(coilLength/totalLength*mNumSpheres)+1;
        Real numSpheresStem = std::ceil(mStemLength/totalLength*mNumSpheres)+1;
        Real numSpheresBar  = std::ceil(mCoilRadius/totalLength*mNumSpheres)+1;

        Real ofrx = std::pow(1.0/mFiberRadius,mFiberExpo);

        Real xci,yci,zci;

        // coils

        for (Integer i=0;i<numSpheresCoil;i++)
        {
            Real tc = 2*MATH_PI*mNumCoils*(i-1)/(numSpheresCoil-1);
            Real tz = (i-1)/(numSpheresCoil-1);

            xci = mFiberXctr + mCoilRadius*std::cos(tc);
            yci = mFiberYctr + mCoilRadius*std::sin(tc);
            zci = mFiberZctr + mPitchLength*tz + mStemLength;

            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }

        // straight section
        for (Integer i=0;i<numSpheresStem;i++)
        {
            Real tz = (i-1)/(numSpheresStem-1);

            xci = mFiberXctr;
            yci = mFiberYctr;
            zci = mFiberZctr + mStemLength*tz;

            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }


        // connector
        for (Integer i=0;i<numSpheresBar;i++)
        {
            Real tz = (i-1)/(numSpheresBar-1);

            xci = mFiberXctr + mCoilRadius*tz;
            yci = mFiberYctr;
            zci = mFiberZctr + mStemLength;

            LSval = std::max(LSval,1.0 - std::pow(  ofrx*std::pow(xcoord-xci, mFiberExpo)
            + ofrx*std::pow(ycoord-yci, mFiberExpo)
            + ofrx*std::pow(zcoord-zci, mFiberExpo), 1.0/mFiberExpo));
        }

        return LSval;
    }


private:
     Integer const mNumSpheres  ;     // number of spheres
     Real    const mFiberRadius    ;    // radius of fibers
     Real    const mCoilRadius     ;     // coil radius
     Real    const mNumCoils       ;     // number of coils
     Real    const mStemLength     ;     // length of lower stem
     Real    const mPitchLength    ;   // hight of all coils
     Real    const mFiberExpo      ;     // cuboid exponent
     Real    const mFiberXctr      ;
     Real    const mFiberYctr      ;
     Real    const mFiberZctr      ;
};
}


#endif /* SRC_GEOMETRY_CL_SPIRAL_HPP_ */
