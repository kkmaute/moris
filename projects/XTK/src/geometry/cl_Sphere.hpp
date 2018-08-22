/*
 * cl_Function_Levelset_Sphere.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOOLS_CL_FUNCTION_LEVELSET_SPHERE_HPP_
#define SRC_TOOLS_CL_FUNCTION_LEVELSET_SPHERE_HPP_

#include <cmath>

#include "linalg/cl_XTK_Matrix.hpp"
#include "geometry/cl_Geometry.hpp"


namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Sphere : public Geometry<Real, Integer, Real_Matrix, Integer_Matrix>
{
public:
    Sphere(Real const & aRadius, Real const & aXCenter, Real const & aYCenter, Real const & aZCenter) :
            mRadius(aRadius),
            mXCenter(aXCenter),
            mYCenter(aYCenter),
            mZCenter(aZCenter)
    {
    }

    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 4;
    }

    Real evaluate_field_value_with_coordinate(Integer const & aRowIndex,
                                              Mat<Real,Real_Matrix> const & aCoordinates) const
    {

        Real tFunctionValue = (aCoordinates(aRowIndex, 0) - mXCenter) * (aCoordinates(aRowIndex, 0) - mXCenter)
                            + (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter)
                            + (aCoordinates(aRowIndex, 2) - mZCenter) * (aCoordinates(aRowIndex, 2) - mZCenter)
                              - (mRadius * mRadius);

        return tFunctionValue;
    }


    Mat<Real,Real_Matrix> evaluate_sensitivity_dphi_dp_with_coordinate( Integer const & aRowIndex,
                                                                        Mat<Real, Real_Matrix> const & aCoordinates) const
    {
        Mat<Real,Real_Matrix> tSensitivityDxDp(4, 3, 0.0);

        Real sign = 0.0;

        // dx/dr
        Real tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter)
                                      - (aCoordinates(aRowIndex, 2) - mZCenter) * (aCoordinates(aRowIndex, 2) - mZCenter);

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
            XTK_ERROR << "zero denominator detected";
        }

        (tSensitivityDxDp)(0, 0) = sign * mRadius / std::sqrt(std::abs(tSqrt));

        //dy/dr
        tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 0) - mXCenter) * (aCoordinates(aRowIndex, 0) - mXCenter)
                                 - (aCoordinates(aRowIndex, 2) - mZCenter) * (aCoordinates(aRowIndex, 2) - mZCenter);

        if(tSqrt < 0.0)
        {
            sign = -1.0;
        }
        else if(tSqrt > 0.0)
        {
            sign = 1.0;
        }
        else
            XTK_ERROR << "zero denominator detected";

        tSensitivityDxDp(0, 1) = mRadius / std::sqrt(std::abs(tSqrt));

        //dz/dr
        tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 0) - mXCenter) * (aCoordinates(aRowIndex, 0) - mXCenter)
                                 - (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter);
        if(tSqrt < 0.0)
        {
            sign = -1.0;
        }
        else if(tSqrt > 0.0)
        {
            sign = 1.0;
        }
        else
            XTK_ERROR << "zero denominator detected";

        tSensitivityDxDp(0, 2) = sign * mRadius / std::sqrt(std::abs(tSqrt));

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

private:
    Real const mRadius;
    Real const mXCenter;
    Real const mYCenter;
    Real const mZCenter;
};
}

#endif /* SRC_TOOLS_CL_FUNCTION_LEVELSET_SPHERE_HPP_ */
