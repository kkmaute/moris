/*
 * cl_Circle.hpp
 *
 *  Created on: Aug 7, 2019
 *      Author: ryan
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CIRCLE_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CIRCLE_HPP_

#include <cmath>

#include "cl_Matrix.hpp"
#include "../geometry/cl_GEN_Geometry.hpp"

namespace moris
{
namespace ge
{

class Circle : public GEN_Geometry
{
public:
    Circle(){}

    Circle(moris::real const & aRadius,
           moris::real const & aXCenter,
           moris::real const & aYCenter) :
            mRadius(aRadius),
            mXCenter(aXCenter),
            mYCenter(aYCenter)
    {

    }


    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 4;
    }

    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                              moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {

        moris::real tFunctionValue = (aCoordinates(aRowIndex, 0) - mXCenter) * (aCoordinates(aRowIndex, 0) - mXCenter)
                            + (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter)
                              - (mRadius * mRadius);

        return tFunctionValue;
    }


    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                                        moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        moris::Matrix< moris::DDRMat > tSensitivityDxDp(3, 2, 0.0);

        moris::real sign = 0.0;

        // dx/dr
        moris::real tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter);

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

        (tSensitivityDxDp)(0, 0) = sign * mRadius / std::sqrt(std::abs(tSqrt));

        //dy/dr
        tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 0) - mXCenter) * (aCoordinates(aRowIndex, 0) - mXCenter);

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

        tSensitivityDxDp(0, 1) = mRadius / std::sqrt(std::abs(tSqrt));

        // fill remaining values in tSensitivity
        tSensitivityDxDp(1,0) = 1.0; // dx/dxc
        tSensitivityDxDp(1,1) = 0.0; // dy/dxc
        tSensitivityDxDp(2,0) = 0.0; // dx/dyc
        tSensitivityDxDp(2,1) = 1.0; // dy/dyc

        return tSensitivityDxDp;
    }

private:
    moris::real mRadius;
    moris::real mXCenter;
    moris::real mYCenter;
};

}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CIRCLE_HPP_ */


