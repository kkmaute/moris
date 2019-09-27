/*
 * cl_Star.hpp
 *
 *  Created on: Aug 7, 2019
 *      Author: ryan
 */

#ifndef SRC_TOOLS_CL_FUNCTION_LEVELSET_STAR_HPP_
#define SRC_TOOLS_CL_FUNCTION_LEVELSET_STAR_HPP_

#include <cmath>

#include "cl_Matrix.hpp"
#include "cl_Geometry.hpp"


namespace xtk
{

class Star : public Geometry
{
public:
    Star(){}

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

        moris::real tPhi = std::atan2( aCoordinates( 0 ), aCoordinates( 1 ) );

        moris::real tFunctionValue = 0.5 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aCoordinates( 0 ), 2 ) + std::pow( aCoordinates( 1 ), 2 ) );

        return tFunctionValue;
    }


    moris::Matrix< moris::DDRMat > evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                                        moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        moris::Matrix< moris::DDRMat > tSensitivityDxDp(4, 3, 0.0);
/*
        moris::real sign = 0.0;

        // dx/dr
        moris::real tSqrt = mRadius * mRadius - (aCoordinates(aRowIndex, 1) - mYCenter) * (aCoordinates(aRowIndex, 1) - mYCenter)
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
            std::cout << "zero denominator detected";
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
            std::cout << "zero denominator detected";

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
            std::cout << "zero denominator detected";

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
*/
        return tSensitivityDxDp;


    }

private:
    moris::real mRadius;
    moris::real mXCenter;
    moris::real mYCenter;
};

}

#endif /* SRC_TOOLS_CL_FUNCTION_LEVELSET_CSTAR_HPP_ */


