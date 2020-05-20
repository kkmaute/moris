/*
 * cl_GEN_Sphere.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_SPHERE_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_SPHERE_HPP_

#include <cmath>

#include "cl_GEN_Geometry_Analytic.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
namespace ge
{
class Sphere : public Geometry_Analytic, public Field_Analytic
{
public:

    /**
     * Constructor with only constant parameters
     *
     * @param aXCenter x-coordinate of the center of the sphere
     * @param aYCenter y-coordiante of the center of the sphere
     * @param aZCenter z-coordinate of the center of the sphere
     * @param aRadius radius of the circle
     */
    Sphere(real aXCenter, real aYCenter, real aZCenter, real aRadius)
    : Field(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aRadius}}))
    {
        mXCenter = aXCenter;
        mYCenter = aYCenter;
        mZCenter = aZCenter;
        mRadius = aRadius;
    }

    //------------------------------------ these will be deleted/modified later ----------------------------------------
    real evaluate_field_value(const moris::Matrix<moris::DDRMat> &aCoordinates)
    {
        return evaluate_field_value_with_coordinate(0, aCoordinates);
    }

    void evaluate_all_sensitivities(const moris::Matrix<moris::DDRMat> &aCoordinates, Matrix<DDRMat>& aSensitivities)
    {
        return evaluate_sensitivity_dphi_dp_with_coordinate(0, aCoordinates, aSensitivities);
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
                            + (aCoordinates(aRowIndex, 2) - mZCenter) * (aCoordinates(aRowIndex, 2) - mZCenter)
                              - (mRadius * mRadius);

        return tFunctionValue;
    }


    void evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                        moris::Matrix< moris::DDRMat > const & aCoordinates,
                                                       Matrix<DDRMat>& aSensitivities)
    {
        aSensitivities.resize(4, 3);

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

        (aSensitivities)(0, 0) = sign * mRadius / std::sqrt(std::abs(tSqrt));

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
        {
            std::cout << "zero denominator detected";
        }

        aSensitivities(0, 1) = mRadius / std::sqrt(std::abs(tSqrt));

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
        {
            std::cout << "zero denominator detected";
        }

        aSensitivities(0, 2) = sign * mRadius / std::sqrt(std::abs(tSqrt));

        aSensitivities(1, 0) = 1.0; // dx/dxc
        aSensitivities(1, 1) = 0.0; // dy/dxc
        aSensitivities(1, 2) = 0.0; // dz/dxc
        aSensitivities(2, 0) = 0.0; // dx/dyc
        aSensitivities(2, 1) = 1.0; // dy/dyc
        aSensitivities(2, 2) = 0.0; // dz/dyc
        aSensitivities(3, 0) = 0.0; // dx/dzc
        aSensitivities(3, 1) = 0.0; // dy/dzc
        aSensitivities(3, 2) = 1.0; // dz/dzc
    }

private:
    moris::real mRadius;
    moris::real mXCenter;
    moris::real mYCenter;
    moris::real mZCenter;
};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_SPHERE_HPP_ */
