/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Superellipsoid.cpp
 *
 */

#include "cl_GEN_Superellipsoid.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Superellipsoid::Superellipsoid(
            real                      aXCenter,
            real                      aYCenter,
            real                      aZCenter,
            real                      aXSemidiameter,
            real                      aYSemidiameter,
            real                      aZSemidiameter,
            real                      aExponent,
            Level_Set_Parameters aParameters)
            : Field_Analytic(Matrix<DDRMat>({{aXCenter, aYCenter, aZCenter, aXSemidiameter, aYSemidiameter, aZSemidiameter, aExponent}}) )
    {
        MORIS_ERROR(*(mVariables(3)) > 0 and *(mVariables(4)) > 0 and *(mVariables(5)) > 0,
                    "A GEN Superellipsoid must be created with positive semidiameters.");
    }

    //--------------------------------------------------------------------------------------------------------------

    real Superellipsoid::get_field_value(const Matrix<DDRMat>& aCoordinates)
    {
        // Get variables
        real tXCenter = *(mVariables(0));
        real tYCenter = *(mVariables(1));
        real tZCenter = *(mVariables(2));
        real tXSemidiameter = *(mVariables(3));
        real tYSemidiameter = *(mVariables(4));
        real tZSemidiameter = *(mVariables(5));
        real tExponent = *(mVariables(6));

        // Evaluate field
        return pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent)
                 + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent)
                 + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent), 1.0 / tExponent) - 1.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Superellipsoid::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
    {
        // Get variables
        real tXCenter = *(mVariables(0));
        real tYCenter = *(mVariables(1));
        real tZCenter = *(mVariables(2));
        real tXSemidiameter = *(mVariables(3));
        real tYSemidiameter = *(mVariables(4));
        real tZSemidiameter = *(mVariables(5));
        real tExponent = *(mVariables(6));

        // Constant in all calculations
        real tConstant = pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent)
                       + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent)
                       + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent);
        tConstant = tConstant ? pow(tConstant, -1.0 + (1.0 / tExponent)) : 0.0;

        // Calculate sensitivities
        mSensitivities(0) = -tConstant * pow(1.0 / tXSemidiameter, tExponent) * (aCoordinates(0) - tXCenter)
                * pow(std::abs(tXCenter - aCoordinates(0)), tExponent - 2.0);
        mSensitivities(1) = -tConstant * pow(1.0 / tYSemidiameter, tExponent) * (aCoordinates(1) - tYCenter)
                * pow(std::abs(tYCenter - aCoordinates(1)), tExponent - 2.0);
        mSensitivities(2) = -tConstant * pow(1.0 / tZSemidiameter, tExponent) * (aCoordinates(2) - tZCenter)
                * pow(std::abs(tZCenter - aCoordinates(2)), tExponent - 2.0);
        mSensitivities(3) = -tConstant * pow(1.0 / tXSemidiameter, tExponent + 1.0)
                * pow(std::abs(tXCenter - aCoordinates(0)), tExponent);
        mSensitivities(4) = -tConstant * pow(1.0 / tYSemidiameter, tExponent + 1.0)
                * pow(std::abs(tYCenter - aCoordinates(1)), tExponent);
        mSensitivities(5) = -tConstant * pow(1.0 / tZSemidiameter, tExponent + 1.0)
                * pow(std::abs(tZCenter - aCoordinates(2)), tExponent);

        // TODO? this uses FD only because the analytical function is super complicated for the exponent derivative
        mSensitivities(6) = (pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent + mEpsilon)
                               + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent + mEpsilon)
                               + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent + mEpsilon), 1.0 / (tExponent + mEpsilon))
                           - pow(pow(std::abs(aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent - mEpsilon)
                               + pow(std::abs(aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent - mEpsilon)
                               + pow(std::abs(aCoordinates(2) - tZCenter) / tZSemidiameter, tExponent - mEpsilon), 1.0 / (tExponent - mEpsilon)))
                           / (2.0 * mEpsilon);

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Superellipsoid::get_dfield_dcoordinates(
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for superellipsoid geometry.");
    }

    //--------------------------------------------------------------------------------------------------------------

}
