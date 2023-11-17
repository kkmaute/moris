/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Superellipse.cpp
 *
 */

#include "cl_GEN_Superellipse.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Superellipse::Superellipse(
            real aXCenter,
            real aYCenter,
            real aXSemidiameter,
            real aYSemidiameter,
            real aExponent,
            real aScaling,
            real aRegularization,
            real aShift )
            : Field_Analytic< 2 >( Matrix<DDRMat>({{aXCenter, aYCenter, aXSemidiameter, aYSemidiameter, aExponent, aScaling, aRegularization, aShift}}) )
    {
        MORIS_ERROR( mADVManager.get_variable( 2 ) > 0 and mADVManager.get_variable( 3 ) > 0,
                "A GEN Super-ellipse must be created with positive semi-diameters.");

        MORIS_ERROR(std::abs(std::fmod( mADVManager.get_variable( 4 ),2.0)) < 1e-12,
                "A GEN Super-ellipse must be created with an even exponent.");
    }

    //--------------------------------------------------------------------------------------------------------------

    real Superellipse::get_field_value(const Matrix<DDRMat>& aCoordinates)
    {
        // Get variables
        real tXCenter        = mADVManager.get_variable( 0 );
        real tYCenter        = mADVManager.get_variable( 1 );
        real tXSemidiameter  = mADVManager.get_variable( 2 );
        real tYSemidiameter  = mADVManager.get_variable( 3 );
        real tExponent       = mADVManager.get_variable( 4 );
        real tScaling        = mADVManager.get_variable( 5 );
        real tRegularization = mADVManager.get_variable( 6 );
        real tShift          = mADVManager.get_variable( 7 );

        // Evaluate field
        real tConstant = pow(
                pow((aCoordinates(0) - tXCenter) / tXSemidiameter, tExponent) +
                pow((aCoordinates(1) - tYCenter) / tYSemidiameter, tExponent) +
                pow( tRegularization                             , tExponent), 1.0 / tExponent) - tRegularization;

        real tLevelset = tScaling * (tConstant - 1.0);

        // Ensure that level set value is not approx. zero at evaluation point
        if ( std::abs(tLevelset) < tShift)
        {
            tLevelset += tLevelset < 0.0 ? -tShift : tShift;
        }

        return tLevelset;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix<DDRMat>& Superellipse::get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates)
    {
        // Get variables
        real tXCenter        = mADVManager.get_variable( 0 );
        real tYCenter        = mADVManager.get_variable( 1 );
        real tXSemidiameter  = mADVManager.get_variable( 2 );
        real tYSemidiameter  = mADVManager.get_variable( 3 );
        real tExponent       = mADVManager.get_variable( 4 );
        real tScaling        = mADVManager.get_variable( 5 );
        real tRegularization = mADVManager.get_variable( 6 );

        // Constant in all calculations
        real tConstant0 = pow(
                pow((aCoordinates(0) - tXCenter)/tXSemidiameter, tExponent) +
                pow((aCoordinates(1) - tYCenter)/tYSemidiameter, tExponent) +
                pow( tRegularization                           , tExponent), (1.0 / tExponent - 1.0));

        real tConstant1 = pow((aCoordinates(0) - tXCenter)/tXSemidiameter,tExponent - 1.0);
        real tConstant2 = pow((aCoordinates(1) - tYCenter)/tYSemidiameter,tExponent - 1.0);

        // Calculate sensitivities
        mSensitivities(0) = -tScaling*tConstant1*tConstant0/tXSemidiameter;
        mSensitivities(1) = -tScaling*tConstant2*tConstant0/tYSemidiameter;

        mSensitivities(2) = -tScaling*(aCoordinates(0) - tXCenter)*tConstant1*tConstant0/pow(tXSemidiameter,2.0);
        mSensitivities(3) = -tScaling*(aCoordinates(1) - tYCenter)*tConstant2*tConstant0/pow(tYSemidiameter,2.0);

        // the reminder sensitivities are typically not used and therefore not calculated
        mSensitivities(4) = MORIS_REAL_MAX;
        mSensitivities(5) = MORIS_REAL_MAX;
        mSensitivities(6) = MORIS_REAL_MAX;
        mSensitivities(7) = MORIS_REAL_MAX;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Superellipse::get_dfield_dcoordinates(
            const Matrix<DDRMat>& aCoordinates,
            Matrix<DDRMat>&       aSensitivities)
    {
        MORIS_ERROR(false, "get_dfield_dcoordinates not implemented for superellipse geometry.");
    }

    //--------------------------------------------------------------------------------------------------------------

}
