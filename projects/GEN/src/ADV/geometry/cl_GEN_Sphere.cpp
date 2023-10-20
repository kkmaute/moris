/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Sphere.cpp
 *
 */

#include "cl_GEN_Sphere.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Sphere::Sphere(
            real                      aXCenter,
            real                      aYCenter,
            real                      aZCenter,
            real                      aRadius,
            Level_Set_Parameters aParameters )
            : Field_Analytic( Matrix< DDRMat >( { { aXCenter, aYCenter, aZCenter, aRadius } } ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Sphere::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVManager.get_variable( 0 );
        real tYCenter = mADVManager.get_variable( 1 );
        real tZCenter = mADVManager.get_variable( 2 );
        real tRadius  = mADVManager.get_variable( 3 );

        // Evaluate field
        return sqrt( pow( aCoordinates( 0 ) - tXCenter, 2 )
                       + pow( aCoordinates( 1 ) - tYCenter, 2 )
                       + pow( aCoordinates( 2 ) - tZCenter, 2 ) )
             - tRadius;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Sphere::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVManager.get_variable( 0 );
        real tYCenter = mADVManager.get_variable( 1 );
        real tZCenter = mADVManager.get_variable( 2 );

        // Calculate sensitivities
        real tConstant = sqrt( pow( aCoordinates( 0 ) - tXCenter, 2 )
                               + pow( aCoordinates( 1 ) - tYCenter, 2 )
                               + pow( aCoordinates( 2 ) - tZCenter, 2 ) );

        tConstant = tConstant ? 1.0 / tConstant : 0.0;

        mSensitivities( 0 ) = tConstant * ( tXCenter - aCoordinates( 0 ) );
        mSensitivities( 1 ) = tConstant * ( tYCenter - aCoordinates( 1 ) );
        mSensitivities( 2 ) = tConstant * ( tZCenter - aCoordinates( 2 ) );
        mSensitivities( 3 ) = -1;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Sphere::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        MORIS_ERROR( false, "get_dfield_dcoordinates not implemented for sphere geometry." );
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace ge
