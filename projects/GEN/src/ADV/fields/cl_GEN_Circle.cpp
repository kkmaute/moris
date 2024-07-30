/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Circle.cpp
 *
 */

#include "cl_GEN_Circle.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Circle::Circle(
            const ADV&  aXCenter,
            const ADV&  aYCenter,
            const ADV&  aRadius,
            std::string aName )
            : Field_Analytic< 2 >( { aXCenter, aYCenter, aRadius }, std::move( aName ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Circle::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVHandler.get_variable( 0 );
        real tYCenter = mADVHandler.get_variable( 1 );
        real tRadius  = mADVHandler.get_variable( 2 );

        // Evaluate field
        return sqrt( pow( aCoordinates( 0 ) - tXCenter, 2 ) + pow( aCoordinates( 1 ) - tYCenter, 2 ) ) - tRadius;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Circle::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVHandler.get_variable( 0 );
        real tYCenter = mADVHandler.get_variable( 1 );

        // Calculate sensitivities
        real tConstant = sqrt( pow( aCoordinates( 0 ) - tXCenter, 2 ) + pow( aCoordinates( 1 ) - tYCenter, 2 ) );

        tConstant = tConstant ? 1.0 / tConstant : 0.0;

        mSensitivities( 0 ) = tConstant * ( tXCenter - aCoordinates( 0 ) );
        mSensitivities( 1 ) = tConstant * ( tYCenter - aCoordinates( 1 ) );
        mSensitivities( 2 ) = -1.0;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Circle::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        // Get variables
        real tXCenter = mADVHandler.get_variable( 0 );
        real tYCenter = mADVHandler.get_variable( 1 );

        // Compute level set value
        real tLevelSet = sqrt( pow( aCoordinates( 0 ) - tXCenter, 2 ) + pow( aCoordinates( 1 ) - tYCenter, 2 ) );

        if ( std::abs( tLevelSet ) < MORIS_REAL_EPS )
        {
            tLevelSet = tLevelSet < 0.0 ? -MORIS_REAL_EPS : MORIS_REAL_EPS;
        }

        aSensitivities( 0 ) = ( aCoordinates( 0 ) - tXCenter ) / tLevelSet;
        aSensitivities( 1 ) = ( aCoordinates( 1 ) - tYCenter ) / tLevelSet;
    }

    //--------------------------------------------------------------------------------------------------------------

}
