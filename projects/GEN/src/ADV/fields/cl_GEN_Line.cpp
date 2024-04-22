/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Plane.cpp
 *
 */

#include "cl_GEN_Line.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Line::Line(
            real aXCenter,
            real aYCenter,
            real aXNormal,
            real aYNormal )
            : Field_Analytic< 2 >( { aXCenter, aYCenter, aXNormal, aYNormal } )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Line::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVManager.get_variable( 0 );
        real tYCenter = mADVManager.get_variable( 1 );
        real tXNormal = mADVManager.get_variable( 2 );
        real tYNormal = mADVManager.get_variable( 3 );
        
        // Evaluate field value
        return tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Line::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        real tXCenter = mADVManager.get_variable( 0 );
        real tYCenter = mADVManager.get_variable( 1 );
        real tXNormal = mADVManager.get_variable( 2 );
        real tYNormal = mADVManager.get_variable( 3 );
        
        // Evaluate sensitivities
        mSensitivities( 0 ) = -tXNormal;
        mSensitivities( 1 ) = -tYNormal;
        mSensitivities( 2 ) = aCoordinates( 0 ) - tXCenter;
        mSensitivities( 3 ) = aCoordinates( 1 ) - tYCenter;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Line::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        aSensitivities( 0 ) = mADVManager.get_variable( 2 );
        aSensitivities( 1 ) = mADVManager.get_variable( 3 );
    }

    //--------------------------------------------------------------------------------------------------------------
}
