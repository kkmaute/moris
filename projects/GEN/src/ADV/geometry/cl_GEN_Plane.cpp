/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Plane.cpp
 *
 */

#include "cl_GEN_Plane.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Plane::Plane(
            real                      aXCenter,
            real                      aYCenter,
            real                      aZCenter,
            real                      aXNormal,
            real                      aYNormal,
            real                      aZNormal,
            Level_Set_Parameters aParameters )
            : Field_Analytic( Matrix< DDRMat >( { { aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal } } ) )
    {
        m_eval_field       = &Plane::eval_field_3d;
        m_eval_sensitivity = &Plane::eval_sensitivity_3d;
    }

    //--------------------------------------------------------------------------------------------------------------

    Plane::Plane(
            real                      aXCenter,
            real                      aYCenter,
            real                      aXNormal,
            real                      aYNormal,
            Level_Set_Parameters aParameters )
            : Field_Analytic( Matrix< DDRMat >( { { aXCenter, aYCenter, aXNormal, aYNormal } } ) )
    {
        m_eval_field       = &Plane::eval_field_2d;
        m_eval_sensitivity = &Plane::eval_sensitivity_2d;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Plane::get_field_value( const Matrix< DDRMat >& aCoordinates )
    {
        return ( this->*m_eval_field )( aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Plane::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
    {
        return ( this->*m_eval_sensitivity )( aCoordinates );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Plane::get_dfield_dcoordinates(
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        if ( aCoordinates.numel() == 2 )
        {
            aSensitivities( 0 ) = *( mVariables( 2 ) );
            aSensitivities( 1 ) = *( mVariables( 3 ) );
        }
        else
        {
            aSensitivities( 0 ) = *( mVariables( 3 ) );
            aSensitivities( 1 ) = *( mVariables( 4 ) );
            aSensitivities( 2 ) = *( mVariables( 5 ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Plane::eval_field_2d( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        const real tXCenter = *( mVariables( 0 ) );
        const real tYCenter = *( mVariables( 1 ) );
        const real tXNormal = *( mVariables( 2 ) );
        const real tYNormal = *( mVariables( 3 ) );

        // Evaluate field value
        return tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Plane::eval_field_3d( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        const real tXCenter = *( mVariables( 0 ) );
        const real tYCenter = *( mVariables( 1 ) );
        const real tZCenter = *( mVariables( 2 ) );
        const real tXNormal = *( mVariables( 3 ) );
        const real tYNormal = *( mVariables( 4 ) );
        const real tZNormal = *( mVariables( 5 ) );

        // Evaluate field value
        return tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter ) + tZNormal * ( aCoordinates( 2 ) - tZCenter );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Plane::eval_sensitivity_2d( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        const real tXCenter = *( mVariables( 0 ) );
        const real tYCenter = *( mVariables( 1 ) );
        const real tXNormal = *( mVariables( 2 ) );
        const real tYNormal = *( mVariables( 3 ) );

        // Evaluate sensitivities
        mSensitivities( 0 ) = -tXNormal;
        mSensitivities( 1 ) = -tYNormal;
        mSensitivities( 2 ) = aCoordinates( 0 ) - tXCenter;
        mSensitivities( 3 ) = aCoordinates( 1 ) - tYCenter;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Plane::eval_sensitivity_3d( const Matrix< DDRMat >& aCoordinates )
    {
        // Get variables
        const real tXCenter = *( mVariables( 0 ) );
        const real tYCenter = *( mVariables( 1 ) );
        const real tZCenter = *( mVariables( 2 ) );
        const real tXNormal = *( mVariables( 3 ) );
        const real tYNormal = *( mVariables( 4 ) );
        const real tZNormal = *( mVariables( 5 ) );

        // Evaluate sensitivities
        mSensitivities( 0 ) = -tXNormal;
        mSensitivities( 1 ) = -tYNormal;
        mSensitivities( 2 ) = -tZNormal;
        mSensitivities( 3 ) = aCoordinates( 0 ) - tXCenter;
        mSensitivities( 4 ) = aCoordinates( 1 ) - tYCenter;
        mSensitivities( 5 ) = aCoordinates( 2 ) - tZCenter;

        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

}
