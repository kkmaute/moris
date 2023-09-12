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

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Plane::Plane(
                real                      aXCenter,
                real                      aYCenter,
                real                      aZCenter,
                real                      aXNormal,
                real                      aYNormal,
                real                      aZNormal,
                Geometry_Field_Parameters aParameters )
                : Field( Matrix< DDRMat >( { { aXCenter, aYCenter, aZCenter, aXNormal, aYNormal, aZNormal } } ), aParameters )
                , Geometry( aParameters )
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
                Geometry_Field_Parameters aParameters )
                : Field( Matrix< DDRMat >( { { aXCenter, aYCenter, aXNormal, aYNormal } } ), aParameters )
                , Geometry( aParameters )
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
                aSensitivities( 0 ) = *( mFieldVariables( 2 ) );
                aSensitivities( 1 ) = *( mFieldVariables( 3 ) );
            }
            else
            {
                aSensitivities( 0 ) = *( mFieldVariables( 3 ) );
                aSensitivities( 1 ) = *( mFieldVariables( 4 ) );
                aSensitivities( 2 ) = *( mFieldVariables( 5 ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Plane::eval_field_2d( const Matrix< DDRMat >& aCoordinates )
        {
            // Get variables
            const real tXCenter = *( mFieldVariables( 0 ) );
            const real tYCenter = *( mFieldVariables( 1 ) );
            const real tXNormal = *( mFieldVariables( 2 ) );
            const real tYNormal = *( mFieldVariables( 3 ) );

            // Evaluate field value
            return tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Plane::eval_field_3d( const Matrix< DDRMat >& aCoordinates )
        {
            // Get variables
            const real tXCenter = *( mFieldVariables( 0 ) );
            const real tYCenter = *( mFieldVariables( 1 ) );
            const real tZCenter = *( mFieldVariables( 2 ) );
            const real tXNormal = *( mFieldVariables( 3 ) );
            const real tYNormal = *( mFieldVariables( 4 ) );
            const real tZNormal = *( mFieldVariables( 5 ) );

            // Evaluate field value
            return tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter ) + tZNormal * ( aCoordinates( 2 ) - tZCenter );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Plane::eval_sensitivity_2d( const Matrix< DDRMat >& aCoordinates )
        {
            // Get variables
            const real tXCenter = *( mFieldVariables( 0 ) );
            const real tYCenter = *( mFieldVariables( 1 ) );
            const real tXNormal = *( mFieldVariables( 2 ) );
            const real tYNormal = *( mFieldVariables( 3 ) );

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
            const real tXCenter = *( mFieldVariables( 0 ) );
            const real tYCenter = *( mFieldVariables( 1 ) );
            const real tZCenter = *( mFieldVariables( 2 ) );
            const real tXNormal = *( mFieldVariables( 3 ) );
            const real tYNormal = *( mFieldVariables( 4 ) );
            const real tZNormal = *( mFieldVariables( 5 ) );

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

    }    // namespace ge
}    // namespace moris
