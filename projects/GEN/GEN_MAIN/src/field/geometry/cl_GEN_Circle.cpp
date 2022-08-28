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

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Circle::Circle(
                real                      aXCenter,
                real                      aYCenter,
                real                      aRadius,
                Geometry_Field_Parameters aParameters )
                : Field( Matrix< DDRMat >( { { aXCenter, aYCenter, aRadius } } ), aParameters )
                , Geometry( aParameters )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Circle::get_field_value( const Matrix< DDRMat >& aCoordinates )
        {
            // Get variables
            real tXCenter = *( mFieldVariables( 0 ) );
            real tYCenter = *( mFieldVariables( 1 ) );
            real tRadius  = *( mFieldVariables( 2 ) );

            // Evaluate field
            return sqrt( pow( aCoordinates( 0 ) - tXCenter, 2 ) + pow( aCoordinates( 1 ) - tYCenter, 2 ) ) - tRadius;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Circle::get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates )
        {
            // Get variables
            real tXCenter = *( mFieldVariables( 0 ) );
            real tYCenter = *( mFieldVariables( 1 ) );

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
            MORIS_ERROR( false, "get_dfield_dcoordinates not implemented for circle geometry." );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

