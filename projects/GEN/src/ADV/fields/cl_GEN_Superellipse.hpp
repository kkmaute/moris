/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Superellipse.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"

namespace moris::gen
{
    class Superellipse : public Field_Analytic< 2 >
    {
      public:
        // Constructor to allow this field to be created with ADVs
        ANALYTIC_FIELD_ADV_CONSTRUCTOR( Superellipse, 2, 7, {
            MORIS_ERROR( mADVHandler.get_variable( 2 ) > 0 and mADVHandler.get_variable( 3 ) > 0,
                    "A GEN Superellipse must be created with positive semi-diameters." );

            MORIS_ERROR( std::abs( std::fmod( mADVHandler.get_variable( 4 ), 2.0 ) ) < 1e-12,
                    "A GEN Superellipse must be created with an even exponent." );
        } )

        /**
         * Constructor
         *
         * @param aXCenter x-coordinate of the center of the superellipse
         * @param aYCenter y-coordiante of the center of the superellipse
         * @param aXSemidiameter Superellipse semi-diameter in the x direction
         * @param aYSemidiameter Superellipse semi-diameter in the y direction
         * @param aExponent Superellipse exponent
         * @param aName Name of this field
         */
        Superellipse(
                const ADV&  aXCenter,
                const ADV&  aYCenter,
                const ADV&  aXSemidiameter,
                const ADV&  aYSemidiameter,
                real        aExponent,
                real        aScaling,
                real        aRegularization,
                std::string aName = "" );

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Field value
         */
        real get_field_value( const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given a node coordinate, evaluates the sensitivity of the field with respect to all of the
         * field variables.
         *
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) override;
    };
}    // namespace moris::gen
