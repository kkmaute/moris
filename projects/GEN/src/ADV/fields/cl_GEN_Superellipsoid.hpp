/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Superellipsoid.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"

namespace moris::gen
{
    class Superellipsoid : public Field_Analytic< 3 >
    {
    private:
        real mEpsilon = 1E-8;

    public:

        // Constructor to allow this field to be created with ADVs
        ANALYTIC_FIELD_ADV_CONSTRUCTOR( Superellipsoid, 3, 7,
        {
            MORIS_ERROR( mADVHandler.get_variable( 3 ) > 0 and mADVHandler.get_variable( 4 ) > 0,
                    "A GEN Superellipsoid must be created with positive semi-diameters.");

            MORIS_ERROR( std::abs( std::fmod( mADVHandler.get_variable( 6 ), 2.0) ) < 1e-12,
                    "A GEN Superellipsoid must be created with an even exponent.");
        } )

        /**
         * Constructor with only constant parameters
         *
         * @param aXCenter x-coordinate of the center of the superellipsoid
         * @param aYCenter y-coordiante of the center of the superellipsoid
         * @param aZCenter z-coordinate of the center of the superellipsoid
         * @param aXSemidiameter Superellipsoid semi-diameter in the x direction
         * @param aYSemidiameter Superellipsoid semi-diameter in the y direction
         * @param aZSemidiameter Superellipsoid semi-diameter in the z direction
         * @param aExponent Superellipsoid exponent
         * @param aName Name of this field
         */
        Superellipsoid(
                const ADV&  aXCenter,
                const ADV&  aYCenter,
                const ADV&  aZCenter,
                const ADV&  aXSemidiameter,
                const ADV&  aYSemidiameter,
                const ADV&  aZSemidiameter,
                real        aExponent,
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
}
