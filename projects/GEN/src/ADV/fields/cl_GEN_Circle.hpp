/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Circle.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"

namespace moris::gen
{
    class Circle : public Field_Analytic< 2 >
    {
    public:

        // Constructor to allow this field to be created with ADVs
        ANALYTIC_FIELD_ADV_CONSTRUCTOR( Circle, 2, 3, {} )

        /**
         * Constructor
         *
         * @param aXCenter x-coordinate of the center of the circle
         * @param aYCenter y-coordiante of the center of the circle
         * @param aRadius Radius of the circle
         * @param aName Name of this field
         */
        Circle( const ADV&  aXCenter,
                const ADV&  aYCenter,
                const ADV&  aRadius,
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
