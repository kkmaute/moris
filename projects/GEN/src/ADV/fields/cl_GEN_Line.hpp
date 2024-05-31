/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Plane.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"

namespace moris::gen
{
    class Line : public Field_Analytic< 2 >
    {
      public:

        // Constructor to allow this field to be created with ADVs
        ANALYTIC_FIELD_ADV_CONSTRUCTOR( Line, 2, 4, {} )

        /**
         * Constructor
         *
         * @param aXCenter x-coordinate of the center of the line
         * @param aYCenter y-coordinate of the center of the line
         * @param aXNormal x normal for the line
         * @param aYNormal y normal for the line
         * @param aName Name of this field
         */
        Line(
                const ADV&  aXCenter,
                const ADV&  aYCenter,
                const ADV&  aXNormal,
                const ADV&  aYNormal,
                std::string aName = "" );

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Field value
         */
        real get_field_value( const Matrix< DDRMat >& aCoordinates );

        /**
         * Given a node coordinate, evaluates the sensitivity of the field with respect to all of the
         * field variables.
         *
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates );

        /**
         * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities );
    };
}
