/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Constant_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"

namespace moris::ge
{
    class Constant_Field : public Field_Analytic< 0 >
    {
      public:

        // Constructor to allow this field to be created with ADVs
        ANALYTIC_FIELD_ADV_CONSTRUCTOR( Constant_Field, 0, 1, { mSensitivities = { { 1.0 } }; } )

        /**
         * Given a node index, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Field value
         */
        real get_field_value( const Matrix< DDRMat >& aCoordinates );

        /**
         * Given a node index, evaluates the sensitivity of the field with respect to all of the
         * field variables.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
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
