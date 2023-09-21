/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Constant_Property.hpp
 *
 */

#pragma once

#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris::ge
{
    class Constant_Property : public Field_Analytic
    {

      public:
        /**
         * Constructor
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aPropertyVariableIndices Indices of property variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the property variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aFieldDependencies Other created fields that this property depends on
         * @param aParameters Additional parameters
         */
        template< typename Vector_Type >
        Constant_Property(
                Vector_Type&              aADVs,
                Matrix< DDUMat >          aPropertyVariableIndices,
                Matrix< DDUMat >          aADVIndices,
                Matrix< DDRMat >          aConstants,
                Property_Parameters aParameters = {} )
                : Field_Analytic( aADVs, aPropertyVariableIndices, aADVIndices, aConstants )
        {
            MORIS_ERROR( mVariables.size() == 1, "A constant property has only one variable." );
            mSensitivities = { { 1.0 } };
        }

        /**
         * default constructor
         */
        ~Constant_Property() {}

        /**
         * Given a node index, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Property value
         */
        real get_field_value( const Matrix< DDRMat >& aCoordinates );

        /**
         * Given a node index, evaluates the sensitivity of the property field with respect to all of the
         * property variables.
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
