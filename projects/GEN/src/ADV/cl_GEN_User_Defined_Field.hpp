/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_User_Defined_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Library_IO.hpp"

namespace moris::ge
{
    // User-defined field functions
    typedef real ( *Field_Function ) (
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aParameters );
    typedef void ( *Sensitivity_Function ) (
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aParameters,
            moris::Matrix< DDRMat >&           aReturnValue);

    class User_Defined_Field : public Field_Analytic
    {

    private:
        Field_Function get_field_value_user_defined;
        Sensitivity_Function get_dfield_dadvs_user_defined;

    public:

        /**
         * Constructor, sets the pointers to advs and constant parameters for evaluations.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aFieldFunction User-defined function for evaluating the field field
         * @param aSensitivityFunction User-defined function for evaluating the field sensitivities
         */
        template <typename Vector_Type>
        User_Defined_Field(
                Vector_Type&         aADVs,
                Matrix<DDUMat>       aFieldVariableIndices,
                Matrix<DDUMat>       aADVIndices,
                Matrix<DDRMat>       aConstants,
                Field_Function       aFieldFunction,
                Sensitivity_Function aSensitivityFunction )
                : Field_Analytic( aADVs, aFieldVariableIndices, aADVIndices, aConstants )
        {
            this->set_user_defined_functions(aFieldFunction, aSensitivityFunction);
        }

        /**
         * Constructor with only constants and no sensitivities.
         *
         * @param aConstants The constant field variables not filled by ADVs
         * @param aFieldFunction User-defined function for evaluating the field field
         */
        User_Defined_Field(
                Matrix<DDRMat>   aConstants,
                Field_Function   aFieldFunction );

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Field value
         */
        real get_field_value(const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node coordinate, evaluates the sensitivity of the field field with respect to all of the
         * field variables.
         *
         * @param aCoordinates Coordinate values
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates);

        /**
         * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities);

    private:

        /**
         * Sets the user-defined functions. Eliminates redundant code since it's the same logic for all constructors.
         *
         * @param aFieldFunction User-defined function for evaluating the field field
         * @param aSensitivityFunction User-defined function for evaluating the field sensitivities
         */
        void set_user_defined_functions(
                Field_Function       aFieldFunction,
                Sensitivity_Function aSensitivityFunction);

        /**
         * Used internally to automatically error out if no sensitivities were provided
         */
        static void no_sensitivities(
                const Matrix<DDRMat>& aCoordinates,
                const Cell<real*>&    aParameters,
                Matrix<DDRMat>&       aSensitivities);

    };
}
