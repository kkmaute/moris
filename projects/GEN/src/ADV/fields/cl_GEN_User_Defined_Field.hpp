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
#include "cl_GEN_Field_Analytic.hpp"
#include "cl_Library_IO.hpp"

namespace moris::gen
{
    // User-defined field functions
    typedef real ( *Field_Function )(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aParameters );
    typedef void ( *Sensitivity_Function )(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aParameters,
            Matrix< DDRMat >&       aSensitivities );

    class User_Defined_Field : public Field_Analytic< 0 >    // TODO consider templating user-defined field against reference dimensions as well
    {

      private:
        Vector< real >       mFieldVariables;
        Field_Function       get_field_value_user_defined;
        Sensitivity_Function get_dfield_dadvs_user_defined;

      public:
        /**
         * Constructor, sets the pointers to advs and constant parameters for evaluations.
         *
         * @param aFieldFunction User-defined function for evaluating the field
         * @param aSensitivityFunction User-defined function for evaluating the field sensitivities
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         */
        User_Defined_Field(
                Field_Function       aFieldFunction,
                Sensitivity_Function aSensitivityFunction,
                ADV_ARG_TYPES )
                : Field_Analytic< 0 >( ADV_ARGS )
                , mFieldVariables( aFieldVariableIndices.size() + aConstants.size() )
                , get_field_value_user_defined( aFieldFunction )
                , get_dfield_dadvs_user_defined( aSensitivityFunction )
        {
            this->import_advs( nullptr );
            this->validate_user_defined_functions();
        }

        /**
         * Constructor with both a user-defined field and user-defined sensitivity function.
         *
         * @param aFieldFunction User-defined function for evaluating the field
         * @param aSensitivityFunction User-defined function for evaluating the field sensitivities
         * @param aADVs The parameters that define this field and may be changed as a part of a design
         * @param aName Name of this field
         */
        User_Defined_Field(
                Field_Function       aFieldFunction,
                Sensitivity_Function aSensitivityFunction,
                const Vector< ADV >& aADVs = {},
                std::string          aName = "" );

        /**
         * Constructor with only a user-defined field function.
         *
         * @param aFieldFunction User-defined function for evaluating the field
         * @param aADVs The parameters that define this field and may be changed as a part of a design
         * @param aName Name of this field
         */
        explicit User_Defined_Field(
                Field_Function       aFieldFunction,
                const Vector< ADV >& aADVs = {},
                std::string          aName = "" );

        /**
         * For the specific case of a user-defined field, this function indicates that new field variables must be set for the user-defined function calls.
         *
         * @param aOwnedADVs Full owned distributed ADV vector (not used for this field)
         */
        void import_advs( sol::Dist_Vector* aOwnedADVs ) final;

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Field value
         */
        real get_field_value( const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given a node coordinate, evaluates the sensitivity of the field field with respect to all of the
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

      private:
        /**
         * Validates the user-defined functions. Eliminates redundant code since it's the same logic for all constructors.
         */
        void validate_user_defined_functions();

        /**
         * Used internally to automatically error out if no sensitivities were provided
         */
        static void no_sensitivities(
                const Matrix< DDRMat >& aCoordinates,
                const Vector< real >&     aParameters,
                Matrix< DDRMat >&       aSensitivities );
    };
}
