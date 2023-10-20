/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Analytic.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"

/**
 * \def ANALYTIC_FIELD_ADV_CONSTRUCTOR( class_name, num_variables, ... )
 * Automatically creates a constructor that has ADV arguments for general field creation.
 *
 * @param class_name Name of the specific field constructor to create
 * @param num_variables Number of variables this field takes in
 * @param ... Additional code to be inserted into the constructor body via __VA_ARGS__
 */
#define ANALYTIC_FIELD_ADV_CONSTRUCTOR( class_name, num_variables, ... ) class_name( ADV_ARG_TYPES ) : Field_Analytic( ADV_ARGS ) { \
        VARIABLE_CHECK( num_variables ); __VA_ARGS__ }

namespace moris::ge
{
    class Field_Analytic : public Field
    {
      public:
        /**
         * Constructor using pointers to ADVs for variable evaluations.
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         */
        Field_Analytic( ADV_ARG_TYPES )
                : Field( ADV_ARGS )
        {
        }

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         * @param aParameters Additional parameters
         */
        explicit Field_Analytic( Matrix< DDRMat > aConstants );

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given a node coordinate, returns the field value
         *
         * @param aCoordinates vector of coordinate values
         * @return Field value
         */
        virtual real get_field_value( const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given a node coordinate, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        virtual const Matrix< DDRMat >& get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) override;

        /**
         * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        virtual void get_dfield_dcoordinates(
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) = 0;
    };
}
