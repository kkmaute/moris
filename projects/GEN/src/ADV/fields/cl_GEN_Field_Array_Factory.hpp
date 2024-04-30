/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Array_Factory.hpp
 *
 */

#pragma once

#include "cl_GEN_Combined_Fields.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{
    /**
     * Class containing parameters for creating a field array. Functions verify and set inputs.
     */
    class Field_Array_Factory
    {
      public:

        // Amount that the reference points will be shifted over for each subsequent row in the specified direction. No restrictions on these parameters.
        real mOffsetPerRowX;
        real mOffsetPerRowY;
        real mOffsetPerRowZ;

      private:
        // Starting points for creating the array in each direction
        real mStartingPointX;
        real mStartingPointY;
        real mStartingPointZ;

        // Number of fields to create in each direction
        uint mNumberOfFieldsX;
        uint mNumberOfFieldsY;
        uint mNumberOfFieldsZ;

        // Spacing in each direction
        real mSpacingX;
        real mSpacingY;
        real mSpacingZ;

        // Additional arguments to pass to the field copy constructor
        Vector< uint > mVariableIndices;
        Vector< real > mCoordinates;

        // Copy field and created fields
        std::shared_ptr< Field > mCopyField;
        Vector< std::shared_ptr< Field > > mCreatedFields;

      public:

        /**
         * Constructor, gets default parameters from PRM
         *
         * @param aFieldArrayParameters Field array parameter list
         */
        explicit Field_Array_Factory( const Parameter_List& aFieldArrayParameters = prm::create_field_array_parameter_list() );

        /**
         * Sets the array parameters in the x direction
         *
         * @param aLowerBound Lower bound for reference x coordinate
         * @param aUpperBound Upper bound for reference x coordinate
         * @param aNumberOfFields Number of fields in the x direction
         * @param aMinimumSpacing Minimum spacing in the x direction
         */
        void set_x_parameters(
                real aLowerBound,
                real aUpperBound,
                sint aNumberOfFields,
                real aMinimumSpacing = 0.0 );

        /**
         * Sets the array parameters in the y direction
         *
         * @param aLowerBound Lower bound for reference y coordinate
         * @param aUpperBound Upper bound for reference y coordinate
         * @param aNumberOfFields Number of fields in the y direction
         * @param aMinimumSpacing Minimum spacing in the y direction
         */
        void set_y_parameters(
                real aLowerBound,
                real aUpperBound,
                sint aNumberOfFields,
                real aMinimumSpacing = 0.0 );

        /**
         * Sets the array parameters in the z direction
         *
         * @param aLowerBound Lower bound for reference z coordinate
         * @param aUpperBound Upper bound for reference z coordinate
         * @param aNumberOfFields Number of fields in the z direction
         * @param aMinimumSpacing Minimum spacing in the z direction
         */
        void set_z_parameters(
                real aLowerBound,
                real aUpperBound,
                sint aNumberOfFields,
                real aMinimumSpacing = 0.0 );

        /**
         * Creates a field array with the previously set parameters, copying the given field
         *
         * @param aCopyField Field to copy into the array
         * @param aMinimum Whether to use the minimum or maximum field value
         * @return Array of fields as a Combined_Field
         */
        std::shared_ptr< Combined_Fields > create_field_array(
                std::shared_ptr< Field > aCopyField,
                bool aMinimum = true );

      private:

        /**
         * Calculates the final number of fields for the array in a given direction
         *
         * @param aNumberOfFields Number of fields specified by the user
         * @param aMinimumSpacing Minimum spacing value between reference coordinates
         * @param aTotalSpacing Total space available in the given direction
         * @return Number of fields to create
         */
        uint calculate_final_number_of_fields(
                sint aNumberOfFields,
                real aMinimumSpacing,
                real aTotalSpacing );

        /**
         * Creates fields in a loop over the x direction; internally calls y
         */
        void create_fields_loop_x();

        /**
         * Creates fields in a loop over the y direction; internally calls z
         *
         * @param aUsedDimensionIndex Previously used dimensions
         */
        void create_fields_loop_y( uint aUsedDimensionIndex );

        /**
         * Creates fields in a loop over the z direction
         *
         * @param aUsedDimensionIndex Previously used dimensions
         */
        void create_fields_loop_z( uint aUsedDimensionIndex );
    };
}
