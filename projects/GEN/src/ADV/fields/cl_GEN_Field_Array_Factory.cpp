/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Array_Factory.cpp
 *
 */

#include "cl_GEN_Field_Array_Factory.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Field_Array_Factory::Field_Array_Factory( const ParameterList& aFieldArrayParameters )
            : mOffsetPerRowX( aFieldArrayParameters.get< real >( "offset_per_row_x") )
            , mOffsetPerRowY( aFieldArrayParameters.get< real >( "offset_per_row_y") )
            , mOffsetPerRowZ( aFieldArrayParameters.get< real >( "offset_per_row_z") )
    {
        this->set_x_parameters(
                aFieldArrayParameters.get< real >( "lower_bound_x" ),
                aFieldArrayParameters.get< real >( "upper_bound_x" ),
                aFieldArrayParameters.get< sint >( "number_of_fields_x" ),
                aFieldArrayParameters.get< real >( "minimum_spacing_x" ) );
        this->set_y_parameters(
                aFieldArrayParameters.get< real >( "lower_bound_y" ),
                aFieldArrayParameters.get< real >( "upper_bound_y" ),
                aFieldArrayParameters.get< sint >( "number_of_fields_y" ),
                aFieldArrayParameters.get< real >( "minimum_spacing_y" ) );
        this->set_z_parameters(
                aFieldArrayParameters.get< real >( "lower_bound_z" ),
                aFieldArrayParameters.get< real >( "upper_bound_z" ),
                aFieldArrayParameters.get< sint >( "number_of_fields_z" ),
                aFieldArrayParameters.get< real >( "minimum_spacing_z" ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Array_Factory::set_x_parameters(
            real aLowerBound,
            real aUpperBound,
            sint aNumberOfFields,
            real aMinimumSpacing )
    {
        // Check upper/lower bounds
        MORIS_ERROR( aUpperBound >= aLowerBound,
                "Upper bounds must be greater than or equal to lower bounds for a GEN field array (x-direction)." );

        // Set number of fields
        mNumberOfFieldsX = this->calculate_final_number_of_fields(
                aNumberOfFields,
                aMinimumSpacing,
                aUpperBound - aLowerBound );

        // Set starting point and spacing from lower bound, or use middle point if there's only one field
        if ( mNumberOfFieldsX > 1 )
        {
            mStartingPointX = aLowerBound;
            mSpacingX = ( aUpperBound - aLowerBound ) / ( mNumberOfFieldsX - 1 );
        }
        else
        {
            mStartingPointX = ( aUpperBound + aLowerBound ) / 2;
            mSpacingX = MORIS_REAL_MAX;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Array_Factory::set_y_parameters(
            real aLowerBound,
            real aUpperBound,
            sint aNumberOfFields,
            real aMinimumSpacing )
    {
        // Check upper/lower bounds
        MORIS_ERROR( aUpperBound >= aLowerBound,
                "Upper bounds must be greater than or equal to lower bounds for a GEN field array (y-direction)." );

        // Set number of fields
        mNumberOfFieldsY = this->calculate_final_number_of_fields(
                aNumberOfFields,
                aMinimumSpacing,
                aUpperBound - aLowerBound );

        // Set starting point and spacing from lower bound, or use middle point if there's only one field
        if ( mNumberOfFieldsY > 1 )
        {
            mStartingPointY = aLowerBound;
            mSpacingY = ( aUpperBound - aLowerBound ) / ( mNumberOfFieldsY - 1 );
        }
        else
        {
            mStartingPointY = ( aUpperBound + aLowerBound ) / 2;
            mSpacingY = MORIS_REAL_MAX;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Array_Factory::set_z_parameters(
            real aLowerBound,
            real aUpperBound,
            sint aNumberOfFields,
            real aMinimumSpacing )
    {
        // Check upper/lower bounds
        MORIS_ERROR( aUpperBound >= aLowerBound,
                "Upper bounds must be greater than or equal to lower bounds for a GEN field array (z-direction)." );

        // Set number of fields
        mNumberOfFieldsZ = this->calculate_final_number_of_fields(
                aNumberOfFields,
                aMinimumSpacing,
                aUpperBound - aLowerBound );

        // Set starting point and spacing from lower bound, or use middle point if there's only one field
        if ( mNumberOfFieldsZ > 1 )
        {
            mStartingPointZ = aLowerBound;
            mSpacingZ = ( aUpperBound - aLowerBound ) / ( mNumberOfFieldsZ - 1 );
        }
        else
        {
            mStartingPointZ = ( aUpperBound + aLowerBound ) / 2;
            mSpacingZ = MORIS_REAL_MAX;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Combined_Fields >
    Field_Array_Factory::create_field_array(
            std::shared_ptr< Field > aCopyField,
            bool aMinimum )
    {
        // Store field for copying
        mCopyField = aCopyField;

        // Reset variable indices
        mVariableIndices.clear();
        mVariableIndices.reserve( 3 );

        // Calculate total number of fields and set variable indices
        uint tTotalNumberOfFields = 1;
        for ( uint iNumberOfFieldsInDirection : { mNumberOfFieldsX, mNumberOfFieldsY, mNumberOfFieldsZ } )
        {
            if ( iNumberOfFieldsInDirection > 0 )
            {
                tTotalNumberOfFields *= iNumberOfFieldsInDirection;
                mVariableIndices.push_back( mVariableIndices.size() );
            }
        }

        // Check number of coordinates
        MORIS_ERROR( mVariableIndices.size() == mCopyField->get_number_of_reference_coordinates(),
                "Number of dimensions used in a field array (%lu) must be the same as the number of coordinates required to define a reference point (%d) for the given field %s",
                mVariableIndices.size(),
                mCopyField->get_number_of_reference_coordinates(),
                mCopyField->get_name().c_str() );

        // Resize new coordinates
        mCoordinates.resize( mVariableIndices.size() );

        // Reset field vector
        mCreatedFields.clear();
        mCreatedFields.reserve( tTotalNumberOfFields );

        // Populate created fields with array of fields, starting with x-direction (calls other directions internally)
        this->create_fields_loop_x();

        // Create combined fields
        return std::make_shared< Combined_Fields >( mCreatedFields, aMinimum, mCopyField->get_name() );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Field_Array_Factory::calculate_final_number_of_fields(
            sint aNumberOfFields,
            real aMinimumSpacing,
            real aTotalSpacing )
    {
        // Determine if either option was specified
        if ( aMinimumSpacing > 0.0 or aNumberOfFields > 0 )
        {
            // Calculate maximum number of fields permitted by spacing
            uint tMaximumNumberOfFields = MORIS_UINT_MAX;
            if ( aMinimumSpacing > 0.0 )
            {
                tMaximumNumberOfFields = std::floor( aTotalSpacing / aMinimumSpacing ) + 1;
            }

            // Just use minimum spacing if desired
            if ( aNumberOfFields <= 0 )
            {
                return tMaximumNumberOfFields;
            }
            // Otherwise, use the minimum of the number of provided fields and the maximum permitted by spacing
            else
            {
                return std::min( static_cast< uint >( aNumberOfFields ), tMaximumNumberOfFields );
            }
        }
        else
        {
            // Dimension is not used
            return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Array_Factory::create_fields_loop_x()
    {
        // Determine if this dimension is used
        if ( mNumberOfFieldsX > 0 )
        {
            // Loop over dimension
            for ( uint iFieldIndex = 0; iFieldIndex < mNumberOfFieldsX; iFieldIndex++ )
            {
                // Set x coordinate
                mCoordinates( 0 ) = mStartingPointX + iFieldIndex * mSpacingX;

                // Offset z coordinate (or y in 2D)
                if ( mNumberOfFieldsZ > 0 )
                {
                    mStartingPointZ += std::fmod( iFieldIndex * mOffsetPerRowZ, mSpacingZ );
                }
                else if ( mNumberOfFieldsY > 0 )
                {
                    mStartingPointY += std::fmod( iFieldIndex * mOffsetPerRowY, mSpacingY );
                }

                // Loop over y
                create_fields_loop_y( 1 );

                // Un-offset z coordinate
                if ( mNumberOfFieldsZ > 0 )
                {
                    mStartingPointZ -= std::fmod( iFieldIndex * mOffsetPerRowZ, mSpacingZ );
                }
                else if ( mNumberOfFieldsY > 0 )
                {
                    mStartingPointY -= std::fmod( iFieldIndex * mOffsetPerRowY, mSpacingY );
                }
            }
        }
        else
        {
            // Proceed to next dimension
            create_fields_loop_y( 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Array_Factory::create_fields_loop_y( uint aUsedDimensionIndex )
    {
        // Determine if this dimension is used
        if ( mNumberOfFieldsY > 0 )
        {
            // Loop over dimension
            for ( uint iFieldIndex = 0; iFieldIndex < mNumberOfFieldsY; iFieldIndex++ )
            {
                // Set y coordinate
                mCoordinates( aUsedDimensionIndex ) = mStartingPointY + iFieldIndex * mSpacingY;

                // Offset x coordinate
                if ( mNumberOfFieldsX > 0 )
                {
                    mCoordinates( 0 ) += std::fmod( iFieldIndex * mOffsetPerRowX, mSpacingX );
                }

                // Loop over z
                create_fields_loop_z( aUsedDimensionIndex + 1 );

                // Un-offset x coordinate
                if ( mNumberOfFieldsX > 0 )
                {
                    mCoordinates( 0 ) -= std::fmod( iFieldIndex * mOffsetPerRowX, mSpacingX );
                }
            }
        }
        else
        {
            // Proceed to next dimension
            create_fields_loop_z( aUsedDimensionIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Array_Factory::create_fields_loop_z( uint aUsedDimensionIndex )
    {
        // Determine if this dimension is used
        if ( mNumberOfFieldsZ > 0 )
        {
            // Loop over dimension
            for ( uint iFieldIndex = 0; iFieldIndex < mNumberOfFieldsZ; iFieldIndex++ )
            {
                // Set z-coordinate
                mCoordinates( aUsedDimensionIndex ) = mStartingPointZ + iFieldIndex * mSpacingZ;

                // Offset y coordinate
                if ( mNumberOfFieldsY > 0 )
                {
                    mCoordinates( aUsedDimensionIndex - 1 ) += std::fmod( iFieldIndex * mOffsetPerRowY, mSpacingY );
                }

                // Create new field
                mCreatedFields.push_back( mCopyField->copy( mVariableIndices, mCoordinates ) );

                // Un-offset y coordinate
                if ( mNumberOfFieldsY > 0 )
                {
                    mCoordinates( aUsedDimensionIndex - 1 ) -= std::fmod( iFieldIndex * mOffsetPerRowY, mSpacingY );
                }
            }
        }
        else
        {
            // Create field without looping
            mCreatedFields.push_back( mCopyField->copy( mVariableIndices, mCoordinates ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
