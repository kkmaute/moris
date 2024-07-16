/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Validator.cpp
 *
 */

#include "cl_Validator.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    static Vector< std::string > gNoSelections;

    //--------------------------------------------------------------------------------------------------------------

    std::string to_string( const std::string& aValue )
    {
        return "\"" + aValue + "\"";
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Vector< std::string >& Validator::get_selection_names()
    {
        return gNoSelections;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    bool Type_Validator< T >::make_valid_parameter( Variant& aVariant )
    {
        return aVariant.index() == variant_index< T >();
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    bool Type_Validator< uint >::make_valid_parameter( Variant& aVariant )
    {
        if ( aVariant.index() == variant_index< uint >() )
        {
            return true;
        }
        else if ( aVariant.index() == variant_index< sint >() )
        {
            // Get signed integer and check to make sure it is greater than or equal to zero
            sint tValue = std::get< sint >( aVariant );
            if ( tValue >= 0 )
            {
                aVariant = static_cast< uint >( tValue );
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    std::string Type_Validator< T >::get_validation_message()
    {
        return get_type_name< T >();
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    Validator* Type_Validator< T >::copy()
    {
        return new Type_Validator< T >();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    template< typename T >
    bool Vector_Validator< T >::make_valid_parameter( Variant& aVariant )
    {
        if ( aVariant.index() == variant_index< Vector< T > >() )
        {
            // Vector value given, okay
            return true;
        }
        else if ( aVariant.index() == variant_index< T >() )
        {
            // Single value given, convert to vector
            T tValue = std::get< T >( aVariant );
            Vector< T > tVector = { tValue };
            aVariant = tVector;
            return true;
        }
        else
        {
            // Other type given
            return false;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    bool Vector_Validator< uint >::make_valid_parameter( Variant& aVariant )
    {
        if ( aVariant.index() == variant_index< Vector< uint > >() )
        {
            // Vector value given, okay
            return true;
        }
        else if ( aVariant.index() == variant_index< uint >() )
        {
            // Single value given, convert to vector
            uint tValue = std::get< uint >( aVariant );
            Vector< uint > tVector = { tValue };
            aVariant = tVector;
            return true;
        }
        else if ( aVariant.index() == variant_index< Vector< sint > >() ) // Additional case for signed integer vector
        {
            // Create vectors
            Vector< sint > tSignedVector = std::get< Vector< sint > >( aVariant );
            Vector< uint > tUnsignedVector( tSignedVector.size() );

            // Check each element
            for ( uint iVectorIndex = 0; iVectorIndex < tSignedVector.size(); iVectorIndex++ )
            {
                // Only copy over if greater than or equal to zero
                if ( tSignedVector( iVectorIndex ) >= 0 )
                {
                    tUnsignedVector( iVectorIndex ) = tSignedVector( iVectorIndex );
                }
                else
                {
                    return false;
                }
            }

            // Success, if we got here
            aVariant = tUnsignedVector;
            return true;
        }
        else if ( aVariant.index() == variant_index< sint >() ) // Additional case for single signed integer
        {
            // Get value
            sint tValue = std::get< sint >( aVariant );

            // Check if greater than or equal to zero
            if ( tValue >= 0 )
            {
                Vector< uint > tVector = { static_cast< uint >( tValue ) };
                aVariant = tVector;
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            // All cases failed
            return false;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    std::string Vector_Validator< T >::get_validation_message()
    {
        return get_type_name< Vector< T > >();
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    Validator* Vector_Validator< T >::copy()
    {
        return new Vector_Validator< T >();
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    bool Range_Validator< T >::make_valid_parameter( Variant& aVariant )
    {
        return aVariant.index() == variant_index< T >()
           and std::get< T >( aVariant ) >= mMinimumValue
           and std::get< T >( aVariant ) <= mMaximumValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    std::string Range_Validator< T >::get_validation_message()
    {
        return get_type_name< T >() + ", [" + std::to_string( mMinimumValue ) + ", " + std::to_string( mMaximumValue ) + "]";
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    Validator* Range_Validator< T >::copy()
    {
        return new Range_Validator( mMinimumValue, mMaximumValue );
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    template< typename T >
    const Vector< std::string >& Selection_Validator< T >::get_selection_names()
    {
        return mValidSelections;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    bool Selection_Validator< T >::make_valid_parameter( Variant& aVariant )
    {
        // Check if index is correct
        if ( aVariant.index() == variant_index< T >() )
        {
            // Check that given variant is in valid selections
            T aParameter = std::get< T >( aVariant );
            for ( const auto& iValidSelection : mValidSelections )
            {
                if ( aParameter == iValidSelection )
                {
                    return true;
                }
            }
        }

        // Index not correct or valid selection not found
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    std::string Selection_Validator< T >::get_validation_message()
    {
        // Use std::to_string for valid arguments
        using namespace std;

        // Create string from the set of valid values
        std::string tValidSelectionString;
        for ( auto iValue : mValidSelections )
        {
            tValidSelectionString += ", " + to_string( iValue );
        }
        tValidSelectionString.erase( 0, 2 );

        // Return string
        return tValidSelectionString;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    Validator* Selection_Validator< T >::copy()
    {
        return new Selection_Validator( mValidSelections );
    }

    //--------------------------------------------------------------------------------------------------------------

    Enum_Validator::Enum_Validator()
            : mEnumStrings() // TODO get from class
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< std::string >& Enum_Validator::get_selection_names()
    {
        return mEnumStrings;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Enum_Validator::make_valid_parameter( Variant& aVariant )
    {
        // Check if uint was given directly
        if ( aVariant.index() == variant_index< uint >() )
        {
            // Check if index lies within vector size
            return std::get< uint >( aVariant ) < mEnumStrings.size();
        }
        else if ( aVariant.index() == variant_index< std::string >() )
        {
            // Check that given variant is a valid enum name
            std::string aParameter = std::get< std::string >( aVariant );
            for ( const auto& iEnumString : mEnumStrings )
            {
                if ( aParameter == iEnumString )
                {
                    return true;
                }
            }
        }

        // Index not correct or valid selection not found
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Enum_Validator::get_validation_message()
    {
        // Create message from enum strings
        std::string tMessage;
        std::string tDelimiter;
        for ( const auto& iEnumString : mEnumStrings )
        {
            tMessage += tDelimiter + iEnumString;
            tDelimiter = ", ";
        }

        // Return string
        return tMessage;
    }

    //--------------------------------------------------------------------------------------------------------------

    Validator* Enum_Validator::copy()
    {
        return new Enum_Validator();
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    bool Type_Validator< Design_Variable >::make_valid_parameter( Variant& aVariant )
    {
        if ( aVariant.index() == variant_index< Design_Variable >() )
        {
            // Design variable given directly, okay
            return true;
        }
        else if ( aVariant.index() == variant_index< real >() )
        {
            // Real value given, convert to constant parameter
            Design_Variable tDesignVariable( std::get< real >( aVariant ) );
            aVariant = tDesignVariable;
            return true;
        }
        else if ( aVariant.index() == variant_index< Vector< real > >() )
        {
            // Vector given, check that there are 3 elements then set to design variable
            Vector< real > tVariableInputs = std::get< Vector< real > >( aVariant );
            if ( tVariableInputs.size() == 3 )
            {
                Design_Variable tDesignVariable( tVariableInputs( 0 ), tVariableInputs( 1 ), tVariableInputs( 2 ) );
                aVariant = tDesignVariable;
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            // Incorrect type
            return false;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    bool Range_Validator< Design_Variable >::make_valid_parameter( Variant& aVariant )
    {
        return Type_Validator< Design_Variable >().make_valid_parameter( aVariant )
           and std::get< Design_Variable >( aVariant ) > mMinimumValue
           and std::get< Design_Variable >( aVariant ) < mMaximumValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    std::string Type_Validator< Design_Variable >::get_validation_message()
    {
        return "Constant value (real) or lower bound, initial value, upper bound (3 reals)";
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    std::string Range_Validator< Design_Variable >::get_validation_message()
    {
        return Type_Validator< Design_Variable >().get_validation_message()
             + ", [" + std::to_string( mMinimumValue.get_value() ) + ", " + std::to_string( mMaximumValue.get_value() ) + "]";
    }

    //--------------------------------------------------------------------------------------------------------------
}
