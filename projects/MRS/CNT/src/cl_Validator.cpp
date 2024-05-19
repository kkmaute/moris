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

    std::string to_string( const std::string& aValue )
    {
        return "\"" + aValue + "\"";
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    bool Type_Validator< T >::make_valid_parameter( Variant& aVariant )
    {
        return aVariant.index() == variant_index< T >();
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    std::string Type_Validator< T >::get_valid_values()
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

    template< typename T >
    std::string Vector_Validator< T >::get_valid_values()
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
    std::string Range_Validator< T >::get_valid_values()
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
    bool Selection_Validator< T >::make_valid_parameter( Variant& aVariant )
    {
        return aVariant.index() == variant_index< T >() and
               mValidSelections.find( std::get< T >( aVariant ) ) != mValidSelections.end();
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename T >
    std::string Selection_Validator< T >::get_valid_values()
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
    std::string Type_Validator< Design_Variable >::get_valid_values()
    {
        return "Constant value (real) or lower bound, initial value, upper bound (3 reals)";
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    std::string Range_Validator< Design_Variable >::get_valid_values()
    {
        return Type_Validator< Design_Variable >().get_valid_values()
             + ", [" + std::to_string( mMinimumValue.get_value() ) + ", " + std::to_string( mMaximumValue.get_value() ) + "]";
    }

    //--------------------------------------------------------------------------------------------------------------
}
