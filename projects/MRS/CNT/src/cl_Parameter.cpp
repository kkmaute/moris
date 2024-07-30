/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Parameter.cpp
 *
 */

#include "cl_Parameter.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Parameter::Parameter( const Vector< std::string >& aEnumStrings )
            : mEntryType( Entry_Type::SELECTION )
    {
        // Set default value without validation
        uint tValue = 0;
        mValue      = make_variant( tValue );

        // Create selection validator
        mValidator = new Enum_Validator( aEnumStrings );
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter::Parameter( const Parameter& aParameter )
            : mValue( aParameter.mValue )
            , mEntryType( aParameter.mEntryType )
            , mNumberOfEntries( aParameter.mNumberOfEntries )
            , mExternalValidator( aParameter.mExternalValidator )
    {
        if ( aParameter.mValidator )
        {
            mValidator = aParameter.mValidator->copy();
        }
        else
        {
            mValidator = nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter::~Parameter()
    {
        delete mValidator;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Variant& Parameter::get_value() const
    {
        return mValue;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string Parameter::get_string() const
    {
        return convert_variant_to_string( mValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parameter::index() const
    {
        return mValue.index();
    }

    //--------------------------------------------------------------------------------------------------------------

    Entry_Type Parameter::get_entry_type() const
    {
        return mEntryType;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parameter::get_number_of_entries() const
    {
        return mNumberOfEntries;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Parameter::is_locked() const
    {
        return mValidator == nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Parameter::needs_linking() const
    {
        return mNeedsLinking;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< std::string >& Parameter::get_selection_names() const
    {
        return mValidator->get_selection_names();
    }

    //--------------------------------------------------------------------------------------------------------------

    const External_Validator& Parameter::get_external_validator() const
    {
        return mExternalValidator;
    }

    //--------------------------------------------------------------------------------------------------------------

    template<>
    Parameter::Parameter(
            const char*         aString,
            Entry_Type          aExternalValidationType,
            std::string         aExternalParameterName,
            Parameter_List_Type aExternalParameterListType,
            uint                aExternalParameterListIndex )
            : Parameter( std::string( aString ), aExternalValidationType, std::move( aExternalParameterName ), aExternalParameterListType, aExternalParameterListIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Parameter::operator==( const Parameter& aOther )
    {
        return this->mValue == aOther.mValue;
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris
