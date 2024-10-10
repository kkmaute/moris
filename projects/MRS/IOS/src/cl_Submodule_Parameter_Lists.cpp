/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Submodule_Parameter_Lists.cpp
 *
 */

#include "cl_Submodule_Parameter_Lists.hpp"

#include <utility>

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Submodule_Parameter_Lists::Submodule_Parameter_Lists( std::string aType )
            : mType( std::move( aType ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::set_type( const std::string& aType )
    {
        mType = aType;
    }

    //--------------------------------------------------------------------------------------------------------------

    const std::string& Submodule_Parameter_Lists::get_type()
    {
        return mType;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Submodule_Parameter_Lists::size() const
    {
        return mParameterLists.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::add_parameter_list( const Parameter_List& aParameterList )
    {
        mParameterLists.push_back( aParameterList );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::append( const Submodule_Parameter_Lists& aSubmoduleToAppend )
    {
        mParameterLists.append( aSubmoduleToAppend.mParameterLists );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Submodule_Parameter_Lists::empty()
    {
        return mParameterLists.empty();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Submodule_Parameter_Lists::insert(
        const std::string&     aName,
        const Design_Variable& aDesignVariable )
    {
        // Insert into parameter list
        mParameterLists.back().insert( aName, aDesignVariable );
    }

    //--------------------------------------------------------------------------------------------------------------

    Parameter_List& Submodule_Parameter_Lists::operator()( uint aParameterListIndex )
    {
        return mParameterLists( aParameterListIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Parameter_List& Submodule_Parameter_Lists::operator()( uint aParameterListIndex ) const
    {
        return mParameterLists( aParameterListIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Submodule_Parameter_Lists::begin()->decltype( mParameterLists.begin() )
    {
        return mParameterLists.begin();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Submodule_Parameter_Lists::end()->decltype( mParameterLists.end() )
    {
        return mParameterLists.end();
    }

    //--------------------------------------------------------------------------------------------------------------
}
