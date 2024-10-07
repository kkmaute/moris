/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Module_Parameter_Lists.cpp
 *
 */

#include "cl_Module_Parameter_Lists.hpp"

namespace moris
{
    //--------------------------------------------------------------------------------------------------------------

    Module_Parameter_Lists::Module_Parameter_Lists( Parameter_List_Type aParameterListType )
            : mParameterListType( aParameterListType )
    {
        // Check for valid parameter list type before populating
        if ( mParameterListType != Parameter_List_Type::END_ENUM )
        {
            mSubModuleParameterLists.resize( get_number_of_sub_parameter_lists_in_module( aParameterListType ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Module_Parameter_Lists::size() const
    {
        return mSubModuleParameterLists.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Module_Parameter_Lists::clear()
    {
        mSubModuleParameterLists.clear();
        mParameterListType = Parameter_List_Type::END_ENUM;
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Module_Parameter_Lists::begin()->decltype( mSubModuleParameterLists.begin() )
    {
        return mSubModuleParameterLists.begin();
    }

    //--------------------------------------------------------------------------------------------------------------

    auto Module_Parameter_Lists::end()->decltype( mSubModuleParameterLists.end() )
    {
        return mSubModuleParameterLists.end();
    }

    //--------------------------------------------------------------------------------------------------------------
}
