/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Edge.cpp
 *
 */

#include "cl_HMR_Background_Edge.hpp"

#include "cl_HMR_Background_Element_Base.hpp"

namespace moris::hmr
{
//-------------------------------------------------------------------------------

    Background_Edge::Background_Edge(
            Background_Element_Base * aElement,
            uint                      aIndex )
    {
        this->insert_element( aElement, aIndex );
    }

//-------------------------------------------------------------------------------
    void Background_Edge::flag()
    {
        mFlag = true;
    }

//-------------------------------------------------------------------------------

    void Background_Edge::unflag()
    {
        mFlag = false;
    }

//-------------------------------------------------------------------------------

    bool Background_Edge::is_flagged() const
    {
        return mFlag;
    }

//-------------------------------------------------------------------------------

    void Background_Edge::insert_element(
            Background_Element_Base * aElement,
            uint                      aIndex )
    {
        MORIS_ASSERT( mElementCounter < 4, "Error in element counter" );
        mIndexInElement[ mElementCounter ] = aIndex;
        mElements[ mElementCounter ] = aElement;
        ++mElementCounter;
    }

//-------------------------------------------------------------------------------

    uint Background_Edge::get_number_of_elements() const
    {
        return mElementCounter;
    }

//-------------------------------------------------------------------------------

    Background_Element_Base * Background_Edge::get_element( uint aIndex )
    {
        return mElements[ aIndex ];
    }

//-------------------------------------------------------------------------------

    uint Background_Edge::get_index_on_element( uint aIndex ) const
    {
        return mIndexInElement[ aIndex ];
    }

//-------------------------------------------------------------------------------
} /* namespace moris */

