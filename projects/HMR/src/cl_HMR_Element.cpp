/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Element.cpp
 *
 */

#include "cl_HMR_Element.hpp" //HMR/src

#include "cl_HMR_Edge.hpp"
#include "cl_HMR_Facet.hpp"

namespace moris::hmr
{
//------------------------------------------------------------------------------

    Element::Element(
            Background_Element_Base  * aElement,
            uint aActivationPattern )
    : mElement( aElement )
    , mActivationPattern( aActivationPattern )
    , mIndex(MORIS_INDEX_MAX)
    {
    }

//------------------------------------------------------------------------------

    Element * Element::get_neighbor(
            moris::Vector< Element * > & aAllElementsOnProc,
            luint aNeighborNumber )
    {
        // get neighbor of background element
        Background_Element_Base* tElement = mElement->get_neighbor( aNeighborNumber );

        // test if neighbor exists
        if ( tElement != nullptr )
        {
            // test if neighbor is on the same level
            if ( tElement->get_level() == mElement->get_level() )
            {
                return aAllElementsOnProc( tElement->get_memory_index() );
            }
            else
            {
                return nullptr;
            }
        }
        else
        {
            return nullptr;
        }
    }

//-------------------------------------------------------------------------------

    Element * Element::get_child(
            moris::Vector< Element * > & aAllElementsOnProc,
            uint aChildIndex )
    {
        if( mElement->has_children() )
        {
            return aAllElementsOnProc( mElement->get_child( aChildIndex )->get_memory_index() );
        }
        else
        {
            return nullptr;
        }
    }

//-------------------------------------------------------------------------------

    // special funciton for HMR
    Facet * Element::get_hmr_facet( uint aIndex )
    {
        MORIS_ERROR( false, "get_lagrange_facet() cannot be called from this element type" );
        return nullptr;
    }

//-------------------------------------------------------------------------------

    void Element::set_hmr_facet( Facet * aFacet, uint aIndex )
    {
        MORIS_ERROR( false, "set_hmr_facet() cannot be called from this element type" );
    }

//-------------------------------------------------------------------------------

    Edge * Element::get_hmr_edge( uint aIndex )
    {
        MORIS_ERROR( false, "get_hmr_edge() cannot be called from this element type" );
        return nullptr;
    }

//-------------------------------------------------------------------------------

    const Edge * Element::get_hmr_edge( uint aIndex ) const
    {
        MORIS_ERROR( false, "get_hmr_edge() cannot be called from this element type" );
        return nullptr;
    }

//-------------------------------------------------------------------------------

    void Element::set_hmr_edge( Edge * aEdge, uint aIndex )
    {
        MORIS_ERROR( false, "set_hmr_edge() cannot be called from this element type" );
    }

//-------------------------------------------------------------------------------
} /* namespace moris */
