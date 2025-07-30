/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Edge.hpp
 *
 */

#pragma once

#include "cl_HMR_Edge.hpp"
namespace moris::hmr
{
// ----------------------------------------------------------------------------

    template< uint D >
    class Lagrange_Edge : public Edge
    {
        Basis * mVertices[ D ] = { nullptr };

// ----------------------------------------------------------------------------
    public :
// ----------------------------------------------------------------------------
        inline
        Lagrange_Edge( Mesh_Base       * aMesh,
                       Background_Edge * aBackgroundEdge) : Edge( aMesh, aBackgroundEdge )
        {
            this->copy_vertex_pointers();
        }

// ----------------------------------------------------------------------------

        ~Lagrange_Edge() override{};

// ----------------------------------------------------------------------------
        inline
        uint get_number_of_vertices() const override
        {
            return D;
        }

// ----------------------------------------------------------------------------

        mtk::Interpolation_Order get_interpolation_order() const override;

// ----------------------------------------------------------------------------

        mtk::Integration_Order get_integration_order() const override;

// ----------------------------------------------------------------------------
        inline
        const Basis * get_basis( const uint aIndex ) const override
        {
            return mVertices[ aIndex ];
        }

// ----------------------------------------------------------------------------
        inline
        Basis * get_basis( const uint aIndex ) override
        {
            return mVertices[ aIndex ];
        }

// ----------------------------------------------------------------------------
    protected:
// ----------------------------------------------------------------------------
        void copy_vertex_pointers() override;

// ----------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------

    template< uint D >
    inline
    mtk::Interpolation_Order Lagrange_Edge< D >::get_interpolation_order() const
    {
        MORIS_ERROR( false, "get_interpolation_order() not implemented for this Lagrange_Edge.");
        return mtk::Interpolation_Order::UNDEFINED;
    }

//------------------------------------------------------------------------------

    template< uint D >
    inline
    mtk::Integration_Order Lagrange_Edge< D >::get_integration_order() const
    {
        MORIS_ERROR( false, "get_integration_order() not implemented for this Lagrange_Edge.");
        return mtk::Integration_Order::UNDEFINED;
    }

//------------------------------------------------------------------------------

    template< uint D >
    inline
    void Lagrange_Edge< D >::copy_vertex_pointers()
    {
        MORIS_ERROR( false, "copy_vertex_pointers() not implemented for this Lagrange_Edge.");
    }

//------------------------------------------------------------------------------
} /* namespace moris */
