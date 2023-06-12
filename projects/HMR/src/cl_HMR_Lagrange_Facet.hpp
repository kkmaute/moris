/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Facet.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_FACET_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_FACET_HPP_

#include "cl_HMR_Facet.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

            /**
            * \brief Lagrange Element templated against
            *
            * uint N: number of dimensions (1, 2, or 3)
            * uint D: number of nodes
            */
           template< uint N, uint D >
           class Lagrange_Facet : public Facet
           {
               //! vertices of this face
               Basis * mVertices[ D ] = { nullptr };

//------------------------------------------------------------------------------
           public:
//------------------------------------------------------------------------------

               inline
               Lagrange_Facet( Mesh_Base        * aMesh,
                               Background_Facet * aBackgroundFacet ) : Facet( aMesh, aBackgroundFacet )
               {
                   // copy vertex pointers
                   this->copy_vertex_pointers( this->get_index_on_leader() );
               }

//------------------------------------------------------------------------------

               ~Lagrange_Facet(){};

//------------------------------------------------------------------------------
               inline
               uint get_number_of_vertices() const
               {
                   return D;
               }

//------------------------------------------------------------------------------
               inline
               Matrix< DDRMat > get_vertex_coords() const
               {
                   // create output matrix
                   Matrix< DDRMat > aCoords( D, N );

                   // loop over all basis
                   for( uint k=0; k<D; ++k )
                   {

                       // fixme: do this in one line
                       Matrix< DDRMat > tNodeCoords = mVertices[ k ]->get_coords();

                       // copy coords from vertex
                       aCoords.set_row( k, tNodeCoords );
                   }

                   return aCoords;
               }

//------------------------------------------------------------------------------

               mtk::Geometry_Type get_geometry_type() const;

//------------------------------------------------------------------------------

               mtk::Interpolation_Order get_interpolation_order() const;

//------------------------------------------------------------------------------

               mtk::Integration_Order get_integration_order() const;

//------------------------------------------------------------------------------
               inline
               const mtk::Vertex * get_vertex( uint aIndex ) const
               {
                   return mVertices[ aIndex ];
               }

// ----------------------------------------------------------------------------
               inline
               const Basis * get_basis( uint aIndex ) const
               {
                   return mVertices[ aIndex ];
               }

// ----------------------------------------------------------------------------
               inline
               Basis * get_basis( uint aIndex )
               {
                   return mVertices[ aIndex ];
               }

//------------------------------------------------------------------------------
           protected:
//------------------------------------------------------------------------------

               /**
                * internal function called by constructor
                */
               void copy_vertex_pointers( uint aIndex );

//------------------------------------------------------------------------------
           };

//------------------------------------------------------------------------------

           template< uint N, uint D >
           inline
           mtk::Geometry_Type Lagrange_Facet< N, D >::get_geometry_type() const
           {
               MORIS_ERROR( false,
                    "get_geometry_type() not implemented for this Lagrange_Facet.");
               return mtk::Geometry_Type::UNDEFINED;
           }

//------------------------------------------------------------------------------

           template< uint N, uint D >
           inline
           mtk::Interpolation_Order Lagrange_Facet< N, D >::get_interpolation_order() const
           {
               MORIS_ERROR( false,
                    "get_interpolation_order() not implemented for this Lagrange_Facet.");
               return mtk::Interpolation_Order::UNDEFINED;
           }

//------------------------------------------------------------------------------

           template< uint N, uint D >
           inline
           mtk::Integration_Order Lagrange_Facet< N, D >::get_integration_order() const
           {
               MORIS_ERROR( false,
                       "Lagrange_Facet::get_integration_order() not implemented.");
               return mtk::Integration_Order::UNDEFINED;
           }

//------------------------------------------------------------------------------

           template< uint N, uint D >
           inline
           void Lagrange_Facet< N, D >::copy_vertex_pointers( uint aIndex )
           {
               MORIS_ERROR( false,
                   "copy_vertex_pointers() not implemented for this Lagrange_Facet.");
           }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_FACET_HPP_ */
