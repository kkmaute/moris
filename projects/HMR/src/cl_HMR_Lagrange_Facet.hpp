/*
 * cl_HMR_Lagrange_Facet.hpp
 *
 *  Created on: Sep 24, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_FACET_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_FACET_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_trans.hpp"
#include "cl_HMR_Facet.hpp"

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

               Lagrange_Facet( Mesh_Base * aMesh,
                               Background_Facet * aBackgroundFacet ) :
                                   Facet( aMesh, aBackgroundFacet )
               {
                   // copy vertex pointers
                   this->copy_vertex_pointers(
                           aBackgroundFacet->get_index_on_master() );
               }

//------------------------------------------------------------------------------

               ~Lagrange_Facet(){};

//------------------------------------------------------------------------------

               uint
               get_number_of_vertices() const
               {
                   return D;
               }

//------------------------------------------------------------------------------

               Matrix< DDRMat >
               get_vertex_coords() const
               {
                   // create output matrix
                   Matrix< DDRMat > aCoords( D, N );

                   // loop over all basis
                   for( uint k=0; k<D; ++k )
                   {
                       // copy coords from vertex
                       aCoords.set_row( k,
                               trans( mVertices[ k ]->get_coords() ) );
                   }

                   return aCoords;
               }

//------------------------------------------------------------------------------

               mtk::Geometry_Type
               get_geometry_type() const;

//------------------------------------------------------------------------------

               mtk::Interpolation_Order
               get_interpolation_order() const;

//------------------------------------------------------------------------------

               const mtk::Vertex *
               get_vertex( const uint & aIndex ) const
               {
                   return mVertices[ aIndex ];
               }



// ----------------------------------------------------------------------------

               const Basis *
               get_basis( const uint & aIndex ) const
               {
                   return mVertices[ aIndex ];
               }

//------------------------------------------------------------------------------
           protected:
//------------------------------------------------------------------------------

               /**
                * internal function called by constructor
                */
               void
               copy_vertex_pointers( const uint & aIndex );

//------------------------------------------------------------------------------
           };

//------------------------------------------------------------------------------

           template< uint N, uint D >
           mtk::Geometry_Type
           Lagrange_Facet< N, D >::get_geometry_type() const
           {
               MORIS_ERROR( false,
                    "get_geometry_type() not implemented for this Lagrange_Facet.");
               return mtk::Geometry_Type::UNDEFINED;
           }

//------------------------------------------------------------------------------

           template< uint N, uint D >
           mtk::Interpolation_Order
           Lagrange_Facet< N, D >::get_interpolation_order() const
           {
               MORIS_ERROR( false,
                    "get_interpolation_order() not implemented for this Lagrange_Facet.");
               return mtk::Interpolation_Order::UNDEFINED;
           }

//------------------------------------------------------------------------------

           template< uint N, uint D >
           void
           Lagrange_Facet< N, D >::copy_vertex_pointers( const uint & aIndex )
           {
               MORIS_ERROR( false,
                   "copy_vertex_pointers() not implemented for this Lagrange_Facet.");
           }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_FACET_HPP_ */
