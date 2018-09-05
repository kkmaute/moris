/*
 * cl_HMR_Lagrange_Element.hpp
 *
 *  Created on: May 24, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HPP_


#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_MTK_Vertex.hpp" //MTK/src

#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Lagrange_Node.hpp" //HMR/src
#include "cl_HMR_Background_Element.hpp" //HMR/src

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
        class Lagrange_Element : public Element
        {

            //! pointer to nodes
            Basis*     mNodes[ D ] = { nullptr };

            //! pointer to twin on B-Spline element
            Element* mTwin = nullptr;

// -----------------------------------------------------------------------------
        public:
// -----------------------------------------------------------------------------

            /**
             * default Lagrange Element constructor
             */
            Lagrange_Element( Background_Element_Base* aElement,
                              const uint & aActivationPattern ) :
                Element( aElement, aActivationPattern )

            {

            }

//------------------------------------------------------------------------------

            /**
             * default destructor
             */
            ~Lagrange_Element()
            {
            }

//------------------------------------------------------------------------------

            /**
             * MTK Interface: returns the number of vertices connected to this
             *                element
             */
            uint
            get_number_of_vertices() const
            {
                return D;
            }

//------------------------------------------------------------------------------

            /**
             * MTK Interface: returns a cell with the vertex pointers of this
             * element
             */
            moris::Cell< mtk::Vertex* >
            get_vertex_pointers()
            {
                moris::Cell< mtk::Vertex* > aVertices( D );
                for( uint k = 0; k<D; ++k )
                {
                    aVertices( k ) = mNodes[ k ];
                }

                return aVertices;
            }

//------------------------------------------------------------------------------

            /**
             * MTK Interface: returns a mat with the vertex IDs
             */
            Mat< luint >
            get_vertex_ids() const
            {
                Mat< luint > aIDs( D, 1 );
                for( uint k = 0; k<D; ++k )
                {
                    // the following line is correct
                    aIDs( k ) =  mNodes[ k ]->get_domain_index();
                }

                return aIDs;
            }

//------------------------------------------------------------------------------
            /**
             * for debugging
             *
             * @return void
             */
            void
            print_connectivity()
            {
                std::fprintf( stdout,
                        "connectivity of element %4lu ( ID %4lu, parent %4lu ):\n",
                        ( long unsigned int ) mElement->get_domain_index( mActivationPattern ),
                        ( long unsigned int ) mElement->get_domain_id(),
                        ( long unsigned int ) mElement->get_parent()->get_domain_id() );
                for( uint k=0; k<D; ++k )
                {
                    // get node
                    Basis* tNode = this->get_basis( k );
                    std::fprintf( stdout,
                            "    %2u :  Node %lu , ID %lu, MEM %lu \n",
                            ( unsigned int ) k,
                            ( long unsigned int ) tNode->get_domain_index(),
                            ( long unsigned int ) tNode->get_domain_id(),
                            ( long unsigned int ) tNode->get_memory_index());
                }
                std::fprintf( stdout, "\n" );
            }


//------------------------------------------------------------------------------

            /**
             * get pointer to node
             *
             * @param[in] aIndex              element local index of node
             *
             * @return    Basis* pointer to Lagrange node
             *
             */
            Basis*
            get_basis( const uint& aIndex )
            {
                return mNodes[ aIndex ];
            }

//------------------------------------------------------------------------------

            /**
             * set pointer of node to specified index and object
             *
             * @param[in]    aIndex  element local index of node
             * @param[in]    aBasis  pointer to Lagrange node
             *
             * @return void
             *
             */
            void
            insert_basis(
                    const  uint  & aIndex,
                    Basis        * aBasis )
            {
                mNodes[ aIndex ] = aBasis;
            }

//------------------------------------------------------------------------------

            /**
             * string needed for gmsh output
             *
             * @return std::string
             *
             */
           std::string
           get_gmsh_string();

//------------------------------------------------------------------------------

           /**
            * VTK ID needed for VTK output
            *
            * @return uint
            */
           uint
           get_vtk_type();

//------------------------------------------------------------------------------

           /**
            * node IDs needed for VTK output
            *
            * @param[out] moris::Mat<luint>
            *
            * @return void
            *
            */
           void
           get_basis_indices_for_vtk( Mat<luint> & aNodes );

//------------------------------------------------------------------------------

           /**
            * Creates all nodes on the coarsest level.
            * Called by Lagrange mesh create_basis_on_level_zero().
            *
            * @param[inout] aAllElementsOnProc   cell containing all Lagrange
            *                                    elements including the aura
            * @param[inout] aNodeCounter         counter to keep track of
            *                                    how many nodes were generated
            * @return void
            */
           void
           create_basis_on_level_zero(
                   moris::Cell< Element * > & aAllElementsOnProc,
                   luint             & aNodeCounter );

//------------------------------------------------------------------------------

           /**
            * Creates nodes for children of refined elements.
            * Called by Lagrange mesh.
            *
            * @param[inout] aAllElementsOnProc   cell containing all Lagrange
            *                                    elements including the aura
            * @param[inout] aBasisCounter        counter to keep track of
            *                                    how many nodes were generated
            * @return void
            */
           void
           create_basis_for_children(
                   moris::Cell< Element * > & aAllElementsOnProc,
                       luint             & aBasisCounter );

//------------------------------------------------------------------------------

           /**
            * returns the ijk position of a given basis
            *
            * @param[in]  aBasisNumber   element local number of basis
            * @param[out] aIJK           proc local ijk position of this basis
            *
            * @return void
            *
            */
           void
           get_ijk_of_basis(
                   const uint & aBasisNumber,
                   luint      * aIJK );

//------------------------------------------------------------------------------

           /**
            * set twin on corresponding B-Spline mesh
            */
           void
           set_twin( Element * aTwin )
           {
               mTwin = aTwin;
           }

//------------------------------------------------------------------------------

           /**
            * returns the mtk geometry type of this element
            */
           mtk::Geometry_Type
           get_geometry_type() const ;

//------------------------------------------------------------------------------

           /**
            * returns a Mat with the node coords
            */
           Mat< real >
           get_vertex_coords() const;

//------------------------------------------------------------------------------

           /**
            * returns the interpolation order of this element
            */
           mtk::Interpolation_Order
           get_interpolation_order() const;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

           /**
            * create new node at position
            *
            * @param[in]    aNodeNumber   element local index of new node
            *
            * @return void
            */
           void
           create_basis( const uint & aBasisNumber )
           {
               // container for basis position
               luint tIJK[ N ];

               // ask element for position
               this->get_ijk_of_basis( aBasisNumber, tIJK );

               // create new Lagrange node
               mNodes[ aBasisNumber ]
                       = new Lagrange_Node< N >(
                            tIJK,
                            mElement->get_level(),
                            mElement->get_owner() );
           }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        template< uint N, uint D >
        void
        Lagrange_Element< N, D >::create_basis_on_level_zero(
                moris::Cell< Element * > & aAllElementsOnProc,
                luint             & aBasisCounter )
        {
            MORIS_ERROR( false, "Don't know how to create Lagrange nodes on level zero.");
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        void
        Lagrange_Element< N, D >::create_basis_for_children(
                moris::Cell< Element * > & aAllElementsOnProc,
                luint             & aBasisCounter )
        {
            MORIS_ERROR( false, "Don't know how to create Lagrange nodes for children.");
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        std::string
        Lagrange_Element< N, D >::get_gmsh_string()
        {
            std::string aString = "GMSH not implemented for this element";
            return aString;
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        uint
        Lagrange_Element< N, D >::get_vtk_type()
        {
            // this element has no VTK id
            return 2;
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        void
        Lagrange_Element< N, D >::get_basis_indices_for_vtk( Mat<luint> & aBasis )
        {
            // do nothing
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        void
        Lagrange_Element< N, D >::get_ijk_of_basis(
                const uint & aBasisNumber,
                luint      * aIJK )
        {
            MORIS_ERROR( false, "Don't know how to get ijk of basis.");
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        Mat< real >
        Lagrange_Element< N, D >::get_vertex_coords() const
        {
            Mat<real> aCoords( D, N );
            for( uint k=0; k<D; ++k )
            {
                const real * tXYZ = mNodes[ k ]->get_xyz();

                for( uint i=0; i<N; ++i )
                {
                    aCoords( k, i ) = tXYZ[ i ];
                }
            }
            return aCoords;
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        mtk::Geometry_Type
        Lagrange_Element< N, D >::get_geometry_type() const
        {
            MORIS_ERROR( false, "get_geometry_type() not available for this element.");
            return mtk::Geometry_Type::UNDEFINED;
        }

//------------------------------------------------------------------------------

        template< uint N, uint D >
        mtk::Interpolation_Order
        Lagrange_Element< N, D >::get_interpolation_order() const
        {
            MORIS_ERROR( false, "get_interpolation_order() not available for this element.");
            return mtk::Interpolation_Order::UNDEFINED;
        }

//------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */



#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HPP_ */
