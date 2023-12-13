/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HPP_

#include "cl_HMR_Background_Element.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"               //HMR/src
#include "cl_HMR_Facet.hpp"                 //HMR/src
#include "cl_HMR_Lagrange_Node.hpp"         //HMR/src
#include "typedefs.hpp"                     //COR/src
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"        //LINALG/src
#include "cl_MTK_Vertex.hpp"    //MTK/src

namespace moris::hmr
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
        // pointer to nodes
        Basis_Function** mNodes;

        // pointers to twin on B-Spline mesh
        moris::Cell< Element* > mTwins;

        // cell with facets
        moris::Cell< Facet* > mFacets;

        // cell with edges
        moris::Cell< Edge* > mEdges;

        // -----------------------------------------------------------------------------

      public:
        // -----------------------------------------------------------------------------

        /**
         * default Lagrange Element constructor
         */
        inline Lagrange_Element(
                Background_Element_Base* aElement,
                uint              aActivationPattern )
                : Element( aElement, aActivationPattern )
        {
            if ( N == 2 )
            {
                mFacets.resize( 4, nullptr );
            }
            else if ( N == 3 )
            {
                mFacets.resize( 6, nullptr );
                mEdges.resize( 12, nullptr );
            }

            this->set_cell_info();
        }

        //------------------------------------------------------------------------------

        /**
         * default destructor
         */
        ~Lagrange_Element()
        {
            this->delete_basis_container();
        }

        //------------------------------------------------------------------------------

        void
        set_cell_info();

        //------------------------------------------------------------------------------

        void
        init_basis_container()
        {
            MORIS_ASSERT( !mHaveBasis,
                    "Basis container of element already initiated" );

            mHaveBasis = true;
            mNodes     = new Basis_Function*[ D ];
            for ( uint k = 0; k < D; ++k )
            {
                mNodes[ k ] = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void
        delete_basis_container()
        {
            if ( mHaveBasis )
            {
                mHaveBasis = false;
                delete[] mNodes;
            }
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
        get_vertex_pointers() const
        {
            moris::Cell< mtk::Vertex* > aVertices( D );
            for ( uint k = 0; k < D; ++k )
            {
                aVertices( k ) = mNodes[ k ];
            }

            return aVertices;
        }

        //------------------------------------------------------------------------------

        // TODO MESHCLEANUP
        void
        remove_vertex_pointer( moris_index aIndex )
        {
            std::cout << "In HMR Lagrange Element" << std::endl;
        }

        //------------------------------------------------------------------------------

        /**
         * MTK Interface: returns a mat with the vertex IDs
         */
        Matrix< IdMat >
        get_vertex_ids() const
        {
            Matrix< IdMat > aIDs( D, 1 );
            for ( uint k = 0; k < D; ++k )
            {
                // the following line is correct
                aIDs( k ) = mNodes[ k ]->get_id();
            }

            return aIDs;
        }

        //------------------------------------------------------------------------------

        /**
         * MTK Interface: returns a mat with the vertex IDs
         */
        Matrix< IndexMat >
        get_vertex_inds() const
        {
            Matrix< IndexMat > aIndices( 1, D );    // FIXME was originally a column vector

            for ( uint k = 0; k < D; ++k )
            {
                // the following line is correct
                aIndices( k ) = mNodes[ k ]->get_index();
            }

            return aIndices;
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
                    (long unsigned int)mElement->get_hmr_index( mActivationPattern ),
                    (long unsigned int)mElement->get_hmr_id(),
                    (long unsigned int)mElement->get_parent()->get_hmr_id() );

            for ( uint k = 0; k < D; ++k )
            {
                // get node
                Basis_Function* tNode = this->get_basis_function( k );
                std::fprintf( stdout,
                        "    %2u :  Node %lu , ID %lu, MEM %lu \n",
                        (unsigned int)k,
                        (long unsigned int)tNode->get_hmr_index(),
                        (long unsigned int)tNode->get_hmr_id(),
                        (long unsigned int)tNode->get_memory_index() );
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
        Basis_Function*
        get_basis_function( uint aIndex )
        {
            if ( mHaveBasis )
            {
                MORIS_ASSERT( aIndex < D, "Try to access mNodes with index being out of bound.\n" );
                return mNodes[ aIndex ];
            }
            else
            {
                return nullptr;
            }
        }

        const Basis_Function*
        get_basis_function( uint aIndex ) const
        {
            if ( mHaveBasis )
            {
                MORIS_ASSERT( aIndex < D, "Try to access mNodes with index being out of bound.\n" );
                return mNodes[ aIndex ];
            }
            else
            {
                return nullptr;
            }
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
        insert_basis_function(
                uint aIndex,
                Basis_Function*      aBasis )
        {
            MORIS_ASSERT( aIndex < D, "Try to insert bases into mNodes with index being out of bound.\n" );
            mNodes[ aIndex ] = aBasis;
        }

        /**
         * Flags all bases in this element
         */
        void flag_all_bases() override
        {
            if ( mHaveBasis )
            {
                for ( uint iNodeIndex = 0; iNodeIndex < D; iNodeIndex++ )
                {
                    mNodes[ iNodeIndex ]->flag();
                }
            }
        }

        //------------------------------------------------------------------------------

        /**
         * string needed for gmsh output
         *
         * @return std::string
         *
         */
        std::string get_gmsh_string();

        //------------------------------------------------------------------------------

        /**
         * VTK ID needed for VTK output
         *
         * @return uint
         */
        uint get_vtk_type();

        //------------------------------------------------------------------------------

        /**
         * node IDs needed for VTK output
         *
         * @param[out] moris::Matrix< DDLUMat >
         *
         * @return void
         *
         */
        void get_basis_indices_for_vtk( Matrix< DDLUMat >& aNodes );

        //------------------------------------------------------------------------------

        /**
         * Creates all bases on the coarsest level.
         *
         * @param aAllElementsOnProc Cell containing all elements including the aura
         * @return Number of created bases
         */
        luint create_basis_on_level_zero(
                moris::Cell< Element* >& aAllElementsOnProc );

        //------------------------------------------------------------------------------

        /**
         * Creates bases for children of refined elements.
         *
         * @param aAllElementsOnProc Cell containing all elements including the aura
         * @return Number of created bases
         */
        luint create_basis_for_children(
                moris::Cell< Element* >& aAllElementsOnProc );

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
        void get_ijk_of_basis_function(
                uint aBasisNumber,
                luint*      aIJK );

        //------------------------------------------------------------------------------

        /**
         * reserve memory for twin container
         */
        void
        allocate_twin_container( const uint aSize )
        {
            mTwins.resize( aSize, nullptr );
        }

        //------------------------------------------------------------------------------

        /**
         * set twin on corresponding B-Spline mesh
         */
        void
        set_twin( const uint aIndex, Element* aTwin )
        {
            mTwins( aIndex ) = aTwin;
        }

        //------------------------------------------------------------------------------

        /**
         * returns a Mat with the node coords
         */
        Matrix< DDRMat > get_vertex_coords() const;

        //------------------------------------------------------------------------------

        Facet*
        get_hmr_facet( uint aIndex )
        {
            return mFacets( aIndex );
        }

        //------------------------------------------------------------------------------

        void
        set_hmr_facet(
                Facet*      aFacet,
                uint aIndex )
        {
            mFacets( aIndex ) = aFacet;
        }

        //------------------------------------------------------------------------------

        Edge*
        get_hmr_edge( uint aIndex )
        {
            return mEdges( aIndex );
        }

        //------------------------------------------------------------------------------

        const Edge*
        get_hmr_edge( uint aIndex ) const
        {
            return mEdges( aIndex );
        }

        //------------------------------------------------------------------------------

        void
        set_hmr_edge(
                Edge*       aEdge,
                uint aIndex )
        {
            mEdges( aIndex ) = aEdge;
        }

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
        create_basis( uint aBasisNumber )
        {
            // container for basis position
            luint tIJK[ N ];

            // ask element for position
            this->get_ijk_of_basis_function( aBasisNumber, tIJK );

            // create new Lagrange node
            mNodes[ aBasisNumber ] = new Lagrange_Node< N >(
                    tIJK,
                    mElement->get_level(),
                    mElement->get_owner() );
        }

        //------------------------------------------------------------------------------

    }; // class Lagrange_Element

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    inline void
    Lagrange_Element< N, D >::set_cell_info()
    {
        MORIS_ERROR( false, "set_cell_info() not available for this element." );
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    luint Lagrange_Element< N, D >::create_basis_on_level_zero(
            moris::Cell< Element* >& aAllElementsOnProc )
    {
        MORIS_ERROR( false, "Don't know how to create Lagrange nodes on level zero." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    luint Lagrange_Element< N, D >::create_basis_for_children(
            moris::Cell< Element* >& aAllElementsOnProc )
    {
        MORIS_ERROR( false, "Don't know how to create Lagrange nodes for children." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    inline std::string
    Lagrange_Element< N, D >::get_gmsh_string()
    {
        std::string aString = "GMSH not implemented for this element";
        return aString;
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    inline uint
    Lagrange_Element< N, D >::get_vtk_type()
    {
        // this element has no VTK id
        return 2;
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    void Lagrange_Element< N, D >::get_basis_indices_for_vtk( Matrix< DDLUMat >& aBasis )
    {
        // Loop over all bases
        for ( uint iBasisIndex = 0; iBasisIndex < D; iBasisIndex++ )
        {
            aBasis( iBasisIndex ) = mNodes[ iBasisIndex ]->get_memory_index();
        }
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    inline void
    Lagrange_Element< N, D >::get_ijk_of_basis_function(
            uint aBasisNumber,
            luint*      aIJK )
    {
        MORIS_ERROR( false, "Don't know how to get ijk of basis." );
    }

    //------------------------------------------------------------------------------

    template< uint N, uint D >
    inline Matrix< DDRMat >
    Lagrange_Element< N, D >::get_vertex_coords() const
    {
        Matrix< DDRMat > aCoords( D, N );
        for ( uint k = 0; k < D; ++k )
        {
            const real* tXYZ = mNodes[ k ]->get_xyz();

            for ( uint i = 0; i < N; ++i )
            {
                aCoords( k, i ) = tXYZ[ i ];
            }
        }
        return aCoords;
    }

    //------------------------------------------------------------------------------

} /* namespace moris */

//------------------------------------------------------------------------------

#include "cl_HMR_Lagrange_Element_Quad4.hpp" 
#include "cl_HMR_Lagrange_Element_Quad9.hpp" 
#include "cl_HMR_Lagrange_Element_Quad16.hpp"
#include "cl_HMR_Lagrange_Element_Hex8.hpp"  
#include "cl_HMR_Lagrange_Element_Hex27.hpp" 
#include "cl_HMR_Lagrange_Element_Hex64.hpp" 

//------------------------------------------------------------------------------

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HPP_ */
