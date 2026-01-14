/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Element.hpp
 *
 */

#pragma once

#include <string>

#include "cl_HMR_Background_Element.hpp"    //HMR/src
#include "moris_typedefs.hpp"               //COR/src

#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"    //MTK/src

namespace moris::hmr
{
    //------------------------------------------------------------------------------

    // forward declaration of node base type
    class Basis;

    // forward declaration of facet type
    class Facet;
    class Edge;

    //------------------------------------------------------------------------------
    class Element : public mtk::Cell
    {
        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------

        // pointer to element on background mesh
        Background_Element_Base* mElement;

        // flag that tells if nodes of children have been processed
        bool mChildrenBasisFlag = false;
        bool mHaveBasis         = false;

        const uint mActivationPattern;

        // Proc local index
        moris_index mIndex;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * constructor for Lagrange element base class
         *
         * @param[in]   aElement   element on background mesh
         */
        Element(
                Background_Element_Base* aElement,
                uint                     aActivationPattern );

        // -----------------------------------------------------------------------------

        /**
         * Virtual destructor. Does nothing.
         */
        ~Element() override{};

        //------------------------------------------------------------------------------

        bool
        have_basis() const
        {
            return mHaveBasis;
        }

        //------------------------------------------------------------------------------

        /**
         * MTK Interface: returns proc owner of element
         *
         * @return uint
         */
        moris_id
        get_owner() const override
        {
            return mElement->get_owner();
        }

        //------------------------------------------------------------------------------

        /**
         * MTK Interface: returns the id of the element
         */
        moris_id
        get_id() const override
        {
            return mElement->get_hmr_index( mActivationPattern ) + 1;    // <-- this is correct
        }

        //------------------------------------------------------------------------------

        /**
         * sets the index of this element
         */
        void
        set_index( const moris_index aIndex ) override
        {
            mIndex = aIndex;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the index of this element
         */
        moris_index
        get_index() const override
        {
            return mIndex;
        }

        //------------------------------------------------------------------------------

        /**
         * returns a pointer to the background element
         *
         * @return Element_Base* Element on background mesh
         */
        Background_Element_Base*
        get_background_element()
        {
            MORIS_ASSERT( mElement != nullptr, "Element::get_background_element(), Lagrange or B-spline element returns nullptr background element" );

            return mElement;
        }
        //------------------------------------------------------------------------------

        /**
         * returns a pointer to the background element ( const version )
         *
         * @return Element_Base* Element on background mesh
         */
        const Background_Element_Base*
        get_background_element() const
        {
            MORIS_ASSERT( mElement != nullptr, "Element::get_background_element(), Lagrange or B-spline element returne nullptr background element" );

            return mElement;
        }

        //------------------------------------------------------------------------------

        /**
         * special function for HMR
         */
        luint
        get_memory_index_of_background_element() const
        {
            return mElement->get_memory_index();
        }

        //------------------------------------------------------------------------------

        virtual Facet* get_hmr_facet( uint aIndex );

        //------------------------------------------------------------------------------

        virtual void set_hmr_facet( Facet* aFacet, uint aIndex );

        //------------------------------------------------------------------------------

        virtual Edge* get_hmr_edge( uint aIndex );

        //------------------------------------------------------------------------------

        virtual const Edge* get_hmr_edge( uint aIndex ) const;

        //------------------------------------------------------------------------------

        virtual void set_hmr_edge( Edge* aEdge, uint aIndex );

        //------------------------------------------------------------------------------

        /**
         * tells if an element is active
         *
         * @param[in]     aPattern   pattern this question refers to
         * @return bool   true if active
         */
        bool
        is_active() const
        {
            return mElement->is_active( mActivationPattern );
        }

        //------------------------------------------------------------------------------

        /**
         * tells if an element is refined
         *
         * @param[in]     aPattern   pattern this question refers to
         * @return bool   true if active
         */
        bool
        is_refined() const
        {
            return mElement->is_refined( mActivationPattern );
        }

        //------------------------------------------------------------------------------

        /**
         * tells if an element is neither active nor refined
         *
         * @param[in]     aPattern   pattern this question refers to
         * @return bool   true if element is irrelevant (i.e. neither active nor refined)
         */
        bool
        is_neither_active_nor_refined() const
        {
            return mElement->is_neither_active_nor_refined( mActivationPattern );
        }

        //------------------------------------------------------------------------------

        // /**
        //  * tells if an element is not active (irrespective of refinement status)
        //  *
        //  * @param[in]     aPattern   pattern this question refers to
        //  * @return bool   true if element is deactivated
        //  */
        // bool
        // is_deactivated() const
        // {
        //     bool tElementIsActive = mElement->is_active( mActivationPattern );
        //     return !tElementIsActive;
        // }

        //------------------------------------------------------------------------------

        /**
         * returns the memory index of the background element
         *
         * @return luint
         */
        auto
        get_memory_index() const -> decltype( mElement->get_memory_index() )
        {
            return mElement->get_memory_index();
        }

        //------------------------------------------------------------------------------

        /**
         * returns the domain index of the background element
         *
         * @return luint
         */
        auto
        get_hmr_index() const
                -> decltype( mElement->get_hmr_index( mActivationPattern ) )
        {
            return mElement->get_hmr_index( mActivationPattern );
        }

        //------------------------------------------------------------------------------

        /**
         * returns a unique system wide ID of the element
         *
         * @return    luint global ID of element
         */
        auto
        get_hmr_id() const -> decltype( mElement->get_hmr_id() )
        {
            return mElement->get_hmr_id();
        }

        //------------------------------------------------------------------------------

        /**
         * returns the level of the background element
         *
         * @return uint
         */
        auto
        get_level() const -> decltype( mElement->get_level() ) override
        {
            return mElement->get_level();
        }

        //------------------------------------------------------------------------------

        /**
         * tells if element is a padding element
         *
         * @return bool
         */
        bool
        is_padding() const
        {
            return mElement->is_padding();
        }

        /**
         * for debugging
         *
         * @return void
         */
        virtual void print_connectivity() = 0;

        //------------------------------------------------------------------------------

        /**
         * get pointer to node
         *
         * @param[in]    aIndex   element local index of node
         *
         * @return Basis* pointer to Lagrange node
         *
         */
        virtual Basis* get_basis( uint aIndex ) = 0;

        virtual const Basis* get_basis( uint aIndex ) const = 0;

        //------------------------------------------------------------------------------

        /**
         * set pointer of node to specified index and object
         *
         * @param[in]    aIndex  element local index of node
         * @param[in]    aNode   pointer to Lagrange node
         *
         * @return void
         *
         */
        virtual void insert_basis(
                uint   aIndex,
                Basis* aBasis ) = 0;

        //------------------------------------------------------------------------------

        /**
         * sets the flag telling that nodes of children have been calculated
         *
         * @return void
         */

        void
        set_children_basis_flag()
        {
            mChildrenBasisFlag = true;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the falue of the children nodes flag
         */
        auto
        children_have_basis() const -> decltype( mChildrenBasisFlag )
        {
            return mChildrenBasisFlag;
        }

        //------------------------------------------------------------------------------

        /**
         * string needed for gmsh output
         *
         * @return std::string
         *
         */
        virtual std::string get_gmsh_string() = 0;
        //------------------------------------------------------------------------------

        /**
         * VTK ID needed for VTK output
         *
         * @return uint
         */
        virtual uint get_vtk_type() = 0;

        //------------------------------------------------------------------------------

        /**
         * node IDs needed for VTK output
         *
         * @param[out] moris::Matrix< DDLUMat >
         *
         * @return void
         *
         */
        virtual void get_basis_indices_for_vtk( Matrix< DDLUMat >& aBasis ) = 0;

        //------------------------------------------------------------------------------

        /**
         * Creates all bases on the coarsest level.
         *
         * @param aAllElementsOnProc Cell containing all elements including aura
         * @return Number of created bases
         */
        virtual luint create_basis_on_level_zero(
                Vector< Element* >& aAllElementsOnProc ) = 0;

        //------------------------------------------------------------------------------

        /**
         * Creates bases for children of refined elements.
         *
         * @param aAllElementsOnProc Cell containing all elements including aura
         * @return Number of created bases
         */
        virtual luint create_basis_for_children(
                Vector< Element* >& aAllElementsOnProc ) = 0;

        //------------------------------------------------------------------------------

        /**
         * returns a neighbor if it exists and is on the same level
         *
         * @param[in] aAllElementsOnProc   cell containing all Lagrange
         *                             elements including the aura
         *
         * @param[in] aNeighborNumber    desired neighbor of element
         */
        Element* get_neighbor(
                Vector< Element* >& aAllElementsOnProc,
                luint               aNeighborNumber );

        //------------------------------------------------------------------------------

        /**
         * returns a child if it exists
         */
        Element* get_child(
                Vector< Element* >& aAllElementsOnProc,
                uint                aChildIndex );

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
        virtual void get_ijk_of_basis(
                uint   aBasisNumber,
                luint* aIJK ) = 0;

        //------------------------------------------------------------------------------

        /**
         * Links each basis of an element with neighbor basis.
         *
         * @param[inout] aAllElementsOnProc   cell containing all Bspline
         *                                elements including the aura
         * @return void
         */
        virtual void
        link_basis_with_neighbors( Vector< Element* >& aAllElementsOnProc )
        {
            MORIS_ERROR( false, "Link basis with neighbors not available for this element." );
        }

        //------------------------------------------------------------------------------

        virtual luint
        refine( Vector< Element* >& aAllElementsOnProc )
        {
            MORIS_ERROR( false, "refine() not available for this element." );
            return 0;
        }

        //------------------------------------------------------------------------------

        virtual void
        allocate_twin_container( const uint aSize )
        {
            MORIS_ERROR( false, "allocate_twin_container() not available for this element." );
        }

        //------------------------------------------------------------------------------

        virtual void
        set_twin( const uint aIndex, Element* aTwin )
        {
            MORIS_ERROR( false, "set_twin() not available for this element." );
        }

        //------------------------------------------------------------------------------

        Vector< mtk::Vertex* >
        get_vertex_pointers() const override
        {
            MORIS_ERROR( false, "get_vertex_pointers() not available for this element." );
            return Vector< mtk::Vertex* >( 0 );
        }

        //------------------------------------------------------------------------------

        // TODO MESHCLEANUP
        void
        remove_vertex_pointer( moris_index aIndex ) override
        {
            std::cout << "In HMR Element" << '\n';
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        get_vertex_coords() const override
        {
            MORIS_ERROR( false, "get_vertex_coords() not available for this element." );
            return Matrix< DDRMat >( 0, 0 );
        }

        //------------------------------------------------------------------------------

        /**
         * tells how many nodes are connected to this element
         */
        uint
        get_number_of_vertices() const override
        {
            MORIS_ERROR( false, " get_number_of_vertices() not available for this element." );
            return 0;
        }

        //------------------------------------------------------------------------------

        /**
         * returns a vector with the ids ( here: domain indices) of the
         * nodes connected to the element
         */
        Matrix< IdMat >
        get_vertex_ids() const override
        {
            MORIS_ERROR( false, "get_vertex_ids() const not available for this element." );
            return Matrix< IdMat >( 0, 0 );
        }

        //------------------------------------------------------------------------------

        /**
         * returns a vector with the ids ( here: domain indices) of the
         * nodes connected to the element
         */
        Matrix< IndexMat >
        get_vertex_inds() const override
        {
            MORIS_ERROR( false, "get_vertex_indices() const not available for this element." );
            return Matrix< IndexMat >( 0, 0 );
        }

        //------------------------------------------------------------------------------

        /**
         * update the min refinement level to this value if it is bigger
         * than the one stored in the memory
         */
        void
        update_min_refinement_level( uint aMinRefinementLevel )
        {
            mElement->update_min_refinement_level( aMinRefinementLevel );
        }

        //------------------------------------------------------------------------------

        /**
         * set min refinement level to this value
         */
        void
        set_min_refinement_level( uint aMinRefinementLevel )
        {
            mElement->set_min_refinement_level( aMinRefinementLevel );
        }

        //------------------------------------------------------------------------------

        /**
         * return the min refinement level to this value
         */
        uint
        get_min_refinement_level() const
        {
            return mElement->get_min_refinement_level();
        }

        //-------------------------------------------------------------------------------

        virtual void
        init_basis_container()
        {
            MORIS_ERROR( false, "init_basis_container() not available for this element." );
        }

        //-------------------------------------------------------------------------------

        virtual void
        delete_basis_container()
        {
            MORIS_ERROR( false, "delete_basis_container() not available for this element." );
        }

        //-------------------------------------------------------------------------------

        const luint*
        get_ijk() const override
        {
            return mElement->get_ijk();
        }

        /**
         * Flags all bases in this element
         */
        virtual void flag_all_bases() = 0;

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------

        /**
         * Creates a new node at specified position.
         * Protected function called by create_basis_on_level_zero.
         *
         * @param[in]    aBasisNumber   element local index of new node
         *
         * @return void
         */
        virtual void create_basis( uint aBasisNumber ) = 0;

        //------------------------------------------------------------------------------

        /**
         * Refines a basis of this element
         *
         * @param aBasisNumber Index of the basis to refine
         * @return Number of created bases
         */
        virtual luint
        refine_basis( uint aBasisNumber )
        {
            MORIS_ERROR( false, "refine_basis() not available for this element." );
            return 0;
        }

        //-------------------------------------------------------------------------------

    };    // class Element

    //------------------------------------------------------------------------------

}    // namespace moris::hmr
