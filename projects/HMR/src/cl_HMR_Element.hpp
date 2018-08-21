/*
 * cl_Element.hpp
 *
 *  Created on: May 24, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_ELEMENT_HPP_
#define SRC_HMR_CL_HMR_ELEMENT_HPP_

#include <string>
#include "typedefs.hpp" //COR/src

#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_HMR_Background_Element.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        // forward declaration of node base type
        class Basis;

//------------------------------------------------------------------------------
        class Element : public mtk::Cell
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! pointer to element on background mesh
            Background_Element_Base* mElement;

            //! flag that tells if nodes of children have been processed
            bool                     mChildrenBasisFlag = false;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor for Lagrange element base class
             *
             * @param[in]   aElement   element on background mesh
             */
            Element(
                    Background_Element_Base*  aElement );

// -----------------------------------------------------------------------------

            /**
             * Virtual destructor. Does nothing.
             */
            virtual ~Element(){};


//------------------------------------------------------------------------------

            /**
             * MTK Interface: returns proc owner of element
             *
             * @return uint
             */
            uint
            get_owner() const
            {
                return mElement->get_owner();
            }


//------------------------------------------------------------------------------

            /**
             * MTK Interface: returns the id of the element
             */
            luint
            get_id() const
            {
                return mElement->get_domain_index(); // <-- this is correct
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
                return mElement;
            }
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

            /**
             * tells if an element is active
             *
             * @return bool   true if active
             */
            bool
            is_active() const
            {
                return mElement->is_active();
            }

//------------------------------------------------------------------------------

            /**
             * returns the memory index of the background element
             *
             * @return luint
             */
            auto
            get_memory_index() const  -> decltype( mElement->get_memory_index() )
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
            get_domain_index() const  -> decltype( mElement->get_domain_index() )
            {
                return mElement->get_domain_index();
            }

//------------------------------------------------------------------------------

            /**
             * returns a unique system wide ID of the element
             *
             * @return    luint global ID of element
             */
            auto
            get_domain_id() const -> decltype( mElement->get_domain_id() )
            {
                return mElement->get_domain_id();
            }

//------------------------------------------------------------------------------

            /**
             * returns the level of the background element
             *
             * @return uint
             */
            auto
            get_level() const  -> decltype( mElement->get_level() )
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
            virtual void
            print_connectivity() = 0;

//------------------------------------------------------------------------------


            /**
             * get pointer to node
             *
             * @param[in]    aIndex   element local index of node
             *
             * @return Basis* pointer to Lagrange node
             *
             */
            virtual Basis*
            get_basis( const uint& aIndex ) = 0;

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
            virtual void
            insert_basis(
                    const uint & aIndex,
                    Basis      * aBasis ) = 0;

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
            children_have_basis() const -> decltype ( mChildrenBasisFlag )
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
            virtual std::string
            get_gmsh_string() = 0;
//------------------------------------------------------------------------------

            /**
             * VTK ID needed for VTK output
             *
             * @return uint
             */
            virtual uint
            get_vtk_type() = 0;

//------------------------------------------------------------------------------

            /**
             * node IDs needed for VTK output
             *
             * @param[out] moris::Mat<luint>
             *
             * @return void
             *
             */
            virtual void
            get_basis_indices_for_vtk( Mat<luint> & aBasis ) = 0;

//------------------------------------------------------------------------------

            /**
             * Creates all nodes on the coarsest level.
             * Called by Lagrange mesh create_basis_on_level_zero().
             *
             * @param[inout] aAllElementsOnProc   cell containing all Lagrange
             *                                    elements including the aura
             * @param[inout] aBasisCounter         counter to keep track of
             *                                    how many nodes were generated
             * @return void
             */
            virtual void
            create_basis_on_level_zero(
                    moris::Cell< Element * > & aAllElementsOnProc,
                    luint                    & aBasisCounter ) = 0;

//------------------------------------------------------------------------------

            /**
             * Creates nodes for children of refined elements.
             * Called by Lagrange mesh create_basiss_on_higher_levels().
             *
             * @param[inout] aAllElementsOnProc   cell containing all Lagrange
             *                                    elements including the aura
             * @param[inout] aNodeCounter         counter to keep track of
             *                                    how many nodes were generated
             * @return void
             */
            virtual void
            create_basis_for_children(
                    moris::Cell< Element * > & aAllElementsOnProc,
                    luint                    & aBasisCounter ) = 0;

//------------------------------------------------------------------------------

            /**
             * returns a neighbor if it exists and is on the same level
             *
             * @param[in] aAllElementsOnProc   cell containing all Lagrange
             *                                 elements including the aura
             *
             * @param[in] aNeighborNumber      desired neighbor of element
             */
            Element * get_neighbor(
                    moris::Cell< Element * > & aAllElementsOnProc,
                    const luint              & aNeighborNumber );

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
            virtual void
            get_ijk_of_basis(
                        const uint & aBasisNumber,
                        luint      * aIJK ) = 0;


//------------------------------------------------------------------------------

            /**
             * Links each basis of an element with neighbor basis.
             *
             * @param[inout] aAllElementsOnProc   cell containing all Bspline
             *                                    elements including the aura
             * @return void
             */
            virtual void
            link_basis_with_neighbors(
                    moris::Cell< Element* > & aAllElementsOnProc )
            {
                MORIS_ERROR( false, "Link basis with neighbors not available for this element.");
            }


//------------------------------------------------------------------------------

            virtual void
            refine(
                    moris::Cell< Element* > & aAllElementsOnProc,
                    luint            & aBasisCounter )
            {
                MORIS_ERROR( false, "refine() not available for this element.");
            }

//------------------------------------------------------------------------------

            virtual void
            set_twin( Element* aTwin )
            {
                MORIS_ERROR( false, "set_twin() not available for this element.");
            }

//------------------------------------------------------------------------------

            virtual moris::Cell< mtk::Vertex* >
            get_vertex_pointers()
            {
                MORIS_ERROR( false, "get_vertex_pointers() not available for this element.");
                return moris::Cell< mtk::Vertex* >(0);
            }

//------------------------------------------------------------------------------

            virtual Mat<real>
            get_vertex_coords() const
            {
                MORIS_ERROR( false, "get_vertex_coords() not available for this element.");
                return Mat<real>(0,0);
            }

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
            virtual void
            create_basis( const uint & aBasisNumber ) = 0;

//------------------------------------------------------------------------------

            /**
             * refines the basis of an element
             *
             * @param[in] aIndex    Index of basis to refine
             * @return void
             */
            virtual void
            refine_basis( const uint & aBasisNumber, luint & aBasisCounter )
            {
                MORIS_ERROR( false, "refine_basis() not available for this element.");
            }

//------------------------------------------------------------------------------

            /**
             * tells how many nodes are connected to this element
             */
            uint
            get_number_of_vertices() const
            {
                MORIS_ERROR( false, " get_number_of_vertices() not available for this element.");
                return 0;
            }

//------------------------------------------------------------------------------

            /**
             * returns a vector with the ids ( here: domain indices) of the
             * nodes connected to the element
             */
            Mat< luint >
            get_vertex_ids() const
            {
                MORIS_ERROR( false, "get_vertex_ids() not available for this element.");
                return Mat< luint > (0,0);
            }

//------------------------------------------------------------------------------

            virtual mtk::Geometry_Type
            get_geometry_type() const
            {
                MORIS_ERROR( false, "get_geometry_type() not available for this element.");
                return mtk::Geometry_Type::UNDEFINED;
            }

//------------------------------------------------------------------------------

            virtual mtk::Interpolation_Order
            get_interpolation_order() const
            {
                MORIS_ERROR( false, "get_interpolation_order() not available for this element.");
                return mtk::Interpolation_Order::UNDEFINED;
            }

//------------------------------------------------------------------------------

            /**
             * set the T-Matrix flag
             */
            void
            set_t_matrix_flag();


//-------------------------------------------------------------------------------

            /**
             * unset the T-Matrix flag
             */
            void
            unset_t_matrix_flag();

//-------------------------------------------------------------------------------

            /**
             * query the T-Matrix flag
             */
            bool
            get_t_matrix_flag() const ;

//-------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_ELEMENT_HPP_ */
