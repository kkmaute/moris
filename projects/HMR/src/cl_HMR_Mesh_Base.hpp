/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh_Base.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_MESH_BASE_HPP_
#define SRC_HMR_CL_HMR_MESH_BASE_HPP_

#include "cl_HMR_Background_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Basis.hpp"                   //HMR/src
#include "cl_HMR_Element.hpp"                 //HMR/src
#include "cl_HMR_Parameters.hpp"              //HMR/src
#include "HMR_Globals.hpp"                    //HMR/src
#include "HMR_Tools.hpp"
#include "typedefs.hpp"    //COR/src

namespace moris
{
    namespace hmr
    {
        //------------------------------------------------------------------------------

        /**
         * \brief Base class for Lagrange and BSpline Mesh
         */
        class Mesh_Base
        {

            // ----------------------------------------------------------------------------

          protected:
            // ----------------------------------------------------------------------------

            //! ref to container of user defined settings
            const Parameters* mParameters;

            //! pointer to background mesh
            Background_Mesh_Base* mBackgroundMesh;

            //! polynomial degree of mesh
            const uint mOrder;

            //! dimensionality from settings
            const uint mNumberOfDimensions;

            //! Quad4, Quad9, Hex8 or Hex27 is supported
            const uint mNumberOfBasisPerElement;

            //! 2, 8 or26
            const uint mNumberOfNeighborsPerElement;

            //! cell of Lagrange elements
            Cell< Element* > mAllElementsOnProc;

            //! cell of coarsest Lagrange elements including aura
            Cell< Element* > mAllCoarsestElementsOnProc;

            //! Cell containing all nodes this proc has
            Cell< Basis* > mAllBasisOnProc;

            //! all elements on proc, including aura, active and refined
            luint mNumberOfAllElementsOnProc = 0;

            //! node counter ( nodes including aura )
            luint mNumberOfAllBasis = 0;

            //! number of active elements used by this proc
            luint mNumberOfElements = 0;

            //! number of active elements used by this proc including Aura
            luint mNumberOfElementsIncludingAura = 0;

            //! refinement pattern this mesh relates to
            const uint mActivationPattern;

            //! index
            moris::uint mMeshIndex = MORIS_UINT_MAX;

            //! maximum level of entities that are active
            uint mMaxLevel = 0;

            // ----------------------------------------------------------------------------

          public:
            // ----------------------------------------------------------------------------
            /**
             * Default Mesh constructor
             *
             * @param[in] aParameters         container of user defined settings
             * @param[in] aBackgroundMesh   pointer to background mesh
             * @param[in] aOrder            polynomial degree of mesh
             */
            Mesh_Base( const Parameters*  aParameters,
                    Background_Mesh_Base* aBackgroundMesh,
                    uint           aOrder,
                    uint           aActivationPattern );

            // ----------------------------------------------------------------------------

            /**
             * Virtual destructor. Does nothing.
             */
            virtual ~Mesh_Base(){};

            // ----------------------------------------------------------------------------

            uint
            get_index() const
            {
                return mMeshIndex;
            }

            // ----------------------------------------------------------------------------

            void
            set_index( uint aIndex )
            {
                mMeshIndex = aIndex;
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to the parameters object
             */
            const Parameters*
            get_parameters() const
            {
                return mParameters;
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to a basis
             *
             * @param[in]    aBasisIndex           number of node in memory
             * @return const Basis*  pointer to node
             */
            Basis*
            get_basis_by_memory_index( luint aBasisIndex )
            {
                return mAllBasisOnProc( aBasisIndex );
            }

            // ----------------------------------------------------------------------------

            /**
             * returns the size of mAllBasisOnProc
             */
            luint
            get_number_of_all_basis_on_proc() const
            {
                return mAllBasisOnProc.size();
            }

            // ----------------------------------------------------------------------------

            /**
             * returns the number of active elements owned by proc
             * ( refers to initialization or last call of update_mesh)
             *
             * @return luint
             */
            luint
            get_number_of_elements() const
            {
                // return mNumberOfElements;
                return mBackgroundMesh->get_number_of_active_elements_on_proc();
            }

            // ----------------------------------------------------------------------------

            /**
             * returns the number of active elements including aura proc
             * ( refers to initialization or last call of update_mesh)
             *
             * @return luint
             */
            luint
            get_number_of_elements_including_aura() const
            {
                // return mNumberOfElements;
                return mBackgroundMesh->get_number_of_active_elements_on_proc_including_aura();
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to an active element
             *
             * @param[in] aElementIndex  number of active element
             *
             * @return const Element*   pointer to node
             */
            Element*
            get_element( const luint aElementIndex )
            {
                return mAllElementsOnProc( mBackgroundMesh->get_element( aElementIndex )
                                                   ->get_memory_index() );
            }
            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to an active element ( const version )
             *
             * @param[in] aElementIndex  number of active element
             *
             * @return const Element*   pointer to node
             */
            const Element*
            get_element( luint aElementIndex ) const
            {
                std::cout << mAllElementsOnProc.size() << ", " << mNumberOfElements << std::endl;
                return mAllElementsOnProc( mBackgroundMesh->get_element( aElementIndex )
                                                   ->get_memory_index() );
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to an active element
             *
             * @param[in] aElementIndex  number of active element
             *
             * @return const Element*   pointer to node
             */
            Element*
            get_element_including_aura( const luint aElementIndex )
            {
                return mAllElementsOnProc( mBackgroundMesh->get_element_including_aura( aElementIndex )
                                                   ->get_memory_index() );
            }
            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to an active element ( const version )
             *
             * @param[in] aElementIndex  number of active element
             *
             * @return const Element*   pointer to node
             */
            const Element*
            get_element_including_aura( luint aElementIndex ) const
            {
                return mAllElementsOnProc( mBackgroundMesh->get_element_including_aura( aElementIndex )
                                                   ->get_memory_index() );
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to an active, refined or deactive element
             *
             * @param[in] aMemoryIndex  number of node in memory
             *
             * @return const Element*   pointer to node
             */
            Element*
            get_element_by_memory_index( luint aMemoryIndex )
            {
                return mAllElementsOnProc( aMemoryIndex );
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to an active, refined or deactive element
             * ( const version )
             *
             * @param[in] aMemoryIndex  number of node in memory
             *
             * @return const Element*   pointer to node
             */
            const Element*
            get_element_by_memory_index( luint aMemoryIndex ) const
            {
                return mAllElementsOnProc( aMemoryIndex );
            }

            // ----------------------------------------------------------------------------

            /**
             * returns the maximum number of elements on this proc, including
             * aura, refined and deactive
             */
            luint
            get_number_of_all_elements_on_proc() const
            {
                return mAllElementsOnProc.size();
            }

            // ----------------------------------------------------------------------------

            /**
             * Creates an Element and links it to corresponding element
             * on the background mesh.
             *
             * @param[in] aElement  pointer to element on background mesh
             *
             * @return Element*  new B-Spline or Lagrange element
             */
            virtual Element* create_element( Background_Element_Base* aElement ) = 0;

            // ----------------------------------------------------------------------------

            /**
             * returns the polynomial degree of the mesh
             */
            uint
            get_order() const
            {
                return mOrder;
            }

            // ----------------------------------------------------------------------------

            /**
             * returns the number of basis per element
             */
            uint
            get_number_of_basis_per_element() const
            {
                return mNumberOfBasisPerElement;
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a pointer to the parent of an element
             */
            Element*
            get_parent_of_element( Element* aElement )
            {
                if ( aElement->get_level() > 0 )
                {
                    return mAllElementsOnProc( aElement->get_background_element()
                                                       ->get_parent()
                                                       ->get_memory_index() );
                }
                else
                {
                    return nullptr;
                }
            }

            // ----------------------------------------------------------------------------

            /**
             * returns a moris::mat of   < number of basis >
             *                         * < number of dimensions >
             *
             * containing the coordinates of the mesh
             */
            void get_basis_coords_of_element(
                    Matrix< DDRMat >& aBasisCoords,
                    luint      aElementIndex );

            // ----------------------------------------------------------------------------

            /**
             * unset the flag of all nodes on proc
             */
            void unflag_all_basis();

            // ----------------------------------------------------------------------------

            /**
             * unset the used flag of all basis on proc
             */
            void unuse_all_basis();

            // ----------------------------------------------------------------------------

            /**
             * set the flag of all nodes on proc
             */
            void flag_all_basis();

            // ----------------------------------------------------------------------------

            /**
             * returns the pattern this mesh refers to
             */
            uint
            get_activation_pattern() const
            {
                return mActivationPattern;
            }

            // ----------------------------------------------------------------------------

            /**
             * tell background to activate pattern of this mesh
             */
            void
            select_activation_pattern()
            {
                mBackgroundMesh->set_activation_pattern( mActivationPattern );
            }

            // ----------------------------------------------------------------------------

            Background_Mesh_Base*
            get_background_mesh()
            {
                return mBackgroundMesh;
            }

            // ----------------------------------------------------------------------------

            uint
            get_max_level() const
            {
                return mMaxLevel;
            }

            // ----------------------------------------------------------------------------

          protected:
            // ----------------------------------------------------------------------------
            /**
             * deletes all nodes and elements
             *
             * @ return void
             */
            void delete_pointers();

            // ----------------------------------------------------------------------------

            /**
             * creates elements (without nodes )
             *
             * @return void
             */
            void create_elements();

            // ----------------------------------------------------------------------------

            /**
             * Determines which nodes are connected to which element.
             * Writes value in mesh data struct
             *
             * @return void
             */
            void determine_elements_connected_to_basis();

            // ----------------------------------------------------------------------------

          public:
            /**
             * Creates a list of basis shared with a neighbor.
             * NOTE: aProcNeighborIndexis the index in which a
             * proc appears in mMyProcNeighbors and not the actual
             * rank of the proc.
             *
             *  @param[in]   aProcNeighborIndex
             *  @param[in]   aUseInverseAura  true:  collect inverse aura instead
             *                                false: collect aura
             *
             *  @param[out]  aNodeList     cell with node pointers on aura
             */
            void collect_basis_from_aura( uint aProcNeighborIndex,
                    uint                       aMode,
                    Cell< Basis* >&                   aNodeList );

          protected:
            // ----------------------------------------------------------------------------

            /**
             * Determines node ownership based on element ownership on
             * background mesh
             *
             * @return void
             */
            void guess_basis_ownership();

            // ----------------------------------------------------------------------------

            /**
             * Make sure that neighbor proc agrees with assumed ownership.
             * Correct otherwise
             */
            void confirm_basis_ownership();

            // ----------------------------------------------------------------------------

            void get_reference_element_of_basis( Basis* aBasis,
                    luint&                              aElementMemoryIndex,
                    uint&                               aElementLocalBasisIndex );
            // ----------------------------------------------------------------------------

            /**
             * this funciton updates the indices for the Elements
             * and is called during the mesh update
             */
            void update_element_indices();

            // ----------------------------------------------------------------------------

          private:
            // ----------------------------------------------------------------------------

            /**
             * collects all coarsest elements on proc including aura
             *
             * @return void
             */
            void collect_coarsest_elements();

            /**
             * returns the number of basis if this element is serendipity
             *
             */
            // uint
            // get_number_of_serendipity_basis();

            // ----------------------------------------------------------------------------

          protected:
            moris::luint get_max_basis_hmr_id();

            //------------------------------------------------------------------------------

            void sanity_check_for_ids_and_ownership();

            //------------------------------------------------------------------------------

            void encode_foreign_basis_path( Cell< Basis* >& aBasis,
                    const moris_id&                         aOwner,
                    Matrix< DDLUMat >&                      aElementAncestors,
                    Matrix< DDUMat >&                       aElementPedigree,
                    Matrix< DDUMat >&                       aElementLocalIndex );
        };

        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_MESH_BASE_HPP_ */
