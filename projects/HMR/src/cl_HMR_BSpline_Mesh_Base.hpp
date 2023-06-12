/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Mesh_Base.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_MESH_BASE_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_MESH_BASE_HPP_

#include "cl_HMR_BSpline.hpp"      //HMR/src
#include "cl_HMR_Element.hpp"      //HMR/src
#include "cl_HMR_Mesh_Base.hpp"    //HMR/src

namespace moris
{
    namespace hmr
    {
        //------------------------------------------------------------------------------

        /**
         * \brief Base class for Bspline Mesh
         */
        class BSpline_Mesh_Base : public Mesh_Base
        {

            // ----------------------------------------------------------------------------

          public:
            // ----------------------------------------------------------------------------

            bool test_sanity();

            // ----------------------------------------------------------------------------

          protected:
            // ----------------------------------------------------------------------------
            //! number of children per basis
            const uint mNumberOfChildrenPerBasis;

            //! max number of elements connected to basis
            const uint mNumberOfElementsPerBasis;

            //! Cell containing all basis this proc knows about
            Cell< Basis* > mAllCoarsestBasisOnProc;

            //! Cell of basis that are assigned an HMR index ( moris ID );
            Cell< Basis* > mIndexedBasis;

            //! number of basis used by this proc
            luint mNumberOfBasis = 0;

            //! number of elements used by this proc
            // luint mNumberOfElements = 0;

            //! number of all basis (including unused on padding)
            luint mNumberOfAllBasis = 0;

            //! number of basis on coarsest level
            luint mNumberOfCoarsestBasisOnProc[ 3 ] = { 0, 0, 0 };

            //! Lookup table containing offset for node IDs
            luint mBasisLevelOffset[ gMaxNumberOfLevels ];

            //! a container that remembers the number of basis per level
            luint mNumberOfBasisPerLevel[ gMaxNumberOfLevels ] = { 0 };

            //! counts the number of basis used and owned
            // luint mNumberOfOwnedBasis = 0;

            luint          mNumberOfActiveBasisOnProc  = 0;
            luint          mNumberOfRefinedBasisOnProc = 0;
            Cell< Basis* > mActiveBasisOnProc;
            Cell< Basis* > mRefinedBasisOnProc;

            Matrix< DDRMat > mChildStencil;

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
            BSpline_Mesh_Base( const Parameters* aParameters,
                    Background_Mesh_Base*        aBackgroundMesh,
                    uint                  aOrder,
                    uint                  aActivationPattern );

            // ----------------------------------------------------------------------------

            /**
             * Virtual destructor. Does nothing.
             */
            virtual ~BSpline_Mesh_Base(){};

            // ----------------------------------------------------------------------------

            /**
             * This function is called by the constructor, but can also be called
             * after the B-Spline mesh is generated, and the background mesh is
             * refined.
             *
             * @return void
             *
             */
            void update_mesh();
            // ----------------------------------------------------------------------------

            /**
             * Saves the basis to a VTK file. Useful for debugging.
             *
             * @param[in] string aFilePath file where mesh is to be stored
             *
             * @return void
             */
            void save_to_vtk( const std::string& aFilePath );

            // ----------------------------------------------------------------------------

            /**
             * returns how many children a basis has
             */
            uint
            get_number_of_children_per_basis() const
            {
                return mNumberOfChildrenPerBasis;
            }

            // ----------------------------------------------------------------------------

            /**
             * returns an active basis by a position in the memory
             */
            Basis*
            get_active_basis( luint aIndex )
            {
                return mActiveBasisOnProc( aIndex );
            }

            // ----------------------------------------------------------------------------

            /**
             * returns an active basis by a position in the memory ( const version )
             */
            const Basis*
            get_active_basis( luint aIndex ) const
            {
                return mActiveBasisOnProc( aIndex );
            }

            // ----------------------------------------------------------------------------

            Basis*
            get_basis_by_index( luint aIndex )
            {
                return mIndexedBasis( aIndex );
            }

            // ----------------------------------------------------------------------------
            uint
            get_number_of_indexed_basis() const
            {
                return mIndexedBasis.size();
            }

            // ----------------------------------------------------------------------------

            /**
             * returns the number of active basis owned
             * and shared by current proc
             */
            auto
            get_number_of_active_basis_on_proc() const
                    -> decltype( mNumberOfActiveBasisOnProc )
            {
                return mNumberOfActiveBasisOnProc;
            }

            // ----------------------------------------------------------------------------

            /*
             * A function that tests if each basis is uniquely generated.
             * Returns false otherwise.
             *
             * @return bool
             */
            bool test_for_double_basis();

            // ----------------------------------------------------------------------------

            /**
             * recalculates the domain indices based on flagged basis
             */
            void calculate_basis_indices( const Matrix< IdMat >& aCommTable );

            // ----------------------------------------------------------------------------

            /**
             * special function for multigrid
             */
            void flag_refined_basis_of_owned_elements();

            // ----------------------------------------------------------------------------

            Matrix< DDSMat > get_children_ind_for_basis( const moris::sint aParentBasind );

            // ----------------------------------------------------------------------------

            Matrix< DDRMat > get_children_weights_for_parent( const moris::sint aParentBasind );

            // ----------------------------------------------------------------------------

            uint get_number_of_basis_connected_to_basis( const moris_index aIndex );

            // ----------------------------------------------------------------------------

            /**
             * calculates XZY coordinates for each basis
             *
             * @return void
             */
            virtual void calculate_basis_coordinates() = 0;

            // ----------------------------------------------------------------------------

          protected:
            // ----------------------------------------------------------------------------

            /**
             * creates a basis depending on polynomial order and dimension
             *
             * @param[in]   aIJK        ijk position of node
             * @param[in]   aLevel      level on which basis exists
             * @param[in]   aOwner      owner of basis
             */
            virtual Basis* create_basis( const luint* aIJK,
                    uint                       aLevel,
                    uint                       aOwner ) = 0;

            // ----------------------------------------------------------------------------
            /**
             * Returns the pointer to a basis on the coarsest level. 2D case.
             */
            Basis* get_coarsest_basis_by_ij( luint aI,
                    luint                          aJ );

            // ----------------------------------------------------------------------------

            /**
             * Returns the pointer to a basis on the coarsest level. 3D case.
             */
            Basis* get_coarsest_basis_by_ijk( luint aI,
                    luint                           aJ,
                    luint                           aK );

            // ----------------------------------------------------------------------------

            /**
             * calculates domain wide unique node ID (1D case)
             * Useful for debugging.
             *
             * @param[in]  aLevel    level of node
             * @param[in]  aI        proc local i-position of node
             *
             * @return uint          domain wide unique ID
             */
            virtual luint calculate_basis_id(
                    uint  aLevel,
                    luint aI ) = 0;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * calculates domain wide unique node ID (2D case)
             * Useful for debugging.
             *
             * @param[in]  aLevel    level of node
             * @param[in]  aI        proc local i-position of node
             * @param[in]  aJ        proc local j-position of node
             * @return uint          domain wide unique ID
             */
            virtual luint calculate_basis_id(
                    uint  aLevel,
                    luint aI,
                    luint aJ ) = 0;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * calculates domain wide unique node ID (3D case)
             * Useful for debugging.
             *
             * @param[in]  aLevel    level of node
             * @param[in]  aI        proc local i-position of node
             * @param[in]  aJ        proc local j-position of node
             * @param[in]  aK        proc local k-position of node
             * @return uint          domain wide unique ID
             */
            virtual luint calculate_basis_id(
                    uint  aLevel,
                    luint aI,
                    luint aJ,
                    luint aK ) = 0;

            // ----------------------------------------------------------------------------

            /**
             * makes sure that if a basis is flagged, it is also flagged
             * on any other proc
             */
            void synchronize_flags( const Matrix< IdMat >& aCommTable );

            // ----------------------------------------------------------------------------

            void collect_active_and_refined_elements_from_level( uint aLevel,
                    Cell< Element* >&                                        aElements );

            // ----------------------------------------------------------------------------

          private:
            // ----------------------------------------------------------------------------

            /**
             * Creates the lookup table needed for basis IDs
             * This table contains the maximal possible offset. E.g. a uniform refined mesh.
             * Function seems to be equal to calculate_lookup_tables()
             */
            void calculate_basis_level_offset();

            // ----------------------------------------------------------------------------
            /**
             * creates B-Splines for this mesh
             *
             * @return void
             */
            void create_basis();

            // ----------------------------------------------------------------------------

            /**
             * creates basis on coarsest level
             *
             * @return void
             */
            void create_basis_on_level_zero();

            // ----------------------------------------------------------------------------

            /**
             * tells elements on coarsest level which basis they have
             *
             * @return void
             */
            void link_basis_to_elements_on_level_zero();

            // ----------------------------------------------------------------------------

            /**
             * Loops over all elements and stores basis in
             * mAllBasisOnProc
             */
            void collect_basis();

            // ----------------------------------------------------------------------------

            /**
             * Provides a Cell of basis that live on a specified level
             *
             * @param[ in    ]  aLevel   level to be investigated
             * @param[ inout ]  aBasis   cell containing found basis
             */
            void collect_basis_from_level( uint aLevel,
                    Cell< Basis* >&                    aBasis );

            // ----------------------------------------------------------------------------

            void process_level( uint aLevel );

            // ----------------------------------------------------------------------------

            /**
             * Provides a cell of all basis on current level.
             * Also:
             *     - resets general purpose flags
             *     - determine if basis is used by this proc
             *     - creates basis to element connectivity
             *     - determines basis ownership
             *     - determines basis neighbors for relevant basis
             */
            void preprocess_basis_from_level( Cell< Element* >& aBackgroundElements,
                    Cell< Basis* >&                             aBasis );

            // ----------------------------------------------------------------------------

            /**
             * identifies basis that are flagged for refinement
             */
            void determine_basis_state( Cell< Basis* >& aBasis );

            // ----------------------------------------------------------------------------

            /**
             * Links B-Splines to parents. Needed for testing.
             */
            void link_basis_to_parents();

            // ----------------------------------------------------------------------------

            /*void
            use_only_basis_in_frame(); */

            // ----------------------------------------------------------------------------

            void calculate_basis_ids();

            // ----------------------------------------------------------------------------
            // void
            // synchronize_active_basis_in_aura();

            // void
            // synchronize_refined_basis_in_aura();

            // ----------------------------------------------------------------------------

            void collect_active_and_refined_basis();

            // ----------------------------------------------------------------------------

            void delete_unused_basis( uint     aLevel,
                    Cell< Background_Element_Base* >& aBackgroundElements,
                    Cell< Basis* >&                   aBasis );

            // ----------------------------------------------------------------------------

            void calculate_child_stencil();
        };
        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_MESH_BASE_HPP_ */
