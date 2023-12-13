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

#include "cl_HMR_Element.hpp"      //HMR/src
#include "cl_HMR_Mesh_Base.hpp"    //HMR/src

namespace moris::hmr
{
    //------------------------------------------------------------------------------

    /**
     * \brief Base class for Bspline Mesh
     */
    class BSpline_Mesh_Base : public Mesh_Base
    {
      protected:
        //! Cell containing all basis functions this proc knows about
        Cell< Basis_Function* > mAllCoarsestBFsOnProc;

        //! Cell of basis functions that are assigned an HMR index ( moris ID );
        Cell< Basis_Function* > mIndexedBasisFunctions;

        luint                   mNumActiveBFsOnProc    = 0;
        luint                   mNumCandidateBFsOnProc = 0;
        luint                   mNumBFsOnProc          = 0;
        Cell< Basis_Function* > mActiveBFsOnProc;
        Cell< Basis_Function* > mRefinedBFsOnProc;
        Cell< Basis_Function* > mCandidateBFsOnProc;

        Matrix< DDRMat > mChildStencil;

      public:
        /**
         * Default mesh constructor
         *
         * @param aParameters Container of user defined settings
         * @param aBackgroundMesh Pointer to background mesh
         * @param aOrder Polynomial degree of mesh
         * @param aActivationPattern Activation pattern of this mesh
         * @param aNumberOfBasesPerElement Number of bases per element on this mesh
         */
        BSpline_Mesh_Base(
                const Parameters*     aParameters,
                Background_Mesh_Base* aBackgroundMesh,
                uint                  aOrder,
                uint                  aActivationPattern,
                uint                  aNumberOfBasesPerElement );

        // ----------------------------------------------------------------------------

        /**
         * Virtual destructor. Does nothing.
         */
        virtual ~BSpline_Mesh_Base(){};

        /**
         * Gets the polynomial order in a specific direction
         * @note Name hiding from base HMR mesh is intentional, a B-spline mesh cannot operate with a singular order
         *
         * @param aDimensionIndex Dimension index (0, 1, or 2)
         * @return Polynomial order
         */
        virtual uint get_order( uint aDimensionIndex ) = 0;

        /**
         * Gets the minimum polynomial order of this mesh
         *
         * @return Minimum polynomial order
         */
        virtual uint get_min_order() = 0;

        /**
         * Gets the maximum polynomial order of this mesh
         *
         * @return Maximum polynomial order
         */
        virtual uint get_max_order() = 0;

        /**
         * Gets the number of bases in a B-spline element on this mesh
         *
         * @return Number of bases
         */
        virtual uint get_number_of_bases_per_element() = 0;

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
         * returns an active basis by a position in the memory
         */
        Basis_Function*
        get_active_basis( luint aIndex )
        {
            return mActiveBFsOnProc( aIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * returns an active basis by a position in the memory ( const version )
         */
        const Basis_Function*
        get_active_basis( luint aIndex ) const
        {
            return mActiveBFsOnProc( aIndex );
        }

        // ----------------------------------------------------------------------------

        Basis_Function*
        get_basis_by_index( luint aIndex )
        {
            return mIndexedBasisFunctions( aIndex );
        }

        // ----------------------------------------------------------------------------
        uint
        get_number_of_indexed_basis() const
        {
            return mIndexedBasisFunctions.size();
        }

        // ----------------------------------------------------------------------------

        /**
         * returns the number of active basis owned
         * and shared by current proc
         */
        auto
        get_number_of_active_basis_on_proc() const
                -> decltype( mNumActiveBFsOnProc )
        {
            return mNumActiveBFsOnProc;
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
        virtual void flag_refined_basis_of_owned_elements() = 0;

        // ----------------------------------------------------------------------------

        Matrix< DDSMat > get_children_ind_for_basis( const moris::sint aParentBfIndex );

        // ----------------------------------------------------------------------------

        Matrix< DDRMat > get_children_weights_for_parent( const moris::sint aParentBfIndex );

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

        /**
         * Calculates domain-wide unique node ID. Useful for debugging.
         *
         * @param aLevel Level of node
         * @param aIJK Processor-local ijk-positions of node
         * @return Domain wide unique ID
         */
        virtual luint calculate_basis_id(
                uint         aLevel,
                const luint* aIJK ) = 0;

        /**
         * Evaluates the child matrices for an element on this mesh
         *
         * @param aBasisIndices Basis indices
         * @param aChildMatrices Child matrices to fill
         */
        virtual void evaluate_child_matrices(
                const Matrix< DDUMat >&   aBasisIndices,
                Cell< Matrix< DDRMat > >& aChildMatrices ) = 0;

        /**
         * Evaluates the truncation weights on this mesh for 2D/3D based on 1D weights
         *
         * @param aTruncationWeights Truncation weights to size and fill
         */
        virtual void evaluate_truncation_weights( Matrix< DDRMat >& aTruncationWeights ) = 0;

        // ----------------------------------------------------------------------------

        /**
         * makes sure that if a basis is flagged, it is also flagged
         * on any other proc
         */
        void synchronize_flags( const Matrix< IdMat >& aCommTable );

        // ----------------------------------------------------------------------------

        void collect_active_and_refined_elements_from_level(
                uint              aLevel,
                Cell< Element* >& aElements );

        // ----------------------------------------------------------------------------

        /**
         * @brief test returns true if the activation pattern of the basis
         * seems correct. The rule is as follows
         *
         *  A basis function is:
         *
         *      - active, if all connected elements are either active,
         *                refined or padding, and at least one element is active
         *
         *      - refined, if all connected elements are refined or padding
         *
         *      - deactivated, if at least one element is deactivated ( or does not exist)
         *
         *
         *  From these rules, some statements ( see comments in source ) have been developed.
         *  Note that fulfilling the checked statements are NECESSARY, BUT NOT SUFFICIENT conditions.
         *
         *  This test is not meant to be run during runtime.
         */
        bool test_sanity();

        // ----------------------------------------------------------------------------

      private:
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
        virtual void create_basis_on_level_zero() = 0;

        // ----------------------------------------------------------------------------

        /**
         * tells elements on coarsest level which basis functions they have
         *
         * @return void
         */
        virtual void link_basis_functions_to_elements_on_level_zero() = 0;

        // ----------------------------------------------------------------------------

        /**
         * Loops over all elements and stores basis functions in
         * mAllBFsOnProc
         */
        void collect_basis_functions();

        // ----------------------------------------------------------------------------

        /**
         * Provides a Cell of basis functions that live on a specified level
         *
         * @param[ in    ]  aLevel   level to be investigated
         * @param[ inout ]  aBasis   cell containing found basis
         */
        virtual void collect_basis_functions_from_level(
                uint                     aLevel,
                Cell< Basis_Function* >& aBasis ) = 0;

        // ----------------------------------------------------------------------------

        void process_level( uint aLevel );

        // ----------------------------------------------------------------------------

        /**
         * Provides a cell of all basis functions on the current level. Also:
         *      - resets general purpose flags
         *      - determines if bases are used by this proc
         *      - creates basis to element connectivity
         *      - determines basis ownership
         *      - determines basis neighbors for relevant bases
         *
         * @param aElements Elements on this mesh
         * @param aBases Bases to fill for the current level
         */
        virtual void preprocess_basis_on_level(
                Cell< Element* >&        aElements,
                Cell< Basis_Function* >& aBases ) = 0;

        // ----------------------------------------------------------------------------

        /**
         * Identifies bases that are flagged for refinement on the current level
         *
         * @param aBases Bases on the current level
         */
        virtual void determine_basis_function_state( Cell< Basis_Function* >& aBases ) = 0;

        // ----------------------------------------------------------------------------

        /**
         * Links B-Splines to parents. Needed for testing.
         */
        virtual void link_basis_functions_to_parents() = 0;

        // ----------------------------------------------------------------------------

        /*void
        use_only_basis_in_frame(); */

        // ----------------------------------------------------------------------------

        void compute_basis_function_ids();

        // ----------------------------------------------------------------------------
        // void
        // synchronize_active_basis_in_aura();

        // void
        // synchronize_refined_basis_in_aura();

        // ----------------------------------------------------------------------------

        void collect_active_refined_and_candidate_basis_functions();

        // ----------------------------------------------------------------------------

        virtual void delete_unused_basis_functions(
                uint                              aLevel,
                Cell< Background_Element_Base* >& aBackgroundElements,
                Cell< Basis_Function* >&          aBasis ) = 0;

        // ----------------------------------------------------------------------------

    };    // class hmr::BSpline_Mesh_Base

    //------------------------------------------------------------------------------

}    // namespace moris::hmr

#endif /* SRC_HMR_CL_HMR_BSPLINE_MESH_BASE_HPP_ */
