/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Database.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_DATABASE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_DATABASE_HPP_

#include <memory>    // <-- database is always a shared pointer, so we need std::memory
#include <string>

#include "cl_HMR_Factory.hpp"          //HMR/src
#include "cl_HMR_Lagrange_Mesh.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"       //HMR/src
#include "cl_HMR_Side_Set.hpp"         //HMR/src
#include "cl_HMR_T_Matrix.hpp"         //HMR/src
#include "cl_Vector.hpp"               //CNT/src
#include "cl_Map.hpp"

#include "cl_MTK_Side_Sets_Info.hpp"

namespace moris::hmr
{
    // -----------------------------------------------------------------------------

    class Field;

    // -----------------------------------------------------------------------------
    class Database : public std::enable_shared_from_this< Database >
    {

      private:
        //! object containing user settings
        Parameters* mParameters = nullptr;

        //! pointer to background mesh
        Background_Mesh_Base* mBackgroundMesh = nullptr;

        //! cell of pointers to B-Spline meshes
        Vector< BSpline_Mesh_Base* > mBSplineMeshes;

        //! cell of pointers to Lagrange meshes
        Vector< Lagrange_Mesh_Base* > mLagrangeMeshes;

        //! cell of pointers to Lagrange meshes.
        // These Lagrange meshes are created on the flight and not in the input file
        Vector< Vector< Lagrange_Mesh_Base* > > mAdditionalLagrangeMeshes;

        //! communication table for this mesh. Created during finalize.
        Matrix< IdMat > mCommunicationTable;

        //! flag telling if parameter pointer is supposed to be deleted on destruction
        bool mDeleteParametersOnDestruction = false;

        //! Side sets for input pattern
        // Vector< Matrix< IdMat > >   mInputSideSets;

        //! Side sets for output pattern
        Vector< Side_Set > mOutputSideSets;

        map< std::string, moris_index > mOutputSideSetMap;

        bool mHaveRefinedAtLeastOneElement = false;

        //! flag telling if T-Matrices for input mesh have been calculated
        bool mHaveInputTMatrix = false;

        // ! checks if finalize was called. Finalize should only be called once at the end of HMR
        bool mFinalizedCalled = false;

      public:
        // -----------------------------------------------------------------------------

        /**
         * Database constructor
         */
        Database( Parameters* aParameters );

        // -----------------------------------------------------------------------------

        /**
         * alternative constructor which loads a mesh from a h5 file
         */
        Database( const std::string& aPath );

        // -----------------------------------------------------------------------------

        /**
         * alternative constructor which loads two patterns
         */
        Database(
                const std::string& aInputPath,
                const std::string& aOutputPath );

        // -----------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Database();

        // -----------------------------------------------------------------------------

        void load_pattern_from_hdf5_file(
                const std::string& aPath,
                const bool         aMode );

        // -----------------------------------------------------------------------------

        void load_refinement_pattern(
                Matrix< DDLUMat >&                aElementCounterPerLevelAndPattern,
                Vector< Matrix< DDLUMat > >& aElementPerPattern,
                Matrix< DDUMat >&                 aPatternListUniqueMat );

        // -----------------------------------------------------------------------------

        /**
         * sets the flag that the parameter object must be deleted
         * by the destructor
         */
        void set_parameter_owning_flag();

        // -----------------------------------------------------------------------------
        /**
         * sets the flag that the parameter object must not be deleted
         * by the destructor
         */
        void unset_parameter_owning_flag();

        // -----------------------------------------------------------------------------

        /**
         * creates a union of two patterns
         */
        void unite_patterns(
                uint aSourceA,
                uint aSourceB,
                uint aTarget );

        // -----------------------------------------------------------------------------

        /**
         * copies a source pattern to a target pattern
         */
        void copy_pattern(
                uint aSource,
                uint aTarget );
        // -----------------------------------------------------------------------------

        /**
         * runs the refinement scheme
         *
         * returns true if at least one element has been refined
         */
        void perform_refinement(
                const uint aActivePattern,
                const bool aResetPattern = true );

        // -----------------------------------------------------------------------------

        /**
         * aTarget must be a refined variant of aSource
         */
        void interpolate_field(
                uint                           aSourcePattern,
                const std::shared_ptr< Field > aSource,
                uint                           aTargetPattern,
                std::shared_ptr< Field >       aTarget );

        // -----------------------------------------------------------------------------

        void change_field_order(
                std::shared_ptr< Field > aSource,
                std::shared_ptr< Field > aTarget );

        // -----------------------------------------------------------------------------

        /**
         * returns the pointer to a Lagrange mesh, needed by interface
         * constructor
         */
        Lagrange_Mesh_Base*
        get_lagrange_mesh_by_index( uint aIndex )
        {
            return mLagrangeMeshes( aIndex );
        }

        // -----------------------------------------------------------------------------
        Lagrange_Mesh_Base*
        get_additional_lagrange_mesh_by_index(
                uint aIndex,
                uint aPattern )
        {
            return mAdditionalLagrangeMeshes( aPattern )( aIndex );
        }

        // -----------------------------------------------------------------------------
        void
        add_lagrange_mesh(
                Lagrange_Mesh_Base* aLagrangeMesh,
                uint                aLagrangePattern )
        {
            mAdditionalLagrangeMeshes( aLagrangePattern ).push_back( aLagrangeMesh );
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the pointer to a Bspline mesh, needed by interface
         * constructor
         */
        BSpline_Mesh_Base*
        get_bspline_mesh_by_index( uint aIndex )
        {
            return mBSplineMeshes( aIndex );
        }

        // -----------------------------------------------------------------------------

        Background_Mesh_Base* get_background_mesh();

        // -----------------------------------------------------------------------------

        Matrix< DDUMat > create_output_pattern_list();

        // -----------------------------------------------------------------------------

        /**
         * returns the number of ( active ) elements on this proc
         */
        auto
        get_number_of_elements_on_proc()
                -> decltype( mBackgroundMesh->get_number_of_active_elements_on_proc() )
        {
            return mBackgroundMesh->get_number_of_active_elements_on_proc();
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the number of ( active ) elements on this proc
         */
        auto
        get_number_of_padding_elements_on_proc()
                -> decltype( mBackgroundMesh->get_number_of_padding_elements_on_proc() )
        {
            return mBackgroundMesh->get_number_of_padding_elements_on_proc();
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the number of dimensions in space
         */
        auto
        get_number_of_dimensions() const
                -> decltype( mParameters->get_number_of_dimensions() )
        {
            return mParameters->get_number_of_dimensions();
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the number of Lagrange meshes
         */
        uint
        get_number_of_lagrange_meshes() const
        {
            return mLagrangeMeshes.size();
        }

        // -----------------------------------------------------------------------------

        uint
        get_number_of_additional_lagrange_meshes( uint aLagrangePattern ) const
        {
            return mAdditionalLagrangeMeshes( aLagrangePattern ).size();
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the number of Bspline meshes
         */
        uint
        get_number_of_bspline_meshes() const
        {
            return mBSplineMeshes.size();
        }

        // -----------------------------------------------------------------------------

        /**
         * set active pattern of background mesh
         */
        void
        set_activation_pattern( uint aPattern )
        {
            mBackgroundMesh->set_activation_pattern( aPattern );
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the active pattern
         */
        auto
        get_activation_pattern() const
                -> decltype( mBackgroundMesh->get_activation_pattern() )
        {
            return mBackgroundMesh->get_activation_pattern();
        }

        // -----------------------------------------------------------------------------

        /**
         * function needed for tests etc
         */
        void flag_element( luint aIndex );

        // -----------------------------------------------------------------------------

        void flag_parent( luint aIndex );

        // -----------------------------------------------------------------------------

        void create_extra_refinement_buffer_for_level( const uint aLevel );

        // -----------------------------------------------------------------------------

        /**
         * returns the communication table that is needed by FEM
         */
        const Matrix< IdMat >&
        get_communication_table() const
        {
            return mCommunicationTable;
        }

        // -----------------------------------------------------------------------------

        /**
         * returns the proc neighbors for this proc
         */
        Matrix< IdMat >
        get_proc_neighbors() const
        {
            Matrix< IdMat > tProcNeighbors = mBackgroundMesh->get_proc_neighbors();

            uint tCounter = 0;
            for ( uint Ik = 0; Ik < tProcNeighbors.numel(); Ik++ )
            {
                if ( tProcNeighbors( Ik ) != gNoProcNeighbor )
                {
                    tCounter++;
                }
            }

            Matrix< IdMat > tProcNeighborsActive( tCounter, 1, gNoProcNeighbor );

            tCounter = 0;
            for ( uint Ik = 0; Ik < tProcNeighbors.numel(); Ik++ )
            {
                if ( tProcNeighbors( Ik ) != gNoProcNeighbor )
                {
                    tProcNeighborsActive( tCounter++ ) = tProcNeighbors( Ik );
                }
            }

            return tProcNeighborsActive;
        }

        // -----------------------------------------------------------------------------

        /**
         * return pointer to parameter object ( const version )
         */
        Parameters*
        get_parameters()
        {
            return mParameters;
        }

        // -----------------------------------------------------------------------------

        /**
         * return pointer to parameter object ( const version )
         */
        const Parameters*
        get_parameters() const
        {
            return mParameters;
        }

        // -----------------------------------------------------------------------------

        /**
         * populates the member variables of the relevant nodes
         * with their T-Matrices
         */
        void finalize();

        // -----------------------------------------------------------------------------

        /**
         * resets Lagrange and B-Spline meshes and the refinement pattern of the background mesh
         * Does not delete Background mesh elements. Resets teh finalize flag.
         */
        void reset_refined_meshes();

        // -----------------------------------------------------------------------------

        /**
         * needed for exodus output of cubic meshes, called by finalize
         */
        //            void add_extra_refinement_step_for_exodus();

        // -----------------------------------------------------------------------------

        /**
         *  this function updates the meshes after an refinement step
         */
        void update_bspline_meshes();

        void update_bspline_meshes( uint aPattern );

        void update_lagrange_meshes();

        void update_lagrange_meshes( uint aPattern );

        // -----------------------------------------------------------------------------

        /**
         *  test that all relevant entitiy IDs are set
         */
        void check_entity_ids();

        // -----------------------------------------------------------------------------

        // tells if at least one element has been refined in this database
        bool
        have_refined_at_least_one_element() const
        {
            return mHaveRefinedAtLeastOneElement;
        }

        // -----------------------------------------------------------------------------

        /**
         * returns a sideset based on its label
         */
        const Side_Set&
        get_output_side_set( const std::string& aLabel ) const
        {
            return mOutputSideSets( mOutputSideSetMap.find( aLabel ) );
        }

        // -----------------------------------------------------------------------------
        /**
         * returns list of all side sets
         */
        const Vector< Side_Set >&
        get_side_sets() const
        {
            return mOutputSideSets;
        }

        // -----------------------------------------------------------------------------

        void calculate_t_matrices_for_input();

        // -----------------------------------------------------------------------------

        /**
         * creates a union mesh of the input and the output patterns
         */
        //            void create_union_pattern()
        //            {
        //                this->unite_patterns( mParameters->get_lagrange_input_pattern(),
        //                                      mParameters->get_lagrange_output_pattern(),
        //                                      mParameters->get_union_pattern() );
        //            }

        // -----------------------------------------------------------------------------

        void
        create_union_pattern(
                const uint aSourceA,
                const uint aSourceB,
                const uint aTarget )
        {
            this->unite_patterns( aSourceA,
                    aSourceB,
                    aTarget );
        }

        // -----------------------------------------------------------------------------

        bool
        is_finalized()
        {
            return mFinalizedCalled;
        }

        // -----------------------------------------------------------------------------

      private:
        // -----------------------------------------------------------------------------

        /**
         * this function initializes the Lagrange and B-Spline Meshes
         * is complete
         *
         * @return void
         */
        void create_meshes();

        // -----------------------------------------------------------------------------

        bool is_lagrange_input_mesh( const uint aMeshIndex );

        bool is_bspline_input_mesh( const uint aMeshIndex );

        // -----------------------------------------------------------------------------

        /**
         * this function deletes the Lagrange and B-Spline meshes
         * the function is called before create_meshes
         */
        void delete_meshes();

        // -----------------------------------------------------------------------------

        void delete_additional_meshes( uint aPattern );

        // -----------------------------------------------------------------------------

        /**
         * creates the communication table and writes it into
         * mCommunicationTable. Must be called after mesh has been finalized.
         */
        void create_communication_table();

        // -----------------------------------------------------------------------------

        /**
         * creates the sidesets
         */
        void create_side_sets();

        // -----------------------------------------------------------------------------

        void create_working_pattern_for_bspline_refinement();

        // -----------------------------------------------------------------------------

        /**
         * free side set memory space
         */
        void delete_side_sets();

        // -----------------------------------------------------------------------------
    };
}    // namespace moris::hmr

#endif /* PROJECTS_HMR_SRC_CL_HMR_DATABASE_HPP_ */
