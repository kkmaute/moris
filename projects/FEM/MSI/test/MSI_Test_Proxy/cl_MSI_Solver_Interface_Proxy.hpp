/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Solver_Interface_Proxy.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src
#include "cl_MSI_Solver_Interface.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris::MSI
{
    class MSI_Solver_Interface_Proxy : public MSI_Solver_Interface
    {
      private:
        moris::uint                            mNumMyDofs;                            // local dimension of the problem
        moris::Matrix< DDSMat >                mMyGlobalElements;                     // local-to-global map
        moris::Matrix< DDSMat >                mMyGlobalElementsOverlapping;          // local-to-global map
        moris::Matrix< IdMat >                 mMyGlobalElementsOverlappingOwners;    // local-to-global map
        moris::uint                            mNumElements;                          // number local elements
        Vector< moris::Matrix< DDSMat > > mEleDofConectivity;                    // element - dof connectivities
        moris::Matrix< DDRMat >                mElementMatrixValues;                  // dense element matrix entries
        moris::uint                            mNumDofsPerElement;                    // dofs per element
        moris::Matrix< DDUMat >                mMyConstraintDofs;                     // constraint dofs
        Vector< moris::Matrix< DDRMat > >        mMyRHSValues;                          // Vector with RHS values
        uint                                   mNumRHS = 1;
        Vector< enum MSI::Dof_Type >      mDummyDofs;
        moris::Matrix< DDRMat >                mElementMassMatrixValues;    // RHS Matrix for Eigen Problem

        Matrix< IdMat > mCommTable;

        bool mUseMatrixMarketFiles;    // determines is matrix and RHS comes from MatrixMarket files

        uint mSwitchToEigenProblem = 0;

        sol::Dist_Vector* mEigVector;

      public:
        // ----------------------------------------------------------------------------------------------

        // default constructor
        MSI_Solver_Interface_Proxy();

        // ----------------------------------------------------------------------------------------------

        // default destructor
        ~MSI_Solver_Interface_Proxy() override {}

        // ----------------------------------------------------------------------------------------------

        // local dimension of the problem
        uint
        get_num_my_dofs() override
        {
            return mNumMyDofs;
        }

        // ----------------------------------------------------------------------------------------------

        uint
        get_num_rhs() override
        {
            return mNumRHS;
        }

        // ----------------------------------------------------------------------------------------------

        uint
        get_num_eigen_vectors() override
        {
            return 0;
        }

        // ----------------------------------------------------------------------------------------------

        uint
        get_max_num_global_dofs() override
        {
            moris::uint tNumMyDofs = mNumMyDofs;

            // sum up all distributed dofs
            moris::uint tMaxNumGlobalDofs = sum_all( tNumMyDofs );

            return tMaxNumGlobalDofs;
        }

        // ----------------------------------------------------------------------------------------------

        void free_block_memory( const uint aBlockInd ) override {};

        // ----------------------------------------------------------------------------------------------

        void set_solution_vector( sol::Dist_Vector* aSolutionVector ) override {};

        // ----------------------------------------------------------------------------------------------
        // local-to-global map

        Matrix< DDSMat >
        get_my_local_global_map() override
        {
            return mMyGlobalElements;
        }

        // ----------------------------------------------------------------------------------------------

        moris::Matrix< DDSMat >
        get_my_local_global_overlapping_map() override
        {
            return mMyGlobalElementsOverlapping;
        }

        // ----------------------------------------------------------------------------------------------

        moris::Matrix< DDSMat >
        get_my_local_global_overlapping_map_owners()
        {
            return mMyGlobalElementsOverlappingOwners;
        }

        // ----------------------------------------------------------------------------------------------
        // number of elements on proc

        uint
        get_num_sets() override
        {
            return 1;
        }

        // ----------------------------------------------------------------------------------------------

        // number of elements on proc
        uint
        get_num_my_elements() override
        {
            return mNumElements;
        }

        // ----------------------------------------------------------------------------------------------

        // number of elements on proc
        uint
        get_num_equation_objects_on_set( uint aBlockInd ) override
        {
            return mNumElements;
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_operator(
                const uint&       aMyElementInd,
                Matrix< DDRMat >& aElementMatrix ) override
        {
            aElementMatrix = mElementMatrixValues;
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_operator(
                const uint&       aMyBlockInd,
                const uint&       aMyElementInd,
                Matrix< DDRMat >& aElementMatrix ) override
        {
            if ( mSwitchToEigenProblem < mNumElements * 1 )    // 1 refers to number of blocks
            {
                aElementMatrix = mElementMatrixValues;
                mSwitchToEigenProblem++;
            }

            else
            {
                aElementMatrix = mElementMassMatrixValues;
            }
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_element_topology(
                const uint&       aMyElementInd,
                Matrix< DDSMat >& aElementTopology ) override
        {
            aElementTopology = mEleDofConectivity( aMyElementInd ).get_column( 0 );
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_element_topology(
                const uint&       aMyBlockInd,
                const uint&       aMyElementInd,
                Matrix< DDSMat >& aElementTopology ) override
        {
            aElementTopology = mEleDofConectivity( aMyElementInd ).get_column( 0 );
        }

        // ----------------------------------------------------------------------------------------------

        Matrix< DDUMat >
        get_constrained_Ids() override
        {
            return mMyConstraintDofs;
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_rhs(
                const uint&               aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS ) override
        {
            aElementRHS.resize( mNumRHS );

            for ( uint Ik = 0; Ik < mNumRHS; Ik++ )
            {
                aElementRHS( Ik ) = mMyRHSValues( Ik ).get_column( aMyElementInd );
            }
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_rhs(
                const uint&               aMyBlockInd,
                const uint&               aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS ) override
        {
            aElementRHS.resize( mNumRHS );

            for ( uint Ik = 0; Ik < mNumRHS; Ik++ )
            {
                aElementRHS( Ik ) = mMyRHSValues( Ik ).get_column( aMyElementInd );
            }
        }

        //------------------------------------------------------------------------------

        void
        get_equation_object_operator_and_rhs(
                const moris::uint&        aMyElementInd,
                Matrix< DDRMat >&         aElementMatrix,
                Vector< Matrix< DDRMat > >& aElementRHS ) override
        {
            aElementMatrix = mElementMatrixValues;

            aElementRHS.resize( mNumRHS );

            for ( uint Ik = 0; Ik < mNumRHS; Ik++ )
            {
                aElementRHS( Ik ) = mMyRHSValues( Ik ).get_column( aMyElementInd );
            }
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_operator_and_rhs(
                const moris::uint&        aMyEquSetInd,
                const moris::uint&        aMyElementInd,
                Matrix< DDRMat >&         aElementMatrix,
                Vector< Matrix< DDRMat > >& aElementRHS ) override
        {
            aElementMatrix = mElementMatrixValues;

            aElementRHS.resize( mNumRHS );

            for ( uint Ik = 0; Ik < mNumRHS; Ik++ )
            {
                aElementRHS( Ik ) = mMyRHSValues( Ik ).get_column( aMyElementInd );
            }
        }

        // ----------------------------------------------------------------------------------------------

        void
        use_matrix_market_files() override
        {
            mUseMatrixMarketFiles = true;
        }

        // ----------------------------------------------------------------------------------------------

        const char*
        get_matrix_market_path() override
        {
            if ( mUseMatrixMarketFiles == true )
            {
                const char* tFilePath = "/home/schmidt/codes/MORIS/test/src/distlinalg/";
                return tFilePath;
            }
            else
            {
                return NULL;
            }
        }

        // ----------------------------------------------------------------------------------------------

        virtual const Vector< enum MSI::Dof_Type >&
        get_requested_dof_types() override
        {
            return mDummyDofs;
        }

        // ---------------------------------------------------------------------------------------------

        sol::Dist_Vector*
        get_eigen_solution_vector() override
        {
            return mEigVector;
        }

        // ---------------------------------------------------------------------------------------------

        // return the dummy value of 5
        virtual uint
        estimate_number_of_nonzero_columns() override
        {
            // return 5, it is a mock number
            return 5;
        }
    };
}    // namespace moris::MSI
