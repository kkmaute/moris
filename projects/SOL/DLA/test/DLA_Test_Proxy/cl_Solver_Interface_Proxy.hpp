/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Solver_Interface_Proxy.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_SOLVER_INPUT_TEST_HPP_
#define SRC_DISTLINALG_CL_SOLVER_INPUT_TEST_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    class Solver_Interface_Proxy : public Solver_Interface
    {
      private:
        moris::uint                       mNumMyDofs;                      // local dimension of the problem
        moris::Matrix< DDSMat >           mMyGlobalElements;               // local-to-global map
        moris::Matrix< DDSMat >           mMyGlobalElementsOverlapping;    // local-to-global map
        moris::uint                       mNumElements;                    // number local elements
        moris::Matrix< DDSMat >           mEleDofConectivity;              // element - dof connectivities
        moris::Matrix< DDRMat >           mElementMatrixValues;            // dense element matrix entries
        moris::uint                       mNumDofsPerElement;              // dofs per element
        moris::Matrix< DDUMat >           mMyConstraintDofs;               // constraint dofs
        Vector< moris::Matrix< DDRMat > >   mMyRHSValues;                    // Vector with RHS values
        uint                              mNumRHS = 1;
        Vector< enum MSI::Dof_Type > mDummyDofs;
        moris::Matrix< DDRMat >           mElementMassMatrixValues;    // RHS Matrix for Eigen Problem

        bool mUseMatrixMarketFiles;    // determines is matrix and RHS comes from MatrixMarket files

        uint mSwitchToEigenProblem = 0;

        sol::Dist_Vector* mEigVector;

      public:
        // ----------------------------------------------------------------------------------------------

        Solver_Interface_Proxy();

        // ----------------------------------------------------------------------------------------------

        Solver_Interface_Proxy( uint aNumRHS );

        // ----------------------------------------------------------------------------------------------

        Solver_Interface_Proxy( std::string aProblem );

        // ----------------------------------------------------------------------------------------------

        ~Solver_Interface_Proxy(){};

        // ----------------------------------------------------------------------------------------------
        // local dimension of the problem

        uint
        get_num_my_dofs()
        {
            return mNumMyDofs;
        }

        // ----------------------------------------------------------------------------------------------

        uint
        get_num_rhs()
        {
            return mNumRHS;
        }

        // ----------------------------------------------------------------------------------------------

        uint
        get_num_eigen_vectors()
        {
            return 0;
        }

        // ----------------------------------------------------------------------------------------------

        uint
        get_max_num_global_dofs()
        {
            moris::uint tNumMyDofs = mNumMyDofs;

            // sum up all distributed dofs
            moris::uint tMaxNumGlobalDofs = sum_all( tNumMyDofs );

            return tMaxNumGlobalDofs;
        };

        // ----------------------------------------------------------------------------------------------

        void free_block_memory( const uint aBlockInd ){};

        // ----------------------------------------------------------------------------------------------

        void set_solution_vector( sol::Dist_Vector* aSolutionVector ){};

        // ----------------------------------------------------------------------------------------------
        // local-to-global map

        Matrix< DDSMat >
        get_my_local_global_map()
        {
            return mMyGlobalElements;
        }

        // ----------------------------------------------------------------------------------------------

        moris::Matrix< DDSMat >
        get_my_local_global_overlapping_map()
        {
            return mMyGlobalElementsOverlapping;
        };

        // ----------------------------------------------------------------------------------------------
        // number of elements on proc

        uint
        get_num_sets()
        {
            return 1;
        }

        // ----------------------------------------------------------------------------------------------

        // number of elements on proc
        uint
        get_num_my_elements()
        {
            return mNumElements;
        }

        // ----------------------------------------------------------------------------------------------

        // number of elements on proc
        uint
        get_num_equation_objects_on_set( uint aBlockInd )
        {
            return mNumElements;
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_operator(
                const uint&       aMyElementInd,
                Matrix< DDRMat >& aElementMatrix )
        {
            aElementMatrix = mElementMatrixValues;
        };

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_operator(
                const uint&       aMyBlockInd,
                const uint&       aMyElementInd,
                Matrix< DDRMat >& aElementMatrix )
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
                Matrix< DDSMat >& aElementTopology )
        {
            aElementTopology = mEleDofConectivity.get_column( aMyElementInd );
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_element_topology(
                const uint&       aMyBlockInd,
                const uint&       aMyElementInd,
                Matrix< DDSMat >& aElementTopology )
        {
            std::cout << aMyElementInd << std::endl;
            aElementTopology = mEleDofConectivity.get_column( aMyElementInd );
        };

        // ----------------------------------------------------------------------------------------------

        Matrix< DDUMat >
        get_constrained_Ids()
        {
            return mMyConstraintDofs;
        }

        // ----------------------------------------------------------------------------------------------

        void
        get_equation_object_rhs(
                const uint&               aMyElementInd,
                Vector< Matrix< DDRMat > >& aElementRHS )
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
                Vector< Matrix< DDRMat > >& aElementRHS )
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
                Vector< Matrix< DDRMat > >& aElementRHS )
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
                Vector< Matrix< DDRMat > >& aElementRHS )
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
        use_matrix_market_files()
        {
            mUseMatrixMarketFiles = true;
        }

        // ----------------------------------------------------------------------------------------------

        const char*
        get_matrix_market_path()
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
        get_requested_dof_types()
        {
            return mDummyDofs;
        };

        // ---------------------------------------------------------------------------------------------

        sol::Dist_Vector*
        get_eigen_solution_vector()
        {
            return mEigVector;
        }
    };
}    // namespace moris
#endif /* SRC_DISTLINALG_CL_SOLVER_INPUT_TEST_HPP_ */
