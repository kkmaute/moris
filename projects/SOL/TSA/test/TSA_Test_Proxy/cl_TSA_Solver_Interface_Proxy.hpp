/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Solver_Interface_Proxy.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_TSA_SOLVER_INTERFACE_PROXY_HPP_
#define SRC_DISTLINALG_CL_TSA_SOLVER_INTERFACE_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
    }
    namespace tsa
    {
        // class Nonlinear_Solver;
        class TSA_Solver_Interface_Proxy : public Solver_Interface
        {
          private:
            moris::uint             mNumMyDofs;                      // local dimension of the problem
            moris::Matrix< DDSMat > mMyGlobalElements;               // local-to-global map
            moris::Matrix< DDSMat > mMyGlobalElementsOverlapping;    // local-to-global map
            moris::uint             mNumElements;                    // number local elements
            moris::Matrix< DDSMat > mEleDofConectivity;              // element - dof conectivities
            moris::Matrix< DDRMat > mElementMatrixValues;            // dense element matrix entries
            moris::Matrix< DDUMat > mMyConstraintDofs;               // constraint dofs
            moris::Matrix< DDRMat > mMyRHSValues;                    // Vector with RHS values

            bool mUseMatrixMarketFiles;    // determines is matrix and RHS comes from MatrixMarket files

            sol::Dist_Vector* mSolutionVector;
            sol::Dist_Vector* mSolutionVectorPrev;
            Matrix< DDRMat >  mMySolVec;
            Matrix< DDRMat >  mMySolVecPrev;

            Vector< enum MSI::Dof_Type > mListOfDofTypes;
            Vector< enum MSI::Dof_Type > mListSecondaryOfDofTypes;

            moris::real      mk = 2;
            Matrix< DDRMat > mT;
            Matrix< DDRMat > mPreviousT;
            moris::real      mDeltaT = 0.0;

            Matrix< DDSMat > mTimeLevelIdsMinus;
            Matrix< DDSMat > mTimeLevelIdsPlus;

          public:
            TSA_Solver_Interface_Proxy();

            // TSA_Solver_Interface_Proxy_( std::shared_ptr< Nonlinear_Solver > aNewtonSolver ){};

            // ----------------------------------------------------------------------------------------------
            ~TSA_Solver_Interface_Proxy() override{};

            // ----------------------------------------------------------------------------------------------

            fem::Element_Type
            get_set_type( uint aMyEquSetInd ) override
            {
                return fem::Element_Type::END_ELEMENT_TYPE;
            }

            void set_solution_vector( sol::Dist_Vector* aSolutionVector ) override;

            void set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector ) override;

            void free_block_memory( const uint aBlockInd ) override{};

            void
            set_time( const Matrix< DDRMat >& aTime ) override
            {
                mT = aTime;
            }

            void compute_IQI() override{};

            void
            set_previous_time( const Matrix< DDRMat >& aTime ) override
            {
                mPreviousT = aTime;
            }
            // ----------------------------------------------------------------------------------------------

            void
            set_requested_dof_types( const Vector< enum MSI::Dof_Type > aListOfDofTypes ) override
            {
                mListOfDofTypes = aListOfDofTypes;
            };

            const Vector< enum MSI::Dof_Type >&
            get_requested_dof_types() override
            {
                return mListOfDofTypes;
            };

            void
            set_secondary_dof_types( const Vector< enum MSI::Dof_Type > aListOfDofTypes ) override
            {
                mListSecondaryOfDofTypes = aListOfDofTypes;
            };

            // local dimension of the problem
            uint
            get_max_num_global_dofs() override
            {
                return 1;
            };

            // ----------------------------------------------------------------------------------------------
            // local dimension of the problem
            uint
            get_num_my_dofs() override
            {
                return mNumMyDofs;
            };

            // ----------------------------------------------------------------------------------------------
            // local-to-global map
            Matrix< DDSMat >
            get_my_local_global_map() override
            {
                mMyGlobalElements.resize( 1, 1 );
                mMyGlobalElements( 0, 0 ) = 0;

                return mMyGlobalElements;
            };

            // ----------------------------------------------------------------------------------------------
            // local-to-global map
            moris::Matrix< DDSMat >
            get_my_local_global_map( const Vector< enum MSI::Dof_Type >& aListOfDofTypes ) override
            {
                mMyGlobalElements.resize( 1, 1 );
                mMyGlobalElements( 0, 0 ) = 0;

                return mMyGlobalElements;
            };

            // ----------------------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_my_local_global_overlapping_map() override
            {
                mMyGlobalElementsOverlapping.resize( 2, 1 );
                mMyGlobalElementsOverlapping( 0, 0 ) = 0;
                mMyGlobalElementsOverlapping( 1, 0 ) = 1;

                return mMyGlobalElementsOverlapping;
            };
            // ----------------------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_my_local_global_overlapping_map( const Vector< enum MSI::Dof_Type >& aListOfDofTypes )
            {
                mMyGlobalElementsOverlapping.resize( 2, 1 );
                mMyGlobalElementsOverlapping( 0, 0 ) = 0;
                mMyGlobalElementsOverlapping( 1, 0 ) = 1;

                return mMyGlobalElementsOverlapping;
            };
            // ----------------------------------------------------------------------------------------------
            // number of elements on processor
            uint
            get_num_my_elements() override
            {
                return mNumElements = 1;
            };

            uint
            get_num_sets() override
            {
                return 1;
            };

            uint
            get_num_equation_objects_on_set( uint aBlockInd ) override
            {
                return mNumElements = 1;
            };

            // ----------------------------------------------------------------------------------------------
            void get_equation_object_operator( const uint& aMyElementInd,
                    Matrix< DDRMat >&                      aElementMatrix ) override;

            void get_equation_object_operator(
                    const uint&       aMyBlockInd,
                    const uint&       aMyElementInd,
                    Matrix< DDRMat >& aElementMatrix ) override;

            // ----------------------------------------------------------------------------------------------

            void
            get_element_topology(
                    const uint&       aMyElementInd,
                    Matrix< DDSMat >& aElementTopology ) override
            {
                aElementTopology.resize( 1, 1 );
                aElementTopology( 0, 0 ) = 0;
            };

            void
            get_element_topology(
                    const uint&       aMyBlockInd,
                    const uint&       aMyElementInd,
                    Matrix< DDSMat >& aElementTopology ) override
            {
                aElementTopology.resize( 1, 1 );
                aElementTopology( 0, 0 ) = 0;
            };

            // ----------------------------------------------------------------------------------------------

            Matrix< DDUMat >
            get_constrained_Ids() override
            {
                return mMyConstraintDofs;
            };

            // ----------------------------------------------------------------------------------------------
            void get_equation_object_rhs(
                    const uint&                 aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override;

            void get_equation_object_rhs(
                    const uint&                 aMyBlockInd,
                    const uint&                 aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override;

            //------------------------------------------------------------------------------

            void get_equation_object_operator_and_rhs(
                    const moris::uint&          aMyElementInd,
                    Matrix< DDRMat >&           aElementMatrix,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override;

            void get_equation_object_operator_and_rhs(
                    const moris::uint&          aMyEquSetInd,
                    const moris::uint&          aMyElementInd,
                    Matrix< DDRMat >&           aElementMatrix,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override;

            // ----------------------------------------------------------------------------------------------

            void
            use_matrix_market_files() override
            {
                mUseMatrixMarketFiles = true;
            };

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
            };

            // ----------------------------------------------------------------------------------------------

            uint
            get_num_rhs() override
            {
                return 1;
            }

            // ----------------------------------------------------------------------------------------------

            uint
            get_num_eigen_vectors() override
            {
                return 0;
            }
        };
    }    // namespace tsa
}    // namespace moris
#endif /* SRC_DISTLINALG_CL_TSA_SOLVER_INTERFACE_PROXY_HPP_ */
