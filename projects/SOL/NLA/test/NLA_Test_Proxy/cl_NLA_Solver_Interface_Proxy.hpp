/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Interface_Proxy.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_NLA_SOLVER_INPUT_TEST_HPP_
#define SRC_DISTLINALG_CL_NLA_SOLVER_INPUT_TEST_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Communication_Tools.hpp"    // COM/src

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
    }
    namespace NLA
    {
        class Nonlinear_Algorithm;
        class NLA_Solver_Interface_Proxy : public Solver_Interface
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

            sol::Dist_Vector* mSolutionVector = nullptr;
            Matrix< DDRMat >  mMySolVec;

            Matrix< DDRMat > ( *mFunctionRes )( const moris::sint aNX, const moris::sint aNY, const moris::real aLambda, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd );
            Matrix< DDRMat > ( *mFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd );
            Vector< Matrix< DDRMat > > ( *mFunctionObj )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd );
            Matrix< DDSMat > ( *mFunctionTopology )( const moris::sint aNX, const moris::sint aNY, const moris::uint aEquationObjectInd );

            Vector< Matrix< DDRMat >> mIQIVal;
            

            moris::sint mNX;
            moris::sint mNY;
            Matrix < DDRMat > mMyValues;

            moris::uint mMyElementInd = 0;

            Matrix< DDRMat > mTime = { { 1.0 }, { 1.0 } };

            Matrix< DDSMat > mTimeLevelIdsMinus;
            Matrix< DDSMat > mTimeLevelIdsPlus;

            sol::Dist_Vector* mSolutionVectorPrev = nullptr;

            Vector< enum MSI::Dof_Type > mListOfDofTypes;

          public:
            // ----------------------------------------------------------------------------

            NLA_Solver_Interface_Proxy();

            // ----------------------------------------------------------------------------

            NLA_Solver_Interface_Proxy(
                    const moris::uint aNumMyDofs,
                    const moris::uint aNumElements,
                    const moris::sint aNX,
                    const moris::sint aNY,
                    Matrix< DDRMat > ( *aFunctionRes )( const moris::sint aNX, const moris::sint aNY, const moris::real aLambda, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd ),
                    Matrix< DDRMat > ( *aFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd ),
                    Matrix< DDSMat > ( *aFunctionTopo )( const moris::sint aNX, const moris::sint aNY, const moris::uint aEquationObjectInd ) );

            // ----------------------------------------------------------------------------

            NLA_Solver_Interface_Proxy(
                    const moris::uint aNumMyDofs,
                    const moris::uint aNumElements,
                    const moris::sint aNX,
                    const moris::sint aNY,
                    const moris::uint aEquationObjectInd,
                    Matrix< DDRMat > ( *aFunctionRes )( const moris::sint aNX, const moris::sint aNY, const moris::real aLambda, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd ),
                    Matrix< DDRMat > ( *aFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd ),
                    Vector< Matrix< DDRMat > > ( *aFunctionObj )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat >& tMyValues, const moris::uint aEquationObjectInd ),
                    Matrix< DDSMat > ( *aFunctionTopo )( const moris::sint aNX, const moris::sint aNY, const moris::uint aEquationObjectInd )
            );


            // ----------------------------------------------------------------------------

            NLA_Solver_Interface_Proxy( const std::shared_ptr< Nonlinear_Algorithm >& aNewtonSolver ) {};

            // ----------------------------------------------------------------------------------------------

            ~NLA_Solver_Interface_Proxy() override {};

            // ----------------------------------------------------------------------------------------------

            void set_time_value(
                    const moris::real& aLambda,
                    moris::uint        aPos = 1 ) override;
            // ----------------------------------------------------------------------------

            void set_time( const Matrix< DDRMat >& aTime ) override;

            // ----------------------------------------------------------------------------

            void set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector ) override;

            // ----------------------------------------------------------------------------------------------

            void set_solution_vector( sol::Dist_Vector* aSolutionVector ) override;

            // ----------------------------------------------------------------------------

            void free_block_memory( const uint aBlockInd ) override {};

            // ----------------------------------------------------------------------------------------------

            void
            set_requested_dof_types( const Vector< enum MSI::Dof_Type >& aListOfDofTypes ) override
            {
                mListOfDofTypes = aListOfDofTypes;
            };

            // ----------------------------------------------------------------------------

            const Vector< enum MSI::Dof_Type >&
            get_requested_dof_types() override
            {
                return mListOfDofTypes;
            };

            // ----------------------------------------------------------------------------

            void set_secondary_dof_types( const Vector< enum MSI::Dof_Type >& aListOfDofTypes ) override {};

            // ----------------------------------------------------------------------------------------------
            // local dimension of the problem

            uint
            get_num_my_dofs() override
            {
                return mNumMyDofs;
            };

            // ----------------------------------------------------------------------------

            void compute_IQI() override;

            // ----------------------------------------------------------------------------
            
            void set_trust_region_IQI_name(Vector< std::string > aIQINames) override;

            // ----------------------------------------------------------------------------

            Vector< moris::Matrix< DDRMat > >& get_IQI() override
            {
                return mIQIVal;

            }

            // ----------------------------------------------------------------------------

            bool is_unit_test() override
            {
                return true;
            }

            // ----------------------------------------------------------------------------

            uint
            get_max_num_global_dofs() override
            {
                moris::uint tNumMyDofs = mNumMyDofs;

                // sum up all distributed dofs
                moris::uint tMaxNumGlobalDofs = sum_all( tNumMyDofs );

                return tMaxNumGlobalDofs;
            };

            // ----------------------------------------------------------------------------------------------

            Matrix< DDSMat >
            get_my_local_global_map() override
            {
                return mMyGlobalElements;
            };
            // ----------------------------------------------------------------------------

            // local-to-global map
            Matrix< DDSMat >
            get_my_local_global_map( const Vector< enum MSI::Dof_Type >& aListOfDofTypes ) override
            {
                return mMyGlobalElements;
            };

            // ----------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_my_local_global_overlapping_map() override
            {
                return mMyGlobalElementsOverlapping;
            };

            // ----------------------------------------------------------------------------------------------

            // number of elements on proc
            uint
            get_num_my_elements() override
            {
                return mNumElements;
            };

            // ----------------------------------------------------------------------------

            uint
            get_num_sets() override
            {
                return 1;
            };

            // ----------------------------------------------------------------------------

            uint
            get_num_equation_objects_on_set( uint aBlockInd ) override
            {
                return mNumElements;
            };
            
            fem::Element_Type
            get_set_type( uint aMyEquSetInd ) override
            {
                return fem::Element_Type::END_ELEMENT_TYPE;
            }

            // ----------------------------------------------------------------------------------------------

            void get_equation_object_operator(
                    const uint&       aMyElementInd,
                    Matrix< DDRMat >& aElementMatrix ) override;

            // ----------------------------------------------------------------------------

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
                aElementTopology = mFunctionTopology( mNX, mNY, aMyElementInd );
            };

            // ----------------------------------------------------------------------------

            void
            get_element_topology(
                    const uint&       aMyBlockInd,
                    const uint&       aMyElementInd,
                    Matrix< DDSMat >& aElementTopology ) override
            {
                aElementTopology = mFunctionTopology( mNX, mNY, aMyElementInd );
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

            // ----------------------------------------------------------------------------

            void get_equation_object_rhs(
                    const uint&                 aMyBlockInd,
                    const uint&                 aMyElementInd,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override;

            //------------------------------------------------------------------------------

            void get_equation_object_operator_and_rhs(
                    const moris::uint&          aMyElementInd,
                    Matrix< DDRMat >&           aElementMatrix,
                    Vector< Matrix< DDRMat > >& aElementRHS ) override;

            // ----------------------------------------------------------------------------

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
                    return nullptr;
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
    }    // namespace NLA
}    // namespace moris
#endif /* SRC_DISTLINALG_CL_NLA_SOLVER_INPUT_TEST_HPP_ */
