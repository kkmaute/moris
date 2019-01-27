/*
 * cl_NLA_Solver_Interface_Proxy.hpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_NLA_SOLVER_INPUT_TEST_HPP_
#define SRC_DISTLINALG_CL_NLA_SOLVER_INPUT_TEST_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

namespace moris
{
class Dist_Vector;
namespace NLA
{
    class Nonlinear_Solver;
    class NLA_Solver_Interface_Proxy : public Solver_Interface
    {
    private:
        moris::uint mNumMyDofs;                           // local dimension of the problem
        moris::Matrix< DDSMat > mMyGlobalElements;        // local-to-global map
        moris::Matrix< DDSMat > mMyGlobalElementsOverlapping;        // local-to-global map
        moris::uint mNumElements;                         // number local elements
        moris::Matrix< DDSMat > mEleDofConectivity;       // element - dof conectivities
        moris::Matrix< DDRMat > mElementMatrixValues;     // dense element matrix entries
        moris::Matrix< DDUMat > mMyConstraintDofs;        // constraint dofs
        moris::Matrix< DDRMat > mMyRHSValues;             // Vector with RHS values

        bool mUseMatrixMarketFiles;                       // determines is matrix and RHS comes from MatrixMarket files

        Dist_Vector * mSolutionVector;
        Matrix< DDRMat > mMySolVec;

        Matrix< DDRMat > ( *mFunctionRes )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat > & tMyValues, const moris::uint aEquationObjectInd );
        Matrix< DDRMat > ( *mFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat > & tMyValues, const moris::uint aEquationObjectInd );
        Matrix< DDSMat > ( *mFunctionTopology )( const moris::sint aNX, const moris::sint aNY, const moris::uint aEquationObjectInd );

        moris::sint mNX;
        moris::sint mNY;

        moris::Cell< enum MSI::Dof_Type > mListOfDofTypes;

    public :
        NLA_Solver_Interface_Proxy();

        NLA_Solver_Interface_Proxy( const moris::uint aNumMyDofs,
                                    const moris::uint aNumElements,
                                    const moris::sint aNX,
                                    const moris::sint aNY,
                                    Matrix< DDRMat > ( *aFunctionRes )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat > & tMyValues, const moris::uint aEquationObjectInd ),
                                    Matrix< DDRMat > ( *aFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat > & tMyValues, const moris::uint aEquationObjectInd ),
                                    Matrix< DDSMat > ( *aFunctionTopo )( const moris::sint aNX, const moris::sint aNY, const moris::uint aEquationObjectInd ) );

        NLA_Solver_Interface_Proxy( std::shared_ptr< Nonlinear_Solver > aNewtonSolver ){};

        // ----------------------------------------------------------------------------------------------
        ~NLA_Solver_Interface_Proxy(){};

        // ----------------------------------------------------------------------------------------------

        void set_solution_vector( Dist_Vector * aSolutionVector );

        // ----------------------------------------------------------------------------------------------

        void set_requested_dof_types( const moris::Cell< enum MSI::Dof_Type > aListOfDofTypes )
        {
           mListOfDofTypes = aListOfDofTypes;
        };


        // ----------------------------------------------------------------------------------------------
        // local dimension of the problem
        uint get_num_my_dofs(){ return mNumMyDofs; };

        uint get_max_num_global_dofs()
        {
            moris::uint tNumMyDofs     = mNumMyDofs;
            moris::uint tMaxNumGlobalDofs = mNumMyDofs;

            // sum up all distributed dofs

            MPI_Allreduce(&tNumMyDofs,&tMaxNumGlobalDofs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

            return tMaxNumGlobalDofs;
        };

        // ----------------------------------------------------------------------------------------------
        // local-to-global map
        Matrix< DDSMat > get_my_local_global_map(){ return mMyGlobalElements; };

        moris::Matrix< DDSMat > get_my_local_global_overlapping_map( ){return mMyGlobalElementsOverlapping; };

        // ----------------------------------------------------------------------------------------------
        // number of elements on proc
        uint get_num_my_elements(){return mNumElements; };

        // ----------------------------------------------------------------------------------------------
        void get_element_matrix(const uint             & aMyElementInd,
                                      Matrix< DDRMat > & aElementMatrix)
        {
            aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
        };

        // ----------------------------------------------------------------------------------------------
        void  get_element_topology(const uint             & aMyElementInd,
                                         Matrix< DDSMat > & aElementTopology)
        {
            aElementTopology = mFunctionTopology( mNX, mNY, aMyElementInd );
        };

        // ----------------------------------------------------------------------------------------------
        Matrix< DDUMat > get_constr_dof(){ return mMyConstraintDofs; };

        // ----------------------------------------------------------------------------------------------
        void get_element_rhs( const uint             & aMyElementInd,
                                    Matrix< DDRMat > & aElementRHS )
        {
            aElementRHS = mFunctionRes( mNX, mNY, mMySolVec, aMyElementInd );
        };

        // ----------------------------------------------------------------------------------------------

        void use_matrix_market_files( )
        {
            mUseMatrixMarketFiles = true;
        };

        // ----------------------------------------------------------------------------------------------

        const char* get_matrix_market_path( )
        {
            if ( mUseMatrixMarketFiles == true )
            {
                const char* tFilePath ="/home/schmidt/codes/MORIS/test/src/distlinalg/";
                return tFilePath;
            }
            else
            {
                return NULL;
            }
        };
    };
}
}
#endif /* SRC_DISTLINALG_CL_NLA_SOLVER_INPUT_TEST_HPP_ */
