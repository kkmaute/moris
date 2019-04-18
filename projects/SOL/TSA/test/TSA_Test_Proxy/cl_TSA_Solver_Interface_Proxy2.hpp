/*
 * cl_TSA_Solver_Interface_Proxy_@.hpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_TSA_SOLVER_INTERFACE_PROXY_HPP_
#define SRC_DISTLINALG_CL_TSA_SOLVER_INTERFACE_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
class Dist_Vector;
namespace tsa
{
    //class Nonlinear_Solver;
    class TSA_Solver_Interface_Proxy_II : public Solver_Interface
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
        Dist_Vector * mSolutionVectorPrev;
        Matrix< DDRMat > mMySolVec;
        Matrix< DDRMat > mMySolVecPrev;

        moris::sint mNX;
        moris::sint mNY;

        moris::Cell< enum MSI::Dof_Type > mListOfDofTypes;

        moris::real mk = 2;
        Matrix< DDRMat> mT;
        moris::real mDeltaT = 0.0;

        Matrix< DDSMat > mTimeLevelIdsMinus;
        Matrix< DDSMat > mTimeLevelIdsPlus;
    public :
        TSA_Solver_Interface_Proxy_II();

        // ----------------------------------------------------------------------------------------------

        ~TSA_Solver_Interface_Proxy_II(){};

        // ----------------------------------------------------------------------------------------------

        void set_solution_vector( Dist_Vector * aSolutionVector );

        void set_solution_vector_prev_time_step( Dist_Vector * aSolutionVector );

        void set_time( const Matrix< DDRMat> & aTime )
        {
            mT = aTime;
        }
        // ----------------------------------------------------------------------------------------------

        void set_requested_dof_types( const moris::Cell< enum MSI::Dof_Type > aListOfDofTypes )
        {
           mListOfDofTypes = aListOfDofTypes;
        };

        // local dimension of the problem
        uint get_max_num_global_dofs(){ return 4; };


        // ----------------------------------------------------------------------------------------------
        // local dimension of the problem
        uint get_num_my_dofs(){ return mNumMyDofs; };

        // ----------------------------------------------------------------------------------------------
        // local-to-global map
        Matrix< DDSMat > get_my_local_global_map()
        {
            if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP)
           {
                   mMyGlobalElements.resize(1,1);
                    mMyGlobalElements(0,0)=0;
           }
           else if( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX)
           {
               mMyGlobalElements.resize(1,1);
               mMyGlobalElements(0,0)=1;
           }

            return mMyGlobalElements;
        };

        // ----------------------------------------------------------------------------------------------

        moris::Matrix< DDSMat > get_my_local_global_overlapping_map( )
        {
            mMyGlobalElementsOverlapping.resize(4,1);
            mMyGlobalElementsOverlapping(0,0)=0;    mMyGlobalElementsOverlapping(1,0)=1;
            mMyGlobalElementsOverlapping(2,0)=2;    mMyGlobalElementsOverlapping(3,0)=3;

            return mMyGlobalElementsOverlapping;
        };

        moris::Matrix< DDSMat > & get_time_level_Ids_minus();

        moris::Matrix< DDSMat > & get_time_level_Ids_plus() ;

        // ----------------------------------------------------------------------------------------------
        // number of elements on proc
        uint get_num_my_elements()
        {
            return mNumElements=1;
        };

        void perform_mapping();

        // ----------------------------------------------------------------------------------------------
        void get_element_matrix(const uint             & aMyElementInd,
                                      Matrix< DDRMat > & aElementMatrix)
        {
            mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
            if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP)
            {
                aElementMatrix.resize(1, 1);
                aElementMatrix(0,0)=( mk + 1/( mDeltaT) );
            }
            else if( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX)
            {
                aElementMatrix.resize(1, 1);
                aElementMatrix(0,0)=( mk + 1/( mDeltaT) );
            }
        };

        // ----------------------------------------------------------------------------------------------
        void  get_element_topology(const uint             & aMyElementInd,
                                         Matrix< DDSMat > & aElementTopology)
        {
//            if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP && mListOfDofTypes( 1 ) == MSI::Dof_Type::UX)
//            {
//                MORIS_ERROR( false, "get_element_topology");
//            }
//            else
            if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
            {
                aElementTopology.resize(1,1);
                aElementTopology(0,0)=0;
            }
            else if( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
            {
                aElementTopology.resize(1,1);
                aElementTopology(0,0)=1;
            }
            else
            {
                MORIS_ERROR( false, "get_element_topology");
            }
        };

        // ----------------------------------------------------------------------------------------------
        Matrix< DDUMat > get_constr_dof(){ return mMyConstraintDofs; };

        // ----------------------------------------------------------------------------------------------
        void get_element_rhs( const uint             & aMyElementInd,
                                    Matrix< DDRMat > & aElementRHS );

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
#endif /* SRC_DISTLINALG_CL_TSA_SOLVER_INTERFACE_PROXY_HPP_ */
