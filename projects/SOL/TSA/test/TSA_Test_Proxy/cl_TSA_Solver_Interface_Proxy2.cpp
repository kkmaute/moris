/*
 * cl_TSA_Solver_Interface_Proxy2.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy2.hpp"
#include "cl_Vector.hpp"
#include "fn_print.hpp"

using namespace moris;
using namespace NLA;
using namespace tsa;

TSA_Solver_Interface_Proxy_II::TSA_Solver_Interface_Proxy_II()
{
}

void TSA_Solver_Interface_Proxy_II::set_solution_vector( Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;

    mSolutionVector->extract_copy( mMySolVec );
}

void TSA_Solver_Interface_Proxy_II::set_solution_vector_prev_time_step( Dist_Vector * aSolutionVector )
{
    mSolutionVectorPrev = aSolutionVector;

    mSolutionVectorPrev->extract_copy( mMySolVecPrev );
}

void TSA_Solver_Interface_Proxy_II::get_element_rhs( const uint             & aMyElementInd,
                                                           Matrix< DDRMat > & aElementRHS )
{
    //std::cout<<*mSolutionVector->get_vector()<<std::endl;
    //print(mMySolVecPrev,"mMySolVecPrev");
//    if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP && mListOfDofTypes( 1 ) == MSI::Dof_Type::UX)
//    {
//        MORIS_ERROR( false, "get_element_rhs");
//    }
//    else
        if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP)
    {
        aElementRHS.resize(1,1);
        aElementRHS(0,0)= ( mk + 1/(  mDeltaT ) ) * mMySolVec( 0,0 ) - mMySolVecPrev( 2, 0 )/( mDeltaT ) - mk * std::cos( mT );
    }
    else if( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX)
    {
        aElementRHS.resize(1,1);
        aElementRHS(0,0)= ( mk + 1/(  mDeltaT ) ) * mMySolVec( 1,0 ) - mMySolVecPrev( 3, 0 )/( mDeltaT ) - mk * std::cos( mT );
    }

}
 moris::Matrix< DDSMat > & TSA_Solver_Interface_Proxy_II::get_time_level_Ids_minus()
{
     if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP)
     {
         mTimeLevelIdsMinus.set_size( 1, 1, 0 );
     }
     else if( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX)
     {
         mTimeLevelIdsMinus.set_size( 1, 1, 1 );
     }

     return mTimeLevelIdsMinus;

}
 moris::Matrix< DDSMat > & TSA_Solver_Interface_Proxy_II::get_time_level_Ids_plus()
{
     if( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP)
     {
         mTimeLevelIdsPlus.set_size( 1, 1, 2 );
     }
     else if( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX)
     {
         mTimeLevelIdsPlus.set_size( 1, 1, 3 );
     }

    return mTimeLevelIdsPlus;
}

