/*
 * cl_TSA_Solver_Interface_Proxy_.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy.hpp"
#include "cl_SOL_Dist_Vector.hpp"

using namespace moris;
using namespace NLA;
using namespace tsa;

TSA_Solver_Interface_Proxy::TSA_Solver_Interface_Proxy()
{
}

void TSA_Solver_Interface_Proxy::set_solution_vector( sol::Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;
}

void TSA_Solver_Interface_Proxy::set_solution_vector_prev_time_step( sol::Dist_Vector * aSolutionVector )
{
    mSolutionVectorPrev = aSolutionVector;

    mSolutionVectorPrev->extract_copy( mMySolVecPrev );
}

void TSA_Solver_Interface_Proxy::get_equation_object_rhs( const uint                     & aMyElementInd,
                                                        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );
    Matrix< DDRMat > tMat;
    Matrix< DDSMat > tMatRows1 = this->get_time_level_Ids_minus();
    Matrix< DDSMat > tMatRows2 = this->get_time_level_Ids_plus();

    //print(mMySolVecPrev,"mMySolVecPrev");
    //print(mMySolVec,"mMySolVec");
    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
    aElementRHS.resize(1);
    aElementRHS(0).resize(1,1);
    aElementRHS(0)(0,0)= ( mk + 1/(  mDeltaT ) ) * mMySolVec( 0,0 ) - mMySolVecPrev( 1, 0 )/( mDeltaT ) - mk * std::cos( mT( 1, 0 ) );
}

void TSA_Solver_Interface_Proxy::get_equation_object_rhs( const uint                     & aMyBlockInd,
                                                  const uint                     & aMyElementInd,
                                                        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    Matrix< DDRMat > tMat;
    Matrix< DDSMat > tMatRows1 = this->get_time_level_Ids_minus();
    Matrix< DDSMat > tMatRows2 = this->get_time_level_Ids_plus();

    //print(mMySolVecPrev,"mMySolVecPrev");
    //print(mMySolVec,"mMySolVec");
    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
    aElementRHS.resize(1);
    aElementRHS(0).resize(1,1);
    aElementRHS(0)(0,0)= ( mk + 1/(  mDeltaT ) ) * mMySolVec( 0,0 ) - mMySolVecPrev( 1, 0 )/( mDeltaT ) - mk * std::cos( mT( 1, 0 ) );
}

 moris::Matrix< DDSMat > & TSA_Solver_Interface_Proxy::get_time_level_Ids_minus()
{
    mTimeLevelIdsMinus.set_size( 1, 1, 0 );
    return mTimeLevelIdsMinus;
}
 moris::Matrix< DDSMat > & TSA_Solver_Interface_Proxy::get_time_level_Ids_plus()
{
    mTimeLevelIdsPlus.set_size( 1, 1 , 1 );
    return mTimeLevelIdsPlus;
}

 void TSA_Solver_Interface_Proxy::perform_mapping()
 {
     moris::Cell< Matrix< DDRMat > > tMat;
     Matrix< DDSMat > tMatRows1 = this->get_time_level_Ids_minus();
     Matrix< DDSMat > tMatRows2 = this->get_time_level_Ids_plus();

     mSolutionVectorPrev->extract_my_values( 1, tMatRows1, 0 , tMat );

     mSolutionVectorPrev->sum_into_global_values( tMatRows2, tMat( 0 ) );

     mSolutionVectorPrev->vector_global_asembly();
 }

 void TSA_Solver_Interface_Proxy::get_equation_object_operator(const uint             & aMyElementInd,
                                         Matrix< DDRMat > & aElementMatrix)
 {
     mSolutionVector->extract_copy( mMySolVec );

     mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
     aElementMatrix.resize(1, 1);
     aElementMatrix(0,0)=( mk + 1/( mDeltaT) );
 }

 void TSA_Solver_Interface_Proxy::get_equation_object_operator( const uint             & aMyBlockInd,
                                    const uint             & aMyElementInd,
                                          Matrix< DDRMat > & aElementMatrix)
 {
     mSolutionVector->extract_copy( mMySolVec );

     mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
     aElementMatrix.resize(1, 1);
     aElementMatrix(0,0)=( mk + 1/( mDeltaT) );
 }

