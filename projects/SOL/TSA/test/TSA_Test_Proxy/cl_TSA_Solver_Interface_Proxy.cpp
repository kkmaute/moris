/*
 * cl_TSA_Solver_Interface_Proxy_.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy.hpp"
#include "cl_Vector.hpp"
#include "fn_print.hpp"

using namespace moris;
using namespace NLA;
using namespace tsa;

TSA_Solver_Interface_Proxy::TSA_Solver_Interface_Proxy()
{
}

void TSA_Solver_Interface_Proxy::set_solution_vector( Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;

    mSolutionVector->extract_copy( mMySolVec );
}

void TSA_Solver_Interface_Proxy::set_solution_vector_prev_time_step( Dist_Vector * aSolutionVector )
{
    mSolutionVectorPrev = aSolutionVector;

    mSolutionVectorPrev->extract_copy( mMySolVecPrev );
}

void TSA_Solver_Interface_Proxy::get_element_rhs( const uint             & aMyElementInd,
                                                     Matrix< DDRMat > & aElementRHS )
{
    //print(mMySolVecPrev,"mMySolVecPrev");
    //print(mMySolVec,"mMySolVec");

        aElementRHS.resize(1,1);
        aElementRHS(0,0)= ( mk + 1/(  mDeltaT ) ) * mMySolVec( 0,0 ) - mMySolVecPrev( 1, 0 )/( mDeltaT ) - mk * std::cos( mT );
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

