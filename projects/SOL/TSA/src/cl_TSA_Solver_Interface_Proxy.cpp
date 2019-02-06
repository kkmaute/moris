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
//    std::cout<<*mSolutionVector->get_vector()<<std::endl;
    //print(mMySolVec,"mMySolVec");

        aElementRHS.resize(1,1);
        aElementRHS(0,0)= ( mk + 1/(  mDeltaT ) ) * mMySolVec( 0,0 ) - mMySolVecPrev( 0, 0 )/( mDeltaT ) - mk * std::cos( mT );

}


