/*
 * cl_NLA_Solver_Interface_Proxy.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Solver_Interface_Proxy2.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Vector.hpp"

using namespace moris;
using namespace NLA;

NLA_Solver_Interface_Proxy_II::NLA_Solver_Interface_Proxy_II()
{
}

void NLA_Solver_Interface_Proxy_II::set_solution_vector( Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;

    mSolutionVector->extract_copy( mMySolVec );
}

void NLA_Solver_Interface_Proxy_II::get_equation_object_rhs( const uint               & aMyElementInd,
                                                     Cell< Matrix< DDRMat > > & aElementRHS )
{
//    std::cout<<*mSolutionVector->get_vector()<<std::endl;
    //print(mMySolVec,"mMySolVec");
	aElementRHS.resize(1);
    if( mListOfDofTypes.size() == 1)
    {
        aElementRHS(0).resize(2,1);
        aElementRHS(0)(0,0)=(0.4 - 10*mMySolVec( 0, 0 ) - 0.4*std::pow(mMySolVec( 1, 0 ),3) + 5*std::pow(mMySolVec( 1, 0 ),2)  -1*mMySolVec( 3, 0 ) -1*mMySolVec( 2, 0 ));
        aElementRHS(0)(1,0)=(0.15 - 0.4*std::pow(mMySolVec( 0, 0 ),3) + 3*std::pow(mMySolVec( 0, 0 ),2) - 10*mMySolVec( 1, 0 )  -1*mMySolVec( 2, 0 ));
    }
    else if( mListOfDofTypes.size() == 2)
    {
        aElementRHS(0).resize(2,1);
        aElementRHS(0)(0,0)=(0.4 - 10*mMySolVec( 2, 0 ) - 0.4*std::pow(mMySolVec( 3, 0 ),3) + 5*std::pow(mMySolVec( 3, 0 ),2) -1*mMySolVec( 0, 0 ));
        aElementRHS(0)(1,0)=(0.15 - 0.4*std::pow(mMySolVec( 2, 0 ),3) + 3*std::pow(mMySolVec( 2, 0 ),2) - 10*mMySolVec( 3, 0 ) -2*mMySolVec( 1, 0 ));
    }
    else if( mListOfDofTypes.size() == 3)
    {
        MORIS_ERROR(false, "NLA_Solver_Interface_Proxy_II::get_equation_object_rhs");
    }
}

void NLA_Solver_Interface_Proxy_II::get_equation_object_rhs( const uint               & aMyBlockInd,
                                                     const uint               & aMyElementInd,
                                                     Cell< Matrix< DDRMat > > & aElementRHS )
{
//    std::cout<<*mSolutionVector->get_vector()<<std::endl;
    //print(mMySolVec,"mMySolVec");
	aElementRHS.resize(1);
    if( mListOfDofTypes.size() == 1)
    {
        aElementRHS(0).resize(2,1);
        aElementRHS(0)(0,0)=(0.4 - 10*mMySolVec( 0, 0 ) - 0.4*std::pow(mMySolVec( 1, 0 ),3) + 5*std::pow(mMySolVec( 1, 0 ),2)  -1*mMySolVec( 3, 0 ) -1*mMySolVec( 2, 0 ));
        aElementRHS(0)(1,0)=(0.15 - 0.4*std::pow(mMySolVec( 0, 0 ),3) + 3*std::pow(mMySolVec( 0, 0 ),2) - 10*mMySolVec( 1, 0 )  -1*mMySolVec( 2, 0 ));
    }
    else if( mListOfDofTypes.size() == 2)
    {
        aElementRHS(0).resize(2,1);
        aElementRHS(0)(0,0)=(0.4 - 10*mMySolVec( 2, 0 ) - 0.4*std::pow(mMySolVec( 3, 0 ),3) + 5*std::pow(mMySolVec( 3, 0 ),2) -1*mMySolVec( 0, 0 ));
        aElementRHS(0)(1,0)=(0.15 - 0.4*std::pow(mMySolVec( 2, 0 ),3) + 3*std::pow(mMySolVec( 2, 0 ),2) - 10*mMySolVec( 3, 0 ) -2*mMySolVec( 1, 0 ));
    }
    else if( mListOfDofTypes.size() == 3)
    {
        MORIS_ERROR(false, "NLA_Solver_Interface_Proxy_II::get_equation_object_rhs");
    }
}


