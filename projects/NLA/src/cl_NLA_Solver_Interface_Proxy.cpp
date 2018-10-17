/*
 * cl_NLA_Solver_Interface_Proxy.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Solver_Interface_Proxy.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_Vector.hpp"

using namespace moris;
using namespace NLA;

NLA_Solver_Interface_Proxy::NLA_Solver_Interface_Proxy()
{
    mUseMatrixMarketFiles = false;

    // Set input values
    mNumMyDofs = 2;
    mNumDofsPerElement = 2;

    mEleDofConectivity.resize( 2, 1 );
    mEleDofConectivity( 0, 0) = 0;   mEleDofConectivity( 1, 0) = 1;

    mMyGlobalElements.resize( mNumMyDofs, 1 );
    mMyGlobalElements(0,0) = 0;    mMyGlobalElements(1,0) = 1;
    mNumElements = 1;
}

void NLA_Solver_Interface_Proxy::set_solution_vector( Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;

    this->set_test_problem();
}

void NLA_Solver_Interface_Proxy::set_test_problem()
{
    Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
    tGlobalIndExtract( 1, 0 ) = 1;
    Matrix< DDRMat > tMyValues;

    mSolutionVector->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues);

//    std::cout<<tMyValues(0,0)<<" Input"<<std::endl;
//    std::cout<<tMyValues(1,0)<<" Input"<<std::endl;

    mElementMatrixValues.resize( 4, 1 );
    mElementMatrixValues( 0, 0 ) = 10;
    mElementMatrixValues( 1, 0 ) = 1.2*std::pow(tMyValues( 0, 0 ),2)-6*tMyValues( 0, 0 );
    mElementMatrixValues( 2, 0 ) = 1.2*std::pow(tMyValues( 1, 0 ),2)-10*tMyValues( 1, 0 );
    mElementMatrixValues( 3, 0 ) = 10;

    mMyRHSValues.resize( 2, 1 );
    mMyRHSValues( 0, 0 ) = 0.4 - 10*tMyValues( 0, 0 ) - 0.4*std::pow(tMyValues( 1, 0 ),3) + 5*std::pow(tMyValues( 1, 0 ),2);
    mMyRHSValues( 1, 0 ) = 0.15 - 0.4*std::pow(tMyValues( 0, 0 ),3) + 3*std::pow(tMyValues( 0, 0 ),2) - 10*tMyValues( 1, 0 );

}

