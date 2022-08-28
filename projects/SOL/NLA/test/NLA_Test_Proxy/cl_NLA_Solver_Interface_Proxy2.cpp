/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Interface_Proxy2.cpp
 *
 */

#include "cl_NLA_Solver_Interface_Proxy2.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_SOL_Dist_Vector.hpp"

using namespace moris;
using namespace NLA;

// ----------------------------------------------------------------------------------------------

NLA_Solver_Interface_Proxy_II::NLA_Solver_Interface_Proxy_II()
{
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy_II::set_solution_vector( sol::Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy_II::get_equation_object_rhs(
        const uint               & aMyElementInd,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

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

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy_II::get_equation_object_rhs(
        const uint               & aMyBlockInd,
        const uint               & aMyElementInd,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

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

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy_II::get_equation_object_operator(
        const uint       & aMyElementInd,
        Matrix< DDRMat > & aElementMatrix)
{
    mSolutionVector->extract_copy( mMySolVec );

    if( mListOfDofTypes.size() == 1)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 0, 0 ),2)+6*mMySolVec( 0, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 1, 0 ),2)+10*mMySolVec( 1, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 2)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 2, 0 ),2)+6*mMySolVec( 2, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 3, 0 ),2)+10*mMySolVec( 3, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 3)
    {
        MORIS_ERROR(false,"NLA_Node_Proxy_II::get_equation_object_operator: not defined");
    }
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy_II::get_equation_object_operator(
        const uint             & aMyBlockInd,
        const uint             & aMyElementInd,
        Matrix< DDRMat > & aElementMatrix)
{
    mSolutionVector->extract_copy( mMySolVec );

    if( mListOfDofTypes.size() == 1)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 0, 0 ),2)+6*mMySolVec( 0, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 1, 0 ),2)+10*mMySolVec( 1, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 2)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 2, 0 ),2)+6*mMySolVec( 2, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 3, 0 ),2)+10*mMySolVec( 3, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 3)
    {
        MORIS_ERROR(false,"NLA_Node_Proxy_II::get_equation_object_operator: not defined");
    }
}

void NLA_Solver_Interface_Proxy_II::get_equation_object_operator_and_rhs(
        const moris::uint        & aMyElementInd,
        Matrix< DDRMat >         & aElementMatrix,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

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

    if( mListOfDofTypes.size() == 1)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 0, 0 ),2)+6*mMySolVec( 0, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 1, 0 ),2)+10*mMySolVec( 1, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 2)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 2, 0 ),2)+6*mMySolVec( 2, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 3, 0 ),2)+10*mMySolVec( 3, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 3)
    {
        MORIS_ERROR(false,"NLA_Node_Proxy_II::get_equation_object_operator: not defined");
    }
}

void NLA_Solver_Interface_Proxy_II::get_equation_object_operator_and_rhs(
        const moris::uint        & aMyEquSetInd,
        const moris::uint        & aMyElementInd,
        Matrix< DDRMat >         & aElementMatrix,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

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

    if( mListOfDofTypes.size() == 1)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 0, 0 ),2)+6*mMySolVec( 0, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 1, 0 ),2)+10*mMySolVec( 1, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 2)
    {
        aElementMatrix.resize(2,2);
        aElementMatrix(0,0)=-10;
        aElementMatrix(0,1)=-1.2*std::pow(mMySolVec( 2, 0 ),2)+6*mMySolVec( 2, 0 );
        aElementMatrix(1,0)=-1.2*std::pow(mMySolVec( 3, 0 ),2)+10*mMySolVec( 3, 0 );
        aElementMatrix(1,1)=-10;
    }
    else if( mListOfDofTypes.size() == 3)
    {
        MORIS_ERROR(false,"NLA_Node_Proxy_II::get_equation_object_operator: not defined");
    }
}

