/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Interface_Proxy.cpp
 *
 */

#include "cl_NLA_Solver_Interface_Proxy.hpp"
#include "cl_Communication_Tools.hpp"    // COM/src
#include "cl_SOL_Dist_Vector.hpp"

using namespace moris;
using namespace NLA;

// ----------------------------------------------------------------------------------------------

NLA_Solver_Interface_Proxy::NLA_Solver_Interface_Proxy()
{
}

// ----------------------------------------------------------------------------------------------

NLA_Solver_Interface_Proxy::NLA_Solver_Interface_Proxy(
        const uint        aNumMyDofs,
        const uint        aNumElements,
        const moris::sint aNX,
        const moris::sint aNY,
        Matrix< DDRMat > ( *aFunctionRes )( const moris::sint aNX, const moris::sint aNY, const real aLambda, const Matrix< DDRMat > &tMyValues, const uint aEquationObjectInd ),
        Matrix< DDRMat > ( *aFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat > &tMyValues, const uint aEquationObjectInd ),
        Matrix< DDSMat > ( *aFunctionTopo )( const moris::sint aNX, const moris::sint aNY, const uint aEquationObjectInd ) )
{
    mUseMatrixMarketFiles = false;

    mFunctionRes      = aFunctionRes;
    mFunctionJac      = aFunctionJac;
    mFunctionTopology = aFunctionTopo;

    mNX = aNX;
    mNY = aNY;

    mNumMyDofs   = aNumMyDofs / par_size();
    mNumElements = aNumElements;

    // mMyGlobalElements.resize( mNumMyDofs, 1 );
    mMyGlobalElements.resize( mNumMyDofs, 1 );

    moris::sint tRank = par_rank();

    for ( uint Ik = ( mNumMyDofs * tRank ); Ik < mNumMyDofs * ( tRank + 1 ); Ik++ )
    {
        mMyGlobalElements( Ik - ( mNumMyDofs * tRank ), 0 ) = Ik;
    }

    mMyGlobalElementsOverlapping.resize( aNumMyDofs, 1 );

    for ( uint Ik = 0; Ik < aNumMyDofs; Ik++ )
    {
        mMyGlobalElementsOverlapping( Ik, 0 ) = Ik;
    }
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_solution_vector( sol::Dist_Vector *aSolutionVector )
{
    mSolutionVector = aSolutionVector;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_time( const Matrix< DDRMat > &aTime )
{
    mTime = aTime;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_solution_vector_prev_time_step( sol::Dist_Vector *aSolutionVector )
{
    mSolutionVectorPrev = aSolutionVector;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_rhs(
        const uint                 &aMyElementInd,
        Vector< Matrix< DDRMat > > &aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementRHS = { mFunctionRes( mNX, mNY, mTime( 1 ), mMySolVec, aMyElementInd ) };
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_rhs(
        const uint                 &aMyBlockInd,
        const uint                 &aMyElementInd,
        Vector< Matrix< DDRMat > > &aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementRHS = { mFunctionRes( mNX, mNY, mTime( 1 ), mMySolVec, aMyElementInd ) };
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator(
        const uint       &aMyElementInd,
        Matrix< DDRMat > &aElementMatrix )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator(
        const uint       &aMyBlockInd,
        const uint       &aMyElementInd,
        Matrix< DDRMat > &aElementMatrix )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator_and_rhs(
        const uint                 &aMyElementInd,
        Matrix< DDRMat >           &aElementMatrix,
        Vector< Matrix< DDRMat > > &aElementRHS )
{
    aElementRHS = { mFunctionRes( mNX, mNY, mTime( 1 ), mMySolVec, aMyElementInd ) };

    mSolutionVector->extract_copy( mMySolVec );

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator_and_rhs(
        const uint                 &aMyEquSetInd,
        const uint                 &aMyElementInd,
        Matrix< DDRMat >           &aElementMatrix,
        Vector< Matrix< DDRMat > > &aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementRHS = { mFunctionRes( mNX, mNY, mTime( 1 ), mMySolVec, aMyElementInd ) };

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}
