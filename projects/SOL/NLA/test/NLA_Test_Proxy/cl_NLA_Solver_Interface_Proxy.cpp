/*
 * cl_NLA_Solver_Interface_Proxy.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Solver_Interface_Proxy.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_SOL_Dist_Vector.hpp"

using namespace moris;
using namespace NLA;

// ----------------------------------------------------------------------------------------------

NLA_Solver_Interface_Proxy::NLA_Solver_Interface_Proxy()
{
}

// ----------------------------------------------------------------------------------------------

NLA_Solver_Interface_Proxy::NLA_Solver_Interface_Proxy(
        const moris::uint aNumMyDofs,
        const moris::uint aNumElements,
        const moris::sint aNX,
        const moris::sint aNY,
        Matrix< DDRMat > ( *aFunctionRes )( const moris::sint aNX, const moris::sint aNY, const moris::real aLambda, const Matrix< DDRMat > & tMyValues, const moris::uint aEquationObjectInd ),
        Matrix< DDRMat > ( *aFunctionJac )( const moris::sint aNX, const moris::sint aNY, const Matrix< DDRMat > & tMyValues, const moris::uint aEquationObjectInd ),
        Matrix< DDSMat > ( *aFunctionTopo )( const moris::sint aNX, const moris::sint aNY, const moris::uint aEquationObjectInd ) )
{
    mUseMatrixMarketFiles = false;

    mFunctionRes = aFunctionRes;
    mFunctionJac = aFunctionJac;
    mFunctionTopology = aFunctionTopo;

    mNX = aNX;
    mNY = aNY;

    mNumMyDofs = aNumMyDofs/par_size();
    mNumElements = aNumElements;

    //mMyGlobalElements.resize( mNumMyDofs, 1 );
    mMyGlobalElements.resize( mNumMyDofs, 1 );

    moris::sint tRank = par_rank();

    for ( moris::uint Ik = ( mNumMyDofs * tRank ); Ik < mNumMyDofs * (tRank+1); Ik++ )
    {
        mMyGlobalElements( Ik-( mNumMyDofs*tRank ), 0 ) = Ik;
    }

    mMyGlobalElementsOverlapping.resize( aNumMyDofs, 1 );

    for ( moris::uint Ik = 0; Ik < aNumMyDofs; Ik++ )
    {
        mMyGlobalElementsOverlapping( Ik, 0 ) = Ik;
    }
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_solution_vector( sol::Dist_Vector * aSolutionVector )
{
    mSolutionVector = aSolutionVector;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_time_value(
        const moris::real & aLambda,
        moris::uint   aPos )
{
    mTime(aPos) = aLambda;
}
// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_time( const Matrix< DDRMat> & aTime )
{
    mTime = aTime;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::set_solution_vector_prev_time_step( sol::Dist_Vector * aSolutionVector )
{
    mSolutionVectorPrev = aSolutionVector;
}

// ----------------------------------------------------------------------------------------------

moris::Matrix< DDSMat > & NLA_Solver_Interface_Proxy::get_time_level_Ids_minus()
{
    mTimeLevelIdsMinus.set_size( 1, 1, 0 );
    return mTimeLevelIdsMinus;
}

// ----------------------------------------------------------------------------------------------

moris::Matrix< DDSMat > & NLA_Solver_Interface_Proxy::get_time_level_Ids_plus()
{
    mTimeLevelIdsPlus.set_size( 1, 1 , 1 );
    return mTimeLevelIdsPlus;
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::perform_mapping()
{
    moris::Cell< Matrix< DDRMat > > tMat;
    Matrix< DDSMat > tMatRows1 = this->get_time_level_Ids_minus();
    Matrix< DDSMat > tMatRows2 = this->get_time_level_Ids_plus();

    mSolutionVectorPrev->extract_my_values( 1, tMatRows1, 0 , tMat );

    mSolutionVectorPrev->sum_into_global_values( tMatRows2, tMat( 0 ) );

    mSolutionVectorPrev->vector_global_asembly();
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_rhs(
        const uint                     & aMyElementInd,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementRHS = { mFunctionRes( mNX, mNY, mTime(1), mMySolVec, aMyElementInd ) };
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_rhs(
        const uint                     & aMyBlockInd,
        const uint                     & aMyElementInd,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementRHS = { mFunctionRes( mNX, mNY, mTime(1), mMySolVec, aMyElementInd ) };
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator(
        const uint             & aMyElementInd,
        Matrix< DDRMat > & aElementMatrix)
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator(
        const uint             & aMyBlockInd,
        const uint             & aMyElementInd,
        Matrix< DDRMat > & aElementMatrix)
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator_and_rhs(
        const moris::uint        & aMyElementInd,
        Matrix< DDRMat >         & aElementMatrix,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    aElementRHS = { mFunctionRes( mNX, mNY, mTime(1), mMySolVec, aMyElementInd ) };

    mSolutionVector->extract_copy( mMySolVec );

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}

// ----------------------------------------------------------------------------------------------

void NLA_Solver_Interface_Proxy::get_equation_object_operator_and_rhs(
        const moris::uint        & aMyEquSetInd,
        const moris::uint        & aMyElementInd,
        Matrix< DDRMat >         & aElementMatrix,
        Cell< Matrix< DDRMat > > & aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    aElementRHS = { mFunctionRes( mNX, mNY, mTime(1), mMySolVec, aMyElementInd ) };

    aElementMatrix = mFunctionJac( mNX, mNY, mMySolVec, aMyElementInd );
}
