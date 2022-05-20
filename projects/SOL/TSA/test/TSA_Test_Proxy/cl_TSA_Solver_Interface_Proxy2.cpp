/*
 * cl_TSA_Solver_Interface_Proxy2.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_Communication_Tools.hpp"    // COM/src
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_TSA_Solver_Interface_Proxy2.hpp"
#include "cl_SOL_Dist_Vector.hpp"

using namespace moris;
using namespace NLA;
using namespace tsa;

// ----------------------------------------------------------------------------------------------

TSA_Solver_Interface_Proxy_II::TSA_Solver_Interface_Proxy_II()
{
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::set_solution_vector( sol::Dist_Vector* aSolutionVector )
{
    mSolutionVector = aSolutionVector;
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector )
{
    mSolutionVectorPrev = aSolutionVector;

    mSolutionVectorPrev->extract_copy( mMySolVecPrev );
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::get_equation_object_rhs( const uint& aMyElementInd,
        Cell< Matrix< DDRMat > >&                                   aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
    aElementRHS.resize( 1 );
    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 0, 0 ) - mMySolVecPrev( 2, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 1, 0 ) - mMySolVecPrev( 3, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::get_equation_object_rhs(
        const uint&               aMyBlockInd,
        const uint&               aMyElementInd,
        Cell< Matrix< DDRMat > >& aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
    aElementRHS.resize( 1 );
    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 0, 0 ) - mMySolVecPrev( 2, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 1, 0 ) - mMySolVecPrev( 3, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::get_equation_object_operator(
        const uint&       aMyElementInd,
        Matrix< DDRMat >& aElementMatrix )
{
    mSolutionVector->extract_copy( mMySolVec );

    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::get_equation_object_operator(
        const uint&       aMyBlockInd,
        const uint&       aMyElementInd,
        Matrix< DDRMat >& aElementMatrix )
{
    mSolutionVector->extract_copy( mMySolVec );

    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );
    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::get_equation_object_operator_and_rhs(
        const moris::uint&        aMyElementInd,
        Matrix< DDRMat >&         aElementMatrix,
        Cell< Matrix< DDRMat > >& aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );

    aElementRHS.resize( 1 );

    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 0, 0 ) - mMySolVecPrev( 2, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 1, 0 ) - mMySolVecPrev( 3, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }

    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
}

// ----------------------------------------------------------------------------------------------

void
TSA_Solver_Interface_Proxy_II::get_equation_object_operator_and_rhs(
        const moris::uint&        aMyEquSetInd,
        const moris::uint&        aMyElementInd,
        Matrix< DDRMat >&         aElementMatrix,
        Cell< Matrix< DDRMat > >& aElementRHS )
{
    mSolutionVector->extract_copy( mMySolVec );

    mDeltaT = mT( 1, 0 ) - mT( 0, 0 );

    aElementRHS.resize( 1 );

    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 0, 0 ) - mMySolVecPrev( 2, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementRHS( 0 ).resize( 1, 1 );
        aElementRHS( 0 )( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) ) * mMySolVec( 1, 0 ) - mMySolVecPrev( 3, 0 ) / (mDeltaT)-mk * std::cos( mT( 1, 0 ) );
    }

    if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::TEMP )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
    else if ( mListOfDofTypes( 0 ) == MSI::Dof_Type::UX )
    {
        aElementMatrix.resize( 1, 1 );
        aElementMatrix( 0, 0 ) = ( mk + 1.0 / ( mDeltaT ) );
    }
}
