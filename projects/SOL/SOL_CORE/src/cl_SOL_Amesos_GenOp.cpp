/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *
 */
/* ------------------------------------------------------------------------- */

#include "cl_SOL_Amesos_GenOp.hpp"
#include "cl_Communication_Tools.hpp"

// C++ system files
#ifdef PARALLEL
#include "mpi.h"
#endif

// ----------------------------------------------------------------------------

Amesos_GenOp::Amesos_GenOp(
        const Teuchos::RCP< Amesos_BaseSolver >& aSolver,
        const Teuchos::RCP< Epetra_Operator >&   aMassMatrix,
        const bool                               tUseTranspose )
        : mSolver( aSolver )
        , mMassMatrix( aMassMatrix )
        , mEpetraProblem( nullptr )
        , mUseTranspose( tUseTranspose )
{
    if ( mSolver.is_null() )
    {
        MORIS_ERROR( false, "Amesos_GenOp constructor: mSolver is null" );
    }
    if ( mMassMatrix.is_null() )
    {
        MORIS_ERROR( false, "Amesos_GenOp constructor: mMassMatrix is null" );
    }
    
    mEpetraProblem = const_cast<Epetra_LinearProblem*> (mSolver->GetProblem ());

    if ( mEpetraProblem == nullptr )
    {
        MORIS_ERROR( false, "Amesos_GenOp constructor: mEpetraProblem is nulll" );
    }

    if ( mSolver->UseTranspose() )
    {
        mSolver->SetUseTranspose( !tUseTranspose );
    }
    else
    {
        mSolver->SetUseTranspose( tUseTranspose );
    }

    if ( mMassMatrix->UseTranspose() )
    {
        mMassMatrix->SetUseTranspose( !tUseTranspose );
    }
    else
    {
        mMassMatrix->SetUseTranspose( tUseTranspose );
    }
}

// ----------------------------------------------------------------------------

int
Amesos_GenOp::SetUseTranspose( bool tUseTranspose )
{
    int tErr = 0;

    if ( mEpetraProblem == nullptr )
    {
        MORIS_ERROR( false, "Amesos_GenOp::set_use_transpose: mEpetraProblem is NULL" );
    }
    if ( mMassMatrix.is_null() )
    {
        MORIS_ERROR( false, "Amesos_GenOp::set_use_transpose: mMassMatrix is NULL" );
    }
    if ( mSolver.is_null() )
    {
        MORIS_ERROR( false, "Amesos_GenOp::set_use_transpose: mSolver is NULL" );
    }

    const bool tSolverUseTranspose = mSolver->UseTranspose();

    if ( tSolverUseTranspose )
    {
        tErr = mSolver->SetUseTranspose( !tUseTranspose );
    }
    else
    {
        tErr = mSolver->SetUseTranspose( tUseTranspose );
    }

    // If set_use_transpose returned zero above, then the Amesos solver
    // doesn't know how to change the transpose state.
    if ( tErr != 0 )
    {
        return tErr;
    }

    if ( mMassMatrix->UseTranspose() )
    {
        tErr = mMassMatrix->SetUseTranspose( !tUseTranspose );
    }
    else
    {
        tErr = mMassMatrix->SetUseTranspose( tUseTranspose );
    }

    // If SetUseTranspose returned zero above, then the mass matrix
    // doesn't know how to change the transpose state.
    if ( tErr != 0 )
    {
        // Put the solver back like we found it.
        (void)mSolver->SetUseTranspose( tSolverUseTranspose );
        return tErr;
    }

    mUseTranspose = tUseTranspose;
    return 0;    // the method completed correctly
}

// ----------------------------------------------------------------------------

int
Amesos_GenOp::Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    if ( mEpetraProblem == nullptr )
    {
        MORIS_ERROR( false, "Amesos_GenOp::apply: mEpetraProblem is NULL" );
    }
    if ( mMassMatrix.is_null() )
    {
        MORIS_ERROR( false, "Amesos_GenOp::apply: mMassMatrix is NULL" );
    }
    if ( mSolver.is_null() )
    {
        MORIS_ERROR( false, "Amesos_GenOp::apply: mSolver is NULL" );
    }

    if ( !mUseTranspose )
    {
        // Storage for M*X
        Epetra_MultiVector MX( X.Map(), X.NumVectors() );

        // Apply M*X
        mMassMatrix->Apply( X, MX );
        Y.PutScalar( 0.0 );

        // Set the LHS and RHS
        mEpetraProblem->SetRHS( &MX );
        mEpetraProblem->SetLHS( &Y );

        // Solve the linear system A*Y = MX
        mSolver->Solve();
    }
    else    // apply the transposed operator
    {
        // Storage for A^{-T}*X
        Epetra_MultiVector ATX( X.Map(), X.NumVectors() );
        Epetra_MultiVector tmpX = const_cast< Epetra_MultiVector& >( X );

        // Set the LHS and RHS
        mEpetraProblem->SetRHS( &tmpX );
        mEpetraProblem->SetLHS( &ATX );

        // Solve the linear system A^T*Y = X
        mSolver->Solve();

        // Apply M*ATX
        mMassMatrix->Apply( ATX, Y );
    }

    return 0;    // the method completed correctly
}
