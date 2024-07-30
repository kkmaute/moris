/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * The Amesos GenOp class.
 *
 */

#ifndef SOLVERS_AMESOS_GEN_OP_HPP_
#define SOLVERS_AMESOS_GEN_OP_HPP_

// C system files

// C++ system files
#include <cstddef>
#include <cassert>

// TPL header files
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// Include header for Amesos solver.
#include "cl_DLA_Linear_Solver_Amesos.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"

// Include selected communicator class required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

class Amesos_GenOp : public virtual Epetra_Operator
{
  private:
    Teuchos::RCP< Amesos_BaseSolver > mSolver;

    Teuchos::RCP< Epetra_Operator > mMassMatrix;

    Epetra_LinearProblem* mEpetraProblem;

    bool mUseTranspose;

  public:
    // Constructor
    Amesos_GenOp(
            const Teuchos::RCP< Amesos_BaseSolver >& aSolver,
            const Teuchos::RCP< Epetra_Operator >&   aMassMatrix,
            const bool                               tUseTranspose = false );

    // Virtual Destructor
    virtual ~Amesos_GenOp() {}

    int Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    // The Operator's (human-readable) label.
    const char*
    Label() const
    {
        return "Operator that applies K^{-1} M or M K^{-T}";
    }

    // Whether to apply \f$K^{-1} M\f$ (false) or \f$K^{-T} M\f$ (true).
    bool
    UseTranspose() const
    {
        return mUseTranspose;
    };

    int SetUseTranspose( bool tUseTranspose );

    // The Operator's communicator.
    const Epetra_Comm&
    Comm() const
    {
        return mSolver->Comm();
    }

    // The Operator's domain Map.
    const Epetra_Map&
    OperatorDomainMap() const
    {
        return mMassMatrix->OperatorDomainMap();
    }

    // The Operator's range Map.
    const Epetra_Map&
    OperatorRangeMap() const
    {
        return mMassMatrix->OperatorRangeMap();
    }

    int
    ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        return -1;
    };

    // NOT IMPLEMENTED: Whether this Operator can compute its infinity norm.
    bool
    HasNormInf() const
    {
        return false;
    }

    double
    NormInf() const
    {
        return -1.0;
    }
};

#endif    // SOLVERS_AMESOS_GEN_OP_HPP_
