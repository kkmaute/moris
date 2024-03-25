/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_ML.cpp
 *
 */
#include "cl_DLA_Linear_Solver_ML.hpp"

// TPL header files
#include "Epetra_ConfigDefs.h"

#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_DLA_Linear_Problem.hpp"

// Teuchos
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// detailed logging
#include "cl_Tracer.hpp"

#include "ml_epetra_preconditioner.h"
#include "ml_struct.h"

using namespace moris;
using namespace dla;

//----------------------------------------------------------------------------------------

Linear_Solver_ML::Linear_Solver_ML( const Parameter_List& aParameterList )
        : Linear_Solver_Algorithm_Trilinos( aParameterList )
{
}

//----------------------------------------------------------------------------------------

bool
Linear_Solver_ML::build_external_preconditioner( const sint& aIter )
{
    // if it is not null
    if ( mPreconditioner )
    {
        // build preconditioner
        mPreconditioner->build( mLinearSystem, aIter );

        // check if preconditioner exists
        return mPreconditioner->exists();
    }

    return false;
}

//----------------------------------------------------------------------------------------

moris::sint
Linear_Solver_ML::solve_linear_system(
        Linear_Problem*   aLinearSystem,
        const moris::sint aIter )
{
    Tracer tTracer( "LinearSolver", "ML", "Solve" );

    // set linear system
    mLinearSystem = aLinearSystem;

    // Get LHS and RHS vectors
    sol::Dist_Vector* tRHS = mLinearSystem->get_solver_RHS();
    sol::Dist_Vector* tLHS = mLinearSystem->get_free_solver_LHS();

    // Determine the number of RHS and LHS
    uint tNumRHS = tRHS->get_num_vectors();
    uint tNumLHS = tLHS->get_num_vectors();

    MORIS_ERROR( 1 == tNumRHS and 1 == tNumLHS, "Number of LHS does not match number of RHS should be 1" );

    // Get underlying Eptra vectors of RHS and LHS
    Epetra_MultiVector* tRHSepetra = static_cast< Vector_Epetra* >( tRHS )->get_epetra_vector();
    Epetra_MultiVector* tLHSepetra = static_cast< Vector_Epetra* >( tLHS )->get_epetra_vector();

    // initialize and build external preconditioner based on input parameters
    bool tPreconditionerIsBuilt =  this->build_external_preconditioner( aIter );

    if ( !tPreconditionerIsBuilt )
    {
        MORIS_ERROR( false, "Preconditioner is not built" );
    }

    // get the ml preconditioner operator which is a epetra operator
    Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > tMLPrec = mPreconditioner->get_ml_prec();
    
    // get the ml object from the ml preconditioner
    ML* tMLObject = const_cast<ML*> ( tMLPrec->GetML() );

    // iterate over the ML 
    ML_Set_Tolerance(tMLObject, mParameterList.get<real>("Convergence_Tolerance"));
    ML_Set_MaxIterations(tMLObject, mParameterList.get<sint>("Max_Iter"));
    ML_Iterate(tMLObject, (*tLHSepetra)[0], (*tRHSepetra)[0]);

    return 0; 
}

//----------------------------------------------------------------------------------------
moris::sint
Linear_Solver_ML::solve_linear_system()
{
    MORIS_ERROR( false, "Linear_Solver_ML::solve_linear_system not implemented yet ");
    return 0;
}

