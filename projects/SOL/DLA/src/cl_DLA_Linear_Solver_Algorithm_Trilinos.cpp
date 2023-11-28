#include "cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp"

#include "fn_convert_epetra_operator_to_matrix.hpp"
#include "fn_cond.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Matrix.hpp"

using namespace moris::dla;

// cast the preconditioner to the correct type(trillons and assign it to the object)
void Linear_Solver_Algorithm_Trilinos::set_preconditioner( Preconditioner* aPreconditioner )
{
    mPreconditioner = dynamic_cast< Preconditioner_Trilinos* >( aPreconditioner );
}

//-----------------------------------------------------------------------------------

void Linear_Solver_Algorithm_Trilinos::compute_operator_condition_number_with_moris()
{
    // // get Epetra matrix
    Epetra_FECrsMatrix*     tOperator         = mLinearSystem->get_matrix()->get_matrix();
    moris::Matrix< DDRMat > tMorisDenseMatrix = convert_epetra_operator_to_arma_sp_mat< moris::Matrix, DDRMat >( *tOperator );

    // compute the condition number
    real tConditionNumber = cond( tMorisDenseMatrix );
    
    // output the condition number
    MORIS_LOG_INFO( "Condition number of the operator is: %f ", tConditionNumber );
}

//-----------------------------------------------------------------------------------

void Linear_Solver_Algorithm_Trilinos::compute_preconditioned_operator_condition_number_with_moris()
{
    // get Epetra matrix
    Epetra_FECrsMatrix* tOperator = mLinearSystem->get_matrix()->get_matrix();

    // this class is responsible for constructing the preconditioner*operator
    class LocalEpetraOperator
    {
      public:
        Epetra_FECrsMatrix* mOperator;
        Epetra_Operator*    mPreconditioner;
        LocalEpetraOperator( Epetra_FECrsMatrix* aOperator, Epetra_Operator* aPrec )
                : mOperator( aOperator )
                , mPreconditioner( aPrec )
        {
            // Create the Epetra_Operator class for the preconditioner with the
            
        };
        void Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y )
        {
            Epetra_MultiVector temp_Y = Epetra_MultiVector( mOperator->OperatorDomainMap(), mOperator->OperatorDomainMap().NumGlobalElements() );

#ifdef MORIS_HAVE_DEBUG
            // Apply M
            int info = mOperator->Apply( X, temp_Y );
            MORIS_ASSERT( !info, "problem with multiplication in preconditioner" );

            // Apply A or A^{-1}
            info = mPreconditioner->ApplyInverse( temp_Y, Y );
            MORIS_ASSERT( !info, "problem with multiplication in preconditioner" );
#else
            // Apply M
            mOperator->Apply( X, temp_Y );

            // Apply A or A^{-1}
            mPreconditioner->ApplyInverse( temp_Y, Y );
#endif
        };
    };
    
    // get the raw pointer to the preconditioner
    Epetra_Operator* tPrec = mPreconditioner->get_operator().get();
    
    // construct the preconditioned operator
    LocalEpetraOperator tLocalOpeator( tOperator, tPrec );

    // convert the preconditioned operator to a matrix
    Epetra_CrsMatrix tPrecOperator = convert_epetra_operator_to_matrix<LocalEpetraOperator>( &tLocalOpeator, tOperator->OperatorDomainMap() );
    
    // convert the preconditioned operator to a moris matrix
    moris::Matrix< DDRMat > tMorisDenseMatrix = convert_epetra_operator_to_arma_sp_mat< moris::Matrix, DDRMat >( tPrecOperator );

    /// compute the condition number
    real tConditionNumber = cond( tMorisDenseMatrix );

    MORIS_LOG_INFO( "Condition number of the preconditioned operator is: %f ", tConditionNumber );
}
