/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Eigen_Solver.cpp
 *
 */

// Project header files
#include "cl_DLA_Eigen_Solver.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_DLA_Solver_Factory.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "cl_SOL_Enums.hpp"

// C++ system files
#include <cstddef>
#include "moris_typedefs.hpp"
#include "cl_Communication_Tools.hpp"
#include <fstream>
#include <iostream>

// Trilinos includes
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"


#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsMatrix.h"


// Include header for AztecOO iterative linear solver, and
// AztecOO_Operator.  The latter wraps an AztecOO solver in an
// Epetra_Operator.
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// TPL header files
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_OperatorOut.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_FEVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"

// Includes for the generalized Davidson method
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "Teuchos_LAPACK.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"

// Includes for the Block Krylov Schur method
#include "AnasaziBlockKrylovSchurSolMgr.hpp"

#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Epetra_InvOperator.h"

// ML
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"

#include "Amesos.h"

using namespace moris;
using namespace dla;

// ----------------------------------------------------------------------------

Eigen_Solver::Eigen_Solver( const Parameter_List* aParameterList )
        : Linear_Solver_Algorithm_Trilinos( *aParameterList )
        , mMat( NULL )
        , mMassMat( NULL )
        , mSolTime( 0.0 )
        , mFreeSolVec( NULL )
        , mNumReturnedEigVals( 0 )
{
    if ( mParameterList.get< bool >( "Verbosity" ) )
    {
        mPrinter.stream( Anasazi::Errors ) << std::endl
                                           << Anasazi::Anasazi_Version() << std::endl;
    }

    if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_BLOCK_DAVIDSON" )
    {
        MORIS_LOG_INFO( "Eigen Algorithm is Block Davidson" );
    }
    else if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_GENERALIZED_DAVIDSON" )
    {
        MORIS_LOG_INFO( "Eigen Algorithm is Generalized Davidson" );
    }
    else if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_BLOCK_KRYLOV_SCHUR" )
    {
        MORIS_LOG_INFO( "Eigen Algorithm is Block Krylov Schur" );
    }
    else if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS" )
    {
        MORIS_LOG_INFO( "Eigen Algorithm is Block Krylov Schur Amesos " );
    }
    else
    {
        MORIS_ERROR( false, "Wrong Eigensolver algorithm specified!\n" );
    }
}

// ----------------------------------------------------------------------------

Eigen_Solver::~Eigen_Solver()
{
    if ( mFreeSolVec )
    {
        delete mFreeSolVec;
        mFreeSolVec = nullptr;
    }
}

// ----------------------------------------------------------------------------

moris::sint
Eigen_Solver::solve_linear_system(
        Linear_Problem*   aLinearSystem,
        const moris::sint aIter )
{
    // get stiffness matrix of class Epetra_FECrs
    mMat = aLinearSystem->get_matrix()->get_matrix();

    // assemble the right hand side of the matrix as it is needed for generalized eigenvalue problem
    aLinearSystem->assemble_rhs_matrix();

    // get mass matrix of class Epetra_FECrs
    mMassMat = aLinearSystem->get_mass_matrix()->get_matrix();

    // get stiffness matrix of distributed matrix type
    mNewMat = aLinearSystem->get_matrix();

    // request type of eigen solver algorithm from parameterlist
    std::string tEigAlgType =
            mParameterList.get< std::string >( "Eigen_Algorithm" );

    if ( tEigAlgType == "EIGALG_BLOCK_DAVIDSON" )
    {
        // set block davidson method as eigen solver algorithm
        return this->solve_block_davidson_system( aLinearSystem );
    }
    else if ( tEigAlgType == "EIGALG_GENERALIZED_DAVIDSON" )
    {
        // set generalized davidson method as eigen solver algorithm
        return this->solve_generalized_davidson_system( aLinearSystem );
    }
    else if ( tEigAlgType == "EIGALG_BLOCK_KRYLOV_SCHUR" )
    {
        // set block krylov schur method as eigen solver algorithm
        return this->solve_block_krylov_schur_system( aLinearSystem );
    }
    else if ( tEigAlgType == "EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS" )
    {
        // set block krylov schur (amesos) method as eigen solver algorithm
        return this->solve_block_krylov_schur_amesos_system( aLinearSystem );
    }
    else
    {
        MORIS_ERROR( false, "Wrong Eigensolver algorithm specified!\n" );
        return 0;
    }
}

// ----------------------------------------------------------------------------

void Eigen_Solver::build_linearized_system( Linear_Problem* aLinearSystem )
{
    mSPmat = Teuchos::RCP( mMat, false );

    // sparse-mass matrix of Teuchos::RCP type
    mSPmassmat = Teuchos::RCP( mMassMat, false );

    // request number of DOFs from parameterlist
    moris::sint tNumFreeDofs = mParameterList.get< moris::sint >( "NumFreeDofs" );

    // request block-size from parameterlist
    moris::sint tBlockSize = mParameterList.get< moris::sint >( "Block_Size" );

    // request number of eigenvalues from parameterlist
    moris::sint tNumEigVals = mParameterList.get< moris::sint >( "Num_Eig_Vals" );

    // stiffness and mass matrix needs to be symmetric for block davidson and generalized davidson
    bool tAssumeSymmetric = true;

    // get map of dist_map class
    mMap = mNewMat->get_map();

    // create eigen solver vector of Vector_Epetra class
    mFreeSolVec = new Vector_Epetra( mMap, 1, false, false );

    // Make sure that the number of blocks and eigenvalues does not exceed NumFreeDofs
    if ( tBlockSize > tNumFreeDofs || tNumEigVals > tNumFreeDofs )
    {
        MORIS_ERROR( false, "EigenSolver::BuildLinearSystem: Number of blocks and/or number of eigenvalues can not exceed NumFreeDofs\n" );
    }

    // Create initial vector for the solver
    mIvec = Teuchos::RCP( new Epetra_MultiVector( mSPmat->OperatorDomainMap(), tBlockSize ) );
    mIvec->Random();

    // if mLeftPreconditioner exists
    if ( mLeftPreconditioner )
    {
        mLeftPreconditioner->build( aLinearSystem, 0 );

        // get the ml precondinioer from the tPrec Object
        Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > tPrec = mLeftPreconditioner->get_ml_prec();

        // Teuchos::RCP<Epetra_Operator> tPrec = mLeftPreconditioner->get_operator();

        // Create an Operator that computes y = M^{-1} K x.
        mSPmat = Teuchos::rcp( new Anasazi::EpetraGenOp( tPrec, mSPmat ) );
    }

    // Create the eigenproblem, except if the algorithm used is BLOCK_KRYLOV_SCHUR or EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS
    // as they build their eigenproblems later locally
    if ( ( mParameterList.get< std::string >( "Eigen_Algorithm" ) != "EIGALG_BLOCK_KRYLOV_SCHUR" ) || ( mParameterList.get< std::string >( "Eigen_Algorithm" ) != "EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS" ) )
    {
        // Create the eigen problem object
        mMyEigProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem< double, MV, OP >( mSPmat, mSPmassmat, mIvec ) );

        mMyEigProblem->setHermitian( tAssumeSymmetric );

        // Set the number of eigenvalues to be computed
        mMyEigProblem->setNEV( tNumEigVals );
    }
}

// -----------------------------------------------------------------------------

void Eigen_Solver::set_eigen_solver_manager_parameters()
{
    // create eigen solver parameterlist
    mParameterList = prm::create_eigen_algorithm_parameter_list();
}

// -----------------------------------------------------------------------------

int Eigen_Solver::solve_block_davidson_system( Linear_Problem* aLinearSystem )
{
    int MyPID = par_rank();

    // build linear system
    this->build_linearized_system( aLinearSystem );
    Teuchos::RCP< Epetra_Operator > tPrecOp;
    Teuchos::RCP< Epetra_Operator > tIfpackPrec;

    // TODO: avoid rebuilding the preconditioner and use it if it exists
    mPreconditioner->build( aLinearSystem, 0 );

    // create an
    tIfpackPrec = mPreconditioner->get_operator();

    tPrecOp = Teuchos::RCP( new Epetra_InvOperator( tIfpackPrec.get() ) );

    // set preconditions to eigenproblem
    mMyEigProblem->setPrec( tPrecOp );

    // check if eigen problem is set
    bool tboolret = mMyEigProblem->setProblem();

    // Check if an error was returned
    if ( tboolret != true && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockDavidsonSystem: Error was returned when setting up the problem!\n" );
        return -1;
    }

    // Set verbosity level
    bool tverbose   = mParameterList.get< bool >( "Verbosity" );
    int  tverbosity = Anasazi::Errors + Anasazi::Warnings;
    if ( tverbose )
    {
        tverbosity += Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
    }

    // Set the parameters and pass to solver manager
    Teuchos::ParameterList MyPL;
    MyPL.set( "Verbosity", tverbosity );
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );
    MyPL.set( "Num Blocks", mParameterList.get< moris::sint >( "Num_Blocks" ) );
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );
    MyPL.set( "Convergence Tolerance", mParameterList.get< moris::real >( "Convergence_Tolerance" ) );
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

    // Create the solver manager
    Anasazi::BlockDavidsonSolMgr< double, MV, OP > MySolverMan( mMyEigProblem, MyPL );

    // Solve the problem
    Anasazi::ReturnType tReturnCode = MySolverMan.solve();

    // Check if the problem solve converged
    if ( tReturnCode != Anasazi::Converged and MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockDavidsonSystem: Solver returned UNconverged.\n" );
        return -1;
    }

    // Get the eigenvalues and eigenvectors from the eigen problem
    mSol = mMyEigProblem->getSolution();

    // store eigenvalues and eigenvectors
    std::vector< Anasazi::Value< double > > evals = mSol.Evals;
    Teuchos::RCP< MV >                      evecs = mSol.Evecs;

    // store eigenvalues into member variable
    mevals = evals;

    // Set eigenvector to solver interface
    bool tUpdateFlag = mParameterList.get< bool >( "Update_Flag" );
    if ( tUpdateFlag )
    {
        sol::Dist_Vector* tDistVec   = aLinearSystem->get_solver_input()->get_eigen_solution_vector();
        MV*               tEpetraVec = ( dynamic_cast< Vector_Epetra* >( tDistVec ) )->get_epetra_vector();

        //  Check for size of both vectors
        // TODO: this does not work in parallel
        if ( evecs->NumVectors() == tEpetraVec->NumVectors() )
        {
            // updates mEigenSolVector with evecs and returns Multivector
            tEpetraVec->Update( 1.0, *evecs, 0.0 );
        }
        else
        {
            MORIS_ERROR( false, "Number of eigen vector and eigen value in parameterlist should be the same." );
        }
    }

    // Compute residuals.
    std::vector< double > normR( mSol.numVecs );
    if ( mSol.numVecs > 0 )
    {
        Teuchos::SerialDenseMatrix< int, double > T( mSol.numVecs, mSol.numVecs );

        Epetra_MultiVector Kevec( mSPmat->OperatorDomainMap(), evecs->NumVectors() );
        Epetra_MultiVector Mevec( mSPmassmat->OperatorDomainMap(), evecs->NumVectors() );

        // store real eigenvalues into serial dense matrix
        T.putScalar( 0.0 );

        for ( int i = 0; i < mSol.numVecs; i++ )
        {
            T( i, i ) = evals[ i ].realpart;
        }

        // Operator application on stiffness and mass matrix
        mSPmat->Apply( *evecs, Kevec );
        mSPmassmat->Apply( *evecs, Mevec );

        MVT::MvTimesMatAddMv( -1.0, Mevec, T, 1.0, Kevec );
        MVT::MvNorm( Kevec, normR );
    }
    mNumReturnedEigVals = mSol.numVecs;

    // Initialize variables needed for the computation of Eigvec * MassMat * Eigvec
    int                      i = 0;
    std::vector< int >       curind( 1 );
    std::vector< double >    massnormsqr( 1 );
    Teuchos::RCP< MV >       tempKevec, Mevecs;
    Teuchos::RCP< const MV > tempeveci, tempMevec, tempEigVec;
    mNormEvecsMassMatEvecs.assign( mNumReturnedEigVals, 0.0 );

    // Compute M * evecs = Mevecs
    Mevecs = Teuchos::rcp( new Epetra_MultiVector( mSPmassmat->OperatorDomainMap(), mNumReturnedEigVals ) );
    OPT::Apply( *mSPmassmat, *evecs, *Mevecs );

    // Loop over all modes and scale eigenvectors
    while ( i < mNumReturnedEigVals )
    {
        // Get a view of the M*evecr
        curind[ 0 ] = i;
        tempMevec   = MVT::CloneView( *Mevecs, curind );

        // Get current eigenvector
        tempEigVec = MVT::CloneView( *evecs, curind );

        // Compute the norm of Eigvec * MassMat * Eigvec
        MVT::MvDot( *tempEigVec, *tempMevec, massnormsqr );
        mNormEvecsMassMatEvecs[ i ] = std::sqrt( massnormsqr[ 0 ] );

        // Increment counter
        i++;
    }

    // Output computed eigenvalues and their direct residuals
    {
        MORIS_LOG_INFO( "------------------------------------------------" );
        for ( int i = 0; i < mSol.numVecs; i++ )
        {
            MORIS_LOG_INFO( "EigenValue: %16f ", evals[ i ].realpart );
            MORIS_LOG_INFO( "Direct Residual: %18e ", normR[ i ] / evals[ i ].realpart );
        }
    }

    // print eigenVector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set leader solution vector to free solution vector
        Vector_Epetra* aLeaderSolVec = mFreeSolVec;

        // get solution for eigen vectors
        this->get_solution( m, mFreeSolVec, aLeaderSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    return 0;
}

// ----------------------------------------------------------------------------

int Eigen_Solver::solve_generalized_davidson_system( Linear_Problem* aLinearSystem )
{
    using std::cout;
    using std::endl;

    int MyPID = par_rank();

    // build linear system
    this->build_linearized_system( aLinearSystem );

    // Set verbosity level
    bool tverbose   = mParameterList.get< bool >( "Verbosity" );
    int  tverbosity = Anasazi::Errors + Anasazi::Warnings;
    if ( tverbose )
    {
        tverbosity += Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
    }

    // Create parameter list to pass into solver
    Teuchos::ParameterList MyPL;
    MyPL.set( "Verbosity", tverbosity );
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );
    MyPL.set( "Convergence Tolerance", mParameterList.get< moris::real >( "Convergence_Tolerance" ) );
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

    mPreconditioner->build( aLinearSystem, 0 );

    Teuchos::RCP< Epetra_Operator > tPrecOp;
    Teuchos::RCP< Epetra_Operator > tIfpackPrec;

    tIfpackPrec = mPreconditioner->get_operator();

    tPrecOp = Teuchos::RCP( new Epetra_InvOperator( tIfpackPrec.get() ) );

    // set preconditions to eigenproblem
    mMyEigProblem->setPrec( tPrecOp );

    // Inform the eigenproblem that you are finished passing it information
    bool tboolret = mMyEigProblem->setProblem();

    // Check if an error was returned
    if ( tboolret != true && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveGeneralizedDavidsonSystem: Error was returned when setting up the problem!\n" );
        return -1;
    }

    // loop over increasing convergence tolerances until converged solution found
    real tTolConv = mParameterList.get< moris::real >( "Convergence_Tolerance" );

    for ( uint iconv = 0; iconv < 5; ++iconv )
    {
        // Set required relative convergence
        MyPL.set( "Convergence Tolerance", tTolConv );

        if ( MyPID == 0 )
        {
            fprintf( stdout, "\n ... Eigensolver Step %d (%d) - set tolerance for eigensolver to %e\n\n", iconv + 1, 5, tTolConv );
        }

        // Initialize the Block Arnoldi solver
        Anasazi::GeneralizedDavidsonSolMgr< double, MV, OP > MySolverMgr( mMyEigProblem, MyPL );

        // Solve the problem to the specified tolerances or length
        Anasazi::ReturnType tReturnCode = MySolverMgr.solve();

        // Get number of converged eigen vectors
        mSol = mMyEigProblem->getSolution();

        // Check if the problem solve converged
        if ( tReturnCode != Anasazi::Converged || mSol.numVecs < mParameterList.get< sint >( "Num_Eig_Vals" ) )
        {
            if ( MyPID == 0 )
            {
                MORIS_ERROR( false, "EigenSolver::SolveGeneralizedDavidsonSystem: Solver returned UNconverged. Increase tolerance." );
            }
        }
        else
        {
            break;
        }

        // Increase tolerance
        tTolConv *= 10.0;
    }

    // Get the eigenvalues and eigenvectors from the eigenproblem
    mSol                                          = mMyEigProblem->getSolution();
    std::vector< Anasazi::Value< double > > evals = mSol.Evals;
    Teuchos::RCP< MV >                      evecs = mSol.Evecs;
    std::vector< int >                      index = mSol.index;
    mNumReturnedEigVals                           = mSol.numVecs;

    // store eigenvalues into member variable
    mevals = evals;

    //  Set eigenvector solution to solver interface
    bool tUpdateFlag = mParameterList.get< bool >( "Update_Flag" );
    if ( tUpdateFlag )
    {
        sol::Dist_Vector* tDistVec   = aLinearSystem->get_solver_input()->get_eigen_solution_vector();
        MV*               tEpetraVec = ( dynamic_cast< Vector_Epetra* >( tDistVec ) )->get_epetra_vector();

        //  Check for size of both vectors
        // TODO: this does not work in parallel
        if ( evecs->NumVectors() == tEpetraVec->NumVectors() )
        {
            // update mEigenSolVector with evecs and return Multivector
            tEpetraVec->Update( 1.0, *evecs, 0.0 );
        }
        else
        {
            MORIS_ERROR( false, "Number of eigen vector and eigen value in parameterlist should be the same." );
        }
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, "\n ... received %d eigenvectors from solver\n\n", mNumReturnedEigVals );
    }

    if ( mNumReturnedEigVals > 0 )
    {
        // Compute residuals.
        Teuchos::LAPACK< int, double > lapack;
        std::vector< double >          normR( mNumReturnedEigVals );
        mNormEvecsMassMatEvecs.assign( mNumReturnedEigVals, 0.0 );

        // The problem is non-Hermitian.
        int                      i = 0;
        std::vector< int >       curind( 1 );
        std::vector< double >    resnorm( 1 ), massnormsqr( 1 ), tempnrm( 1 );
        Teuchos::RCP< MV >       tempKevec, Mevecs;
        Teuchos::RCP< const MV > tempeveci, tempMevec, tempEigVec;

        // Compute K * evecs = Kevecs
        Epetra_MultiVector Kevecs( mSPmat->OperatorDomainMap(), mNumReturnedEigVals );
        OPT::Apply( *mSPmat, *evecs, Kevecs );

        // Compute M * evecs = Mevecs
        Mevecs = Teuchos::rcp( new Epetra_MultiVector( mSPmassmat->OperatorDomainMap(), mNumReturnedEigVals ) );
        OPT::Apply( *mSPmassmat, *evecs, *Mevecs );

        Teuchos::SerialDenseMatrix< int, double > Breal( 1, 1 ), Bimag( 1, 1 );

        while ( i < mNumReturnedEigVals )
        {
            if ( index[ i ] == 0 )
            {
                // Get a view of the M*evecr
                curind[ 0 ] = i;
                tempMevec   = MVT::CloneView( *Mevecs, curind );

                // Get a copy of A*evecr
                tempKevec = MVT::CloneCopy( Kevecs, curind );

                // Compute K*evecr - lambda*M*evecr
                Breal( 0, 0 ) = evals[ i ].realpart;
                MVT::MvTimesMatAddMv( -1.0, *tempMevec, Breal, 1.0, *tempKevec );

                // Compute the norm of the residual
                MVT::MvNorm( *tempKevec, resnorm );
                normR[ i ] = resnorm[ 0 ];

                // Get current eigenvector
                tempEigVec = MVT::CloneView( *evecs, curind );

                // Compute the norm of Eigvec * MassMat * Eigvec
                MVT::MvDot( *tempEigVec, *tempMevec, massnormsqr );
                mNormEvecsMassMatEvecs[ i ] = std::sqrt( massnormsqr[ 0 ] );

                // Increment counter
                i++;
            }
            else
            {
                // Get a view of the real part of M*evecr
                curind[ 0 ] = i;
                tempMevec   = MVT::CloneView( *Mevecs, curind );

                // Get a copy of K*evecr
                tempKevec = MVT::CloneCopy( Kevecs, curind );

                // Get a view of the imaginary part of the eigenvector (eveci)
                curind[ 0 ] = i + 1;
                tempeveci   = MVT::CloneView( *Mevecs, curind );

                // Set the eigenvalue into Breal and Bimag
                Breal( 0, 0 ) = evals[ i ].realpart;
                Bimag( 0, 0 ) = evals[ i ].imagpart;

                // Compute K*evecr - M*evecr*lambdar + M*eveci*lambdai
                MVT::MvTimesMatAddMv( -1.0, *tempMevec, Breal, 1.0, *tempKevec );
                MVT::MvTimesMatAddMv( 1.0, *tempeveci, Bimag, 1.0, *tempKevec );
                MVT::MvNorm( *tempKevec, tempnrm );

                // Get a copy of K*eveci
                tempKevec = MVT::CloneCopy( Kevecs, curind );

                // Compute K*eveci - M*eveci*lambdar - M*evecr*lambdai
                MVT::MvTimesMatAddMv( -1.0, *tempMevec, Bimag, 1.0, *tempKevec );
                MVT::MvTimesMatAddMv( -1.0, *tempeveci, Breal, 1.0, *tempKevec );
                MVT::MvNorm( *tempKevec, resnorm );

                // Compute the norms and scale by magnitude of eigenvalue
                normR[ i ]     = lapack.LAPY2( tempnrm[ 0 ], resnorm[ 0 ] );
                normR[ i + 1 ] = normR[ i ];

                // Get current eigenvector
                tempEigVec = MVT::CloneView( *evecs, curind );

                // Compute the norm of Eigvec * MassMat * Eigvec
                MVT::MvDot( *tempEigVec, *tempMevec, massnormsqr );
                mNormEvecsMassMatEvecs[ i ]     = std::sqrt( massnormsqr[ 0 ] );
                mNormEvecsMassMatEvecs[ i + 1 ] = mNormEvecsMassMatEvecs[ i ];

                // Increment counter
                i = i + 2;
            }
        }

        // Output computed eigenvalues and their direct residuals
        MORIS_LOG_INFO( "===============================================" );
        for ( int j = 0; j < mNumReturnedEigVals; j++ )
        {
            MORIS_LOG_INFO( "Real Part: %16f ", evals[ j ].realpart );
            MORIS_LOG_INFO( "Imaginary Part: %16f", evals[ j ].imagpart );
            MORIS_LOG_INFO( "Direct Residual: %18e ", normR[ j ] );
        }
    }

    // print eigenvector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set leader solution vector to free solution vector
        Vector_Epetra* aLeaderSolVec = mFreeSolVec;

        // get solution for eigen vectors
        this->get_solution( m, mFreeSolVec, aLeaderSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, "\n ... Finished EigenValue Solve using a Generalized Davidson Solver\n\n" );
    }

    return 0;
}

// -----------------------------------------------------------------------------

int Eigen_Solver::solve_block_krylov_schur_system( Linear_Problem* aLinearSystem )
{

    this->build_linearized_system( aLinearSystem );

    // NOTE: This implementation follows trilinos/packages/anasazi/epetra/example/BlockKrylovSchur/BlockKrylovSchurEpetraExGenAztecOO.cpp

    int MyPID = par_rank();

    // build preconditioner
    mPreconditioner->build( aLinearSystem, 0 );
    Teuchos::RCP< Epetra_Operator > tPrecOp = mPreconditioner->get_operator();

    // Tell the linear problem about the mass matrix M
    mEpetraProblem.SetOperator( mSPmat.get() );

    // Create AztecOO iterative solver for solving linear systems with K.
    AztecOO precSolver( mEpetraProblem );
    precSolver.SetPrecOperator( tPrecOp.get() );
    precSolver.SetAztecOption( AZ_output, AZ_none );     // Don't print output
    precSolver.SetAztecOption( AZ_solver, AZ_gmres );    // Use GMRES

    // Use the above AztecOO solver to create the AztecOO_Operator.
    Teuchos::RCP< AztecOO_Operator > precOperator = Teuchos::rcp( new AztecOO_Operator( &precSolver, mMat->NumGlobalRows(), 1e-12 ) );

    // Create an Operator that computes y = M^{-1} K x.
    Teuchos::RCP< Anasazi::EpetraGenOp > Aop = Teuchos::rcp( new Anasazi::EpetraGenOp( precOperator, mSPmassmat ) );

    // Create the eigen problem object
    mMyEigProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem< double, MV, OP >( Aop, mSPmassmat, mIvec ) );

    // request number of eigenvalues from parameterlist
    moris::sint tNumEigVals = mParameterList.get< moris::sint >( "Num_Eig_Vals" );

    // set whether linear system is symmetric or not
    bool tAssumeSymmetric = false;

    // Set whether AMat and BMat are Hermitian = symmetric
    mMyEigProblem->setHermitian( tAssumeSymmetric );

    // Set the number of eigenvalues to be computed
    mMyEigProblem->setNEV( tNumEigVals );

    // check if eigenproblem is set or not
    bool tboolret = mMyEigProblem->setProblem();

    // Check if an error was returned
    if ( tboolret != true && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurSystem: Error was returned when setting up the problem!\n" );
        return -1;
    }

    // Set verbosity level
    bool tverbose   = mParameterList.get< bool >( "Verbosity" );
    int  tverbosity = Anasazi::Errors + Anasazi::Warnings;
    if ( tverbose )
    {
        tverbosity += Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
    }

    // Set the parameters and pass to solver manager
    Teuchos::ParameterList MyPL;
    MyPL.set( "Verbosity", tverbosity );
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );
    MyPL.set( "Num Blocks", mParameterList.get< moris::sint >( "Num_Blocks" ) );
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );
    MyPL.set( "Convergence Tolerance", mParameterList.get< real >( "Convergence_Tolerance" ) );
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

    // Create the solver manager
    Anasazi::BlockKrylovSchurSolMgr< double, MV, OP > MySolverMan( mMyEigProblem, MyPL );

    // Solve the problem
    Anasazi::ReturnType tReturnCode = MySolverMan.solve();

    // Check if the problem solve converged
    if ( tReturnCode != Anasazi::Converged && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurSystem: Solver returned UNconverged.\n" );
        return -1;
    }

    // Get the eigenvalues and eigenvectors from the eigen problem
    mSol                                          = mMyEigProblem->getSolution();
    std::vector< Anasazi::Value< double > > evals = mSol.Evals;
    Teuchos::RCP< MV >                      evecs = mSol.Evecs;

    // store eigenvalues into member variable
    mevals = evals;

    // Set eigenvector solution to solver interface
    bool tUpdateFlag = mParameterList.get< bool >( "Update_Flag" );
    if ( tUpdateFlag )
    {
        sol::Dist_Vector* tDistVec   = aLinearSystem->get_solver_input()->get_eigen_solution_vector();
        MV*               tEpetraVec = ( dynamic_cast< Vector_Epetra* >( tDistVec ) )->get_epetra_vector();

        //  Check for size of both vectors
        // TODO: this does not work in parallel
        if ( evecs->NumVectors() == tEpetraVec->NumVectors() )
        {
            // update mEigenSolVector with evecs and return Multivector
            tEpetraVec->Update( 1.0, *evecs, 0.0 );
        }
        else
        {
            MORIS_ERROR( false, "Number of eigen vector and eigen value in parameterlist should be the same." );
        }
    }

    // Compute residuals.
    std::vector< double > normR( mSol.numVecs );
    if ( mSol.numVecs > 0 )
    {
        Teuchos::SerialDenseMatrix< int, double > T( mSol.numVecs, mSol.numVecs );
        Epetra_MultiVector                        Kevec( mSPmat->OperatorDomainMap(), evecs->NumVectors() );
        Epetra_MultiVector                        Mevec( mSPmassmat->OperatorDomainMap(), evecs->NumVectors() );
        T.putScalar( 0.0 );

        // Create diagonal matrix containing all eigenvalues
        for ( int i = 0; i < mSol.numVecs; i++ )
        {
            T( i, i ) = evals[ i ].realpart;
        }

        // Compute Kevec = Kmat * evec
        mSPmat->Apply( *evecs, Kevec );

        // Compute Mevec = Mmat * evec
        mSPmassmat->Apply( *evecs, Mevec );

        // Compute Kevec = Kevec - T * Mevec
        MVT::MvTimesMatAddMv( -1.0, Mevec, T, 1.0, Kevec );

        // Get the residual norm of Kevec
        MVT::MvNorm( Kevec, normR );
    }
    mNumReturnedEigVals = mSol.numVecs;

    // Initialize variables needed for the computation of Eigvec * MassMat * Eigvec
    int                      i = 0;
    std::vector< int >       curind( 1 );
    std::vector< double >    massnormsqr( 1 );
    Teuchos::RCP< MV >       tempKevec, Mevecs;
    Teuchos::RCP< const MV > tempeveci, tempMevec, tempEigVec;
    mNormEvecsMassMatEvecs.assign( mNumReturnedEigVals, 0.0 );

    // Compute M * evecs = Mevecs
    Mevecs = Teuchos::rcp( new Epetra_MultiVector( mSPmassmat->OperatorDomainMap(), mNumReturnedEigVals ) );
    OPT::Apply( *mSPmassmat, *evecs, *Mevecs );

    // Loop over all modes and scale eigenvectors
    while ( i < mNumReturnedEigVals )
    {
        // Get a view of the M*evecr
        curind[ 0 ] = i;
        tempMevec   = MVT::CloneView( *Mevecs, curind );

        // Get current eigenvector
        tempEigVec = MVT::CloneView( *evecs, curind );

        // Compute the norm of Eigvec * MassMat * Eigvec
        MVT::MvDot( *tempEigVec, *tempMevec, massnormsqr );
        mNormEvecsMassMatEvecs[ i ] = std::sqrt( massnormsqr[ 0 ] );

        // Increment counter
        i++;
    }

    // Output computed eigenvalues and their direct residuals

    MORIS_LOG_INFO( "===============================================" );
    for ( int i = 0; i < mSol.numVecs; i++ )
    {
        MORIS_LOG_INFO( "EigenValue: %16f", 1 / evals[ i ].realpart );
        MORIS_LOG_INFO( "Direct Residual: %18e", normR[ i ] / evals[ i ].realpart );
    }


    // print eigen vector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set leader solution vector to free solution vector
        Vector_Epetra* aLeaderSolVec = mFreeSolVec;

        // get solution for eigen vector
        this->get_solution( m, mFreeSolVec, aLeaderSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, "\n ... Finished EigenValue Solve using a Block Krylov Schur Iterative (AZTEC) Solver\n\n" );
    }

    return 0;
}

// -----------------------------------------------------------------------------

int Eigen_Solver::solve_block_krylov_schur_amesos_system( Linear_Problem* aLinearSystem )
{
    this->build_linearized_system( aLinearSystem );

    // NOTE: This implementation follows trilinos/packages/anasazi/epetra/example/BlockKrylovSchur/BlockKrylovSchurEpetraExGenAmesos.cpp
    int MyPID = par_rank();

    // Tell the linear problem about the mass matrix M
    mEpetraProblem.SetOperator( mSPmat.get() );

    // Create Amesos factory and solver for solving linear systems with K
    Amesos tAmesosFactory;

    // Note that the AmesosProblem object "absorbs" M.  Anasazi doesn't see M, M is needed to orthogonalize the eigenvectors.
    Teuchos::RCP< Amesos_BaseSolver > mAmesosSolver;

    std::string tAmesosDirectSolverType = mParameterList.get< std::string >( "Amesos_Direct_Solver_Type" );

    // Determine which direct solver method to use first by ensuring a correct solver name
    if ( tAmesosFactory.Query( tAmesosDirectSolverType ) )
    {
        mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( tAmesosDirectSolverType, mEpetraProblem ) );

        if ( mAmesosSolver.is_null() )
        {
            MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurAmesosSystem: Error was returned when creating the Amesos solver!\n" );
        }
    }
    else
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurAmesosSystem: The specified mEigSolMethod is not yet supported, check spelling!\n" );
    }

    // The AmesosGenOp class assumes that the symbolic and numeric factorizations have already been performed on the linear problem
    mAmesosSolver->SymbolicFactorization();
    mAmesosSolver->NumericFactorization();

    // Set verbosity level
    bool tverbose   = mParameterList.get< bool >( "Verbosity" );
    int  tverbosity = Anasazi::Errors + Anasazi::Warnings;
    if ( tverbose )
    {
        tverbosity += Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
    }

    // Set the parameters and pass to solver manager
    Teuchos::ParameterList MyPL;
    MyPL.set( "Verbosity", tverbosity );
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );
    MyPL.set( "Num Blocks", mParameterList.get< moris::sint >( "Num_Blocks" ) );
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );
    MyPL.set( "Convergence Tolerance", mParameterList.get< real >( "Convergence_Tolerance" ) );
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

    // Create an initial set of vectors to start the eigensolver.  Note:
    // This needs to have the same number of columns as the block size.
    Teuchos::RCP< MV > mIvec = Teuchos::rcp( new MV( mSPmassmat->OperatorDomainMap(), mParameterList.get< moris::sint >( "Num_Blocks" ) ) );
    mIvec->Random();

    // Create the Epetra_Operator for the spectral transformation using the Amesos direct solver
    Teuchos::RCP< Amesos_GenOp > Aop = Teuchos::rcp( new Amesos_GenOp( mAmesosSolver, mSPmassmat ) );

    // Create the eigen problem object
    mMyEigProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem< double, MV, OP >( Aop, mSPmassmat, mIvec ) );

    bool tAssumeSymmetric = true;
    // Set whether AMat and BMat are Hermitian = symmetric
    mMyEigProblem->setHermitian( tAssumeSymmetric );

    // Set the number of eigenvalues to be computed
    mMyEigProblem->setNEV( mParameterList.get< moris::sint >( "Num_Eig_Vals" ) );

    // check if the eigenproblem is set or not
    bool tboolret = mMyEigProblem->setProblem();

    // Check if an error was returned
    if ( tboolret != true && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurAmesosSystem: Error was returned when setting up the problem!\n" );
        return -1;
    }

    // Create the solver manager
    Anasazi::BlockKrylovSchurSolMgr< double, MV, OP > MySolverMan( mMyEigProblem, MyPL );

    // Solve the problem
    Anasazi::ReturnType tReturnCode = MySolverMan.solve();

    // Check if the problem solve converged
    if ( tReturnCode != Anasazi::Converged && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurAmesosSystem: Solver returned UNconverged.\n" );
        return -1;
    }

    // Get the eigenvalues and eigenvectors from the eigen problem
    mSol                                          = mMyEigProblem->getSolution();
    std::vector< Anasazi::Value< double > > evals = mSol.Evals;
    Teuchos::RCP< MV >                      evecs = mSol.Evecs;

    // store eigenvalues into member variable
    mevals = evals;

    // Set eigenvector solution to solver interface
    bool tUpdateFlag = mParameterList.get< bool >( "Update_Flag" );
    if ( tUpdateFlag )
    {
        sol::Dist_Vector* tDistVec   = aLinearSystem->get_solver_input()->get_eigen_solution_vector();
        MV*               tEpetraVec = ( dynamic_cast< Vector_Epetra* >( tDistVec ) )->get_epetra_vector();

        //  Check for size of both vectors
        // TODO: this does not work in parallel
        if ( evecs->NumVectors() == tEpetraVec->NumVectors() )
        {
            // update mEigenSolVector with evecs and return Multivector
            tEpetraVec->Update( 1.0, *evecs, 0.0 );
        }
        else
        {
            MORIS_ERROR( false, "Number of eigen vector and eigen value in parameterlist should be the same." );
        }
    }


    // Compute residuals.
    std::vector< double > normR( mSol.numVecs );
    if ( mSol.numVecs > 0 )
    {
        Teuchos::SerialDenseMatrix< int, double > T( mSol.numVecs, mSol.numVecs );
        Epetra_MultiVector                        Kevec( mSPmassmat->OperatorDomainMap(), evecs->NumVectors() );
        Epetra_MultiVector                        Mevec( mSPmassmat->OperatorDomainMap(), evecs->NumVectors() );
        T.putScalar( 0.0 );

        // Create diagonal matrix containing all eigenvalues
        for ( int i = 0; i < mSol.numVecs; i++ )
        {
            T( i, i ) = evals[ i ].realpart;
        }

        // Compute Kevec = Kmat * evec
        mSPmat->Apply( *evecs, Kevec );

        // Compute Mevec = Mmat * evec
        mSPmassmat->Apply( *evecs, Mevec );

        // Compute Kevec = Kevec - T * Mevec
        MVT::MvTimesMatAddMv( -1.0, Mevec, T, 1.0, Kevec );

        // Get the residual norm of Kevec
        MVT::MvNorm( Kevec, normR );
    }
    mNumReturnedEigVals = mSol.numVecs;

    // Initialize variables needed for the computation of Eigvec * MassMat * Eigvec
    int                      i = 0;
    std::vector< int >       curind( 1 );
    std::vector< double >    massnormsqr( 1 );
    Teuchos::RCP< MV >       tempKevec, Mevecs;
    Teuchos::RCP< const MV > tempeveci, tempMevec, tempEigVec;
    mNormEvecsMassMatEvecs.assign( mNumReturnedEigVals, 0.0 );

    // Compute M * evecs = Mevecs
    Mevecs = Teuchos::rcp( new Epetra_MultiVector( mSPmassmat->OperatorDomainMap(), mNumReturnedEigVals ) );
    OPT::Apply( *mSPmassmat, *evecs, *Mevecs );

    // Loop over all modes and scale eigenvectors
    while ( i < mNumReturnedEigVals )
    {
        // Get a view of the M*evecr
        curind[ 0 ] = i;
        tempMevec   = MVT::CloneView( *Mevecs, curind );

        // Get current eigenvector
        tempEigVec = MVT::CloneView( *evecs, curind );

        // Compute the norm of Eigvec * MassMat * Eigvec
        MVT::MvDot( *tempEigVec, *tempMevec, massnormsqr );
        mNormEvecsMassMatEvecs[ i ] = std::sqrt( massnormsqr[ 0 ] );

        // Increment counter
        i++;
    }

    // Output computed eigenvalues and their direct residuals
    MORIS_LOG_INFO( "------------------------------------------------" );
    for ( int i = 0; i < mSol.numVecs; i++ )
    {
        MORIS_LOG_INFO( "EigenValue: %16f", 1.0 / evals[ i ].realpart );
        MORIS_LOG_INFO( "Direct Residual: %18e", normR[ i ] * evals[ i ].realpart );
    }

    // print eigen vector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set leader solution vector to free solution vector
        Vector_Epetra* aLeaderSolVec = mFreeSolVec;

        // get solution for eigen vector
        this->get_solution( m, mFreeSolVec, aLeaderSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, " ... Finished EigenValue Solve using a Block Krylov Schur Direct (AMESOS) Solver\n\n" );
    }

    return 0;
}

// ----------------------------------------------------------------------------

int Eigen_Solver::get_solution(
        uint           aEigValIndex,
        Vector_Epetra* aSolVec,
        Vector_Epetra* aLeaderSolVec,
        real&          aEigValReal,
        real&          aEigValImag )
{

    // Get the eigenvalues and eigenvectors from the eigenproblem
    std::vector< Anasazi::Value< double > > evals = mSol.Evals;
    Teuchos::RCP< MV >                      evecs = mSol.Evecs;

    // get number of eigenvalues to be computed from parameterlist
    moris::sint tNumEigVals = mParameterList.get< moris::sint >( "Num_Eig_Vals" );

    // Make sure that the number of computed Eigenpairs and the number of requested Eigenvalues match
    if ( aEigValIndex >= (uint)mSol.numVecs )
    {
        fprintf( stdout, "numVecs = %i, tNumEigVals = %i, eigValIndex = %i\n", mSol.numVecs, tNumEigVals, aEigValIndex );
        fprintf( stdout, "Possible solution to fix this error, set mBlockSize = tNumEigVals in your input file\n" );
        MORIS_ERROR( false, "EigenSolver::GetSolution: Number of computed Eigenpairs is smaller than the number of requested Eigenvalues! Will use largest computed eigenvector \n" );
        aEigValIndex = mSol.numVecs - 1;
    }

    if ( (sint)aEigValIndex < 0 )
    {
        MORIS_ERROR( false, "EigvalIndex is smaller zero; likely no eigen vector computed; check previous warnings\n" );
    }

    // Get the real and imaginary part of the current eigenvalue
    aEigValReal = evals[ aEigValIndex ].realpart;
    aEigValImag = evals[ aEigValIndex ].imagpart;

    // Extract the current eigenvector from evecs
    Epetra_MultiVector* CurrentEigenVector = evecs->operator()( (int)aEigValIndex );

    // Get the current normalization factor  = ||Eigvec*MassMat*Eigvec||
    real NormFactor = 1.0 / mNormEvecsMassMatEvecs.at( aEigValIndex );

    // Normalize eigenvector with NormFactor
    mFreeSolVec->get_epetra_vector()->Update( NormFactor, *CurrentEigenVector, 0.0 );

    // Add index of eigen vector
    char tIndex[ 12 ];
    std::sprintf( tIndex, "_%01u", aEigValIndex + 1 );

    // Assemble name of the normalized eigen vector
    char tEigVector[ 100 ];
    std::strcpy( tEigVector, "NormalizedEigenVector" );
    std::strcat( tEigVector, tIndex );
    std::strcat( tEigVector, ".h5\0" );

    // save current eigen vector to HDF5 file
    mFreeSolVec->save_vector_to_HDF5( tEigVector );

    return 0;
}

//-----------------------------------------------------------------------------

void
Eigen_Solver::compute_preconditioned_operator()
{
    // 
//     const Epetra_Map& tBlockMap = mSPmat->OperatorDomainMap();
//     int tMatrixSize = tBlockMap.MaxAllGID(); 
    
//     Epetra_MultiVector tIn = Epetra_MultiVector(tBlockMap,tMatrixSize);
//     tIn.
}
