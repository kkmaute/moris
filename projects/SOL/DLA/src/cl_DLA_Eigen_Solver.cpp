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
#include "typedefs.hpp"
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

using namespace moris;
using namespace dla;

// ----------------------------------------------------------------------------

Eigen_Solver::Eigen_Solver()
{
    // set eigen solver manager parameters
    this->set_eigen_solver_manager_parameters();
}

// ----------------------------------------------------------------------------

Eigen_Solver::Eigen_Solver( const ParameterList* aParameterList )
        : Linear_Solver_Algorithm( *aParameterList )
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
        std::cout << "Eigen Algorithm is Block Davidson" << std::endl;
    }
    else if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_GENERALIZED_DAVIDSON" )
    {
        std::cout << "Eigen Algorithm is Generalized Davidson" << std::endl;
    }
    else if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_BLOCK_KRYLOV_SCHUR" )
    {
        std::cout << "Eigen Algorithm is Block Krylov Schur" << std::endl;
    }
    else if ( mParameterList.get< std::string >( "Eigen_Algorithm" ) == "EIGALG_BLOCK_KRYLOV_SCHUR_AMESOS" )
    {
        std::cout << "Eigen Algorithm is Block Krylov Schur Amesos " << std::endl;
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

void
Eigen_Solver::build_linearized_system()
{
    // sparse-stiffness matrix of Teuchos::RCP type
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
    mFreeSolVec = new Vector_Epetra( mMap, 1, true, false );

    // Make sure that the number of blocks and eigenvalues does not exceed NumFreeDofs
    if ( tBlockSize > tNumFreeDofs || tNumEigVals > tNumFreeDofs )
    {
        MORIS_ERROR( false, "EigenSolver::BuildLinearSystem: Number of blocks and/or number of eigenvalues can not exceed NumFreeDofs\n" );
    }

    // Create initial vector for the solver
    mIvec = Teuchos::RCP( new Epetra_MultiVector( mSPmat->OperatorDomainMap(), tBlockSize ) );
    mIvec->Random();

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

void
Eigen_Solver::set_eigen_solver_manager_parameters()
{
    // create eigen solver parameterlist
    mParameterList = prm::create_eigen_algorithm_parameter_list();
}

// -----------------------------------------------------------------------------

int
Eigen_Solver::solve_block_davidson_system( Linear_Problem* aLinearSystem )
{
    int MyPID = par_rank();

    // build linear system
    this->build_linearized_system();

    // Set the pre-conditioner
    Teuchos::RCP< Epetra_Operator >                     tPrecOp;
    Teuchos::RCP< Ifpack_Preconditioner >               tIfpackPrec;
    Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > tMlPrec;

    // initialize preconditioner
    mPrec.initialize( mParameterList, aLinearSystem );

    // request ifpack type preconditioner as a string
    std::string tIfpackPrectype = mParameterList.get< std::string >( "ifpack_prec_type" );

    // request solver type for ifpack preconditioner as a string
    std::string tSolvtype = mParameterList.get< std::string >( "amesos: solver type" );

    // request multilvel preconditioner type as a string
    std::string tmlprectype = mParameterList.get< std::string >( "ml_prec_type" );

    // build preconditioner
    if ( !tIfpackPrectype.empty() && !tSolvtype.empty() )
    {
        MORIS_LOG_INFO( "Constructing %s type preconditioner", tSolvtype.c_str() );
    }
    else if ( !tmlprectype.empty() )
    {
        MORIS_LOG_INFO( "Constructing %s type Multilevel Preconditioner", tmlprectype.c_str() );
    }
    else
    {
        MORIS_ERROR( false, "Incorrect preconditioner type" );
    }
    mPrec.build();

    // get ifpack preconditioner if it exists
    tIfpackPrec = mPrec.get_ifpack_prec();

    // get ml preconditioner if it exists
    tMlPrec = mPrec.get_ml_prec();

    // get preconditioned operator
    if ( tIfpackPrec.get() != NULL )
    {
        tPrecOp = Teuchos::RCP( new Epetra_InvOperator( tIfpackPrec.get() ) );
    }
    else if ( tMlPrec.get() != NULL )
    {
        tPrecOp = Teuchos::RCP( new Epetra_InvOperator( tMlPrec.get() ) );
    }

    // Set the pre-conditions to the eigenproblem
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

    // set verbosity
    MyPL.set( "Verbosity", tverbosity );

    // set sorting type for eigenvalues
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );

    // set block size
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );

    // set number of blocks
    MyPL.set( "Num Blocks", mParameterList.get< moris::sint >( "Num_Blocks" ) );

    // set maximum subspace dimensions
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );

    // set maximum restarts
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );

    // set convergence tolerance
    MyPL.set( "Convergence Tolerance", mParameterList.get< moris::real >( "Convergence_Tolerance" ) );

    // set relative convergence tolerance
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

    // Create the solver manager
    Anasazi::BlockDavidsonSolMgr< double, MV, OP > MySolverMan( mMyEigProblem, MyPL );

    // Solve the problem
    Anasazi::ReturnType tReturnCode = MySolverMan.solve();

    // Check if the problem solve converged
    if ( tReturnCode != Anasazi::Converged && MyPID == 0 )
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
        if ( evecs->GlobalLength() == tEpetraVec->GlobalLength() && evecs->NumVectors() == tEpetraVec->NumVectors() )
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
    if ( tverbose == 0 && MyPID == 0 )
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
        // set master solution vector to free solution vector
        Vector_Epetra* aMasterSolVec = mFreeSolVec;

        // get solution for eigen vectors
        this->get_solution( m, mFreeSolVec, aMasterSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, "\n ... Finished EigenValue Solve using a Block Davidson Solver\n" );
    }

    return 0;
}

// ----------------------------------------------------------------------------

int
Eigen_Solver::solve_generalized_davidson_system( Linear_Problem* aLinearSystem )
{
    using std::cout;
    using std::endl;

    int MyPID = par_rank();

    // build linear system
    this->build_linearized_system();

    // Set verbosity level
    bool tverbose   = mParameterList.get< bool >( "Verbosity" );
    int  tverbosity = Anasazi::Errors + Anasazi::Warnings;
    if ( tverbose )
    {
        tverbosity += Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
    }

    // Create parameter list to pass into solver
    Teuchos::ParameterList MyPL;

    // set verbosity
    MyPL.set( "Verbosity", tverbosity );

    // set sorting type for eigen values
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );

    // set block size
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );

    // set maximum subspace dimension
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );

    // set maximum restarts
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );

    // set convergence tolerance
    MyPL.set( "Convergence Tolerance", mParameterList.get< moris::real >( "Convergence_Tolerance" ) );

    // set relative convergence tolerance
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

    // Set the pre-conditioner
    Teuchos::RCP< Epetra_Operator >                     tPrecOp;
    Teuchos::RCP< Ifpack_Preconditioner >               tIfpackPrec;
    Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > tMlPrec;

    // initialize preconditioner
    mPrec.initialize( mParameterList, aLinearSystem );

    // request ifpack type preconditioner as a string
    std::string tIfpackPrectype = mParameterList.get< std::string >( "ifpack_prec_type" );

    // request solver type for ifpack preconditioner as a string
    std::string tSolvtype = mParameterList.get< std::string >( "amesos: solver type" );

    // request multilvel preconditioner type as a string
    std::string tmlprectype = mParameterList.get< std::string >( "ml_prec_type" );

    // build preconditioner
    if ( !tIfpackPrectype.empty() && !tSolvtype.empty() )
    {
        MORIS_LOG_INFO( "Constructing %s type preconditioner", tSolvtype.c_str() );
    }
    else if ( !tmlprectype.empty() )
    {
        MORIS_LOG_INFO( "Constructing %s type Multilevel Preconditioner", tmlprectype.c_str() );
    }
    else
    {
        MORIS_ERROR( false, "Incorrect preconditioner type" );
    }
    mPrec.build();

    // get ifpack preconditioner if it exists
    tIfpackPrec = mPrec.get_ifpack_prec();

    // get ml preconditioner if it exists
    tMlPrec = mPrec.get_ml_prec();

    // get preconditioned operator
    if ( tIfpackPrec.get() != NULL )
    {
        tPrecOp = Teuchos::RCP( new Epetra_InvOperator( mPrec.get_ifpack_prec().get() ) );
    }
    else if ( tMlPrec.get() != NULL )
    {
        tPrecOp = Teuchos::RCP( new Epetra_InvOperator( mPrec.get_ml_prec().get() ) );
    }
    else
    {
        MORIS_ERROR( false, "Incorrect preconditioner type" );
    }

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
        if ( evecs->GlobalLength() == tEpetraVec->GlobalLength() && evecs->NumVectors() == tEpetraVec->NumVectors() )
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
        if ( tverbose == 0 && MyPID == 0 )
        {
            MORIS_LOG_INFO( "===============================================" );
            for ( int j = 0; j < mNumReturnedEigVals; j++ )
            {
                MORIS_LOG_INFO( "Real Part: %16f ", evals[ j ].realpart );
                MORIS_LOG_INFO( "Imaginary Part: %16f", evals[ j ].imagpart );
                MORIS_LOG_INFO( "Direct Residual: %18e ", normR[ j ] );
            }
        }    // end of if verbose
    }

    // print eigenvector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set master solution vector to free solution vector
        Vector_Epetra* aMasterSolVec = mFreeSolVec;

        // get solution for eigen vectors
        this->get_solution( m, mFreeSolVec, aMasterSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, "\n ... Finished EigenValue Solve using a Generalized Davidson Solver\n\n" );
    }

    return 0;
}

// -----------------------------------------------------------------------------

int
Eigen_Solver::solve_block_krylov_schur_system( Linear_Problem* aLinearSystem )
{
    // NOTE: This implementation follows trilinos/packages/anasazi/epetra/example/BlockKrylovSchur/BlockKrylovSchurEpetraExGenAztecOO.cpp

    int MyPID = par_rank();

    // Set the preconditioner
    Teuchos::RCP< Epetra_Operator >                     tPrecOp;
    Teuchos::RCP< Ifpack_Preconditioner >               tIfpackPrec;
    Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > tMlPrec;

    // initialize preconditioner
    mPrec.initialize( mParameterList, aLinearSystem );

    // request ifpack type preconditioner as a string
    std::string tIfpackPrectype = mParameterList.get< std::string >( "ifpack_prec_type" );

    // request solver type for ifpack preconditioner as a string
    std::string tSolvtype = mParameterList.get< std::string >( "amesos: solver type" );

    // request multilvel preconditioner type as a string
    std::string tmlprectype = mParameterList.get< std::string >( "ml_prec_type" );

    // build preconditioner
    if ( !tIfpackPrectype.empty() && !tSolvtype.empty() )
    {
        MORIS_LOG_INFO( "Constructing %s type preconditioner", tSolvtype.c_str() );
    }
    else if ( !tmlprectype.empty() )
    {
        MORIS_LOG_INFO( "Constructing %s type Multilevel Preconditioner", tmlprectype.c_str() );
    }
    else
    {
        MORIS_ERROR( false, "Incorrect preconditioner type" );
    }

    // build preconditioner
    mPrec.build();

    // get ifpack preconditioner if it exists
    tIfpackPrec = mPrec.get_ifpack_prec();

    // get Multilevel preconditioner if it exists
    tMlPrec = mPrec.get_ml_prec();

    // Tell the linear problem about the mass matrix M
    mEpetraProblem.SetOperator( mSPmassmat.getRawPtr() );

    // Create AztecOO iterative solver for solving linear systems with K.
    AztecOO precSolver( mEpetraProblem );

    // Tell the solver to use the Ifpack preconditioner we created above.
    if ( tIfpackPrec.get() != NULL )
    {
        precSolver.SetPrecOperator( tIfpackPrec.get() );
    }
    else if ( tMlPrec.get() != NULL )
    {
        precSolver.SetPrecOperator( tMlPrec.get() );
    }
    else
    {
        MORIS_ERROR( false, "No valid preconditioner set" );
    };

    // Set AztecOO solver options.
    precSolver.SetAztecOption( AZ_output, AZ_none );     // Don't print output
    precSolver.SetAztecOption( AZ_solver, AZ_gmres );    // Use GMRES

    // Use the above AztecOO solver to create the AztecOO_Operator.
    Teuchos::RCP< AztecOO_Operator > precOperator = Teuchos::rcp( new AztecOO_Operator( &precSolver, mSPmat->NumGlobalRows(), 1e-12 ) );

    // Create an Operator that computes y = M^{-1} K x.
    Teuchos::RCP< Anasazi::EpetraGenOp > Aop = Teuchos::rcp( new Anasazi::EpetraGenOp( precOperator, mSPmat ) );

    // Create the eigen problem object
    mMyEigProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem< double, MV, OP >( Aop, mSPmat, mIvec ) );

    // request number of eigenvalues from parameterlist
    moris::sint tNumEigVals = mParameterList.get< moris::sint >( "Num_Eig_Vals" );

    // set whether linear system is symmetric or not
    bool tAssumeSymmetric = false;

    // Set whether AMat and BMat are Hermitian = symmetric
    mMyEigProblem->setHermitian( tAssumeSymmetric );

    // Set the number of eigenvalues to be computed
    mMyEigProblem->setNEV( tNumEigVals );

    // Set the pre-conditions to the eigenproblem
    mMyEigProblem->setPrec( tPrecOp );

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

    // set verbosity
    MyPL.set( "Verbosity", tverbosity );

    // set sorting type for eigenvalues
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );

    // set block-size
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );

    // set number of blocks
    MyPL.set( "Num Blocks", mParameterList.get< moris::sint >( "NumBlocks" ) );

    // set maximum subspace dimensions
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );

    // set maximum restarts
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );

    // set convergence tolerance
    MyPL.set( "Convergence Tolerance", mParameterList.get< real >( "Covergence_Tolerance" ) );

    // set relative convergence tolerance
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
        if ( evecs->GlobalLength() == tEpetraVec->GlobalLength() && evecs->NumVectors() == tEpetraVec->NumVectors() )
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
    if ( tverbose && MyPID == 0 )
    {
        MORIS_LOG_INFO( "===============================================" );
        for ( int i = 0; i < mSol.numVecs; i++ )
        {
            MORIS_LOG_INFO( "EigenValue: %16f", evals[ i ].realpart );
            MORIS_LOG_INFO( "Direct Residual: %18e", normR[ i ] / evals[ i ].realpart );
        }
    }

    // print eigen vector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set master solution vector to free solution vector
        Vector_Epetra* aMasterSolVec = mFreeSolVec;

        // get solution for eigen vector
        this->get_solution( m, mFreeSolVec, aMasterSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, "\n ... Finished EigenValue Solve using a Block Krylov Schur Iterative (AZTEC) Solver\n\n" );
    }

    return 0;
}

// -----------------------------------------------------------------------------

int
Eigen_Solver::solve_block_krylov_schur_amesos_system( Linear_Problem* aLinearSystem )
{
    // NOTE: This implementation follows trilinos/packages/anasazi/epetra/example/BlockKrylovSchur/BlockKrylovSchurEpetraExGenAmesos.cpp

    int MyPID = par_rank();

    // Tell the linear problem about the mass matrix M
    mEpetraProblem.SetOperator( mSPmassmat.getRawPtr() );

    // Create Amesos factory and solver for solving linear systems with K
    Amesos tAmesosFactory;

    // Note that the AmesosProblem object "absorbs" M.  Anasazi doesn't see M, just the operator that implements K^{-1} M
    Teuchos::RCP< Amesos_BaseSolver > mAmesosSolver;

    // Determine which direct solver method to use
    switch ( (uint)( mEigSolMethod ) )
    {
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_KLU ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Klu", mEpetraProblem ) );
            break;
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_UMFPACK ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Umfpack", mEpetraProblem ) );
            break;
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_DSCPACK ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Dscpack", mEpetraProblem ) );
            break;
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_MUMPS ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Mumps", mEpetraProblem ) );
            break;
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_LAPACK ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Lapack", mEpetraProblem ) );
            break;
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_SCALAPACK ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Scalapack", mEpetraProblem ) );
            break;
        case (uint)( sol::EigSolMethod::LINSOL_AMESOS_PARDISO ):
            mAmesosSolver = Teuchos::rcp( tAmesosFactory.Create( "Amesos_Pardiso", mEpetraProblem ) );
            break;
        default:
        {
            MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurAmesosSystem: The specified mEigSolMethod is not yet supported!\n" );
        }
    }

    // The AmesosGenOp class assumes that the symbolic and numeric factorizations have already been performed on the linear problem
    mAmesosSolver->SymbolicFactorization();
    mAmesosSolver->NumericFactorization();

    // Create the Epetra_Operator for the spectral transformation using the Amesos direct solver
    Teuchos::RCP< Amesos_GenOp > Aop = Teuchos::rcp( new Amesos_GenOp( mAmesosSolver, mSPmassmat ) );

    // Create the eigen problem object
    mMyEigProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem< double, MV, OP >( Aop, mSPmat, mIvec ) );

    // request number of eigenvalues from parameterlist
    moris::sint tNumEigVals = mParameterList.get< moris::sint >( "Num_Eig_Vals" );

    // set linear system to be symmetric or not
    bool tAssumeSymmetric = false;

    // Set whether AMat and BMat are Hermitian = symmetric
    mMyEigProblem->setHermitian( tAssumeSymmetric );

    // Set the number of eigenvalues to be computed
    mMyEigProblem->setNEV( tNumEigVals );

    // check if the eigenproblem is set or not
    bool tboolret = mMyEigProblem->setProblem();

    // Check if an error was returned
    if ( tboolret != true && MyPID == 0 )
    {
        MORIS_ERROR( false, "EigenSolver::SolveBlockKrylovSchurAmesosSystem: Error was returned when setting up the problem!\n" );
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

    // set verbosity
    MyPL.set( "Verbosity", tverbosity );

    // set sorting type for eigenvalues
    MyPL.set( "Which", mParameterList.get< std::string >( "Which" ) );

    // set block-size
    MyPL.set( "Block Size", mParameterList.get< moris::sint >( "Block_Size" ) );

    // set number of blocks
    MyPL.set( "Num Blocks", mParameterList.get< moris::sint >( "NumBlocks" ) );

    // set maximum subspace dimension
    MyPL.set( "Maximum SubSpace Dimension", mParameterList.get< moris::sint >( "MaxSubSpaceDims" ) );

    // set maximum restart
    MyPL.set( "Maximum Restarts", mParameterList.get< moris::sint >( "MaxRestarts" ) );

    // set convergence tolerance
    MyPL.set( "Convergence Tolerance", mParameterList.get< real >( "Covergence_Tolerance" ) );

    // set relative convergence tolerance
    MyPL.set( "Relative Convergence Tolerance", mParameterList.get< bool >( "Relative_Convergence_Tolerance" ) );

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
        if ( evecs->GlobalLength() == tEpetraVec->GlobalLength() && evecs->NumVectors() == tEpetraVec->NumVectors() )
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
    if ( tverbose && MyPID == 0 )
    {
        MORIS_LOG_INFO( "------------------------------------------------" );
        for ( int i = 0; i < mSol.numVecs; i++ )
        {
            MORIS_LOG_INFO( "EigenValue: %16f", evals[ i ].realpart );
            MORIS_LOG_INFO( "Direct Residual: %18e", normR[ i ] / evals[ i ].realpart );
        }
    }

    // print eigen vector
    for ( int m = 0; m < mNumReturnedEigVals; m++ )
    {
        // set master solution vector to free solution vector
        Vector_Epetra* aMasterSolVec = mFreeSolVec;

        // get solution for eigen vector
        this->get_solution( m, mFreeSolVec, aMasterSolVec, evals[ m ].realpart, evals[ m ].imagpart );
    }

    if ( MyPID == 0 )
    {
        fprintf( stdout, " ... Finished EigenValue Solve using a Block Krylov Schur Direct (AMESOS) Solver\n\n" );
    }

    return 0;
}

// ----------------------------------------------------------------------------

int
Eigen_Solver::get_solution(
        uint           aEigValIndex,
        Vector_Epetra* aSolVec,
        Vector_Epetra* aMasterSolVec,
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
    char tIndex[ 3 ];
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