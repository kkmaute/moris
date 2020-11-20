/*
 * cl_DLA_Linear_Solver_Belos.cpp
 *
 *  Created on: Feb 06, 2020
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Solver_Belos.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"

#include "cl_DLA_Preconditioner_Trilinos.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include <BelosSolverManager.hpp>
#include "BelosSolverFactory.hpp"

#include "Teuchos_ParameterList.hpp"

#include "BelosEpetraAdapter.hpp"

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

using namespace moris;
using namespace dla;

Linear_Solver_Belos::Linear_Solver_Belos( const moris::ParameterList aParameterlist ) : Linear_Solver_Algorithm( aParameterlist )
{
}

Linear_Solver_Belos::Linear_Solver_Belos( Linear_Problem * aLinearSystem )
{
    // Set chosen solver options
    this->set_solver_parameters();
}

Linear_Solver_Belos::Linear_Solver_Belos()
{
    // Set chosen solver options
    this->set_solver_parameters();
}


Linear_Solver_Belos::~Linear_Solver_Belos()
{
    //    delete mAmesosSolver;
}

void Linear_Solver_Belos::set_solver_parameters()
{
//    // ASSIGN DEFAULT PARAMETER VALUES
//    // https://docs.trilinos.org/dev/packages/belos/doc/html/classBelos_1_1SolverFactory.html#ad86e61fb180a73c6dd5dbf458df6a86f
//
//    // Determine which solver is used by string
//    // options are: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
//    //              Pseudoblock TFQMR, Seed GMRES, Seed CG
//    mParameterList.insert( "Solver Type" ,  "GMRES" );
//
//    mParameterList.insert( "Verbosity" ,  INT_MAX );
//
//    // Allowable Aztec solver iterations
//    mParameterList.insert( "Block Size", INT_MAX   );
//
//    // Allowable Belos solver iterations
//    mParameterList.insert( "Maximum Iterations" , INT_MAX );
//
//    // set Az_conv -convergence criteria
//    // options are
//    mParameterList.insert( "Convergence Tolerance" ,  1e-08 );

    mParameterList = prm::create_linear_algorithm_parameter_list_belos();

    mParameterList.set( "ifpack_prec_type", "ILU");
}


moris::sint Linear_Solver_Belos::solve_linear_system(
        Linear_Problem *  aLinearSystem,
        const moris::sint aIter )
{
    this->set_solver_internal_parameters();
    mLinearSystem = aLinearSystem;

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Belos::SolverFactory;


    RCP<Belos::EpetraPrecOp> belosPrec;

    if( !mParameterList.get< std::string >( "ifpack_prec_type" ).empty() )
    {
        Preconditioner_Trilinos tPreconditioiner( mParameterList, mLinearSystem );
        tPreconditioiner.build_ifpack_preconditioner();
        belosPrec = rcp ( new Belos::EpetraPrecOp ( tPreconditioiner.get_ifpack_prec() ) );
    }
    else if( !mParameterList.get< std::string >( "ml_prec_type" ).empty() )
    {
        Preconditioner_Trilinos tPreconditioiner( mParameterList, mLinearSystem );
        tPreconditioiner.build_ml_preconditioner();
        belosPrec = rcp ( new Belos::EpetraPrecOp ( tPreconditioiner.get_ml_prec() ) );
    }
    else
    {
        MORIS_ERROR( false, "Linear_Solver_Belos::solve_linear_system(), no preconditioner specified");
    }

//    // FIXME move preconditioners in own class
//    // =============================================================== //
//    // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
//    // =============================================================== //
//
//    // Allocate an IFPACK factory.  The object contains no data, only
//    // the Create() method for creating preconditioners.
//    Ifpack Factory;
//
//    // Hard-code options - should be done through input parameters
//    std::string PrecType; // 0: incomplete LU 1: ILUT (not working)
//    int tPrecTypeInd = 0;
//    int OverlapLevel = 1;
//
//    ParameterList List;
//
//    // Specify parameters for ILU.  ILU is local to each MPI process.
//    switch (tPrecTypeInd)
//    {
//        case 0:
//        {
//            PrecType = "ILU";
//            List.set (" fact: drop tolerance"    , 1e-9 );
//            List.set (" fact: level-of-fill"     , 1    );
//            break;
//        }
//        case 1:
//        {
//            PrecType = "ILUT";
//            List.set ( "fact: drop tolerance"    , 1.e-9 );
//            List.set ( "fact: ilut level-of-fill", 12.0  );
//            break;
//        }
//        default:
//        {
//            MORIS_ERROR(false,"Incorrect preconditioner type.\n");
//        }
//    }
//
//    List.set ( "schwarz: combine mode"     , "Add");
    //List.set ( "schwarz: combine mode"     , "Zero");
    //List.set ( "schwarz: compute condest"  , false );
    //List.set ( "schwarz: filter singletons", false );
    //List.set ( "schwarz: reordering type"  , "rcm" );

    // Get pointer to operator
    //Epetra_RowMatrix & tOperator = *(aLinearSystem->get_matrix()->get_matrix());

    //EpetraExt::RowMatrixToMatlabFile( "mat.mat", tOperator );

    // Create the preconditioner.  For the list of PrecType values that
    //RCP< Ifpack_Preconditioner > Prec = rcp (Factory.Create (PrecType, tOperator, OverlapLevel));

//    RCP< Ifpack_Preconditioner > Prec =  rcp (Factory.Create (PrecType, &*aLinearSystem->get_matrix()->get_matrix(), OverlapLevel));
//
//    // Set the parameters.
//    IFPACK_CHK_ERR(Prec->SetParameters(List));
//
//    // Initialize the preconditioner.
//    IFPACK_CHK_ERR(Prec->Initialize());
//
//    // Build the preconditioner, by looking at the values of the matrix.
//    IFPACK_CHK_ERR(Prec->Compute());


    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    //    RCP<Belos::EpetraPrecOp> belosPrec = rcp ( new Belos::EpetraPrecOp ( Prec ) );


    // =================================================== //
    // E N D   O F   I F P A C K   C O N S T R U C T I O N //
    // =================================================== //

    //      aLinearSystem->get_free_solver_LHS()->print();
    //      aLinearSystem->get_solver_RHS()->print();

    RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > problem =
            rcp (new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>(
                    rcp( dynamic_cast< Epetra_CrsMatrix* > ( aLinearSystem->get_matrix()->get_matrix() ), false ),
                    rcp( dynamic_cast<Vector_Epetra*>(aLinearSystem->get_free_solver_LHS())->get_epetra_vector(), false ),
                    rcp( dynamic_cast<Vector_Epetra*>(aLinearSystem->get_solver_RHS())->get_epetra_vector(), false ) ) );

    problem->setLeftPrec( belosPrec );

    bool set = problem->setProblem();
    if (set == false)
    {
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    }

    //    RCP< Belos::SolverManager<double,MV,OP> > newSolver = rcp( new Belos::GCRODRSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));
    //    //
    //    // **********Print out information about problem*******************
    //    //
    ////    if (proc_verbose) {
    ////      std::cout << std::endl << std::endl;
    ////      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    ////      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    ////      std::cout << "Max number of restarts allowed: " << maxrestarts << std::endl;
    ////      std::cout << "Max number of iterations per restart cycle: " << maxiters << std::endl;
    ////      std::cout << "Relative residual tolerance: " << tol << std::endl;
    ////      std::cout << std::endl;
    ////    }


    SolverFactory<double, Epetra_MultiVector,Epetra_Operator> factory;
    // Create the GMRES solver.
    std::string tSolverType = mParameterList.get< std::string >( "Solver Type" );
    RCP<Belos::SolverManager<double, Epetra_MultiVector,Epetra_Operator> > solver =  factory.create ( tSolverType, mMyPl );

    // Tell the solver what problem you want to solve.
    solver->setProblem (problem);

    //Belos::ReturnType result = solver->solve();
    solver->solve();
    // Ask the solver how many iterations the last solve() took.
    //const int numIters = solver->getNumIters();
    //
    //std::cout<<"iter : "<<numIters<<std::endl;

    return 0;
}

void Linear_Solver_Belos::set_solver_internal_parameters()
{
    mMyPl = Teuchos::parameterList();

    if (mParameterList.get< moris::sint >( "Verbosity" ) != INT_MAX)
    {
        mMyPl->set( "Verbosity", mParameterList.get< moris::sint >( "Verbosity" ) );
    }

    if (mParameterList.get< moris::sint >( "Block Size" ) != INT_MAX)
    {
        mMyPl->set ( "Block Size", mParameterList.get< moris::sint >( "Block Size" ) );
    }

    if (mParameterList.get< moris::sint >( "Maximum Iterations" ) != INT_MAX)
    {
        mMyPl->set ( "Maximum Iterations", mParameterList.get< moris::sint >( "Maximum Iterations" ) );
    }

    if (mParameterList.get< moris::real >( "Convergence Tolerance" ) != 1e-08)
    {
        mMyPl->set ( "Convergence Tolerance", mParameterList.get< moris::real >( "Convergence Tolerance" ) );
    }
}
