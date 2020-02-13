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

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include <BelosSolverManager.hpp>
#include "BelosSolverFactory.hpp"

#include "Teuchos_ParameterList.hpp"

#include "BelosEpetraAdapter.hpp"


using namespace moris;
using namespace dla;

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
    // ASSIGN DEFAULT PARAMETER VALUES
    // https://docs.trilinos.org/dev/packages/belos/doc/html/classBelos_1_1SolverFactory.html#ad86e61fb180a73c6dd5dbf458df6a86f
	
	

    // Determine which solver is used by string
    // options are: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
	//              Pseudoblock TFQMR, Seed GMRES, Seed CG
    mParameterList.insert( "Solver Type" ,  "GMRES" );
		
    mParameterList.insert( "Verbosity" ,  INT_MAX );

    // Allowable Aztec solver iterations
    mParameterList.insert( "Block Size", INT_MAX   );

    // Allowable Belos solver iterations
    mParameterList.insert( "Maximum Iterations" , INT_MAX );

    // set Az_conv -convergence criteria
    // options are 
    mParameterList.insert( "Convergence Tolerance" ,  1e-08 );
}


moris::sint Linear_Solver_Belos::solve_linear_system(       Linear_Problem * aLinearSystem,
                                                      const moris::sint      aIter )
{
	this->set_solver_internal_parameters();
    mLinearSystem = aLinearSystem;

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Belos::SolverFactory;


    RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > problem
                   = rcp (new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> (rcp( dynamic_cast< Epetra_CrsMatrix* > ( aLinearSystem->get_matrix()->get_matrix() ), false ),
                                                                                               rcp( aLinearSystem->get_free_solver_LHS()->get_epetra_vector(), false ),
                                                                                               rcp( aLinearSystem->get_solver_RHS()->get_epetra_vector(), false ) ) );

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
//	std::string tSolverType = mParameterList.get< std::string >( "Solver Type" );
    RCP<Belos::SolverManager<double, Epetra_MultiVector,Epetra_Operator> > solver =  factory.create ( "GMRES", mMyPl );

    // Tell the solver what problem you want to solve.
    solver->setProblem (problem);

//	Belos::ReturnType result = solver->solve();
	solver->solve();
	// Ask the solver how many iterations the last solve() took.
	const int numIters = solver->getNumIters();

	std::cout<<"iter : "<<numIters<<std::endl;

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
