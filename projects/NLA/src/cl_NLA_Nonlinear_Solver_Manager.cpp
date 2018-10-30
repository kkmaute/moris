
#include "cl_NLA_Nonlinear_Solver_Manager.hpp"

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;

Nonlinear_Solver_Manager::Nonlinear_Solver_Manager()
{
//    // create solver factory
//    Solver_Factory  tSolFactory;
//
//    // create solver object
//    std::shared_ptr< Linear_Solver > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//    std::shared_ptr< Linear_Solver > tLinSolver2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//    std::shared_ptr< Linear_Solver > tLinSolver3 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//    std::shared_ptr< Linear_Solver > tLinSolver4 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//    std::shared_ptr< Linear_Solver > tLinSolver5 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//    tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
//    tLinSolver1->set_param("AZ_output") = AZ_none;
////    tLinSolver1->set_param("AZ_solver") = AZ_cg;
//
//    mLinearSolverList.clear();
//
//    mLinearSolverList.push_back( tLinSolver1 );
//    mLinearSolverList.push_back( tLinSolver2 );
//    mLinearSolverList.push_back( tLinSolver3 );
//    mLinearSolverList.push_back( tLinSolver4 );
//    mLinearSolverList.push_back( tLinSolver5 );

    this->set_nonlinear_solver_manager_parameters();
}

Nonlinear_Solver_Manager::~Nonlinear_Solver_Manager()
{}

//--------------------------------------------------------------------------------------------------
void Nonlinear_Solver_Manager::set_nonlinear_solver( std::shared_ptr< Nonlinear_Solver > aNonLinSolver )
{
    if( mCallCounter == 0 )
    {
        mNonLinearSolverList.clear();

        mNonLinearSolverList.push_back( aNonLinSolver );
    }
    else
    {
        mNonLinearSolverList.push_back( aNonLinSolver );
    }

    mCallCounter = mCallCounter + 1;
}

//-------------------------------------------------------------------------------------------------------
void Nonlinear_Solver_Manager::set_nonlinear_solver( const moris::uint aListEntry,
                                                     std::shared_ptr< Nonlinear_Solver > aNonLinSolver )
{
    mNonLinearSolverList( aListEntry ) = aNonLinSolver;
}

//-------------------------------------------------------------------------------------------------------
void Nonlinear_Solver_Manager::solver_nonlinear_system( std::shared_ptr< Nonlinear_Problem > aLinearProblem, const moris::sint aIter )
{
//    moris::sint tErrorStatus = 0;
//    moris::sint tMaxNumLinRestarts  = mParameterListLinearSolver.get< moris::sint >( "DLA_max_lin_solver_restarts" );
//    moris::sint tTryRestartOnFailIt = 1;
//
//    tErrorStatus = mLinearSolverList( 0 )->solve_linear_system( aLinearProblem, aIter );
//
//    // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed
//    while ( tErrorStatus !=0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && ( moris::sint )mLinearSolverList.size() <= tMaxNumLinRestarts )
//    {
//        if ( par_rank() == 0 )
//        {
//            // Compute current solution vector norm
//            moris::real tSolVecNorm = aLinearProblem->get_free_solver_LHS()->vec_norm2();
//
//            fprintf( stdout, " ... Previous linear solve failed. Trying restart %i of %i, using current solution with SolVecNorm = %5.15e as an initial guess. \n",
//                                   tTryRestartOnFailIt, tMaxNumLinRestarts, tSolVecNorm);
//        }
//
//        // Re-solve scaled linear system with current solution as an initial guess
//        tErrorStatus = mLinearSolverList( tTryRestartOnFailIt )->solve_linear_system( aLinearProblem, aIter );
//
//        // Iterate TryRestartOnFailIt counter
//        tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
//    }
//
//
//
//    if( ( tErrorStatus != 0 ) )
//    {
//        if( par_rank() == 0)
//        {
//            fprintf( stdout, "\n Linear Solver status absolute value = %i\n", tErrorStatus );
//            fprintf( stdout, "Linear Solver did not exit with status 0!\n" );
//        }
//    }
}

//--------------------------------------------------------------------------------------------------------------------------
    void Nonlinear_Solver_Manager::set_nonlinear_solver_manager_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListNonLinearSolver.insert( "NLA_max_non_lin_solver_restarts" , 0 );
    }

