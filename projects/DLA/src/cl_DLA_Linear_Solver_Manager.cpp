
#include "cl_DLA_Linear_Solver_Manager.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace dla;

Linear_Solver_Manager::Linear_Solver_Manager()
{
    // create solver factory
    Solver_Factory  tSolFactory;

    // create solver object
    std::shared_ptr< Linear_Solver > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
    std::shared_ptr< Linear_Solver > tLinSolver2 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
    std::shared_ptr< Linear_Solver > tLinSolver3 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
    std::shared_ptr< Linear_Solver > tLinSolver4 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
    std::shared_ptr< Linear_Solver > tLinSolver5 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

    tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
    tLinSolver1->set_param("AZ_output") = AZ_none;
//    tLinSolver1->set_param("AZ_solver") = AZ_cg;

    mLinearSolverList.clear();

    mLinearSolverList.push_back( tLinSolver1 );
    mLinearSolverList.push_back( tLinSolver2 );
    mLinearSolverList.push_back( tLinSolver3 );
    mLinearSolverList.push_back( tLinSolver4 );
    mLinearSolverList.push_back( tLinSolver5 );

    this->set_linear_solver_manager_parameters();
}

Linear_Solver_Manager::~Linear_Solver_Manager()
{}

//--------------------------------------------------------------------------------------------------
void Linear_Solver_Manager::set_linear_solver( std::shared_ptr< Linear_Solver > aLinSolver )
{
    if( mCallCounter == 0 )
    {
        mLinearSolverList.clear();

        mLinearSolverList.push_back( aLinSolver );
    }
    else
    {
        mLinearSolverList.push_back( aLinSolver );
    }

    mCallCounter = mCallCounter + 1;
}

//-------------------------------------------------------------------------------------------------------
void Linear_Solver_Manager::set_linear_solver( const moris::uint aListEntry,
                                                     std::shared_ptr< Linear_Solver > aLinSolver )
{
    mLinearSolverList( aListEntry ) = aLinSolver;
}

//-------------------------------------------------------------------------------------------------------
void Linear_Solver_Manager::solver_linear_system( std::shared_ptr< dla::Linear_Problem > aLinearProblem, const moris::sint aIter )
{
    moris::sint tErrorStatus = 0;
    moris::sint tMaxNumLinRestarts  = mParameterListLinearSolver.get< moris::sint >( "DLA_max_lin_solver_restarts" );
    moris::sint tTryRestartOnFailIt = 1;

    tErrorStatus = mLinearSolverList( 0 )->solve_linear_system( aLinearProblem, aIter );

    // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed
    while ( tErrorStatus !=0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && ( moris::sint )mLinearSolverList.size() <= tMaxNumLinRestarts )
    {
        if ( par_rank() == 0 )
        {
            // Compute current solution vector norm
            moris::real tSolVecNorm = aLinearProblem->get_free_solver_LHS()->vec_norm2();

            fprintf( stdout, " ... Previous linear solve failed. Trying restart %i of %i, using current solution with SolVecNorm = %5.15e as an initial guess. \n",
                                   tTryRestartOnFailIt, tMaxNumLinRestarts, tSolVecNorm);
        }

        // Re-solve scaled linear system with current solution as an initial guess
        tErrorStatus = mLinearSolverList( tTryRestartOnFailIt )->solve_linear_system( aLinearProblem, aIter );

        // Iterate TryRestartOnFailIt counter
        tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
    }

//    if ( ( tErrorStatus != 0 && mParameterListLinearSolver.get< bool >( "DLA_hard_break" ) ) && !mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_on_fail" ) )
//    {
//        if( par_rank() == 0 )
//        {
//            fprintf( stdout, "\n Linear Solver status absolute value = %i\n", tErrorStatus );
//            MORIS_ERROR( false, "Linear Solver did not exit with status 0!\n" );
//        }
//    }

    if( ( tErrorStatus != 0 ) )
    {
        if( par_rank() == 0)
        {
            fprintf( stdout, "\n Linear Solver status absolute value = %i\n", tErrorStatus );
            fprintf( stdout, "Linear Solver did not exit with status 0!\n" );
        }
    }

//    if ( tErrorStatus != 0 && mParameterListLinearSolver.get< bool >( "DLA_hard_break" ) )
//    {
//        aIter = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_iter" );
//        aHardBreak = true;
//    }
}

//--------------------------------------------------------------------------------------------------------------------------
    void Linear_Solver_Manager::set_linear_solver_manager_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListLinearSolver.insert( "DLA_max_lin_solver_restarts" , 0 );

        // Maximal number of linear solver restarts on fail
        mParameterListLinearSolver.insert( "DLA_hard_break" , true );

        // Determines if lin solve should restart on fail
        mParameterListLinearSolver.insert( "DLA_rebuild_lin_solver_on_fail" , false );
    }

