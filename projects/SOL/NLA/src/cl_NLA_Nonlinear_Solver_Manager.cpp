
#include "cl_NLA_Nonlinear_Solver_Manager.hpp"
#include "cl_NLA_Nonlinear_Manager.hpp"

#include "cl_DLA_Solver_Interface.hpp"

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;

Nonlinear_Solver_Manager::Nonlinear_Solver_Manager( const enum NonlinearSolverType aNonLinSolverType )
{
    mNonLinSolverType = aNonLinSolverType;

    // create solver factory
    Nonlinear_Solver_Factory  tSolFactory;

    mLinSolManager = new dla::Linear_Solver_Manager();

    // create solver object
    std::shared_ptr< Nonlinear_Solver > tNonLinSolver1 = tSolFactory.create_nonlinear_solver( aNonLinSolverType );
    std::shared_ptr< Nonlinear_Solver > tNonLinSolver2 = tSolFactory.create_nonlinear_solver( aNonLinSolverType );

    tNonLinSolver1->set_linear_solvers( mLinSolManager );
    tNonLinSolver2->set_linear_solvers( mLinSolManager );

    mNonLinearSolverList.clear();

    mNonLinearSolverList.push_back( tNonLinSolver1 );
    mNonLinearSolverList.push_back( tNonLinSolver2 );

    mStaggeredDofTypeList.resize( 0 );

    //--------------------------------------------------------------------------------------------------
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

    void Nonlinear_Solver_Manager::set_nonlinear_manager( Nonlinear_Manager * aNonlinearManager )
    {
        mNonlinearManager = aNonlinearManager;

         mSolverInput = mNonlinearManager->get_solver_interface() ;
    }

//-------------------------------------------------------------------------------------------------------

    void Nonlinear_Solver_Manager::solve()
    {
        //MORIS_ERROR( mNonLinSolverType != NonlinearSolverType::NLBGS_SOLVER, "Nonlinear_Solver_Manager::solve(); Nonlinear Solver is not NLBGS" );

        moris::Cell< enum MSI::Dof_Type > tDofTypeUnion = this->get_dof_type_union();

        mSolverInput->set_requested_dof_types( tDofTypeUnion );

        if ( mNonLinSolverType == NonlinearSolverType::NLBGS_SOLVER)
        {
            mNonlinerProblem = new Nonlinear_Problem( mSolverInput, false );
        }
        else
        {
            mNonlinerProblem = new Nonlinear_Problem( mSolverInput );
        }

        mNonLinearSolverList( 0 )->set_nonlinear_solver_manager( this );

        mNonLinearSolverList( 0 )->solver_nonlinear_system( mNonlinerProblem );
    }

    void Nonlinear_Solver_Manager::solve( Nonlinear_Problem * aNonlinearProblem )
    {

        MORIS_ERROR( mNonLinSolverType != NonlinearSolverType::NLBGS_SOLVER, "Nonlinear_Solver_Manager::solve(); Nonliner Solver is NLBGS" );

        mNonLinearSolverList( 0 )->solver_nonlinear_system( aNonlinearProblem );

        // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed
//        while ( tErrorStatus !=0 && tTryRestartOnFailIt <= tMaxNumLinRestarts && ( moris::sint )mLinearSolverList.size() <= tMaxNumLinRestarts )
//        {
//            if ( par_rank() == 0 )
//            {
//                // Compute current solution vector norm
//                moris::real tSolVecNorm = aLinearProblem->get_free_solver_LHS()->vec_norm2();
//
//                fprintf( stdout, " ... Previous linear solve failed. Trying restart %i of %i, using current solution with SolVecNorm = %5.15e as an initial guess. \n",
//                                       tTryRestartOnFailIt, tMaxNumLinRestarts, tSolVecNorm);
//            }
//
//            // Re-solve scaled linear system with current solution as an initial guess
//            tErrorStatus = mLinearSolverList( tTryRestartOnFailIt )->solve_linear_system( aLinearProblem, aIter );
//
//            // Iterate TryRestartOnFailIt counter
//            tTryRestartOnFailIt = tTryRestartOnFailIt + 1;
//        }
    }

//--------------------------------------------------------------------------------------------------------------------------
    void Nonlinear_Solver_Manager::set_nonlinear_solver_manager_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListNonLinearSolver.insert( "NLA_max_non_lin_solver_restarts" , 0 );
    }

