
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

//--------------------------------------------------------------------------------------------------------------------------
    void Nonlinear_Solver_Manager::set_nonlinear_solver_manager_parameters()
    {
        // Maximal number of linear solver restarts on fail
        mParameterListNonLinearSolver.insert( "NLA_max_non_lin_solver_restarts" , 0 );
    }

