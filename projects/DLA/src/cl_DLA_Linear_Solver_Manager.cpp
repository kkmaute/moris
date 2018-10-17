
#include "cl_DLA_Linear_Solver_Manager.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Enums.hpp"

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

    mLinearSolverList.push_back( tLinSolver1 );
    mLinearSolverList.push_back( tLinSolver2 );
    mLinearSolverList.push_back( tLinSolver3 );
    mLinearSolverList.push_back( tLinSolver4 );
    mLinearSolverList.push_back( tLinSolver5 );
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


