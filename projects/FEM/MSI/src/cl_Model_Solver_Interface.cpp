/*
 * cl_Dof_Manager.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */

#include "cl_Model_Solver_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_Solver_Factory.hpp" // DLA/src
#include "cl_Solver_Input.hpp"

//using namespace moris;
//using namespace MSI;
//namespace moris
//{
//    namespace MSI
//    {
//    Model_Solver_Interface::Model_Solver_Interface( const moris::uint aNumEquationObj,
//                                                          moris::Cell < Equation_Object* > & aListEqnObj ) : mNumEquationObjects( aNumEquationObj ),
//                                                                                                             mEquationObjectList( aListEqnObj ),
//                                                                                                             mDofMgn( aNumEquationObj, aListEqnObj )
//{
//    //Dof_Manager tDofMgn ( aNumEquationObj, aListEqnObj );
//}
//
//    Model_Solver_Interface::~Model_Solver_Interface()
//    {
//
//    }

    void moris::MSI::Model_Solver_Interface::solve_system()
    {
        // create solver input object
        moris::MSI::MSI_Solver_Interface *  tSolverInput;
        tSolverInput = new moris::MSI::MSI_Solver_Interface( this, &mDofMgn );

        // create solver factory
        moris::Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput );

        tLin->solve_linear_system();

        moris::Mat< moris::real > tSol ( 2, 1 );
        tLin->get_solution( tSol );

        delete tSolverInput;

        std::cout<<tSol(0,0)<<std::endl;
        std::cout<<tSol(1,0)<<std::endl;
    }

    void moris::MSI::Model_Solver_Interface::solve_system( moris::Mat< moris::real > & aSolution )
    {
        // create solver input object
        moris::MSI::MSI_Solver_Interface *  tSolverInput;
        tSolverInput = new moris::MSI::MSI_Solver_Interface( this, &mDofMgn );

        // create solver factory
        moris::Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput );

        tLin->solve_linear_system();

        moris::Mat< moris::real > tSol ( 2, 1 );
        tLin->get_solution( tSol );

        delete tSolverInput;

        aSolution = tSol;
    }
//}
//}
