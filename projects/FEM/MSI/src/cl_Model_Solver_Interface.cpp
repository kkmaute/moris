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
    void moris::MSI::Model_Solver_Interface::solve_system( )
    {
        // create solver input object
        moris::MSI::MSI_Solver_Interface *  tSolverInput;
        tSolverInput = new moris::MSI::MSI_Solver_Interface( this, &mDofMgn );

        // create solver factory
        moris::Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput );

        tLin->solve_linear_system();

        moris::Mat< moris::real > tSol;
        tLin->get_solution( tSol );

        delete tSolverInput;
    }

void moris::MSI::Model_Solver_Interface::solve_system( moris::Cell< moris::MSI::Equation_Object* > & aListEqnObj )
    {
        // create solver input object
        moris::MSI::MSI_Solver_Interface *  tSolverInput;
        tSolverInput = new moris::MSI::MSI_Solver_Interface( this, &mDofMgn );

        // create solver factory
        moris::Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput );

        //FIXME just temporary for testing. will be deleted
        for ( auto tElement : aListEqnObj )
        {
            tElement->set_solver( tLin );
        }

        tLin->solve_linear_system();

        // print result
        moris::Mat< moris::real > tSol;
        tLin->get_solution( tSol );

        delete tSolverInput;
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

        moris::Mat< moris::real > tSol;
        tLin->get_solution( tSol );

        //=====================================================
//        tLin->import();
//
//        moris::uint tA=3;                       // number of values
//        moris::Mat< moris::sint > tB(tA,1,1);   // global id of dof
//        tB(1,0)=1;
//        tB(1,0)=2;
//        moris::uint tC = 0;                     // offset
//        moris::Mat< moris::real> tD(tA,1);      // resulting moris:Mat
//
//        tLin->extract_my_values(tA, tB, tC, tD);
//
//        std::cout<<tD(0,0)<<std::endl;
//        std::cout<<tD(1,0)<<std::endl;
//        std::cout<<tD(2,0)<<std::endl;

        //==========================================================
        delete tSolverInput;

        aSolution = tSol;
    }
//}
//}
