/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Newton_Solver_Test.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

#define protected public
#define private public
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Newton_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Solver_Interface_Proxy.hpp"
#undef protected
#undef private

namespace moris
{
    Matrix< DDRMat >
    test_residualt1(
            const moris::sint       aNX,
            const moris::sint       aNY,
            const moris::real       aLambda,
            const Matrix< DDRMat >& tMyValues,
            const moris::uint       aEquationObjectInd )
    {
        Matrix< DDRMat > tResidual( 3, 1, 0.0 );

        std::cout<<tMyValues( 0 )<<'\n';
        std::cout<<tMyValues( 1 )<<'\n';
        std::cout<<tMyValues( 2 )<<'\n';
        tResidual( 0, 0 ) = 2.0*tMyValues( 0 , 0 ) * (std::pow(tMyValues( 1 , 0 ) , 2)) + tMyValues( 1 , 0 ) + 2.0 * tMyValues( 0 ,0 ) + std::cos( tMyValues( 0 ,0 ) ) + 1.0;
        tResidual( 1, 0 ) = 2.0*tMyValues( 1 , 0 ) * (std::pow(tMyValues( 0 , 0 ) , 2)) + tMyValues( 0 , 0 ) + ((3.0 * tMyValues( 1 , 0 )) / 5.0) - (3.0/50.0);
        tResidual( 2, 0 ) = ((2.0*tMyValues( 2 , 0 ))/5.0) - (1.0/10.0);

        std::cout<< "Residual Norm in func  " << std::sqrt(tResidual( 0, 0 )*tResidual( 0, 0 ) + tResidual( 1, 0 )*tResidual( 1, 0 ) + tResidual( 2, 0 )*tResidual( 2, 0 ))<<'\n';

        return tResidual;
    }

    Matrix< DDRMat >
    test_jacobiant1(
            const moris::sint       aNX,
            const moris::sint       aNY,
            const Matrix< DDRMat >& tMyValues,
            const moris::uint       aEquationObjectInd )
    {
        Matrix< DDRMat > tJacobian( 3, 3, 0.0 );

        tJacobian( 0, 0 ) =  2.0 * std::pow( tMyValues( 1, 0 ), 2 ) - std::sin(tMyValues( 0, 0 )) + 2.0;
        tJacobian( 0, 1 ) =  4.0 * tMyValues( 0, 0 ) * tMyValues( 1, 0 ) + 1.0;
        tJacobian( 1, 0 ) =  4.0 * tMyValues( 0, 0 ) * tMyValues( 1, 0 ) + 1.0;
        tJacobian( 1, 1 ) =  2.0 * std::pow( tMyValues( 0, 0 ), 2) + (3.0/5.0);
        tJacobian( 2, 2 ) =  2.0/5.0;

        std::cout<<tJacobian( 0, 0 )<<'\n';
        std::cout<<tJacobian( 0, 1 )<<'\n';
        std::cout<<tJacobian( 0, 2 )<<'\n';
        std::cout<<tJacobian( 1, 0 )<<'\n';
        std::cout<<tJacobian( 1, 1 )<<'\n';
        std::cout<<tJacobian( 1, 2 )<<'\n';
        std::cout<<tJacobian( 2, 0 )<<'\n';
        std::cout<<tJacobian( 2, 1 )<<'\n';
        std::cout<<tJacobian( 2, 2 )<<'\n';

        return tJacobian;
    }

    Vector< Matrix< DDRMat > >
    test_objt1(
            const moris::sint       aNX,
            const moris::sint       aNY,
            const Matrix< DDRMat >& tMyValues,
            const moris::uint       aEquationObjectInd )
    {
        Vector< Matrix< DDRMat > > tObjval;

        Matrix< DDRMat > tObj( 1, 1, 0.0 );

        tObj( 0 ) = tMyValues( 0 , 0 )*(tMyValues( 0 , 0 ) + 1.0) + 0.3*tMyValues( 1 , 0 )*(tMyValues( 1 , 0 )-0.2) + 0.2*tMyValues( 2 , 0 )*(tMyValues( 2 , 0 )-0.5) + std::pow(tMyValues( 0 , 0 ), 2)*std::pow(tMyValues( 1 , 0 ), 2) + tMyValues( 0 , 0 )*tMyValues( 1 , 0 ) + std::sin(tMyValues( 0 , 0 ));

        tObjval.push_back( tObj );

        return tObjval;
    }

    Matrix< DDSMat >
    test_topot1(
            const moris::sint aNX,
            const moris::sint aNY,
            const moris::uint aEquationObjectInd )
    {
        moris::Matrix< moris::DDSMat > tTopo( 3, 1, -1 );

        for ( moris::sint Ik = 0; Ik < 3; Ik++ )
        {
            tTopo( Ik, 0 ) = Ik;
        }

        return tTopo;
    }

    
    

    //------------------------------------------------------------------------------

    namespace NLA
    {
        //#ifdef MORIS_HAVE_PETSC
        TEST_CASE( "Trust Region Solver Test Petsc", "[NLA],[NLA_Test_Trust_Region_Petsc]" )
        {
            if ( moris::par_size() == 1 )
            {
                moris::Solver_Interface* tSolverInput = new NLA_Solver_Interface_Proxy( 3, 1, 1, 1, 0, test_residualt1, test_jacobiant1, test_objt1, test_topot1 );

                dla::Linear_Solver tLinSolManager;
                Nonlinear_Solver   tNonLinSolManager;
                tNonLinSolManager.set_solver_interface( tSolverInput );

                Nonlinear_Problem tNonlinearProblem( tSolverInput, 0, true, sol::MapType::Petsc );

                Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
                tNonlinearSolverParameterList.set( "NLA_max_iter", 1 );
                tNonlinearSolverParameterList.set( "NLA_hard_break", false );
                tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
                tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
                tNonlinearSolverParameterList.set( "NLA_Solver_Implementation", NLA::NonlinearSolverType::TRUST_REGION_SOLVER );
                tNonlinearSolverParameterList.set( "NLA_max_trust_region_iter" , 30);
                Nonlinear_Solver_Factory               tNonlinFactory;
                std::shared_ptr< Nonlinear_Algorithm > tNonlLinSolverAlgorithm =
                        tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

                tNonlLinSolverAlgorithm->set_linear_solver( &tLinSolManager );

                tNonLinSolManager.set_nonlinear_algorithm( tNonlLinSolverAlgorithm, 0 );

                // Set nasty initial guess
                // sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Petsc );
                // sol::Dist_Vector* tInitGuess = tMatFactory.create_vector(tNonlinearProblem.get_solver_interface(),tNonlinearProblem.get_full_vector()->get_map(), 1);
                // sint tLength = tInitGuess->vec_global_length();
                // std::cout<<"Length"<< tLength;
                // // declare IDs and vec
                // Matrix< DDSMat > tIDs(3,1,0);
                // Matrix< DDRMat > tVals(3,1,0.0);
                // tIDs( 1, 0 ) = 1;
                // tIDs( 2, 0 ) = 2;

                // tVals( 0, 0 ) = (2.0);
                // tVals( 1, 0 ) = (7.0);
                // tVals( 2, 0 ) = (-1.0);
                // tInitGuess->sum_into_global_values(tIDs, tVals);

                // // // Replace original initial guess with current initial guess
                // tSolverInput->set_solution_vector( tInitGuess );


                dla::Solver_Factory tSolFactory;
                Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_petsc();
                Parameter_List tPreconditionerParameterList = prm::create_preconditioner_parameter_list(sol::PreconditionerType::PETSC);
                tLinearSolverParameterList.set( "KSPType", std::string( "stcg" ) );
                tLinearSolverParameterList.set( "PCType", std::string( "cholesky" ));
                dla::Preconditioner* tPrec1 = tSolFactory.create_preconditioner( tPreconditionerParameterList );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver1 = tSolFactory.create_solver( tLinearSolverParameterList );
                //std::shared_ptr< dla::Linear_Solver_Algorithm > tLinSolver2 = tSolFactory.create_solver( tLinearSolverParameterList );
                tLinSolver1->set_preconditioner( tPrec1 );

                tLinSolManager.set_linear_algorithm( 0, tLinSolver1 );
                
                //tLinSolManager.set_linear_algorithm( 1, tLinSolver2 );

                tNonLinSolManager.solve( &tNonlinearProblem );

                Matrix< DDSMat > tGlobalIndExtract( 3, 1, 0 );
                tGlobalIndExtract( 2, 0 ) = 2;
                tGlobalIndExtract( 1, 0 ) = 1;
                Vector< Matrix< DDRMat > > tMyValues;

                tNonlLinSolverAlgorithm->extract_my_values( 3, tGlobalIndExtract, 0, tMyValues );

                std::cout<< "Val 1" << tMyValues( 0 )<<"\n";


                CHECK( equal_to( tMyValues( 0 )( 0, 0 ), -0.87264434, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 1, 0 ), 0.43930155, 1.0e+08 ) );
                CHECK( equal_to( tMyValues( 0 )( 2, 0 ), 0.25, 1.0e+08 ) );

                //        delete( tNonlinearProblem );
                delete ( tSolverInput );
            }
        }
//#endif
        

    }
}