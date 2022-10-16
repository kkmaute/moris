/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * opt_test.cpp
 *
 */

#include <catch.hpp>

#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_OPT_Manager.hpp"
#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Problem_User_Defined.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "fn_OPT_Rosenbrock.hpp"
#include "fn_OPT_Test_Interface.hpp"

namespace moris
{
    namespace opt
    {
        TEST_CASE( "[optimization]" )
        {

            // ---------------------------------------------------------------------------------------------------------

            SECTION( "GCMMA" )
            {
                // Set up default parameter lists
                ParameterList tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
                ParameterList tAlgorithmParameterList = moris::prm::create_gcmma_parameter_list();

                tAlgorithmParameterList.set( "version", 1 );
                tAlgorithmParameterList.set( "norm_drop", 2.5e-5 );
                tAlgorithmParameterList.set( "asymp_adaptc", 1.1 );

                // Create interface
                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        &get_dcriteria_dadv_rosenbrock );

                // Create Problem
                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        &compute_dobjective_dadv_rosenbrock,
                        &compute_dobjective_dcriteria_rosenbrock,
                        &compute_dconstraint_dadv_rosenbrock,
                        &compute_dconstraint_dcriteria_rosenbrock );

                // Create manager
                moris::opt::Manager tManager( { tAlgorithmParameterList }, tProblem );

                // Solve optimization problem
                tManager.perform();

                // Check Solution
                if ( par_rank() == 0 )
                {
                    REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                    REQUIRE( norm( tManager.get_advs() - 1.0 ) < 1E-4 );             // check value of design variables
                }
            }

            // ---------------------------------------------------------------------------------------------------------

            SECTION( "GCMMA 07" )
            {
                // Set up default parameter lists
                ParameterList tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
                ParameterList tAlgorithmParameterList = moris::prm::create_gcmma_parameter_list();

                tAlgorithmParameterList.set( "version", 2 );
                tAlgorithmParameterList.set( "norm_drop", 1e-6 );
                tAlgorithmParameterList.set( "max_inner_its", 3 );
                tAlgorithmParameterList.set( "max_its", 11 );

                // Create interface
                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        &get_dcriteria_dadv_rosenbrock );

                // Create Problem
                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        &compute_dobjective_dadv_rosenbrock,
                        &compute_dobjective_dcriteria_rosenbrock,
                        &compute_dconstraint_dadv_rosenbrock,
                        &compute_dconstraint_dcriteria_rosenbrock );

                // Create manager
                moris::opt::Manager tManager( { tAlgorithmParameterList }, tProblem );

                // Solve optimization problem
                tManager.perform();

                // Check Solution
                if ( par_rank() == 0 )
                {
                    REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                    REQUIRE( norm( tManager.get_advs() - 1.0 ) < 1E-4 );             // check value of design variables
                }
            }

            // ---------------------------------------------------------------------------------------------------------

            SECTION( "SQP" )
            {
                // Set up default parameter lists
                ParameterList tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
                ParameterList tAlgorithmParameterList = moris::prm::create_sqp_parameter_list();

                // Create interface
                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        &get_dcriteria_dadv_rosenbrock );

                // Create Problem
                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        &compute_dobjective_dadv_rosenbrock,
                        &compute_dobjective_dcriteria_rosenbrock,
                        &compute_dconstraint_dadv_rosenbrock,
                        &compute_dconstraint_dcriteria_rosenbrock );

                // Create manager
                moris::opt::Manager tManager( { tAlgorithmParameterList }, tProblem );

                // Solve optimization problem
                tManager.perform();

                // Check Solution
                if ( par_rank() == 0 )
                {
                    REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 5E-4 );    // check value of objective
                    REQUIRE( norm( tManager.get_advs() - 1.0 ) < 1E-3 );             // check value of design variables
                }
            }

            // ---------------------------------------------------------------------------------------------------------

            SECTION( "LBFGS" )
            {
                // This optimization problem does not use the constraints!
                // Set up default parameter lists
                ParameterList tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
                ParameterList tAlgorithmParameterList = moris::prm::create_lbfgs_parameter_list();

                // Create interface
                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        &get_dcriteria_dadv_rosenbrock );

                // Create Problem
                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        &compute_dobjective_dadv_rosenbrock,
                        &compute_dobjective_dcriteria_rosenbrock,
                        &compute_dconstraint_dadv_rosenbrock,
                        &compute_dconstraint_dcriteria_rosenbrock );

                // Create manager
                moris::opt::Manager tManager( { tAlgorithmParameterList }, tProblem );

                // Solve optimization problem
                tManager.perform();

                // Check Solution
                if ( par_rank() == 0 )
                {
                    REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                    REQUIRE( norm( tManager.get_advs() - 1.0 ) < 1E-4 );             // check value of design variables
                }
            }

            // ---------------------------------------------------------------------------------------------------------

            SECTION( "Sweep" )
            {
                // MORIS output
                std::string tMorisOutput = std::getenv( "MORISOUTPUT" );

                MORIS_ERROR( tMorisOutput.size() > 0,
                        "Environment variable MORISOUTPUT not set." );

                // Set up default parameter lists
                ParameterList tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
                ParameterList tAlgorithmParameterList = moris::prm::create_sweep_parameter_list();

                // Change parameters
                tAlgorithmParameterList.set( "hdf5_path", tMorisOutput + "sweep.hdf5" );
                tAlgorithmParameterList.set( "num_evaluations_per_adv", "3, 2" );
                tAlgorithmParameterList.set( "finite_difference_type", "all" );
                tAlgorithmParameterList.set( "finite_difference_epsilons", "0.001, 0.01; 0.00001, 0.00001" );
                tAlgorithmParameterList.set( "finite_difference_adv_indices", "0,1" );

                // Create interface
                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        &get_dcriteria_dadv_rosenbrock );

                // Create Problem
                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        &compute_dobjective_dadv_rosenbrock,
                        &compute_dobjective_dcriteria_rosenbrock,
                        &compute_dconstraint_dadv_rosenbrock,
                        &compute_dconstraint_dcriteria_rosenbrock );

                // Create manager
                moris::opt::Manager tManager( { tAlgorithmParameterList }, tProblem );

                // Solve optimization problem
                tManager.perform();

                // Sweep without sensitivities
                tAlgorithmParameterList.set( "evaluate_objective_gradients", false );
                tAlgorithmParameterList.set( "evaluate_constraint_gradients", false );

                // Create interface
                tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        nullptr );

                // Create Problem
                tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        nullptr,
                        nullptr,
                        nullptr,
                        nullptr );

                // Create manager
                tManager = moris::opt::Manager( { tAlgorithmParameterList }, tProblem );

                // Solve optimization problem
                tManager.perform();
            }

            // ---------------------------------------------------------------------------------------------------------

            SECTION( "Interface" )
            {
                if ( par_size() == 4 or par_size() == 8 )
                {
                    // Create interfaces
                    Cell< std::shared_ptr< Criteria_Interface > > tInterfaces( 4 );
                    tInterfaces( 0 ) = std::make_shared< Interface_User_Defined >(
                            &initialize_test_1,
                            &get_criteria_test,
                            &get_dcriteria_dadv_test );
                    tInterfaces( 1 ) = std::make_shared< Interface_User_Defined >(
                            &initialize_test_2,
                            &get_criteria_test,
                            &get_dcriteria_dadv_test );
                    tInterfaces( 2 ) = std::make_shared< Interface_User_Defined >(
                            &initialize_test_3,
                            &get_criteria_test,
                            &get_dcriteria_dadv_test );
                    tInterfaces( 3 ) = std::make_shared< Interface_User_Defined >(
                            &initialize_test_4,
                            &get_criteria_test,
                            &get_dcriteria_dadv_test );

                    // Interface manager parameter list
                    ParameterList tInterfaceManagerParameterList = moris::prm::create_opt_interface_manager_parameter_list();

                    // Set manager parameters
                    tInterfaceManagerParameterList.set( "parallel", true );
                    tInterfaceManagerParameterList.set( "shared_advs", false );
                    if ( par_size() == 4 )
                    {
                        tInterfaceManagerParameterList.set( "num_processors_per_interface", "1, 1, 1, 1" );
                    }
                    else
                    {
                        tInterfaceManagerParameterList.set( "num_processors_per_interface", "1, 2, 3, 2" );
                    }

                    // Create manager
                    std::shared_ptr< Criteria_Interface > tInterface = create_interface(
                            { tInterfaceManagerParameterList },
                            tInterfaces );

                    // Test manager in parallel
                    Matrix< DDRMat > tADVs;
                    Matrix< DDRMat > tLowerBounds;
                    Matrix< DDRMat > tUpperBounds;
                    Matrix< IdMat >  tDummy;
                    tInterface->initialize( tADVs, tLowerBounds, tUpperBounds, tDummy );
                    for ( uint tADVIndex = 0; tADVIndex < 8; tADVIndex++ )
                    {
                        REQUIRE( tADVs( tADVIndex ) == tADVIndex + 1 );
                    }
                    tInterface->get_criteria( tADVs );
                    Matrix< DDRMat > tCriteriaGradients = tInterface->get_dcriteria_dadv();
                    for ( uint tADVIndex = 0; tADVIndex < 8; tADVIndex++ )
                    {
                        REQUIRE( tCriteriaGradients( tADVIndex, tADVIndex ) == tADVIndex + 1 );
                    }
                }
            }

            // ---------------------------------------------------------------------------------------------------------
        }
    }    // namespace opt
}    // namespace moris
