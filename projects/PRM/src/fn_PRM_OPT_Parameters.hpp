/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_OPT_Parameters.hpp
 *
 */

#ifndef MORIS_FN_PRM_OPT_PARAMETERS_HPP
#define MORIS_FN_PRM_OPT_PARAMETERS_HPP

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {
        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_opt_problem_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "is_optimization_problem", false );        // Whether or not to use OPT
            tParameterList.insert( "workflow", "HMR_XTK" );                   // Workflow to use, HMR_XTK - standard workflow, STK_XTK
            tParameterList.insert( "problem", "user_defined" );               // OPT Problem class type
            tParameterList.insert( "restart_file", "" );                      // Name of restart file
            tParameterList.insert( "finite_difference_type", "none" );        // Type of finite differencing for gradients;
                                                                              // central, forward, backward, or none
            tParameterList.insert( "finite_difference_epsilons", "1E-8" );    // Epsilon(s) to use per ADV for finite differencing
            tParameterList.insert( "library", "" );                           // Path to a shared object file for user-defined functions

            tParameterList.insert( "reinitialize_interface_iter", INT_MAX );          // number of iterations until the interface will be reinitialized
            tParameterList.insert( "first_reinitialize_interface_iter", INT_MAX );    // number of iterations until the interface will be reinitialized

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_opt_interface_manager_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "shared_advs", false );                  // If all of the ADVs are shared between criteria interfaces
            tParameterList.insert( "parallel", false );                     // If to execute criteria evaluations in parallel
            tParameterList.insert( "num_processors_per_interface", "" );    // Matrix of processors to use per interface

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_opt_interface_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "type", "user_defined" );    // OPT Interface class type
            tParameterList.insert( "library", "" );             // Path to a shared object file for user-defined functions

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_gcmma_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "algorithm", "gcmma" );    // Algorithm name, don't change
            tParameterList.insert( "restart_index", 0 );      // Restart iteration index
            tParameterList.insert( "max_its", 100 );          // Maximum number of iterations
            tParameterList.insert( "max_inner_its", 0 );      // Maximum inner iterations per every optimization iteration
            tParameterList.insert( "norm_drop", 1e-4 );       // Relative change in objective convergence criteria
            tParameterList.insert( "asymp_adapt0", 0.5 );     // Initial asymptote adaptation factor
            tParameterList.insert( "asymp_adaptb", 0.7 );     // Shrinking asymptote adaptation factor
            tParameterList.insert( "asymp_adaptc", 1.2 );     // Expanding asymptote adaptation factor
            tParameterList.insert( "step_size", 0.01 );       // GCMMA step size
            tParameterList.insert( "penalty", 100.0 );        // GCMMA constraint penalty
            tParameterList.insert( "version", 1 );            // GCMMA version

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_lbfgs_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "algorithm", "lbfgs" );                   // Algorithm name, don't change
            tParameterList.insert( "restart_index", 0 );                     // Restart iteration index
            tParameterList.insert( "max_its", 10 );                          // maximum optimization iterations allowed
            tParameterList.insert( "num_corr", 5 );                          // number of limited memory corrections used in the BFGS update
            tParameterList.insert( "num_function_evaluations", 1000 );       // maximum number of function calls in inner iteration
            tParameterList.insert( "norm_drop", 1.0e-8 );                    // convergence criterion (converted internally to criterion used by LBFGS algorithm)
            tParameterList.insert( "grad_tol", 0.0 );                        // convergence criterion based on projected gradients
            tParameterList.insert( "internal_lbfgs_print_severity", -1 );    // -1 supresses all printing

            tParameterList.insert( "step_size", "1.0" );                 // default step size, 1.0 is for no modification
            tParameterList.insert( "step_size_index", "0" );             // iterations to start the lbfgs
            tParameterList.insert( "number_inner_iterations", "10" );    // number of inner iterations

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_sqp_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "algorithm", "sqp" );    // Algorithm name, don't change

            tParameterList.insert( "restart_index", 0 );    // Restart iteration index

            // Printing
            tParameterList.insert( "Major print level", 1 );      // Controls the amount of output to print each major iteration.
            tParameterList.insert( "Minor print level", 1 );      // Controls the amount of output during the QP subproblems.
            tParameterList.insert( "Print file", 0 );             // Change to >0 to print output to file.
            tParameterList.insert( "Summary file", 0 );           // Change to >0 to print summary to file.
            tParameterList.insert( "Print frequency", 100 );      // Every nth minor iteration we print output to file.
            tParameterList.insert( "Log frequency", 100 );        // Related to print frequency.
            tParameterList.insert( "Summary frequency", 100 );    // Every nth minor iteration we print output to file.
            tParameterList.insert( "Timing level", 3 );           // prints CPU times

            // SQP method
            tParameterList.insert( "Major iterations limit", 1000 );                // number of allowed major iterations
            tParameterList.insert( "Minor iterations limit", 500 );                 // number of allowed minor iterations
            tParameterList.insert( "Iterations limit", 10000 );                     // number of total minor iterations allowed over all major iterations
            tParameterList.insert( "Major step limit", 2.0 );                       // limits the change in variable during linesearch
            tParameterList.insert( "Superbasics limit", 500 );                      // places a limit on the storage of superbasic variables
            tParameterList.insert( "New superbasics limit", 99 );                   // controls early termination of QPs
            tParameterList.insert( "linesearch_type", "Derivative linesearch" );    // other options are:
                                                                                    // "Nonderivative linesearch"
            tParameterList.insert( "Linesearch tolerance", 0.9 );                   // controls accuracy of linesearch
            tParameterList.insert( "Function precision", 3e-13 );                   // relative accuracy with which nonlinear functions are computed
            tParameterList.insert( "Difference interval", 5.5e-7 );                 // sets the interval for forward differencing
            tParameterList.insert( "Central difference interval", 6.7e-5 );         // sets the interval for central differencing
            tParameterList.insert( "Proximal point method", 1 );                    // satisfies linear constraints near initial guess
            tParameterList.insert( "Violation limit", 10.0 );                       // limit on maximum constraint violation after linesearch
            tParameterList.insert( "Unbounded step size", 1.0e18 );                 // determines unboundedness of linesearch step size
            tParameterList.insert( "Unbounded objective value", 1.0e15 );           // determines unboundedness of objective
            tParameterList.insert( "Infinite bound size", 1.0e+20 );                // any upper bound greater than this value is regarded as infinity

            // QP subproblems
            tParameterList.insert( "Elastic weight", 2.0e+4 );      // weighting of infeasibilities in the objective of the QP subproblem
            tParameterList.insert( "Partial price", 1 );            // reduces the work required for each "pricing" operation
            tParameterList.insert( "Pivot tolerance", 3.7e-11 );    // guards the basis matrix from becoming singular

            // Hessian approximation
            tParameterList.insert( "hessian_type", "Hessian Full memory" );    // Method for storing and updating the Hessian.
                                                                               // Set to "Hessian Limited memory" for variables > 75.
            tParameterList.insert( "Hessian frequency", 999999 );              // for full memory Hessian
            tParameterList.insert( "Hessian updates", 20 );                    // for limited memory Hessian

            // Frequencies
            tParameterList.insert( "Expand frequency", 10000 );        // for anti-cycling procedure
            tParameterList.insert( "Factorization frequency", 50 );    // for basis updates

            // LU options
            tParameterList.insert( "LU factor tolerance", 10.0 );                  // limits size of multipliers in L
            tParameterList.insert( "LU update tolerance", 10.0 );                  // limits size of multipliers in L during updates
            tParameterList.insert( "LU density tolerance", 0.6 );                  // handles sparsity of LU factorization
            tParameterList.insert( "LU singularity tolerance", 2e-6 );             // handles guard against singularity during factorization
            tParameterList.insert( "lu_pivoting_type", "LU Partial Pivoting" );    // Related to LU factorization. Other options are
                                                                                   // "LU Rook Pivoting" - more costly and stable
                                                                                   // "LU Complete Pivoting" - more costly and stable

            // Convergence Tolerances
            tParameterList.insert( "Major optimality tolerance", 1e-6 );     // target accuracy of the dual variable
            tParameterList.insert( "Minor optimality tolerance", 5e-7 );     // Also related to target accuracy of the dual variable
            tParameterList.insert( "Major feasibility tolerance", 1e-6 );    // feasibility with respect to nonlinear constraints
            tParameterList.insert( "Feasibility tolerance", 1e-6 );          // See minor feasibility tolerance (deprecated)
            tParameterList.insert( "Minor feasibility tolerance", 1e-6 );    // feasibility with respect to linear constraints

            // Derivative checking
            tParameterList.insert( "Verify level", 0 );    // Finite difference check on derivatives computed by user-provided routines

            // Scaling
            tParameterList.insert( "Scale option", 1 );         // flag for scaling of constraints and variables
            tParameterList.insert( "Scale tolerance", 0.9 );    // affects effectiveness with which constraint matrix is scaled

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        inline ParameterList
        create_sweep_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "algorithm", "sweep" );                     // Algorithm name, don't change
            tParameterList.insert( "num_evaluations_per_adv", "10" );          // Uniformly sweep each adv with this many evaluation points per adv
                                                                               // Can specify different number per adv, or one value (applies to all advs)
            tParameterList.insert( "custom_adv_evaluations", "" );             // Evaluate with ADVs at specified values, overrides num_evaluations_per_adv
            tParameterList.insert( "include_bounds", true );                   // Allow evaluations with ADVs at the lower and upper bounds
            tParameterList.insert( "evaluate_objectives", true );              // Calculate and output the objective at each point
            tParameterList.insert( "evaluate_constraints", true );             // Calculate and output the constraints at each point
            tParameterList.insert( "evaluate_objective_gradients", true );     // Calculate and output the objective gradients at each point
            tParameterList.insert( "evaluate_constraint_gradients", true );    // Calculate and output the constraint gradients at each point
            tParameterList.insert( "finite_difference_type", "none" );         // Type of finite differencing for gradients;
                                                                               // central, forward, backward, all, or none
            tParameterList.insert( "finite_difference_epsilons", "1E-8" );     // Use finite differencing to obtain gradients with these epsilons
            tParameterList.insert( "finite_difference_adv_indices", "" );      // Indices of ADVs with respect to which sensitivities are computed by FD
            tParameterList.insert( "save", true );                             // Save the sweep evaluations in "hdf5_path"
            tParameterList.insert( "print", false );                           // Print the sweep evaluations to the screen with moris::print
            tParameterList.insert( "hdf5_path", "" );                          // Path and file name for saving if "save" is set to true

            return tParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif    // MORIS_FN_PRM_OPT_PARAMETERS_HPP
