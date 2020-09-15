/*
 * fn_PRM_SOL_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: schmidt
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_SOL_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_SOL_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"

#include "cl_Param_List.hpp"

#include "cl_SOL_Enums.hpp"
#include "cl_NLA_Nonlinear_Solver_Enums.hpp"       //CON/src
#include "cl_TSA_Time_Solver_Enums.hpp"            //CON/src


namespace moris
{
    namespace prm
    {

        //------------------------------------------------------------------------------

        ParameterList create_solver_warehouse_parameterlist( )
        {
            ParameterList tSolverWarehouseList;

            tSolverWarehouseList.insert( "SOL_TPL_Type" , static_cast< uint >( sol::MapType::Epetra ) );

            return tSolverWarehouseList;
        }

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        ParameterList create_linear_algorithm_parameter_list_aztec( )
        {
            ParameterList tLinAlgorithmParameterList;

            enum moris::sol::SolverType tType = moris::sol::SolverType::AZTEC_IMPL;

            tLinAlgorithmParameterList.insert( "Solver_Implementation" , static_cast< uint >( tType ) );

            // ASSIGN DEFAULT PARAMETER VALUES
            // AztecOO User Guide, SAND REPORT, SAND2004-3796, https://trilinos.org/oldsite/packages/aztecoo/AztecOOUserGuide.pdf

            // Determine which solver is used
            //options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
            tLinAlgorithmParameterList.insert( "AZ_solver" ,  INT_MAX );

            // Allowable Aztec solver iterations
            tLinAlgorithmParameterList.insert( "AZ_max_iter", INT_MAX   );

            // Allowable Aztec iterative residual
            tLinAlgorithmParameterList.insert( "rel_residual" , 1e-08 );

            // set Az_conv -convergence criteria
            // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
            tLinAlgorithmParameterList.insert( "AZ_conv" ,  INT_MAX );

            // set Az_diagnostic parameters
            // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
            tLinAlgorithmParameterList.insert( "AZ_diagnostics" ,  INT_MAX );

            // set AZ_output options
            // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
            tLinAlgorithmParameterList.insert( "AZ_output" ,  INT_MAX );

            // Determines the submatrices factored with the domain decomposition algorithms
            // Option to specify with how many rows from other processors each processor’s local submatrix is augmented.
            tLinAlgorithmParameterList.insert( "AZ_overlap" , INT_MAX );

            // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
            // Options are AZ_standard, AZ_symmetric
            tLinAlgorithmParameterList.insert( "AZ_type_overlap" , INT_MAX );

            // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
            // Option to enable (=1) or disable (=0) the Reverse Cuthill–McKee (RCM) algorithm to reorder system equations for smaller bandwidth
            tLinAlgorithmParameterList.insert( "AZ_reorder" , INT_MAX );

            // Use preconditioner from a previous Iterate() call
            // Option are AZ_calc, AZ_recalc, AZ_reuse
            tLinAlgorithmParameterList.insert( "AZ_pre_calc" , INT_MAX );

            // Determines  whether  matrix  factorization  information will be kept after this solve
            // for example for preconditioner_recalculation
            tLinAlgorithmParameterList.insert( "AZ_keep_info" , INT_MAX );

            //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
            // Set AZ_kspace
            // Krylov subspace size for restarted GMRES
            // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption. For very difficult problems, set it equal to the maximum number of iterations.
            tLinAlgorithmParameterList.insert( "AZ_kspace" ,INT_MAX );

            // Set AZ_orthog
            //AZ_classic or AZ_modified
            tLinAlgorithmParameterList.insert( "AZ_orthog" , INT_MAX );

            // Set AZ_rthresh
            // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
            tLinAlgorithmParameterList.insert( "AZ_rthresh" , -1.0 );

            // Set AZ_athresh
            //Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
            tLinAlgorithmParameterList.insert( "AZ_athresh" , -1.0 );

            //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
            // Determine which preconditioner is used
            // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
            tLinAlgorithmParameterList.insert( "AZ_precond" ,  INT_MAX );

            // Set preconditioner subdomain solve - direct solve or incomplete
            // Options are AZ_lu, AZ_ilut, , AZ_rilu, AZ_bilu, AZ_icc
            tLinAlgorithmParameterList.insert( "AZ_subdomain_solve" ,  INT_MAX );

            // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
            tLinAlgorithmParameterList.insert( "AZ_poly_ord" ,  INT_MAX );

            // Set drop tolerance - for LU, ILUT
            tLinAlgorithmParameterList.insert( "AZ_drop" ,  -1.0 );

            // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
            tLinAlgorithmParameterList.insert( "AZ_graph_fill" ,  INT_MAX );

            // Set ilut fill
            tLinAlgorithmParameterList.insert( "AZ_ilut_fill" ,  -1.0 );

            // Set Damping or relaxation parameter used for RILU
            tLinAlgorithmParameterList.insert( "AZ_omega" ,  -1.0 );

            // Set Damping or relaxation parameter used for RILU
            tLinAlgorithmParameterList.insert( "Use_ML_Prec" ,  false );

            // Set Damping or relaxation parameter used for RILU
            tLinAlgorithmParameterList.insert( "ML_reuse" ,  false );

            return tLinAlgorithmParameterList;
        }

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        ParameterList create_linear_algorithm_parameter_list_amesos( )
        {
            ParameterList tLinAlgorithmParameterList;

            enum moris::sol::SolverType tType = moris::sol::SolverType::AMESOS_IMPL;

            tLinAlgorithmParameterList.insert( "Solver_Implementation" , static_cast< uint >( tType ) );

            return tLinAlgorithmParameterList;
        }

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        ParameterList create_linear_algorithm_parameter_list_belos( )
        {
            ParameterList tLinAlgorithmParameterList;

            enum moris::sol::SolverType tType = moris::sol::SolverType::BELOS_IMPL;

            tLinAlgorithmParameterList.insert( "Solver_Implementation" , static_cast< uint >( tType ) );

            // ASSIGN DEFAULT PARAMETER VALUES
            // https://docs.trilinos.org/dev/packages/belos/doc/html/classBelos_1_1SolverFactory.html#ad86e61fb180a73c6dd5dbf458df6a86f

            // Determine which solver is used by string
            // options are: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
            //              Pseudoblock TFQMR, Seed GMRES, Seed CG
            tLinAlgorithmParameterList.insert( "Solver Type" ,  "GMRES" );

            tLinAlgorithmParameterList.insert( "Verbosity" ,  INT_MAX );

            tLinAlgorithmParameterList.insert( "Num Blocks", INT_MAX   );

            tLinAlgorithmParameterList.insert( "Block Size", INT_MAX   );

            // Allowable Belos solver iterations
            tLinAlgorithmParameterList.insert( "Maximum Iterations" , INT_MAX );

            // Allowable Belos solver iterations
            tLinAlgorithmParameterList.insert( "Maximum Restarts" , INT_MAX );

            // set convergence criteria
            tLinAlgorithmParameterList.insert( "Convergence Tolerance" ,  1e-08 );

            return tLinAlgorithmParameterList;
        }

        //------------------------------------------------------------------------------
        // creates a parameter list with default inputs
        ParameterList create_linear_algorithm_parameter_list_petsc( )
        {
            ParameterList tLinAlgorithmParameterList;

            enum moris::sol::SolverType tType = moris::sol::SolverType::PETSC;

            tLinAlgorithmParameterList.insert( "Solver_Implementation" , static_cast< uint >( tType ) );

            // Set KSP type
            tLinAlgorithmParameterList.insert( "KSPType", std::string( "gmres" ) );

            // Set default preconditioner
            tLinAlgorithmParameterList.insert( "PCType", std::string( "ilu" ) );

            // Sets maximal iters for KSP
            tLinAlgorithmParameterList.insert( "KSPMaxits", 1000 );

            // Sets KSP gmres restart
            tLinAlgorithmParameterList.insert( "KSPMGMRESRestart", 500 );

            // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES
            tLinAlgorithmParameterList.insert( "KSPGMRESHapTol", 1e-10 );

            // Sets tolerance for KSP
            tLinAlgorithmParameterList.insert( "KSPTol", 1e-10 );

            // Sets the number of levels of fill to use for ILU
            tLinAlgorithmParameterList.insert( "ILUFill", 0 );

            // Sets drop tolerance for ilu
            tLinAlgorithmParameterList.insert( "ILUTol", 1e-6 );

            // Set multigrid levels
            tLinAlgorithmParameterList.insert( "MultigridLevels", 3 );

            // Schwarz preconditioner volume fraction threshold
            tLinAlgorithmParameterList.insert( "MG_use_schwarz_smoother", false );

            // Schwarz smoothing iterations
            tLinAlgorithmParameterList.insert( "MG_schwarz_smoothing_iters", 1 );

            // Schwarz preconditioner volume fraction threshold
            tLinAlgorithmParameterList.insert( "ASM_volume_fraction_threshold", 0.1 );


            return tLinAlgorithmParameterList;
        }

        //------------------------------------------------------------------------------

        ParameterList create_linear_solver_parameter_list()
        {
            ParameterList tLinSolverParameterList;

            tLinSolverParameterList.insert( "DLA_Linear_solver_algorithms" , std::string("0") );

            // Maximal number of linear solver restarts on fail
            tLinSolverParameterList.insert( "DLA_max_lin_solver_restarts" , 0 );

            // Maximal number of linear solver restarts on fail
            tLinSolverParameterList.insert( "DLA_hard_break" , true );

            // Determines if linear solve should restart on fail
            tLinSolverParameterList.insert( "DLA_rebuild_lin_solver_on_fail" , false );

            return tLinSolverParameterList;
        }

        //------------------------------------------------------------------------------

        ParameterList create_nonlinear_algorithm_parameter_list()
        {
            ParameterList tNonLinAlgorithmParameterList;

            enum moris::NLA::NonlinearSolverType NonlinearSolverType = moris::NLA::NonlinearSolverType::NEWTON_SOLVER;

            tNonLinAlgorithmParameterList.insert( "NLA_Solver_Implementation" , static_cast< uint >( NonlinearSolverType ) );

            tNonLinAlgorithmParameterList.insert( "NLA_Linear_solver" , 0 );

            tNonLinAlgorithmParameterList.insert( "NLA_linear_solver_for_adjoint_solve" , -1 );

            // Allowable Newton solver iterations
            tNonLinAlgorithmParameterList.insert( "NLA_max_iter", 10 );

            // Allowable Newton solver iterations
            tNonLinAlgorithmParameterList.insert( "NLA_restart", 0 );

            // Desired total residual norm drop
            tNonLinAlgorithmParameterList.insert( "NLA_rel_res_norm_drop" , 1e-08 );

            // Desired total residual norm
            tNonLinAlgorithmParameterList.insert( "NLA_tot_res_norm" , 1e-12 );

            // Maximal residual norm
            tNonLinAlgorithmParameterList.insert( "NLA_max_rel_res_norm" , 1e12 );

            // Maximal number of linear solver restarts on fail
            tNonLinAlgorithmParameterList.insert( "NLA_max_lin_solver_restarts" , 0 );

            // Maximal number of linear solver restarts on fail
            tNonLinAlgorithmParameterList.insert( "NLA_relaxation_parameter" , 1.0 );

            // Maximal number of linear solver restarts on fail
            tNonLinAlgorithmParameterList.insert( "NLA_hard_break" , false );

            // Determines if linear solve should restart on fail
            tNonLinAlgorithmParameterList.insert( "NLA_rebuild_lin_solv_on_fail" , false );

            // Determines if jacobian is rebuild for every nonlinear iteration
            tNonLinAlgorithmParameterList.insert( "NLA_rebuild_jacobian" , true );

            // Determines if linear solve should restart on fail
            tNonLinAlgorithmParameterList.insert( "NLA_combined_res_jac_assembly" , true );

            // Determines if Newton should restart on fail
            tNonLinAlgorithmParameterList.insert( "NLA_rebuild_nonlin_solv_on_fail" , false );

            // Specifying the number of Newton retries
            tNonLinAlgorithmParameterList.insert( "NLA_num_nonlin_rebuild_iterations" , 1 );

            // Determines relaxation multiplier
            tNonLinAlgorithmParameterList.insert( "NLA_relaxation_multiplier_on_fail" , 0.5 );

            // Determines Newton maxits multiplier
            tNonLinAlgorithmParameterList.insert( "NLA_maxits_multiplier_on_fail" , 2 );

            return tNonLinAlgorithmParameterList;
        }

        //------------------------------------------------------------------------------

        ParameterList create_nonlinear_solver_parameter_list()
        {
            ParameterList tNonLinSolverParameterList;

            enum moris::NLA::NonlinearSolverType NonlinearSolverType = moris::NLA::NonlinearSolverType::NEWTON_SOLVER;

            tNonLinSolverParameterList.insert( "NLA_Solver_Implementation" , static_cast< uint >( NonlinearSolverType ) );

            tNonLinSolverParameterList.insert( "NLA_DofTypes" , std::string("UNDEFINED") );

            tNonLinSolverParameterList.insert( "NLA_Sub_Nonlinear_Solver" , std::string("") );

            tNonLinSolverParameterList.insert( "NLA_Nonlinear_solver_algorithms" , std::string("0") );

            // Maximal number of linear solver restarts on fail
            tNonLinSolverParameterList.insert( "NLA_max_non_lin_solver_restarts" , 0 );

            return tNonLinSolverParameterList;
        }

        //------------------------------------------------------------------------------

        ParameterList create_time_solver_algorithm_parameter_list()
        {
            ParameterList tTimeAlgorithmParameterList;

            enum moris::tsa::TimeSolverType tType = tsa::TimeSolverType::MONOLITHIC;

            tTimeAlgorithmParameterList.insert( "TSA_Solver_Implementation" , static_cast< uint >( tType ) );

            tTimeAlgorithmParameterList.insert( "TSA_Nonlinear_solver" , 0 );

            tTimeAlgorithmParameterList.insert( "TSA_nonlinear_solver_for_adjoint_solve" , -1 );

            // Number of time steps
            tTimeAlgorithmParameterList.insert( "TSA_Num_Time_Steps", 1 );

            // Time Frame
            tTimeAlgorithmParameterList.insert( "TSA_Time_Frame", 1.0 );

            return tTimeAlgorithmParameterList;
        }

        //------------------------------------------------------------------------------

        ParameterList create_time_solver_parameter_list()
        {
            ParameterList tTimeParameterList;

            tTimeParameterList.insert( "TSA_TPL_Type" , static_cast< uint >( sol::MapType::Epetra ) );

            tTimeParameterList.insert( "TSA_Solver_algorithms" , std::string("0") );

            tTimeParameterList.insert( "TSA_DofTypes" , std::string("UNDEFINED") );

            // Maximal number of linear solver restarts on fail
            tTimeParameterList.insert( "TSA_Max_Time_Solver_Restarts" , 0 );

            tTimeParameterList.insert( "TSA_Output_Indices" , std::string("") );

            tTimeParameterList.insert( "TSA_Output_Crteria" , std::string("") );

            tTimeParameterList.insert( "TSA_Initialize_Sol_Vec" , std::string("") );

            tTimeParameterList.insert( "TSA_time_level_per_type" , std::string("") );

            return tTimeParameterList;
        }

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        ParameterList create_linear_algorithm_parameter_list( const enum moris::sol::SolverType aSolverType,
                const uint                        aIndex = 0 )
        {
            ParameterList tParameterList;

            switch( aSolverType )
            {
                case ( sol::SolverType::AZTEC_IMPL ):
                                return create_linear_algorithm_parameter_list_aztec( );
                break;
                case ( sol::SolverType::AMESOS_IMPL ):
                                return create_linear_algorithm_parameter_list_amesos();
                break;
                case ( sol::SolverType::BELOS_IMPL ):
                                return create_linear_algorithm_parameter_list_belos();
                break;
                case ( sol::SolverType::PETSC ):
                                return create_linear_algorithm_parameter_list_petsc( );
                break;
                default:
                    MORIS_ERROR( false, "Parameter list for this solver not implemented yet" );
                    break;
            }
            return tParameterList;
        }

        //------------------------------------------------------------------------------

        //    // creates a parameter list with default inputs
        //    ParameterList create_nonlinear_algorithm_parameter_list( const enum moris::NLA::NonlinearSolverType aSolverType,
        //                                                             const uint                        aIndex = 0 )
        //    {
        //        ParameterList tParameterList;
        //
        //        switch( aSolverType )
        //        {
        //        case ( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ):
        //            return create_linear_algorithm_parameter_list_aztec( );
        //            break;
        //        case ( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ):
        //    		MORIS_ERROR( false, "No implemented yet" );
        //            break;
        //
        //        default:
        //            MORIS_ERROR( false, "Parameterlist for this solver not implemented yet" );
        //            break;
        //        }
        //
        //    return tParameterList;
        //}

        //------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_FN_PRM_SOL_PARAMETERS_HPP_ */
