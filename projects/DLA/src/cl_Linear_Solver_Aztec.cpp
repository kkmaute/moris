/*
 * cl_Linear_Solver_Aztec.cpp
 *
 *  Created on: May 14, 2018
 *      Author: schmidt
 */

#include "cl_Linear_Solver_Aztec.hpp"

// TPL header files
#include "Epetra_ConfigDefs.h"
#include "AztecOO_ConfigDefs.h"
#include "AztecOO_ConditionNumber.h"

// Ifpack
#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"

//#include "Ifpack_CrsIct.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_LocalFilter.h"

// ML
//#include "ml_include.h"
//#include "ml_epetra_utils.h"
//#include "ml_epetra_preconditioner.h"

// Teuchos
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

using namespace moris;

Linear_Solver_Aztec::Linear_Solver_Aztec( Solver_Input * aInput ) : Linear_Solver_Trilinos ( aInput ),
                                                                    mAztecSolver ( mEpetraProblem ),
                                                                    mMlPrec ( NULL )
{
    this->set_solver_parameters();
}

 //----------------------------------------------------------------------------------------

Linear_Solver_Aztec::~Linear_Solver_Aztec()
{
    delete( mMlPrec );
}

void Linear_Solver_Aztec::set_solver_parameters()
{
    // ASSIGN DEFAULT PARAMETER VALUES
    // AztecOO User Guide, SAND REPORT, SAND2004-3796, https://trilinos.org/oldsite/packages/aztecoo/AztecOOUserGuide.pdf

    // Determine which solver is used
    //options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
    mParameterList.insert( "AZ_solver" ,  INT_MAX );

    // Allowable Aztec solver iterations
    mParameterList.insert( "AZ_max_iter"      , INT_MAX   );

    // Allowable Aztec irelative residual
    mParameterList.insert( "rel_residual" , 1e-08 );

    // set Az_conv -convergence criteria
    // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
    mParameterList.insert( "AZ_conv" ,  INT_MAX );

    // set Az_diagnostic parameters
    // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
    mParameterList.insert( "AZ_diagnostics" ,  INT_MAX );

    // set AZ_output options
    // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
    mParameterList.insert( "AZ_output" ,  INT_MAX );

    // Determines the submatrices factored with the domain decomposition algorithms
    // Option to specify with how many rows from other processors each processor’s local submatrix is augmented.
    mParameterList.insert( "AZ_overlap" , INT_MAX );

    // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
    // Options are AZ_standard, AZ_symmetric
    mParameterList.insert( "AZ_type_overlap" , INT_MAX );

    // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
    // Option to enable (=1) or disable (=0) the Reverse Cuthill–McKee (RCM) algorithm to reorder system equations for smaller bandwidth
    mParameterList.insert( "AZ_reorder" , INT_MAX );

    // Use preconditioner from a previous Iterate() call
    // Option are AZ_calc, AZ_recalc, AZ_reuse
    mParameterList.insert( "AZ_pre_calc" , INT_MAX );

    // Determines  whether  matrix  factorization  information will be kept after this solve
    // for example for preconditioner_recalculation
    mParameterList.insert( "AZ_keep_info" , INT_MAX );

    //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
    // Set AZ_kspace
    // Krylov subspace size for restarted GMRES
    // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption. For very difficult problems, set it equal to the maximum number of iterations.
    mParameterList.insert( "AZ_kspace" ,INT_MAX );

    // Set AZ_orthog
    //AZ_classic or AZ_modified
    mParameterList.insert( "AZ_orthog" , INT_MAX );

    // Set AZ_rthresh
    // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
    mParameterList.insert( "AZ_rthresh" , -1.0 );

    // Set AZ_athresh
    //Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
    mParameterList.insert( "AZ_athresh" , -1.0 );

    //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
    // Determine which preconditioner is used
    // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
    mParameterList.insert( "AZ_precond" ,  INT_MAX );

    // Set preconditioner subdomain solve - direct solve or incomplete
    // Options are AZ_lu, AZ_ilut, AZ_ilu, AZ_rilu, AZ_bilu, AZ_icc
    mParameterList.insert( "AZ_subdomain_solve" ,  INT_MAX );

    // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
    mParameterList.insert( "AZ_poly_ord" ,  INT_MAX );

    // Set drop tolerance - for LU, ILUT
    mParameterList.insert( "AZ_drop" ,  -1.0 );

    // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
    mParameterList.insert( "AZ_graph_fill" ,  INT_MAX );

    // Set ilut fill
    mParameterList.insert( "AZ_ilut_fill" ,  -1.0 );

    // Set Damping or relaxation parameter used for RILU
    mParameterList.insert( "AZ_omega" ,  -1.0 );

//==============================================================================================================
//                             ML Preconditioner settings
//==============================================================================================================

    mParameterList.insert( "ML output"                  ,  -1.0 );
    mParameterList.insert( "print unused"               ,  -1.0 );
    mParameterList.insert( "ML print initial list"      ,  -1.0 );
    mParameterList.insert( "ML print final list"        ,  -1.0 );
    mParameterList.insert( "eigen-analysis: type"       ,  -1.0 );
    mParameterList.insert( "eigen-analysis: iterations" ,  -1.0 );

    mParameterList.insert( "cycle applications"       ,  -1.0 );
    mParameterList.insert( "max levels"               ,  -1.0 );
    mParameterList.insert( "increasing or decreasing" ,  -1.0 );
    mParameterList.insert( "prec type"                ,  -1.0 );

    mParameterList.insert( "aggregation: type"                      ,  -1.0 );
    mParameterList.insert( "aggregation: threshold"                 ,  -1.0 );
    mParameterList.insert( "aggregation: damping factor"            ,  -1.0 );
    mParameterList.insert( "aggregation: smoothing sweeps"          ,  -1.0 );
    mParameterList.insert( "aggregation: use tentative restriction" ,  -1.0 );
    mParameterList.insert( "aggregation: symmetrize"                ,  -1.0 );
    mParameterList.insert( "aggregation: local aggregates"          ,  -1.0 );
    mParameterList.insert( "aggregation: local aggregates"          ,  -1.0 );
    mParameterList.insert( "aggregation: nodes per aggregate"       ,  -1.0 );

    mParameterList.insert( "energy minimization: enable"  ,  -1.0 );
    mParameterList.insert( "energy minimization: type"    ,  -1.0 );
    mParameterList.insert( "energy minimization: droptol" ,  -1.0 );
    mParameterList.insert( "energy minimization: cheap"   ,  -1.0 );

    mParameterList.insert( "smoother: type"                      ,  -1.0 );
    mParameterList.insert( "smoother: sweeps"                    ,  -1.0 );
    mParameterList.insert( "smoother: damping factor"            ,  -1.0 );
    mParameterList.insert( "smoother: pre or post"               ,  -1.0 );
    mParameterList.insert( "smoother: Aztec as solver"           ,  -1.0 );
    mParameterList.insert( "smoother: Aztec options"             ,  -1.0 );
    mParameterList.insert( "smoother: Aztec params"              ,  -1.0 );
    mParameterList.insert( "smoother: ifpack level-of-fill"      ,  -1.0 );
    mParameterList.insert( "smoother: ifpack overlap"            ,  -1.0 );
    mParameterList.insert( "smoother: ifpack absolute threshold" ,  -1.0 );
    mParameterList.insert( "smoother: ifpack relative threshold" ,  -1.0 );

    mParameterList.insert( "coarse: type"            ,  -1.0 );
    mParameterList.insert( "coarse: max size"        ,  -1.0 );
    mParameterList.insert( "coarse: pre or post"     ,  -1.0 );
    mParameterList.insert( "coarse: sweeps"          ,  -1.0 );
    mParameterList.insert( "coarse: damping factor"  ,  -1.0 );
    mParameterList.insert( "coarse: Chebyshev alpha" ,  -1.0 );
    mParameterList.insert( "coarse: max processes"   ,  -1.0 );

    mParameterList.insert( "repartition: enable"      ,  -1.0 );
    mParameterList.insert( "repartition: partitioner" ,  -1.0 );

    mParameterList.insert( "analyze memory"               ,  -1.0 );
    mParameterList.insert( "viz: enable"                  ,  -1.0 );
    mParameterList.insert( "viz: output format"           ,  -1.0 );
    mParameterList.insert( "viz: print starting solution" ,  -1.0 );

    mParameterList.insert( "null space: type"                ,  -1.0 );
    mParameterList.insert( "null space: dimension"           ,  -1.0 );
    mParameterList.insert( "null space: vectors"             ,  -1.0 );
    mParameterList.insert( "null space: vectors to compute"  ,  -1.0 );
    mParameterList.insert( "null space: add default vectors" ,  -1.0 );
}

moris::sint Linear_Solver_Aztec::solve_linear_system()
{
    moris::sint error = 0;
    // Set all Aztec options
    this->set_solver_internal_parameters();

    moris::sint tMaxIt  = mParameterList.get< moris::sint >( "AZ_max_iter" );
    moris::real tRelRes = mParameterList.get< moris::real >( "rel_residual" );

    // M L   Preconditioning
//    if ( mMlPrec != NULL  )
//    {
//        clock_t startPrecTime = clock();
//        {
//            mMlPrec->ComputePreconditioner();
//
//            mAztecSolver.SetPrecOperator ( mMlPrec );
//            //mIsPastFirstSolve = true;
//        }
//        mPreCondTime = moris::real ( clock() - startPrecTime ) / CLOCKS_PER_SEC;
//    }

    // Solve the linear system
    error = mAztecSolver.Iterate( tMaxIt, tRelRes );

    MORIS_ERROR( error==0, "Error in solving linear system with Aztec" );

    // Get linear solution info
    mSolNumIters       = mAztecSolver.NumIters();
    mSolTrueResidual   = mAztecSolver.TrueResidual();
    mSolScaledResidual = mAztecSolver.ScaledResidual();
    mSolTime           = mAztecSolver.SolveTime();

    return error;
}

void Linear_Solver_Aztec::set_solver_internal_parameters()
{

    // Generic iterative solver parameters

    // Solver Type
    if (mParameterList.get< moris::sint >( "AZ_solver" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption( AZ_solver, mParameterList.get< moris::sint >( "AZ_solver" ) );
    }

    // Set AZ_overlap
    // Determines the submatrices factored with the domain decomposition algorithms
    // Option to specify with how many rows from other processors each processor’s local submatrix is augmented.
    if (mParameterList.get< moris::sint >( "AZ_overlap" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_overlap, mParameterList.get< moris::sint >( "AZ_overlap" ) );
    }

    // Set AZ_type_overlap
    // AZ_standard = The resulting value of an unknown is determined by the processor owning that unknown. Information from other processors about that unknown is discarded.
    if (mParameterList.get< moris::sint >( "AZ_type_overlap" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_type_overlap, mParameterList.get< moris::sint >( "AZ_type_overlap" ) );
    }

    // Set AZ_reorder
    // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
    // Option to enable (=1) or disable (=0) the Reverse Cuthill–McKee (RCM) algorithm to reorder system equations for smaller bandwidth
    if (mParameterList.get< moris::sint >( "AZ_reorder" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_reorder, mParameterList.get< moris::sint >( "AZ_reorder" ) );
    }

    // Set AZ_aux_vec
    // AZ_resid = r_tilde is set to the initial residual vector
    //mAztecSolver.SetAztecOption ( AZ_aux_vec, AZ_resid );

    // GMRES specific solver parameters
    // Set AZ_kspace
    // Krylov subspace size for restarted GMRES
    // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption. For very difficult problems, set it equal to the maximum number of iterations.
    if (mParameterList.get< moris::sint >( "AZ_kspace" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_kspace,  mParameterList.get< moris::sint >( "AZ_kspace" ));
    }

    // Set AZ_orthog
    if (mParameterList.get< moris::sint >( "AZ_orthog" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_orthog, mParameterList.get< moris::sint >( "AZ_orthog" ) );
    }

    // Set AZ_rthresh
    // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
    if (mParameterList.get< moris::real >( "AZ_rthresh" ) != -1.0)
    {
        mAztecSolver.SetAztecParam ( AZ_rthresh, mParameterList.get< moris::real >( "AZ_rthresh" ) );
    }

    // Set AZ_athresh
    //Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
    if (mParameterList.get< moris::real >( "AZ_athresh" ) != -1.0)
    {
        mAztecSolver.SetAztecParam ( AZ_athresh, mParameterList.get< moris::real >( "AZ_athresh" ));
    }

    //---------------------------------------------------------------------------------------------------------------
    // Set AZ_conv criteria
    if (mParameterList.get< moris::sint >( "AZ_conv" ) != INT_MAX)
    {
        mAztecSolver.SetAztecParam ( AZ_conv, mParameterList.get< moris::sint >( "AZ_conv" ));
    }

    // Set AZ_diagnostics
    if (mParameterList.get< moris::sint >( "AZ_diagnostics" ) != INT_MAX)
    {
        mAztecSolver.SetAztecParam ( AZ_diagnostics, mParameterList.get< moris::sint >( "AZ_diagnostics" ));
    }

    // Set AZ_output
    if (mParameterList.get< moris::sint >( "AZ_output" ) != INT_MAX)
    {
        mAztecSolver.SetAztecParam ( AZ_output, mParameterList.get< int >( "AZ_output" ));
    }

    // Set if preconditioner is recalculated
    if (mParameterList.get< moris::sint >( "AZ_pre_calc" ) != INT_MAX)
    {
        mAztecSolver.SetAztecParam ( AZ_pre_calc, mParameterList.get< moris::sint >( "AZ_pre_calc" ));
    }

    // Set if preconditioner is recalculated
    if (mParameterList.get< moris::sint >( "AZ_keep_info" ) != INT_MAX)
    {
        mAztecSolver.SetAztecParam ( AZ_keep_info, mParameterList.get< moris::sint >( "AZ_keep_info" ));
    }

    // Determine which preconditioner is used
    if (mParameterList.get< moris::sint >( "AZ_precond" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_precond,  mParameterList.get< moris::sint >( "AZ_precond" ) );
    }

    // Set preconditioner subdomain solve - direct solve or incomplete
    if (mParameterList.get< moris::sint >( "AZ_subdomain_solve" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_subdomain_solve, mParameterList.get< moris::sint >( "AZ_subdomain_solve" ) );
    }

    // Set preconditioner polynomial order
    if (mParameterList.get< moris::sint >( "AZ_poly_ord" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_poly_ord, mParameterList.get< moris::sint >( "AZ_poly_ord" ) );
    }

    // Set drop tolerance - for LU, ILUT
    if (mParameterList.get< moris::real >( "AZ_drop" ) != -1.0 )
    {
        mAztecSolver.SetAztecParam ( AZ_drop, mParameterList.get< moris::real >( "AZ_drop" ) );
    }

    // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
    if (mParameterList.get< moris::sint >( "AZ_graph_fill" ) != INT_MAX)
    {
        mAztecSolver.SetAztecOption ( AZ_graph_fill, mParameterList.get< moris::sint >( "AZ_graph_fill" ) );
    }

    // Set Damping or relaxation parameter used for RILU
    if (mParameterList.get< moris::real >( "AZ_omega" ) != -1.0 )
    {
        mAztecSolver.SetAztecParam ( AZ_omega, mParameterList.get< moris::real >( "AZ_omega" ) );
    }

    // Set ilut fill
    if (mParameterList.get< moris::real >( "AZ_ilut_fill" ) != -1.0 )
    {
        mAztecSolver.SetAztecParam ( AZ_ilut_fill, mParameterList.get< moris::real >( "AZ_ilut_fill" ) );
    }

//==============================================================================================================
//                             ML Preconditioner settings
//==============================================================================================================

    //ML_Epetra::SetDefaults ( mLinearSolverData->mMl->Defaults,mlParams );

    //mlParams.set ( "ML output",mLinearSolverData->mMl->General.MlOutput );
}

