/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Problem.cpp
 *
 */
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_Stopwatch.hpp"    //CHR/src

// detailed logging package
#include "cl_Tracer.hpp"

namespace moris::dla
{
    //----------------------------------------------------------------------------------------
    sol::Dist_Vector*
    Linear_Problem::get_full_solver_LHS()
    {
        // zero out full LHS vec
        mFullVectorLHS->vec_put_scalar( 0.0 );

        mFreeVectorLHS->vec_plus_vec( 1.0, *mPointVectorLHS, 0.0 );

        // Import free LHS to full LHS
        mFullVectorLHS->import_local_to_global( *mFreeVectorLHS );

        return mFullVectorLHS;
    }

    //----------------------------------------------------------------------------------------
    void
    Linear_Problem::set_free_solver_LHS( sol::Dist_Vector* aFullSolVector )
    {
        mFreeVectorLHS->import_local_to_global( *aFullSolVector );

        mPointVectorLHS->vec_plus_vec( 1.0, *mFreeVectorLHS, 0.0 );
    }

    //----------------------------------------------------------------------------------------
    void
    Linear_Problem::assemble_residual()
    {
        Tracer tTracer( "LinearProblem", "AssembleResidual" );

        // Zero out RHS
        mPointVectorRHS->vec_put_scalar( 0.0 );

        // assemble RHS
        mSolverInterface->assemble_RHS( mPointVectorRHS );

        // additional contributions for adjoint sensitivities
        if ( !mSolverInterface->is_forward_analysis() )
        {
            std::cout << "need fix in Linear_Problem::assemble_residual() \n";

            mSolverInterface->assemble_additional_DqDs_RHS_contribution( mPointVectorRHS );

            this->compute_residual_for_adjoint_solve();
        }

        if ( mSolverWarehouse )
        {
            if ( !mSolverWarehouse->get_output_to_matlab_string().empty() )
            {
                // Get the nonlinear system index
                uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

                // construct string for file name
                std::string tResFileName = mSolverWarehouse->get_output_to_matlab_string() + "." + std::to_string( tNonlinearSystemIndex ) + ".res.dat";

                // save to file
                mPointVectorRHS->save_vector_to_matlab_file( tResFileName.c_str() );

                // log that output was successful
                MORIS_LOG_INFO( "Saved Residual to Matlab File: %s ", tResFileName.c_str() );
            }
        }
    }

    //----------------------------------------------------------------------------------------

    void
    Linear_Problem::assemble_staggered_residual_contribution()
    {
        Tracer tTracer( "LinearProblem", "AssembleStaggeredResidualContribution" );

        mSolverInterface->assemble_staggered_RHS_contribution( mPointVectorRHS );
    }

    //----------------------------------------------------------------------------------------

    void
    Linear_Problem::assemble_jacobian()
    {
        Tracer tTracer( "LinearProblem", "AssembleJacobian" );

        mMat->mat_put_scalar( 0.0 );

        // assemble Jacobian
        mSolverInterface->assemble_jacobian( mMat );

        if ( mSolverWarehouse )
        {
            if ( !mSolverWarehouse->get_output_to_matlab_string().empty() )
            {
                // Get the nonlinear system index
                uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

                // construct string for file name
                std::string tJacFileName = mSolverWarehouse->get_output_to_matlab_string() + "." + std::to_string( tNonlinearSystemIndex ) + ".jac.dat";

                // save to file
                mMat->save_matrix_to_matlab_file( tJacFileName.c_str() );

                // log that output was successful
                MORIS_LOG_INFO( "Saved Jacobian to Matlab File: %s ", tJacFileName.c_str() );
            }
        }
    }

    //----------------------------------------------------------------------------------------
    void
    Linear_Problem::assemble_residual_and_jacobian()
    {
        Tracer tTracer( "LinearProblem", "AssembleResidualAndJacobian" );

        mPointVectorRHS->vec_put_scalar( 0.0 );
        mMat->mat_put_scalar( 0.0 );

        mSolverInterface->fill_matrix_and_RHS( mMat, mPointVectorRHS );

        // additional contributions for adjoint sensitivities
        if ( !mSolverInterface->is_forward_analysis() )
        {
            std::cout << "need fix in Linear_Problem::assemble_residual_and_jacobian \n";

            mSolverInterface->assemble_additional_DqDs_RHS_contribution( mPointVectorRHS );

            this->compute_residual_for_adjoint_solve();
        }

        if ( mSolverWarehouse )
        {
            if ( !mSolverWarehouse->get_output_to_matlab_string().empty() )
            {
                // Get the nonlinear system index
                uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

                // get pseudo time index of most recent Newton iteration (last parameter set to true)
                uint tNewtonIndex = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );

                MORIS_ERROR( tNewtonIndex > 0,
                        "Linear_Problem::assemble_residual_and_jacobian - Newton iteration index smaller equal zero" );

                // construct strings for file names
                std::string tJacFileName = mSolverWarehouse->get_output_to_matlab_string()
                                         + "." + std::to_string( tNonlinearSystemIndex )
                                         + "." + std::to_string( tNewtonIndex - 1 )
                                         + ".jac.dat";

                std::string tResFileName = mSolverWarehouse->get_output_to_matlab_string()
                                         + "." + std::to_string( tNonlinearSystemIndex )
                                         + "." + std::to_string( tNewtonIndex - 1 )
                                         + ".res.dat";

                if ( !mSolverInterface->is_forward_analysis() )
                {
                    tJacFileName = "SensitivityAnalysis." + tJacFileName;
                    tResFileName = "SensitivityAnalysis." + tResFileName;
                }

                // output to matlab .dat file
                mMat->save_matrix_to_matlab_file( tJacFileName.c_str() );
                mPointVectorRHS->save_vector_to_matlab_file( tResFileName.c_str() );

                // log that output was successful
                MORIS_LOG_INFO( "Saved Jacobian and Residual to Matlab File: %s %s %s",
                        tJacFileName.c_str(),
                        " and ",
                        tResFileName.c_str() );
            }
        }
    }

    //----------------------------------------------------------------------------------------

    real
    Linear_Problem::compute_static_residual_norm()
    {
        Tracer tTracer( "LinearProblem", "ComputeStaticResidual" );

        // create factor to create distributed vector
        sol::Matrix_Vector_Factory tVecFactory( mTplType );

        // create auxiliary vector for jacobian times sol vec
        sol::Dist_Vector* tDynRes = tVecFactory.create_vector(
                mSolverInterface,
                mPointVectorRHS->get_map(),
                mSolverInterface->get_num_rhs() );

        // initialize dynamic residual
        tDynRes->vec_put_scalar( 0.0 );

        // assemble dynamic residual
        mSolverInterface->assemble_RHS( tDynRes, fem::Time_Continuity_Flag::TIME_CONTINUITY_ONLY );

        MORIS_LOG_INFO( "Norm of dynamic residual: %e", tDynRes->vec_norm2()( 0 ) );
        MORIS_LOG_INFO( "Norm of total residual  : %e", mPointVectorRHS->vec_norm2()( 0 ) );

        // subtract dynamic from total residual to obtain static residual
        tDynRes->vec_plus_vec( -1.0, *mPointVectorRHS, 1.0 );

        // compute norm of static residual
        real tStaticResNorm = tDynRes->vec_norm2()( 0 );

        MORIS_LOG_INFO( "Norm of static residual : %e", tStaticResNorm );
        ;

        // delete auxiliary vector
        delete tDynRes;

        // return norm of static residual
        return tStaticResNorm;
    }

    //----------------------------------------------------------------------------------------

    void
    Linear_Problem::compute_residual_for_adjoint_solve()
    {
        std::cout << "need fix in Linear_Problem::compute_residual_for_adjoint_solve \n";

        // create factor to create distributed vector
        sol::Matrix_Vector_Factory tMatFactory( mTplType );

        // create auxiliary vector for jacobian times sol vec
        sol::Dist_Vector* tMatTimesSolVec = tMatFactory.create_vector(
                mSolverInterface,
                mPointVectorLHS->get_map(),
                mSolverInterface->get_num_rhs() );

        // multiply jacobian with previous solution vector
        mMat->mat_vec_product( *mPointVectorLHS, *tMatTimesSolVec, false );

        // add contribution to RHS
        mPointVectorRHS->vec_plus_vec( 1.0, *tMatTimesSolVec, 1.0 );

        //        std::cout << "Linear_Problem::compute_residual_for_adjoint_solve";
        //        mPointVectorRHS->print();

        // delete auxiliary vector
        delete tMatTimesSolVec;
    }

    //----------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Linear_Problem::compute_residual_of_linear_system()
    {
        sol::Matrix_Vector_Factory tMatFactory( mTplType );

        // get number of RHS
        uint tNumberOfRHS = mSolverInterface->get_num_rhs();

        // create vector for jacobian times sol vec
        sol::Dist_Vector* tResVec = tMatFactory.create_vector(
                mSolverInterface,
                mPointVectorLHS->get_map(),
                tNumberOfRHS );

        // multiply jacobian with previous solution vector
        mMat->mat_vec_product( *mPointVectorLHS, *tResVec, false );

        // add contribution to RHS
        tResVec->vec_plus_vec( -1.0, *mPointVectorRHS, 1.0 );

        // norm of residual and RHS
        Vector< real > tResNorm = tResVec->vec_norm2();
        Vector< real > tRhsNorm = mPointVectorRHS->vec_norm2();

        // allocate vector for relative residuals
        Matrix< DDRMat > tRelativeResidualNorm( tNumberOfRHS, 1 );

        // compute and store relative norms
        for ( uint tRhsIndex = 0; tRhsIndex < tNumberOfRHS; ++tRhsIndex )
        {
            tRelativeResidualNorm( tRhsIndex ) =
                    tResNorm( tRhsIndex ) / std::max( tRhsNorm( tRhsIndex ), MORIS_REAL_EPS );
        }

        // delete temporary vector
        delete tResVec;

        // return vector with relative residuals
        return tRelativeResidualNorm;
    }

    //----------------------------------------------------------------------------------------

    void
    Linear_Problem::construct_rhs_matrix()
    {
        MORIS_ERROR( false, "It is not implemented yet" );
    }

    //----------------------------------------------------------------------------------------

    void
    Linear_Problem::assemble_rhs_matrix()
    {
        this->construct_rhs_matrix();

        if ( mRHSMatType != "" )
        {
            if ( mRHSMatType == "MassMat" )
            {
                mMassMat->mat_put_scalar( 0.0 );

                // assemble jacobian
                mSolverInterface->assemble_jacobian( mMassMat, moris::fem::Time_Continuity_Flag::TIME_CONTINUITY_ONLY );
            }
            else if ( mRHSMatType == "GeomStiffMat" )
            {
                mMassMat->mat_put_scalar( 0.0 );

                // assemble jacobian
                mSolverInterface->assemble_jacobian( mMassMat, moris::fem::Time_Continuity_Flag::GEOMETRIC_STIFFNESS_ONLY );
            }
            else if ( mRHSMatType == "IdentityMat" )
            {
                mMassMat->mat_put_scalar( 0.0 );

                mPointVectorRHS->vec_put_scalar( 1.0 );

                // create identity matrix by replacing diagonal values to 1.0
                mMassMat->replace_diagonal_values( *mPointVectorRHS );
            }
            else
            {
                MORIS_ERROR( false, "RHS Matrix Type not correct" );
            }
        }

        if ( !mSolverWarehouse->get_output_to_matlab_string().empty() )
        {
            // Get the nonlinear system index
            uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver", LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

            // construct string for file name
            std::string tJacFileName = mSolverWarehouse->get_output_to_matlab_string() + "." + std::to_string( tNonlinearSystemIndex ) + ".jac.dat";

            std::string tMassFileName = "Mass_" + tJacFileName;

            mMassMat->save_matrix_to_matlab_file( tMassFileName.c_str() );

            // log that output was successful
            MORIS_LOG_INFO( "Saved Jacobian to Matlab File: %s ", tMassFileName.c_str() );
        }
    }
}    // namespace moris::dla
