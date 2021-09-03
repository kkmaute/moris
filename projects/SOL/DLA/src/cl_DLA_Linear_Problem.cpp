/*
 * cl_DLA_Linear_Problem.cpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_Stopwatch.hpp" //CHR/src

// detailed logging package
#include "cl_Tracer.hpp"

namespace moris
{
    namespace dla
    {
        sol::Dist_Vector * Linear_Problem::get_full_solver_LHS()
        {
            // zero out full LHS vec
            mFullVectorLHS->vec_put_scalar( 0.0 );

            mFreeVectorLHS->vec_plus_vec( 1.0, *mPointVectorLHS, 0.0 );

            // Import free LHS to full LHS
            mFullVectorLHS->import_local_to_global( *mFreeVectorLHS );

            return mFullVectorLHS;
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::set_free_solver_LHS( sol::Dist_Vector * aFullSolVector)
        {
            mFreeVectorLHS->import_local_to_global( *aFullSolVector );

            mPointVectorLHS->vec_plus_vec( 1.0, *mFreeVectorLHS, 0.0 );
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_residual_and_jacobian( sol::Dist_Vector * aFullSolutionVector )
        {
            // zero out RHS
            mVectorRHS->vec_put_scalar( 0.0 );
            mPointVectorRHS->vec_put_scalar( 0.0 );

            // zero out matrix
            mMat->mat_put_scalar( 0.0 );

            //
            mSolverInterface->fill_matrix_and_RHS( mMat, mPointVectorRHS, aFullSolutionVector );
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_residual()
        {
            Tracer tTracer( "LinearProblem", "AssembleResidual" );

            // Zero out RHS
            mVectorRHS->vec_put_scalar( 0.0 );
            mPointVectorRHS->vec_put_scalar( 0.0 );

            // start timer
            tic tTimer;

            // assemble RHS
            mSolverInterface->assemble_RHS( mPointVectorRHS );

            mSolverInterface->assemble_additional_DqDs_RHS_contribution( mPointVectorRHS );

            if( !mSolverInterface->get_is_forward_analysis() )
            {
                this->compute_residual_for_adjoint_solve();
            }

            if( mSolverWarehouse )
            {
                if( !mSolverWarehouse->get_output_to_matlab_string().empty() )
                {
                    // Get the nonlinear system index
                    uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver" , LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

                    // construct string for file name
                    std::string tResFileName = mSolverWarehouse->get_output_to_matlab_string() + 
                            "." + std::to_string(tNonlinearSystemIndex) + ".res.dat";
                            
                    // save to file
                    mPointVectorRHS->save_vector_to_matlab_file( tResFileName.c_str() );

                    // log that output was successful
                    MORIS_LOG_INFO( "Saved Residual to Matlab File: ", tResFileName.c_str() );
                }
            }

            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
            MORIS_LOG_INFO( " Assembly of residual on processor %u took %5.3f seconds.", ( uint ) par_rank(), ( double ) tElapsedTime / 1000);

            //std::cout<<"first RHS"<<std::endl;
            //mVectorRHS->print();
        }

        //----------------------------------------------------------------------------------------

        void Linear_Problem::assemble_staggered_residual_contribution()
        {
            mSolverInterface->assemble_staggerd_RHS_contribution( mPointVectorRHS );
            //std::cout<<"second RHS"<<std::endl;
            //mVectorRHS->print();
        }

        //----------------------------------------------------------------------------------------

        void Linear_Problem::assemble_jacobian()
        {
            Tracer tTracer( "LinearProblem","AssembleJacobian" );

            mMat->mat_put_scalar( 0.0 );

            // start timer
            tic tTimer;

            // assemble Jacobian
            mSolverInterface->assemble_jacobian( mMat );

            if( mSolverWarehouse )
            {
                if( !mSolverWarehouse->get_output_to_matlab_string().empty() )
                {
                    // Get the nonlinear system index
                    uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver" , LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

                    // construct string for file name
                    std::string tJacFileName = mSolverWarehouse->get_output_to_matlab_string() + 
                            "." + std::to_string(tNonlinearSystemIndex) + ".jac.dat";

                    // save to file
                    mMat->save_matrix_to_matlab_file( tJacFileName.c_str() );

                    // log that output was successful
                    MORIS_LOG_INFO( "Saved Jacobian to Matlab File: ", tJacFileName.c_str() );
                }
            }

            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            MORIS_LOG_INFO( " Assembly of Jacobian on processor %u took %5.3f seconds.", ( uint ) par_rank(), ( double ) tElapsedTime / 1000);
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_residual_and_jacobian( )
        {
            mVectorRHS->vec_put_scalar( 0.0 );
            mPointVectorRHS->vec_put_scalar( 0.0 );
            mMat->mat_put_scalar( 0.0 );

            // start timer
            tic tTimer;

            mSolverInterface->fill_matrix_and_RHS( mMat, mPointVectorRHS );

            mSolverInterface->assemble_additional_DqDs_RHS_contribution( mPointVectorRHS );

            if( !mSolverInterface->get_is_forward_analysis() )
            {
                this->compute_residual_for_adjoint_solve();
            }

            if( mSolverWarehouse )
            {
                if( !mSolverWarehouse->get_output_to_matlab_string().empty() )
                {
                    // Get the nonlinear system index
                    uint tNonlinearSystemIndex = gLogger.get_iteration( "NonLinearSolver" , LOGGER_ARBITRARY_DESCRIPTOR, LOGGER_ARBITRARY_DESCRIPTOR );

                    // construct strings for file names
                    std::string tJacFileName = mSolverWarehouse->get_output_to_matlab_string() + 
                            "." + std::to_string(tNonlinearSystemIndex) + ".jac.dat";
                    std::string tResFileName = mSolverWarehouse->get_output_to_matlab_string() + 
                            "." + std::to_string(tNonlinearSystemIndex) + ".res.dat";

                    // output to matlab .dat file
                    mMat->save_matrix_to_matlab_file( tJacFileName.c_str() );
                    mPointVectorRHS->save_vector_to_matlab_file( tResFileName.c_str() );

                    // log that output was successful
                    MORIS_LOG_INFO( "Saved Jacobian and Residual to Matlab File: ", tJacFileName.c_str(), " and ", tResFileName.c_str() );
                
                }
            }

            //mMat->print();
            // stop timer
            //real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
            //MORIS_LOG_INFO( "Assembly of Residual and Jacobian on processor %u took %5.3f seconds.", ( uint ) par_rank(), ( double ) tElapsedTime / 1000);
        }

        //----------------------------------------------------------------------------------------

        void Linear_Problem::compute_residual_for_adjoint_solve()
        {
            sol::Matrix_Vector_Factory tMatFactory( mTplType );

            // create vector for jacobian times sol vec
            sol::Dist_Vector * tMatTimesSolVec = tMatFactory.create_vector(
                    mSolverInterface,
                    mVectorRHS->get_map(),
                    mSolverInterface->get_num_rhs() );

            // multiply jacobian with previous solution vector
            mMat->mat_vec_product( *mPointVectorLHS, *tMatTimesSolVec, false );

            // add contribution to RHS
            mPointVectorRHS->vec_plus_vec( 1.0, *tMatTimesSolVec, 1.0 );

            delete tMatTimesSolVec;
        }

        //----------------------------------------------------------------------------------------

        Matrix<DDRMat> Linear_Problem::compute_residual_of_linear_system()
        {
            sol::Matrix_Vector_Factory tMatFactory( mTplType );

            // get number of RHS
            uint tNumberOfRHS = mSolverInterface->get_num_rhs();

            // create vector for jacobian times sol vec
            sol::Dist_Vector * tResVec = tMatFactory.create_vector(
                    mSolverInterface,
                    mVectorRHS->get_map(),
                    tNumberOfRHS );

            // multiply jacobian with previous solution vector
            mMat->mat_vec_product( *mPointVectorLHS, *tResVec, false );

            // add contribution to RHS
            tResVec->vec_plus_vec( -1.0, *mPointVectorRHS, 1.0 );

            // norm of residual and RHS
            Cell<real> tResNorm = tResVec->vec_norm2();
            Cell<real> tRhsNorm = mPointVectorRHS->vec_norm2();

            // allocate vector for relative residuals
            Matrix<DDRMat> tRelativeResidualNorm(tNumberOfRHS,1);

            // compute and store relative norms
            for ( uint tRhsIndex=0; tRhsIndex<tNumberOfRHS;++tRhsIndex )
            {
                tRelativeResidualNorm( tRhsIndex ) =
                        tResNorm( tRhsIndex )/ std::max(tRhsNorm(tRhsIndex), MORIS_REAL_EPS);
            }

            // delete temporary vector
            delete tResVec;

            // return vector with relative residuals
            return tRelativeResidualNorm;
        }
    }
}
