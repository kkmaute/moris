/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Trust_Region_Optimizer.cpp
 *
 */

#include "cl_Communication_Tools.hpp"
//------------------------------------------------------------------------------
// nla includes
#include "cl_NLA_Trust_Region_Solver.hpp"
#include "cl_NLA_Convergence.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
//------------------------------------------------------------------------------
// dla includes
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Linear_Problem.hpp"
//------------------------------------------------------------------------------
// linalg includes
#include "cl_Matrix.hpp"
#include "fn_linspace.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
//------------------------------------------------------------------------------
// std includes
#include <ctime>
#include <cmath>
#include <tuple>
//------------------------------------------------------------------------------
// logging includes
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
//------------------------------------------------------------------------------
using namespace moris;
using namespace NLA;
using namespace dla;

//--------------------------------------------------------------------------------------------------------------------------

Trust_Region_Solver::Trust_Region_Solver( const Parameter_List& aParameterlist )
        : Nonlinear_Algorithm( aParameterlist )
{

}
//--------------------------------------------------------------------------------------------------------------------------
Trust_Region_Solver::Trust_Region_Solver(  )
        : Nonlinear_Algorithm()
{
    mLinSolverManager = new dla::Linear_Solver();
    
  
}

//--------------------------------------------------------------------------------------------------------------------------

Trust_Region_Solver::Trust_Region_Solver( dla::Linear_Solver * aLinSolver )
        : Nonlinear_Algorithm()
{
    mLinSolverManager = aLinSolver;
    
}

//--------------------------------------------------------------------------------------------------------------------------

Trust_Region_Solver::~Trust_Region_Solver()
{
}

//--------------------------------------------------------------------------------------------------------------------------
void Trust_Region_Solver::solver_nonlinear_system( Nonlinear_Problem *aNonlinearProblem )
{
    Tracer tTracer( "NonLinearAlgorithm", "Trust Region Solver", "Solve" );

    // Define flag to exit out of trust region solver loop
    bool tConverged = false;
    
    sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Petsc );

    std::string tStepType;

    // set nonlinear system
    mNonlinearProblem = aNonlinearProblem;

    // Zero out the initial guess
    //mNonlinearProblem->get_full_vector()->vec_put_scalar( 0.0 );

    // Set trust region flag to true for solver interface
    mNonlinearProblem->get_solver_interface()->set_trust_region_flag( true );

    // Set IQI (objective functional for the trust region solver to minimize)
    //std::string tIQIs = mParameterListNonlinearSolver.get< Vector< std::string > >( "NLA_trust_region_objective" );
    //mNonlinearProblem->get_solver_interface()->set_requested_IQI_names( tIQIs );

    //const char* tFileName = "Initial_Sol_Vec.hdf5";
    //tInitGuess->read_vector_from_HDF5(tFileName, "Res_Vec");
    //tInitGuess->vector_global_assembly();
    //declare IDs and vec
    // Matrix< DDSMat > tIDs(3,1,0);
    // Matrix< DDRMat > tVals(3,1,0.0);
    // tIDs( 1, 0 ) = 1;
    // tIDs( 2, 0 ) = 2;

    // tVals( 0, 0 ) = (2.0);
    // tVals( 1, 0 ) = (7.0);
    // tVals( 2, 0 ) = (-1.0);
    // mNonlinearProblem->get_full_vector()->sum_into_global_values( tIDs, tVals );
    // mNonlinearProblem->get_full_vector()->vector_global_assembly();

    Matrix< DDRMat > tUpdateOld;

    // Get trust region size
    real tTrSize = mParameterListNonlinearSolver.get< real >( "NLA_trust_region_size" );

    // get maximum number of iterations
    sint tMaxIts = mParameterListNonlinearSolver.get< sint >( "NLA_max_trust_region_iter" );

    // increase maximum number of iterations by one if static residual is required
    if ( mMyNonLinSolverManager->get_compute_static_residual_flag() )
    {
        tMaxIts++;
    }

    //Compute initial obj val
    mNonlinearProblem->get_solver_interface()->compute_IQI();
    Vector< Matrix< DDRMat > > tInitialIQI = mNonlinearProblem->get_solver_interface()->get_IQI();

    Matrix< DDRMat > tInitObjVal = {{0.0}};
    Matrix< DDRMat > tIQIValue = {{0.0}};

    // Sum it up since the initial obj val is strain energy + dirichlet energy + traction
    for (uint iIQIval = 0; iIQIval < tInitialIQI.size(); iIQIval++ )
    {
        tInitObjVal += tInitialIQI( iIQIval );
        
    }

    // get iteration id when references norm are computed
    sint tRefIts = mParameterListNonlinearSolver.get< sint >( "NLA_ref_iter" );

    // get option for computing residual and jacobian: separate or together
    bool tCombinedResJacAssembly = mParameterListNonlinearSolver.get< bool >( "NLA_combined_res_jac_assembly" );

    // initialize flags
    bool tRebuildJacobian = true;

    // check for convergence
    bool tHardBreak = false;

    // initialize convergence monitoring
    Convergence tConvergence( tRefIts );

    // trust region loop
    for ( sint It = 1; It <= tMaxIts; ++It )
    {
        // log solver iteration
        MORIS_LOG_ITERATION();

        // Set trust region size - to be used by linear solver later
        mLinSolverManager->set_trust_region_size( tTrSize );
        // assemble RHS and Jac
        if ( It > 1 )
        {
            tRebuildJacobian = mParameterListNonlinearSolver.get< bool >( "NLA_rebuild_jacobian" );
        }

        // For sensitivity analysis only: set current solution to LHS of linear system as residual is defined by A x - b
        if ( !mMyNonLinSolverManager->get_solver_interface()->is_forward_analysis() )
        {
            mNonlinearProblem->get_linearized_problem()->set_free_solver_LHS( mNonlinearProblem->get_full_vector() );
        }

        // build residual and jacobian
        mNonlinearProblem->build_linearized_problem( tRebuildJacobian, tCombinedResJacAssembly, It );

        // Create copies of tangent stiffness and residual
        mJac = mNonlinearProblem->get_linearized_problem()->get_matrix();
        mGlobalRHS = (mNonlinearProblem->get_linearized_problem()->get_solver_RHS());

        // Get residual norm
        Vector< real > tNorm = mGlobalRHS->vec_norm2();
    

        // Declare Cauchy point dist vector and its norm
        sol::Dist_Vector* tCauchyPoint = tMatFactory.create_vector(mNonlinearProblem->get_solver_interface(),mNonlinearProblem->get_full_vector()->get_map(), 1);
        tCauchyPoint->vec_put_scalar(0.0);

        Vector< real > tCauchyPointNorm;
        real tCauchyPointNormSquared;

        // Declare Solution update
        sol::Dist_Vector* tD = tMatFactory.create_vector(mNonlinearProblem->get_solver_interface(),mNonlinearProblem->get_full_vector()->get_map(), 1);
        tD->vec_put_scalar( 0.0 );
        Matrix< DDRMat > tmatD;

        // Get hessian - gradient product (dist vec as well as matrix)
        sol::Dist_Vector* tKg = tMatFactory.create_vector(mNonlinearProblem->get_solver_interface(),mNonlinearProblem->get_full_vector()->get_map(), 1);
        tKg->vec_put_scalar( 0.0 );

        mJac->mat_vec_product( *mGlobalRHS, *tKg, false );
        Matrix< DDRMat > tmatKg;
        tKg->extract_copy(tmatKg);

        // Premultiply hes-grad prod with grad
        Matrix< DDRMat > tmatGlobalRHS;
        mGlobalRHS->extract_copy( tmatGlobalRHS );
        Matrix< DDRMat > tgKg = trans(tmatGlobalRHS)*tmatKg;

        // Based on sign of grad-hes-grad product determine the Cauchy solution
        if ( tgKg( 0 ) > 0 )
        {
            real tAlpha = std::pow( norm( tmatGlobalRHS ), 2) / tgKg( 0 );            
            tCauchyPoint->vec_plus_vec( -tAlpha, *mGlobalRHS, 0.0);
            tCauchyPointNorm = tCauchyPoint->vec_norm2();
            tCauchyPointNormSquared = std::pow(tCauchyPointNorm( 0 ),2);
            
        }
        else
        {
            real tScaleFact = tTrSize * norm( tmatGlobalRHS );
            tCauchyPoint->vec_plus_vec( -tScaleFact, *mGlobalRHS, 0.0);
            tCauchyPointNormSquared = std::pow( tTrSize, 2 ); 
            MORIS_LOG("Negative curvature unpreconditioned cauchy point found");

        }

        // Based on cauchypoint norm squared decide whether to solve linear subproblem or not.
        if (tCauchyPointNormSquared >= std::pow(tTrSize,2))
        {
            MORIS_LOG("Unpreconditioned gradient cauchy point outside trust region");
            tCauchyPoint->scale_vector(tTrSize/std::pow(tCauchyPointNormSquared,0.5));
            mNonlinearProblem->get_linearized_problem()->set_free_solver_LHS( tCauchyPoint ); 
            tCauchyPointNormSquared = std::pow( tTrSize, 2 );
            tStepType = "Boundary";
        }
        else
        {
            this->solve_linear_system(It,tHardBreak);
            
        }
        Matrix< DDRMat > testResidual;
        mGlobalRHS->extract_copy(testResidual);
        Matrix< DDRMat > tSol;
        (mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS())->extract_copy(tSol);

        bool tAcceptTrSize = false;

        while (tAcceptTrSize == false)
        {
            tD = this->dogleg_step( (mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS()), 
            tCauchyPoint, tTrSize);
            tD->extract_copy(tUpdateOld);

            // compute Jd, dJd, and modelObjective
            sol::Dist_Vector* tJd = tMatFactory.create_vector(mNonlinearProblem->get_solver_interface(),mNonlinearProblem->get_full_vector()->get_map(), 1);
            mJac->mat_vec_product(*tD, *tJd, false);
        
            Matrix< DDRMat > tmatJd;
        
            tD->extract_copy( tmatD );
            tJd->extract_copy( tmatJd );
            mGlobalRHS->extract_copy(tmatGlobalRHS);
        
            Matrix< DDRMat > tDJd = trans(tmatD) * tmatJd;
            Matrix< DDRMat > tModelObjective = trans(tmatGlobalRHS) * tmatD + 0.5*tDJd;

            // Update solution and Rebuild linear system
             ( mNonlinearProblem->get_full_vector() )->vec_plus_vec(    //
                1.0,
                *tD,
                1.0 );
            mNonlinearProblem->build_linearized_problem( tRebuildJacobian, tCombinedResJacAssembly, It);

            bool tIsConverged = tConvergence.check_for_convergence(
                this,
                It,
                tMaxIts,
                tHardBreak );

            // exit if convergence criterion is met
            if ( tIsConverged )
            {
               MORIS_LOG_INFO( "Number of Iterations (Convergence): %d", It );
               tConverged = tIsConverged;
               break;
            }

            // check if hard break is triggered or maximum iterations are reached
            if ( tHardBreak or ( It == tMaxIts && tMaxIts > 1 ) )
            {
               MORIS_LOG_INFO( "Number of Iterations (Hard Stop): %d", It );
               break;
            }

            // Compute rho value
            mNonlinearProblem->get_solver_interface()->compute_IQI();
            Vector< Matrix< DDRMat > > tIQI = mNonlinearProblem->get_solver_interface()->get_IQI();
            
            // Sum up all IQI values since that's the objective
            for (uint iIQIval = 0; iIQIval < tIQI.size(); iIQIval++ )
            {
                tIQIValue += tIQI( iIQIval );
            }

            // Compute rho
            real tRho = -(tIQIValue( 0 ) - tInitObjVal( 0 ))/(-tModelObjective( 0 ));

            //real tRho = -(tIQI(0)(0) - tInitialIQI(0)(0))/(-tModelObjective(0));

            // Update trust region size
            // Get convergence reason from KSP
            if (mLinSolverManager->get_conv_reason())
            {
                tStepType = "Boundary";
            }
            if ((tRho <= mParameterListNonlinearSolver.get< real >("NLA_trust_region_eta2")))
            {
                tTrSize = mParameterListNonlinearSolver.get< real >("NLA_trust_region_t1") * tTrSize; 
                mLinSolverManager->set_trust_region_size( tTrSize );
            }
            else if ((tRho > mParameterListNonlinearSolver.get< real >("NLA_trust_region_eta3")) && (tStepType == "Boundary"))
            { 
                tTrSize = mParameterListNonlinearSolver.get< real >("NLA_trust_region_t2") * tTrSize;
                mLinSolverManager->set_trust_region_size( tTrSize );
            }
            // Compute updated residual norm
            Vector< real > tNormNew = mGlobalRHS->vec_norm2();

            bool tWillAccept = false;

            if ((tRho >= mParameterListNonlinearSolver.get< real >("NLA_trust_region_eta1") ) || ((tRho >= -0.0) && ( tNormNew( 0 ) < tNorm( 0 ) ))) 
            {
                tInitialIQI = tIQI;
                tWillAccept = true;
                tAcceptTrSize = true;
            }

            if ( tWillAccept == false )
            {
                // revert back to the previous iteration solution
                ( mNonlinearProblem->get_full_vector() )->vec_plus_vec(    //
                -1.0,
                *tD,
                1.0 );

                MORIS_LOG_INFO("Trust region solver iteration not accepted. Rebuilding linear system");
                mNonlinearProblem->build_linearized_problem( tRebuildJacobian, tCombinedResJacAssembly, It );

                tStepType = "Boundary";
            
            }


        }

        if ( tConverged )
        {
            MORIS_LOG( "Trust region solver converged ");
            
            break;
        }



        
        //tD->extract_copy(tUpdateOld);

        //std::cout<<"Print out values"<<tSol(0);

          

    }
    //Matrix< DDRMat > tUpdate;
    //mNonlinearProblem->get_linearized_problem()->get_full_solver_LHS()->extract_copy(tUpdate);
    if (tConverged ==  false)
    {
        MORIS_LOG("Trust region solver failed to converge. Try again");
    }
    mNonlinearProblem->get_solver_interface()->set_trust_region_flag( false );

    //std::cout<<"Print out values"<<tUpdate(0);
    //std::cout<<"Print out values"<<tUpdateOld(0);
    
       
}    // end trust-region algorithm

//--------------------------------------------------------------------------------------------------------------------------
/*Trust_Region_Minimizer_Output Trust_Region_Solver::solver_trust_region_minimizer_system( Matrix< DDRMat > &aX,
                                                                Matrix< DDRMat > &aR,                                                                
                                                                real &aTrSize,
                                                                Nonlinear_Problem *aNonlinearProblem )                 
{
    // Preallocate z vectors
    Matrix< DDRMat > tZ( aX.n_rows(), aX.n_cols() , 0.0 );

    Matrix< DDRMat > tResNormSq = trans(aR) * aR;

    // Exit criteria with preconditioner
    real tCGInexactRelTol = mParameterListNonlinearSolver.get< moris::real >( "CG_INEXACT_REL_TOL" );
    real tCGTol           = mParameterListNonlinearSolver.get< moris::real >( "CG_TOL" );
    real tCGTolSquared    = std::max( tCGInexactRelTol * tCGInexactRelTol * ( tResNormSq(0) ), tCGTol * tCGTol );

    if ( tResNormSq( 0 ) < tCGTolSquared)
    {
        //  we can return
        Trust_Region_Minimizer_Output tOutput;
        tOutput.mZ = tZ;
        tOutput.mQ = tZ;
        tOutput.mStepType = "Interior" ;
        return tOutput;
    }

    // Preconditioner r product
    Matrix< DDRMat > tPr ;//= //this->mult_with_preconditioner( tR );
    
    // compute d
    Matrix< DDRMat > tD = -tPr;

    // Compute Cauchy point
    Matrix< DDRMat > tQ = tD;

    // Compute r , prec, r product
    Matrix< DDRMat > tRPr = trans( aR ) * tPr;

    // Preallocate zz and zd
    real tZz = 0.0;
    real tZd = 0.0;
    real tDd;

    // Get max number of CG iterations
    moris::uint tMaxCGIts = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_cg_iter" );

    // Decide which zz to use based on whether use preconditioned inner products or not
    if ( mParameterListNonlinearSolver.get< bool >( "NLA_use_preconditioned_inner_products" ) )
    {
        real tDd = tRPr( 0 );
    }
    else
    {
        Matrix< DDRMat > tDsq = trans( tD ) * tD;
        real tDd = tDsq( 0 );
        
    }

    // start CG iterations
    for ( uint iCGIters = 0; iCGIters < tMaxCGIts; iCGIters++ )
    {
        // Compute curvature
        Matrix< DDRMat > tCurvature = trans( tD ) * ( mJac->mat_vec_product( tD ) );

        // compute alpha
        real tAlpha = tRPr( 0 ) / tCurvature(0);

        // interim z value
        Matrix< DDRMat > tZNp1 = tZ + tAlpha * tD;
        real tZzNp1 = this->update_step_length_squared( tAlpha, tZz, tZd, tCurvature(0) );

        // If point of negative curvature, we are done
        if ( tCurvature(0) <= 0.0 )
        {
            Matrix< DDRMat > tZOut = this->project_to_boundary_with_coeffs( tZ, tD, aTrSize, tZz, tZd, tDd );
            
            return {tZOut, tQ, "Negative Curvature"};
            MORIS_LOG_INFO( "Number of CG iterations %d", iCGIters + 1 );
        }
        
        // If point is outside trust region, we need to project it back
        else if ( tZzNp1 > aTrSize * aTrSize )
        {
            // Project z back to boundary
            Matrix< DDRMat > tZOut = this->project_to_boundary_with_coeffs( tZ, tD, aTrSize, tZz, tZd, tDd );
            
            // Update step type
            Trust_Region_Minimizer_Output tOutput;
            tOutput.mZ = tZOut;
            tOutput.mQ = tQ;
            tOutput.mStepType = "Boundary" ;
            return tOutput;
            //return {tZOut, tQ, "Boundary"};
            MORIS_LOG_INFO( "Number of CG iterations %d", iCGIters + 1 );
        }

        // Assign new z
        tZ = tZNp1;
        
        // Update value of r
        aR += tAlpha * this->hessian_vector_product( tD );

        // Update preconditioner r product
        tPr = ( aR );

        // Update r, Pr, r product
        Matrix< DDRMat > tRPrNp1 = trans( aR ) * tPr;

        // Exit if r is small enough
        if ( norm( aR )*norm( aR )* < tCGTolSquared )
        {
            //  we can return
            Trust_Region_Minimizer_Output tOutput;
            tOutput.mZ = tZ;
            tOutput.mQ = tQ;
            tOutput.mStepType = "Interior" ;
            return tOutput;
            
            MORIS_LOG_INFO( "Number of CG iterations %d", tMaxCGIts );
        }

        // Update beta
        real tBeta = tRPrNp1( 0 ) / (tRPr)(0);
        tRPr = tRPrNp1;
        
        // Update d
        tD = -tPr + tBeta * tD;

        // Update zz
        tZz = tZzNp1;
        
        // Compute new zd and dd
        if ( mParameterListNonlinearSolver.get< bool >( "NLA_use_preconditioned_inner_products" ) )
        {       
            std::tie( tZd, tDd ) = this->cg_inner_products_preconditioned( tAlpha, tBeta, tZd, tDd, tRPr(0), tZ, tD );
        }
        else
        {
            std::tie( tZd, tDd ) = this->cg_inner_products_unpreconditioned( tZ, tD );
        }

    }// end loop over CG iterations

    // If we reach here, we have not converged, so we return the last z and q
    Trust_Region_Minimizer_Output tOutput;
    tOutput.mZ = tZ;
    tOutput.mQ = tQ;
    tOutput.mStepType = "Interior_" ;
    return tOutput;
    //return {tZ, tQ, "Interior_"};
    MORIS_LOG_INFO( "Number of CG iterations %d", tMaxCGIts );
}    */
//--------------------------------------------------------------------------------------------------------------------------
void Trust_Region_Solver::solve_linear_system(
        sint& aIter,
        bool& aHardBreak )
{
    // Solve linear system
    mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), aIter );
    
}

//--------------------------------------------------------------------------------------------------------------------------
void Trust_Region_Solver::get_full_solution( moris::Matrix< DDRMat > &LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------
void Trust_Region_Solver::get_solution( moris::Matrix< DDRMat > &LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_copy( LHSValues );
}
//--------------------------------------------------------------------------------------------------------------------------
void Trust_Region_Solver::extract_my_values(
        const uint&                 aNumIndices,
        const Matrix< DDSMat >&     aGlobalBlockRows,
        const uint&                 aBlockRowOffsets,
        Vector< Matrix< DDRMat > >& LHSValues )
{
    mNonlinearProblem->get_full_vector()->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}

//--------------------------------------------------------------------------------------------------------------------------
sol::Dist_Vector* Trust_Region_Solver::preconditioned_project_to_boundary( sol::Dist_Vector* aZ,
                                                              sol::Dist_Vector* aD,
                                                              real          &aTrSize,
                                                              real          &aZz)
                                                
{
    // find tau s.t. (z + tau*d)^2 = trSize^2

    //Matrix< DDRMat > tPd = (aD);
    Matrix< DDRMat > tD;
    Matrix< DDRMat > tZ;
    aD->extract_copy(tD);
    aZ->extract_copy(tZ);

    Matrix< DDRMat > tDd = trans( tD ) * tD;
    Matrix< DDRMat > tZd = trans( tZ ) * tD;

    moris::real tTau = ( std::pow( (aTrSize * aTrSize - aZz ) * tDd(0) + tZd(0) * tZd(0), 0.5 ) - tZd(0) ) / tDd(0);

    //return aZ + tTau * aD;
    aZ->vec_plus_vec( tTau, *aD, 1.0 );
    return aZ;
}

//--------------------------------------------------------------------------------------------------------------------------
Matrix< DDRMat > Trust_Region_Solver::project_to_boundary_with_coeffs( Matrix< DDRMat > &aZ,
                                                                       Matrix< DDRMat > &aD,
                                                                       real          &aTrSize,
                                                                       real          &aZz,
                                                                       real          &aZd,
                                                                       real          &aDd    )
                                                
{
    // find tau s.t. (z + tau*d)^2 = trSize^2

    real tTau = ( std::pow( (aTrSize*aTrSize - aZz)*aDd + aZd*aZd, 0.5 ) - aZd ) / aDd;

    return aZ + tTau * aD;
}
//--------------------------------------------------------------------------------------------------------------------------
real Trust_Region_Solver::update_step_length_squared( real          &aAlpha,
                                                      real          &aZz,
                                                      real          &aZd,
                                                      real          &aDd )
                                                
{
    return aZz + 2.0*aAlpha*aZd + aAlpha*aAlpha*aDd;
}
//--------------------------------------------------------------------------------------------------------------------------
std::tuple< real, real > Trust_Region_Solver::cg_inner_products_unpreconditioned( Matrix< DDRMat >         &aZ,
                                                                                  Matrix< DDRMat >         &aD )
                                                
{
    Matrix< DDRMat > tZd = trans(aZ) * aD;
    Matrix< DDRMat > tDd = trans(aD) * aD;

    return std::make_tuple(tZd( 0 ), tDd( 0 ));
}
//--------------------------------------------------------------------------------------------------------------------------
std::tuple< real, real > Trust_Region_Solver::cg_inner_products_preconditioned(   real         &aAlpha,
                                                                                  real         &aBeta,
                                                                                  real         &aZd,
                                                                                  real         &aDd,
                                                                                  real         &aRPr,
                                                                                  Matrix< DDRMat > &aZ,
                                                                                  Matrix< DDRMat > &aD )
                                                
{
    // recurrence formulas from Gould et al. doi:10.1137/S1052623497322735
    aZd = aBeta * ( aZd + aAlpha * aDd );
    aDd = aRPr + ( aBeta * aBeta * aDd ) ;

    return std::make_tuple(aZd, aDd);

}

//--------------------------------------------------------------------------------------------------------------------------
sol::Dist_Vector* Trust_Region_Solver::dogleg_step(sol:: Dist_Vector* aZ,
                                                  sol:: Dist_Vector*  aQ,
                                                  real          &aTrSize)
                                                                                   
                                                
{
    //Matrix< DDRMat > tCc = trans( aQ ) * aQ;
    //Matrix< DDRMat > tNn = trans( aZ ) * aZ;

    sol::Matrix_Vector_Factory tMatFactory(sol::MapType::Petsc);
    sol::Dist_Vector* tD = tMatFactory.create_vector( mNonlinearProblem->get_solver_interface(), aZ->get_map(), 1 );

    Vector< real > tCc = aQ->vec_norm2();
    Vector< real > tNn = aZ->vec_norm2();
    real tTrustRegionSizeSq = aTrSize * aTrSize;


    // If cauchy point is outside trust region, we need to project it back
    if ( std::pow(tCc(0),2) > tTrustRegionSizeSq )
    {
        // Project cauchy point back to boundary
        MORIS_LOG_INFO( "Cauchy point outside trust region at dist = %f", std::sqrt( tCc(0) ) );
        //aQ = ( aTrSize / std::sqrt( tCc(0) ) ) * aQ;
        tD->vec_plus_vec(aTrSize / tCc(0), *aQ, 0.0);
        return tD;
    }
    if ( std::pow(tCc(0),2) > std::pow(tNn(0),2) )
    {
        // Likely due to inaccurate preconditioner
        MORIS_LOG_INFO( "Cauchy point outside newton, check preconditioner" );
        tD->vec_plus_vec(1.0, *aQ, 0.0);
        return tD;
    }
    
    if ( std::pow(tNn(0),2) > tTrustRegionSizeSq )
    {
        real tCCrev = std::pow( tCc(0), 2);
        aZ->vec_plus_vec(-1.0, *aQ, -1.0);
        sol::Dist_Vector* aZRev = this->preconditioned_project_to_boundary( aQ, aZ, aTrSize, tCCrev );
        tD->vec_plus_vec(1.0, *aZRev, 0.0);

        return tD;
    }

    tD->vec_plus_vec(-1.0, *aZ, 0.0);

    return tD;

}
// void Trust_Region_Solver::build_preconditioner( Nonlinear_Problem *aNonlinearProblem )
// {
//     // Get Hessian of objective (tangent stiffness matrix)
//     //sol::Dist_Matrix tJac = this->mJac;
// }


//--------------------------------------------------------------------------------------------------------------------------
void Trust_Region_Solver::initialize_variables( Nonlinear_Problem *aNonlinearProblem )
{
    mNonlinearProblem = aNonlinearProblem;

    std::cout <<" Initializing variables \n";

    bool tRebuildJacobian = true;
    sint tDummy           = 1;

    mNonlinearProblem->build_linearized_problem( tRebuildJacobian, false, tDummy );    // build the linearized problem
                         
    mGlobalRHS = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();    // set pointer to RHS
    mJac       = mNonlinearProblem->get_linearized_problem()->get_matrix();        // set pointer to jacobian matrix

    // Solve to get initial solution
    //mLinSolverManager->solver_linear_system( mNonlinearProblem->get_linearized_problem(), tDummy );   
}
//--------------------------------------------------------------------------------------------------------------------------

