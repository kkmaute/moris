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
    // Set nonlinear problem object
    mNonlinearProblem = aNonlinearProblem;

    // Compute initial objective value
    mNonlinearProblem->get_solver_interface()->compute_IQI();

    // get trust region size
    moris::real tTrSize = mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_size" );

    // Get max trust region iterations
    moris::sint tMaxTrustRegionIts = mParameterListNonlinearSolver.get< moris::sint >( "NLA_max_trust_region_iter" );

    // Declare step type variable
    std::string tStepType;

    // Declare Cauchy point, newton point and cauchy point norm variables
    Matrix< DDRMat > tX;
    Matrix< DDRMat > tQ;
    Matrix< DDRMat > tZ;
    Matrix< DDRMat > tCauchyPointNormSquared;

    // Declare dist vectors for Cauchy point and Newton point
    sol::Dist_Vector* tXd = mNonlinearProblem->get_full_vector();
    tXd->vec_put_scalar(0.0);

    sol::Dist_Vector* tQd = mNonlinearProblem->get_full_vector();
    tQd->vec_put_scalar(0.0);

    sol::Dist_Vector* tZd = mNonlinearProblem->get_full_vector();
    tZd->vec_put_scalar(0.0);


    // Start trust region iterations
    for ( moris::sint iTrustRegionIter = 0; iTrustRegionIter < tMaxTrustRegionIts; iTrustRegionIter++ )
    {
        bool tRebuildJacobian = true;

        mNonlinearProblem->build_linearized_problem( tRebuildJacobian, true, iTrustRegionIter );
        mGlobalRHS = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();
        mJac = mNonlinearProblem->get_linearized_problem()->get_matrix();

        // Compute update in the IQI value
        mNonlinearProblem->get_solver_interface()->compute_IQI();

        // Obtain the IQI values
        Vector< moris::Matrix< DDRMat > > tInitialObjective = mNonlinearProblem->get_solver_interface()->get_IQI();
        
        sol::Dist_Vector* tKg =  mNonlinearProblem->get_full_vector();
        mJac->mat_vec_product( *mGlobalRHS, *tKg, false );
        Matrix< DDRMat > tKG;
        tKg->extract_copy( tKG );
        

        // Assign tr size member variable
        mTrSize = tTrSize;

        // Inform lin solver of the current trust region size
        mLinSolverManager->set_trust_region_size( tTrSize );

        // Get product of g with Kg
        Matrix< DDRMat > tG;
        Matrix< DDRMat > tH;
        mGlobalRHS->extract_copy(tG);

        Matrix< DDRMat > tGKG = trans(tG) * tKG;
        
        bool tHardBreak = false;

        // conditional check for computing newton point and cauchy point
        if ( tGKG( 0 ) > 0.0 )
        {
            // Define alpha
            Matrix< DDRMat > tAlpha = -trans(tG) * tG / ( tGKG(0) );
            
            // Define Cauchy Point
            //tQ = tAlpha( 0 ) * tG;
            tQd = mGlobalRHS;
            tQd->scale_vector( tAlpha( 0 ) );
            // Copy into matrix
            tQd->extract_copy( tQ );

            // Define norm of Cauchy Point
            tCauchyPointNormSquared = trans( tQ ) * tQ;

        }
        else
        { 
            // Norm of grad sq
            Matrix< DDRMat >tGNorm = trans(tG)*tG;
            
            // Define Cauchy Point
            //tQ = -( tTrSize / std::sqrt( tGNorm(0) ) )*tG;
            tQd = mGlobalRHS;
            tQd->scale_vector( -tTrSize / std::sqrt( tGNorm(0) ) );
            // Copy into matrix
            tQd->extract_copy( tQ );

            // Define norm of Cauchy Point
            tCauchyPointNormSquared( 0 ) = tTrSize * tTrSize;

            // Let user know negative curvature was detected
            MORIS_LOG_INFO( "Negative curvature unpreconditioned cauchy point direction found." );

        }

        if ( tCauchyPointNormSquared( 0 )  > tTrSize * tTrSize )
        {
            // Project Cauchy Point back to boundary
            MORIS_LOG_INFO( "unpreconditioned gradient cauchy point outside trust region at dist = %f", std::sqrt( tCauchyPointNormSquared( 0 )  ) );
            //Matrix< DDRMat > tQ;
            //tQ = ( tTrSize / std::sqrt( tCauchyPointNormSquared( 0 )  ) ) * tQ;

            tQd->scale_vector( tTrSize / std::sqrt( tCauchyPointNormSquared( 0 )  ) );

            // Define norm of Cauchy Point
            tCauchyPointNormSquared( 0 )  = tTrSize * tTrSize;

            // Define newton point
            tZd = tQd;

            // Copy into matrix
            tZd->extract_copy(tZ);

            // Define Step type 
            tStepType = "Boundary";
            
            MORIS_LOG_INFO( "CG Iters = %d", 1 );
        }

        else
        {
            // Compute Newton Point using trust region minimizer
            this->solve_linear_system( iTrustRegionIter, tHardBreak );

            // Get the Newton vector
            dla::Linear_Problem* tLinProb = mNonlinearProblem->get_linearized_problem();
            tZd = tLinProb->get_free_solver_LHS();
            // Extract Copy 
            tZd->extract_copy(tZ);
            
        }
        
        // Flag for trust region size update
        bool tTrustRegionSizeUpdateFlag = false;
        
        while ( tTrustRegionSizeUpdateFlag == false )
        {
            sol::Dist_Vector* tDd = this->dogleg_step( tZd, tQd, tTrSize );

            sol::Dist_Vector* tJd = mNonlinearProblem->get_full_vector();
            tJd->vec_put_scalar( 0.0 );
            Matrix< DDRMat > tJdd;
            Matrix< DDRMat > tDdd;
            
            mJac->mat_vec_product(*tDd , *tJd , false);
            tJd->extract_copy( tJdd );
            tDd->extract_copy( tDdd );
            Matrix< DDRMat > tDJd = trans( tDdd ) * tJdd;
            Matrix< DDRMat > tGradUpProd = trans( tG )*( tDdd );
            moris::real tModelObjective = tGradUpProd(0) + 0.5*tDJd(0);

            real tModelNorm = norm( tG + tJdd );

            // Store old solution
            sol::Dist_Vector* tOldSol = mNonlinearProblem->get_full_vector();

            // Update solution
            //sol::Dist_Vector* tY = mNonlinearProblem->get_full_vector();
            //tY->vec_plus_vec( 1.0 , *tDd , 1.0 );

            // Compute new residual and tangent stiffness matrix       
            // Make a copy of the aNonlinearProblem
            // Nonlinear_Problem* tNonlinearProblem = aNonlinearProblem;

            // tNonlinearProblem->get_solver_interface()->set_solution_vector( tY );
            // tNonlinearProblem->build_linearized_problem( tRebuildJacobian, false, iTrustRegionIter ); 
     
            // sol::Dist_Vector* tGy = tNonlinearProblem->get_linearized_problem()->get_solver_RHS();

            // // Get the new residual norm
            // Vector< moris::real > tResNormy = tGy->vec_norm2();

            // If residual converged then exit
            Convergence tConvergence;

            (mNonlinearProblem->get_full_vector())->vec_plus_vec( 1.0 , *tDd , 1.0 );
            mNonlinearProblem->build_linearized_problem( tRebuildJacobian, false, iTrustRegionIter );

            // Get new residual and its norm
            sol::Dist_Vector* tGy = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();

            Vector< real > tResNormUpdate = tGy->vec_norm2();

            // Compute update in the IQI value
            mNonlinearProblem->get_solver_interface()->compute_IQI();

            // Obtain the IQI values
            Vector< moris::Matrix< DDRMat > > tIQIs = mNonlinearProblem->get_solver_interface()->get_IQI();

            // Based on convergence criteria decide whether to update the preconditioner and other things.
            real tRealImprove = -(tIQIs( 0 )( 0 ) - tInitialObjective( 0 )( 0 ));
            

            bool tHardBreak = false;

            bool tIsConverged = tConvergence.check_for_convergence( this , iTrustRegionIter, tMaxTrustRegionIts, tHardBreak );

            if ( tIsConverged )
            {
               if ( tHardBreak )
               {
                   continue;
               }
               
               break;
            }


            real tRho = tRealImprove/(-tModelObjective);

            if ((tRho <= mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_eta2" )) & (tRho != tRho))
            {
                tTrSize = tTrSize*mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_t1" );
            }
            else if( (tRho > mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_eta3" )) & (tStepType == "Boundary"))
            {
                tTrSize = tTrSize*mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_t2" );
            }

            bool tWillAccept = false;

            if ( (tRho >= mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_eta1" )) || ( (tRho >= -0.0) & (tModelNorm <= tResNormUpdate( 0 )) ) )
            {
                 tWillAccept = true;
                 MORIS_LOG("Iteration accepted");
            }

            if ( tWillAccept )
            {
                 // Update gradient value
                 mGlobalRHS = mNonlinearProblem->get_linearized_problem()->get_solver_RHS();

                 // Now happy with trust region size
                 tTrustRegionSizeUpdateFlag = true;

                 // Update the value of mJac
                 mJac = mNonlinearProblem->get_linearized_problem()->get_matrix();
            }
            else
            {
                 tStepType = "Boundary";
                 MORIS_LOG("Trust region iteration not accepted \n");

                // Go back to the old solution - iteration is not accepted
                mNonlinearProblem->get_full_vector()->vec_plus_vec(1.0 , *tOldSol, 0.0);

                // Rebuild Linear problem
                //mNonlinearProblem->build_linearized_problem( true, iTrustRegionIter, 0 );


            }
            
            // Throw error if tr size falls below threshold
            MORIS_ASSERT(tTrSize > mParameterListNonlinearSolver.get< moris::real >( "NLA_trust_region_min_size" ),"Trust region size smaller than threshold. Please debug");
            

            // Update preconditioner if CG iterations exceed threshold
           
        }

    }
    
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

    Vector< real > tCc = aQ->vec_norm2();
    Vector< real > tNn = aZ->vec_norm2();
    real tTrustRegionSizeSq = aTrSize * aTrSize;


    // If cauchy point is outside trust region, we need to project it back
    if ( std::pow(tCc(0),2) > tTrustRegionSizeSq )
    {
        // Project cauchy point back to boundary
        MORIS_LOG_INFO( "Cauchy point outside trust region at dist = %f", std::sqrt( tCc(0) ) );
        //aQ = ( aTrSize / std::sqrt( tCc(0) ) ) * aQ;
        aQ->scale_vector(aTrSize / tCc(0));

        return aQ;
    }
    if ( std::pow(tCc(0),2) > std::pow(tNn(0),2) )
    {
        // Likely due to inaccurate preconditioner
        MORIS_LOG_INFO( "Cauchy point outside newton, check preconditioner" );
        return aQ;
    }
    
    if ( std::pow(tNn(0),2) > tTrustRegionSizeSq )
    {
        real tCCrev = std::pow( tCc(0), 2);
        aZ->vec_plus_vec(-1.0, *aQ, 1.0);
        sol::Dist_Vector* aZRev = this->preconditioned_project_to_boundary( aQ, aZ, aTrSize, tCCrev );

        return aZRev;
    }

    return aZ;

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

