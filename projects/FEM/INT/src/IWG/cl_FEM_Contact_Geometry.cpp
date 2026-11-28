/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include "cl_FEM_Contact_Geometry.hpp"
#include "cl_FEM_Contact_Gap_Data.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_MTK_Ray_Line_Intersection.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include "fn_inv.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"
#include <algorithm>
#include <cmath>

namespace moris::fem
{
    Matrix< DDRMat > Contact_Geometry::remap_nonconformal_rays_undeformed_geometry(
            Field_Interpolator_Manager* aLeaderFieldInterpolatorManager,
            Field_Interpolator_Manager* aFollowerFieldInterpolatorManager )
    {
        // get geometry interpolator for leader and follower IG elements
        Geometry_Interpolator* tLeaderIGGI   = aLeaderFieldInterpolatorManager->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI = aFollowerFieldInterpolatorManager->get_IG_geometry_interpolator();

        //        // reset evaluation flags for the geometry interpolators -
        //        // FIXME: not sure whether this is needed
        //        tLeaderIGGI->reset_eval_flags();
        //        tFollowerIGGI->reset_eval_flags();

        // get the integration point location
        const Matrix< DDRMat >& tIGPoint = tLeaderIGGI->valx();

        // get the normal vector at the integration point
        const Matrix< DDRMat >& tNormal = tLeaderIGGI->get_normal();

        // Get the nodes of the follower element
        const Matrix< DDRMat >& tFollowerCoordinates = tFollowerIGGI->get_space_coeff();

        // check that we are dealing with line elements
        MORIS_ASSERT( tFollowerCoordinates.n_cols() == 2 && tNormal.n_rows() == 2,
                "Contact_Geometry::remap_nonconformal_rays - Nonconformal FD scheme is only implemented for line elements." );

        // Perform the mapping
        mtk::Ray_Line_Intersection tRLI( tFollowerCoordinates.n_cols() );
        tRLI.set_ray_origin( trans( tIGPoint ) );
        tRLI.set_ray_direction( tNormal );
        tRLI.set_target_origin( trans( tFollowerCoordinates.get_row( 0 ) ) );
        tRLI.set_target_span( trans( tFollowerCoordinates.get_row( 1 ) - tFollowerCoordinates.get_row( 0 ) ) );
        tRLI.perform_raytracing();

        Matrix< DDRMat > tFollowerSpaceTime;
        tLeaderIGGI->get_space_time( tFollowerSpaceTime );    // use leader as follower may not be initialized

        if ( tRLI.has_intersection() )
        {
            tFollowerSpaceTime( 0 ) = tRLI.get_intersection_parametric()( 0 );

            // set the space time coordinates of the follower side
            aFollowerFieldInterpolatorManager->set_space_time_from_local_IG_point( tFollowerSpaceTime );
        }
        else
        {
            tFollowerSpaceTime( 0 ) = -1.0;    // set to -1 if no intersection found
        }

        return tFollowerSpaceTime;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Contact_Geometry::remap_nonconformal_rays(
            const bool                     aUseDeformedGeometryForGap,
            const Vector< MSI::Dof_Type >& aDisplDofTypes,
            Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
            Field_Interpolator_Manager*    aLeaderPreviousFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerPreviousFieldInterpolatorManager,
            std::unique_ptr< GapData >&    aGapData,
            bool                           aEvaluateSecondDerivatives )
    {
        if ( aUseDeformedGeometryForGap )
        {
            return Contact_Geometry::remap_nonconformal_rays_deformed_geometry(
                    aDisplDofTypes,
                    aLeaderFieldInterpolatorManager,
                    aFollowerFieldInterpolatorManager,
                    aLeaderPreviousFieldInterpolatorManager,
                    aFollowerPreviousFieldInterpolatorManager,
                    aGapData,
                    aEvaluateSecondDerivatives );
        }
        else
        {
            return Contact_Geometry::remap_nonconformal_rays_undeformed_geometry(
                    aLeaderFieldInterpolatorManager,
                    aFollowerFieldInterpolatorManager );
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > Contact_Geometry::remap_nonconformal_rays_deformed_geometry(
            const Vector< MSI::Dof_Type >& aDisplDofTypes,
            Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
            Field_Interpolator_Manager*    aLeaderPreviousFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerPreviousFieldInterpolatorManager,
            std::unique_ptr< GapData >&    aGapData,
            bool                           aEvaluateSecondDerivatives )
    {
        Field_Interpolator* tLeaderFieldInterpolator =
                aLeaderFieldInterpolatorManager->get_field_interpolators_for_type( aDisplDofTypes( 0 ) );
        Field_Interpolator* tFollowerFieldInterpolator =
                aFollowerFieldInterpolatorManager->get_field_interpolators_for_type( aDisplDofTypes( 0 ) );

        // get the space dimension
        uint const tSpaceDim     = tLeaderFieldInterpolator->get_space_dim();
        uint const tNumTimeBases = tLeaderFieldInterpolator->get_number_of_time_bases();

        // create GapData pointer
        if ( aGapData == nullptr )
        {
            aGapData = std::make_unique< GapData >();
        }

        // check if the gap data need to be evaluated
        if ( !aGapData->mEval
                && ( !aEvaluateSecondDerivatives || aGapData->mHasSecondOrderDerivatives
                        || std::abs( aGapData->mEta( 0 ) ) > 1.0 ) )
        {
            Matrix< DDRMat > tFollowerParamPoint( tSpaceDim, 1, 0.0 );
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = aGapData->mEta.matrix_data();
            return tFollowerParamPoint;
        }

        // An upgrade starts from the converged intersection of this same point.
        // Geometry edits must set mEval, so a cached root is never reused after an edit.
        const bool tUpgradeCachedRay = !aGapData->mEval;

        // Mark all orders invalid during evaluation. Full requests reuse their
        // allocated Hessian storage; first-order requests and failed rays clear it.
        aGapData->mEval                      = true;
        aGapData->mHasSecondOrderDerivatives = false;

        // check that we are dealing with line elements
        MORIS_ERROR( tSpaceDim == 2,
                "Contact_Geometry::remap_nonconformal_rays - only implemented for line elements." );
        MORIS_ERROR( tNumTimeBases == 1,
                "Contact_Geometry::remap_nonconformal_rays - only implemented for constant time slabs." );

        // get IG geometry interpolator for leader and follower
        Geometry_Interpolator* tLeaderIGGI   = aLeaderFieldInterpolatorManager->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI = aFollowerFieldInterpolatorManager->get_IG_geometry_interpolator();

        // get coordinates (global and parametric) of IG nodes
        // const Matrix< DDRMat >& tLeaderIGNodes   = tLeaderIGGI->get_space_coeff();
        const Matrix< DDRMat >& tFollowerIGNodes = tFollowerIGGI->get_space_coeff();

        // const Matrix< DDRMat >& tLeaderIGNodesParam   = tLeaderIGGI->get_space_param_coeff();
        const Matrix< DDRMat >& tFollowerIGNodesParam = tFollowerIGGI->get_space_param_coeff();

        // get IP geometry interpolator for leader and follower
        Geometry_Interpolator* tLeaderIPGI   = tLeaderFieldInterpolator->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIPGI = tFollowerFieldInterpolator->get_IG_geometry_interpolator();

        // get displacements for leader and follower
        const Matrix< DDRMat >& tLeaderUHat   = tLeaderFieldInterpolator->get_coeff();    // degrees of freedom at leader IG nodes (displacements)
        const Matrix< DDRMat >& tFollowerUHat = tFollowerFieldInterpolator->get_coeff();

        // determine number of IP nodes and dofs (assume same for leader and follower)
        const uint tNumNodes = tLeaderUHat.n_rows();
        const uint tNumDofs  = tLeaderUHat.numel();

        MORIS_ASSERT( tNumDofs == tNumNodes * tSpaceDim,
                "Contact_Geometry::remap_nonconformal_rays - Number of dofs does not match number of nodes." );

        // get geometry and displacements for leader quadrature point
        const Matrix< DDRMat >& tLeaderXgp = tLeaderIGGI->valx();
        // const Matrix< DDRMat >& tLeaderNXgp     = tLeaderIGGI->NXi();
        // const Matrix< DDRMat >& tLeaderdNXgpdxi = tLeaderIGGI->dNdXi();

        const Matrix< DDRMat >& tLeaderUgp  = tLeaderFieldInterpolator->val();
        const Matrix< DDRMat >& tLeaderNUgp = tLeaderIPGI->NXi();
        // const Matrix< DDRMat >& tLeaderdNUgpdr = tLeaderIPGI->dNdXi();

        // compute normal vector and its derivatives
        Matrix< DDRMat > tLeaderNormal;
        Matrix< DDRMat > tLeaderdNormaldU( tSpaceDim, tNumDofs );
        Matrix< DDRMat > tLeaderdNormal2dU2;
        if ( aEvaluateSecondDerivatives )
        {
            tLeaderdNormal2dU2.set_size( tSpaceDim, tNumDofs * tNumDofs );
        }
        Matrix< DDRMat > tLeaderRefNormal;

        GapData::compute_outward_normal_at_gp_for_consistent_deformed_geometry(
                tLeaderFieldInterpolator,
                tLeaderIGGI,
                tLeaderNormal,
                tLeaderRefNormal,
                tLeaderdNormaldU,
                tLeaderdNormal2dU2,
                true,
                aEvaluateSecondDerivatives );

        // compute ray intersection with follower element
        Matrix< DDRMat > tSolRay( tSpaceDim, 1, 0.0 );        // FIXME: only implemented for constant time slabs
        Matrix< DDRMat > tSolRayPrev( tSpaceDim, 1, 1.0 );    // FIXME: only implemented for constant time slabs
        Matrix< DDRMat > tJacRayInv;

        if ( tUpgradeCachedRay )
        {
            tSolRay( 0 )                              = aGapData->mGap;
            tSolRay( { 1, tSpaceDim - 1 }, { 0, 0 } ) = aGapData->mEta.matrix_data();
        }

        Matrix< DDRMat > tFollowerParamPoint( tSpaceDim, 1, 0.0 );

        real       tRayResRefNorm      = 0.0;
        real       tRayResAbsNormTol   = 0.0;
        const real tRayResRefNormTol   = 1e-10;
        const real tRayResAbsTolFactor = 1e-12;
        const real tStagnationNormTol  = 1e-12;
        const uint tMaxNewtonSteps     = 100;
        bool       tNewtonConverged    = false;

        for ( uint inew = 0; inew < tMaxNewtonSteps; inew++ )
        {
            // clip eta to be between -1 and 1
            tSolRay( 1 ) = std::min( 1.0, std::max( -1.0, tSolRay( 1 ) ) );    // FIXME: only works in 2D

            // set local parameter point for follower
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = tSolRay( { 1, tSpaceDim - 1 }, { 0, 0 } );

            aFollowerFieldInterpolatorManager->set_space_time_from_local_IG_point( tFollowerParamPoint );

            // get geometry and displacements for follower quadrature point
            const Matrix< DDRMat >& tFollowerYgp     = tFollowerIGGI->valx();
            const Matrix< DDRMat >& tFollowerVgp     = tFollowerFieldInterpolator->val();
            const Matrix< DDRMat >& tFollowerdNVgpdq = tFollowerIPGI->dNdXi();

            const Matrix< DDRMat >& tFollowerdNYgpdeta = tFollowerIGGI->dNdXi();
            Matrix< DDRMat >        tFollowerdYgpdeta  = trans( tFollowerdNYgpdeta * tFollowerIGNodes );
            // xxxx Matrix< DDRMat > tFollowerQgp      = tFollowerNYgp * tFollowerIGNodesParam;
            Matrix< DDRMat > tFollowerdQgpdeta = tFollowerdNYgpdeta * tFollowerIGNodesParam;
            Matrix< DDRMat > tFollowerdVgpdq   = tFollowerdNVgpdq * tFollowerUHat;

            Matrix< DDRMat > tFollowerdVgpdeta = trans( tFollowerdQgpdeta * tFollowerdVgpdq );
            // xxxx Matrix< DDRMat > tFollowerdNVgpdeta = trans( tFollowerdQgpdeta * tFollowerdNVgpdq );

            // compute the residual
            Matrix< DDRMat > tResRay =
                    trans( tLeaderXgp - tFollowerYgp ) + tLeaderUgp - tFollowerVgp + tSolRay( 0 ) * tLeaderNormal;

            // Fix the residual scale at the initial ray iterate. The reference
            // line spans [-1,1], so twice the deformed tangent norm gives a
            // characteristic face length (the exact length for a straight face).
            // Scaling the absolute tolerance by this length preserves units.
            const real tResRayNorm = norm( tResRay );
            if ( inew == 0 )
            {
                tRayResRefNorm         = tResRayNorm;
                const real tFaceLength = 2.0 * norm( tFollowerdYgpdeta + tFollowerdVgpdeta );
                tRayResAbsNormTol      = tRayResAbsTolFactor * tFaceLength;
            }

            // compute jacobian and its inverse
            Matrix< DDRMat > tJacRay( tSpaceDim, tSpaceDim );
            tJacRay( { 0, tSpaceDim - 1 }, { 0, 0 } )             = tLeaderNormal.matrix_data();
            tJacRay( { 0, tSpaceDim - 1 }, { 1, tSpaceDim - 1 } ) = -tFollowerdVgpdeta - tFollowerdYgpdeta;

            // check determination of jacobian
            if ( std::abs( det( tJacRay ) ) < MORIS_REAL_EPS )
            {
                break;
            }

            // compute inverse of jacobian
            tJacRayInv = inv( tJacRay );

            // Accept exact and near-zero residuals before testing stagnation.
            // Keep this after inversion: gap derivatives need the inverse even
            // when the initial ray iterate already solves the intersection.
            if ( tResRayNorm <= std::max( tRayResRefNormTol * tRayResRefNorm, tRayResAbsNormTol ) )
            {
                // MORIS_LOG_INFO( "Remap - Newton converged in %d steps.\n", inew + 1 );
                tNewtonConverged = true;
                break;
            }

            // check for stagnation - stagnation is not considered convergence
            if ( norm( tSolRay - tSolRayPrev ) < tStagnationNormTol )
            {
                // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                //                 fprintf( stdout, "inew = %i, gap = %e, eta = %e, Resnorm = %e\n", inew, tSolRay( 0 ), tSolRay( 1 ), norm( tResRay ) );
                // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                break;
            }
            tSolRayPrev = tSolRay;

            // update solution
            tSolRay -= tJacRayInv * tResRay;
        }

        bool tPlotRays = false;

        // check if Newton converged
        if ( !tNewtonConverged )
        {
            aGapData->clear_second_order_derivatives();

            // set eta to outside bounds
            aGapData->mEta.set_size( tSpaceDim - 1, 1, -2.0 );

            // Cache the failed mapping until the next geometry invalidation.
            aGapData->mEval = false;

            // return a point with eta set to -2.0
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ).fill( -2.0 );

            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            if ( tPlotRays )
            {
                sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
                Matrix< DDRMat > tRayCastPoint = trans( tLeaderXgp ) + tLeaderUgp;
                Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * tLeaderNormal;

                fprintf( stdout, "Niter = %d NotConverged %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        tTgtPoint( 0 ),
                        tTgtPoint( 1 ) );
            }
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            return tFollowerParamPoint;
        }

        // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if ( tPlotRays )
        {
            sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
            Matrix< DDRMat > tRayCastPoint = trans( tLeaderXgp ) + tLeaderUgp;
            Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * tLeaderNormal;

            fprintf( stdout, "Niter = %d Converged %e  %e  %e  %e\n",    //
                    tNiter,
                    tRayCastPoint( 0 ),
                    tRayCastPoint( 1 ),
                    tTgtPoint( 0 ),
                    tTgtPoint( 1 ) );
        }
        // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        MORIS_ASSERT( tNewtonConverged, "Contact_Geometry::remap_nonconformal_rays - Newton did not converge." );

        // check if valid intersection was found
        MORIS_ASSERT( tSolRay( 1 ) >= -1.0 && tSolRay( 1 ) <= 1.0,
                "Contact_Geometry::remap_nonconformal_rays - No valid intersection found (eta(0) out side bounds)." );

        if ( tSpaceDim == 3 )
        {
            MORIS_ASSERT( tSolRay( 2 ) >= -1.0 && tSolRay( 2 ) <= 1.0,
                    "Contact_Geometry::remap_nonconformal_rays - No valid intersection found (eta(1) out side bounds)." );
        }

        // store solution of ray intersection
        const real             tGap = tSolRay( 0 );
        const Matrix< DDRMat > tEta = tSolRay( { 1, tSpaceDim - 1 } );

        const Matrix< DDRMat > tGapVector = tGap * tLeaderNormal;

        // First-order mapping data is also needed by residual assembly.
        const Matrix< DDRMat >& tFollowerNVgp      = tFollowerIPGI->NXi();
        const Matrix< DDRMat >& tFollowerVgp       = tFollowerFieldInterpolator->val();
        const Matrix< DDRMat >& tFollowerdNVgpdq   = tFollowerIPGI->dNdXi();
        const Matrix< DDRMat >  tFollowerdQgpdeta  = tFollowerIGGI->dNdXi() * tFollowerIGNodesParam;
        const Matrix< DDRMat >  tFollowerdYgpdeta  = trans( tFollowerIGGI->dNdXi() * tFollowerIGNodes );
        const Matrix< DDRMat >  tFollowerdNVgpdeta = trans( tFollowerdQgpdeta * tFollowerdNVgpdq );
        const Matrix< DDRMat >  tIdentity          = eye( tSpaceDim, tSpaceDim );
        Matrix< DDRMat >        tFollowerdResRaydv( tSpaceDim, tNumDofs );
        Matrix< DDRMat >        tLeaderdResRaydu( tSpaceDim, tNumDofs );
        for ( uint in = 0; in < tNumNodes; ++in )
        {
            tFollowerdResRaydv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) = tFollowerNVgp( in ) * tIdentity;
            tLeaderdResRaydu( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    -tLeaderNUgp( in ) * tIdentity
                    - tGap * tLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } );
        }

        const Matrix< DDRMat > tdSolRaydv = tJacRayInv * tFollowerdResRaydv;
        const Matrix< DDRMat > tdSolRaydu = tJacRayInv * tLeaderdResRaydu;
        const Matrix< DDRMat > tdGapdv    = tdSolRaydv.get_row( 0 );
        const Matrix< DDRMat > tdGapdu    = tdSolRaydu.get_row( 0 );
        const Matrix< DDRMat > tdEtadv    = tdSolRaydv( { 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } );
        const Matrix< DDRMat > tdEtadu    = tdSolRaydu( { 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } );
        const Matrix< DDRMat > tdGapvecdu = tGap * tLeaderdNormaldU + tLeaderNormal * tdGapdu;
        const Matrix< DDRMat > tdGapvecdv = tLeaderNormal * tdGapdv;

        aGapData->set_matrix_sizes( tSpaceDim, tNumDofs, tNumNodes, aEvaluateSecondDerivatives );

        // Only the Mlika Jacobian and explicit second-order clients need these
        // DOF-pair loops, quadratic buffers, and second geometry derivatives.
        if ( aEvaluateSecondDerivatives )
        {
            // get global coordinates of  follower IP nodes
            const Matrix< DDRMat >& tFollowerIPNodes = tFollowerIPGI->get_space_coeff();

            // get geometry and displacements for follower quadrature point
            // Note: the follower parameter point was already set within Newton loop
            const Matrix< DDRMat >& tFollowerdNYgp2deta2 = tFollowerIGGI->d2NdXi2();

            const Matrix< DDRMat >& tFollowerdNVgp2dq2 = tFollowerIPGI->d2NdXi2();

            const Matrix< DDRMat > tFollowerdYgpdq   = tFollowerdNVgpdq * tFollowerIPNodes;
            const Matrix< DDRMat > tFollowerdYgp2dq2 = tFollowerdNVgp2dq2 * tFollowerIPNodes;

            const Matrix< DDRMat > tFollowerdQgp2deta2 = tFollowerdNYgp2deta2 * tFollowerIGNodesParam;
            const Matrix< DDRMat > tFollowerdVgpdq     = tFollowerdNVgpdq * tFollowerUHat;
            const Matrix< DDRMat > tFollowerdVgp2dq2   = tFollowerdNVgp2dq2 * tFollowerUHat;


            // Define matrix to convert 2nd derivative matrix to vector
            // FIXME: this is only valid for 2D line elements
            const Matrix< DDRMat > tExtractMat = { { 1, 0, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 0, 1, 0 } };

            // compute second derivative of follower velocity with respect to local coordinate on follower side
            const Matrix< DDRMat > tFollowerdVgp2deta2 =
                    trans( tExtractMat * ( tFollowerdYgp2dq2 + tFollowerdVgp2dq2 ) )         //
                            * vectorize( trans( tFollowerdQgpdeta ) * tFollowerdQgpdeta )    //
                    + trans( tFollowerdQgp2deta2 * ( tFollowerdYgpdq + tFollowerdVgpdq ) );

            Matrix< DDRMat > tFollowerdResRay2detadv( tSpaceDim, tNumDofs );
            for ( uint in = 0; in < tNumNodes; ++in )
            {
                tFollowerdResRay2detadv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                        tFollowerdNVgpdeta( in ) * tIdentity;
            }

            // compute the second derivatives of gap and eta wrt follower displacement dofs
            Matrix< DDRMat > tdGap2dv2( tNumDofs, tNumDofs );
            Matrix< DDRMat > tdEta2dv2( tNumDofs, tNumDofs );

            for ( uint idof = 0; idof < tNumDofs; idof++ )
            {
                for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
                {
                    Matrix< DDRMat > tdSolRay2dv2 =
                            tJacRayInv * ( tdEtadv( 0, idof ) * tFollowerdResRay2detadv.get_column( jdof ) +    //
                                           tdEtadv( 0, jdof ) * tFollowerdResRay2detadv.get_column( idof ) +    //
                                           tFollowerdVgp2deta2 * tdEtadv( 0, idof ) * tdEtadv( 0, jdof ) );

                    tdGap2dv2( idof, jdof ) = tdSolRay2dv2( 0 );
                    tdEta2dv2( idof, jdof ) = tdSolRay2dv2( 1 );    // FIXME: this is only valid for 2D line elements
                }
            }

            // compute the second derivatives of gap and eta wrt leader displacement dofs
            Matrix< DDRMat > tdGap2du2( tNumDofs, tNumDofs );
            Matrix< DDRMat > tdEta2du2( tNumDofs, tNumDofs );
            Matrix< DDRMat > tdGap2duv( tNumDofs, tNumDofs );
            Matrix< DDRMat > tdEta2duv( tNumDofs, tNumDofs );

            for ( uint idof = 0; idof < tNumDofs; idof++ )
            {
                for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
                {
                    Matrix< DDRMat > tdSolRay2du2 =
                            tJacRayInv * ( -tGap * tLeaderdNormal2dU2.get_column( idof * tNumDofs + jdof )    //
                                           - tLeaderdNormaldU.get_column( idof ) * tdGapdu( jdof )            //
                                           - tLeaderdNormaldU.get_column( jdof ) * tdGapdu( idof )            //
                                           + tFollowerdVgp2deta2 * tdEtadu( 0, idof ) * tdEtadu( 0, jdof ) );

                    tdGap2du2( idof, jdof ) = tdSolRay2du2( 0 );
                    tdEta2du2( idof, jdof ) = tdSolRay2du2( 1 );    // FIXME: this is only valid for 2D line elements

                    tdSolRay2du2 = tJacRayInv * ( -tLeaderdNormaldU.get_column( idof ) * tdGapdv( jdof )               //
                                                  + tFollowerdResRay2detadv.get_column( jdof ) * tdEtadu( 0, idof )    //
                                                  + tFollowerdVgp2deta2 * tdEtadu( 0, idof ) * tdEtadv( 0, jdof ) );

                    tdGap2duv( idof, jdof ) = tdSolRay2du2( 0 );
                    tdEta2duv( idof, jdof ) = tdSolRay2du2( 1 );    // FIXME: this is only valid for 2D line elements
                }
            }

            Matrix< DDRMat > tdGapvec2du2( tSpaceDim, tNumDofs * tNumDofs );
            Matrix< DDRMat > tdGapvec2duv( tSpaceDim, tNumDofs * tNumDofs );
            Matrix< DDRMat > tdGapvec2dv2 = tLeaderNormal * trans( vectorize( tdGap2dv2 ) );

            uint tCounter = 0;
            for ( uint idof = 0; idof < tNumDofs; idof++ )
            {
                for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
                {
                    tdGapvec2du2.get_column( tCounter ) = tGap * tLeaderdNormal2dU2.get_column( tCounter )         //
                                                        + tLeaderdNormaldU.get_column( idof ) * tdGapdu( jdof )    //
                                                        + tLeaderdNormaldU.get_column( jdof ) * tdGapdu( idof )    //
                                                        + tLeaderNormal * tdGap2du2( idof, jdof );

                    tdGapvec2duv.get_column( tCounter ) = tLeaderNormal * tdGap2duv( idof, jdof )    //
                                                        + tLeaderdNormaldU.get_column( idof ) * tdGapdv( jdof );
                    tCounter++;
                }
            }

            aGapData->set_second_order_derivatives(
                    tSpaceDim,
                    tNumNodes,
                    tNumDofs,
                    tdGap2du2,
                    tdGap2dv2,
                    tdGap2duv,
                    tdEta2du2,
                    tdEta2dv2,
                    tdEta2duv,
                    tLeaderdNormal2dU2,
                    tdGapvec2du2,
                    tdGapvec2dv2,
                    tdGapvec2duv );
        }

        // compute derivative follower shape functions
        const Matrix< DDRMat > tFollowerShapeFunctiondu = tFollowerdNVgpdeta * tdEtadu;
        const Matrix< DDRMat > tFollowerShapeFunctiondv = tFollowerdNVgpdeta * tdEtadv;

        // compute derivatives of follower displacements
        const Matrix< DDRMat > tFollowerdDisplacementsdu = trans( tFollowerUHat ) * tFollowerShapeFunctiondu;
        Matrix< DDRMat >       tFollowerdDisplacementsdv = trans( tFollowerUHat ) * tFollowerShapeFunctiondv;

        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tFollowerdDisplacementsdv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) += tFollowerNVgp( in ) * tIdentity;
        }

        // compute tangential velocity on leader side
        Matrix< DDRMat > tFollowerPreviousVgp;
        Matrix< DDRMat > tLeaderTangentProjection;
        Matrix< DDRMat > tLeaderTangentVelocity;
        Matrix< DDRMat > tLeaderdTangentVelocitydu;
        Matrix< DDRMat > tLeaderdTangentVelocitydv;
        Matrix< DDRMat > tLeaderdTangentVelocityduprev;
        Matrix< DDRMat > tLeaderdTangentVelocitydvprev;
        Matrix< DDRMat > tFollowerPreviousdDisplacementsdu;
        Matrix< DDRMat > tFollowerPreviousdDisplacementsdv;

        if ( aLeaderPreviousFieldInterpolatorManager != nullptr && aFollowerPreviousFieldInterpolatorManager != nullptr )
        {
            Field_Interpolator* tLeaderPreviousFieldInterpolator =
                    aLeaderPreviousFieldInterpolatorManager->get_field_interpolators_for_type( aDisplDofTypes( 0 ) );
            Field_Interpolator* tFollowerPreviousFieldInterpolator =
                    aFollowerPreviousFieldInterpolatorManager->get_field_interpolators_for_type( aDisplDofTypes( 0 ) );

            Geometry_Interpolator* tLeaderPreviousIPGI   = tLeaderPreviousFieldInterpolator->get_IG_geometry_interpolator();
            Geometry_Interpolator* tFollowerPreviousIPGI = tFollowerPreviousFieldInterpolator->get_IG_geometry_interpolator();

            aFollowerPreviousFieldInterpolatorManager->set_space_time_from_local_IG_point( tFollowerParamPoint );

            const Matrix< DDRMat >& tFollowerYgp = tFollowerIGGI->valx();

            const Matrix< DDRMat >& tLeaderPreviousUgp = tLeaderPreviousFieldInterpolator->val();

            tFollowerPreviousVgp = tFollowerPreviousFieldInterpolator->val();

            // get the shape functions for the leader and follower at the previous time step and place values with respect to dofs
            // not dof ordering here is different than in fem using FI->N()
            const Matrix< DDRMat >& tLeaderPreviousNUgplocal   = tLeaderPreviousIPGI->NXi();
            const Matrix< DDRMat >& tFollowerPreviousNVgplocal = tFollowerPreviousIPGI->NXi();

            Matrix< DDRMat > tLeaderPreviousNUgp( tSpaceDim, tNumDofs );
            Matrix< DDRMat > tFollowerPreviousNVgp( tSpaceDim, tNumDofs );

            for ( uint in = 0; in < tNumNodes; in++ )
            {
                tLeaderPreviousNUgp( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =    //
                        tLeaderPreviousNUgplocal( in ) * tIdentity;
                tFollowerPreviousNVgp( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =    //
                        tFollowerPreviousNVgplocal( in ) * tIdentity;
            }

            tLeaderTangentProjection = tIdentity - tLeaderNormal * trans( tLeaderNormal );

            Matrix< DDRMat > tLeaderPreviousNormal;
            Matrix< DDRMat > tLeaderPreviousdNormaldU( tSpaceDim, tNumDofs );
            Matrix< DDRMat > tLeaderPreviousdNormal2dU2;
            Matrix< DDRMat > tLeaderPreviousRefNormal;

            GapData::compute_outward_normal_at_gp_for_consistent_deformed_geometry(
                    tLeaderPreviousFieldInterpolator,
                    tLeaderIGGI,
                    tLeaderPreviousNormal,
                    tLeaderPreviousRefNormal,
                    tLeaderPreviousdNormaldU,
                    tLeaderPreviousdNormal2dU2,
                    true,
                    false );

            const Matrix< DDRMat > tVector = trans( tLeaderXgp ) + tLeaderPreviousUgp - trans( tFollowerYgp ) - tFollowerPreviousVgp + tGap * tLeaderPreviousNormal;

            tLeaderTangentVelocity = tLeaderTangentProjection * tVector;

            // compute derivatives of tangential velocity wrt leader and follower displacement dofs
            const Matrix< DDRMat >& tFollowerYgpdu = tFollowerdYgpdeta * tdEtadu;
            const Matrix< DDRMat >& tFollowerYgpdv = tFollowerdYgpdeta * tdEtadv;

            const Matrix< DDRMat >& tFollowerPreviousUHat      = tFollowerPreviousFieldInterpolator->get_coeff();
            const Matrix< DDRMat >  tFollowerPreviousdNVgpdq   = tFollowerdNVgpdq * tFollowerPreviousUHat;
            const Matrix< DDRMat >  tFollowerPreviousdNVgpdeta = trans( tFollowerdQgpdeta * tFollowerPreviousdNVgpdq );

            const Matrix< DDRMat > tFollowerPreviousVgpdu = tFollowerPreviousdNVgpdeta * tdEtadu;
            const Matrix< DDRMat > tFollowerPreviousVgpdv = tFollowerPreviousdNVgpdeta * tdEtadv;

            tLeaderdTangentVelocitydu.set_size( tSpaceDim, tNumDofs, 0.0 );
            tLeaderdTangentVelocitydv.set_size( tSpaceDim, tNumDofs, 0.0 );
            tLeaderdTangentVelocityduprev.set_size( tSpaceDim, tNumDofs, 0.0 );
            tLeaderdTangentVelocitydvprev.set_size( tSpaceDim, tNumDofs, 0.0 );

            // test multiply_leader_dtangentprojectiondu function
            const Matrix< DDRMat > tMultiplyLeaderTangentProjectionDu = -dot( tVector, tLeaderNormal ) * tLeaderdNormaldU    //
                                                                      - tLeaderNormal * trans( tVector ) * tLeaderdNormaldU;

            //            for ( uint idof = 0; idof < tNumDofs; idof++ )
            //            {
            //                const Matrix< DDRMat > tAlternativeMultiplyLeaderTangentProjectionDu =
            //                        -( tLeaderdNormaldU.get_column( idof ) * trans( tLeaderNormal ) + tLeaderNormal * trans( tLeaderdNormaldU.get_column( idof ) ) )    //
            //                        * ( trans( tLeaderXgp ) + tLeaderPreviousUgp - trans( tFollowerYgp ) - tFollowerPreviousVgp + tGap * tLeaderPreviousNormal );
            //
            //                if ( norm( tMultiplyLeaderTangentProjectionDu.get_column( idof ) - tAlternativeMultiplyLeaderTangentProjectionDu ) > 1e-12 )
            //                {
            //                    MORIS_ERROR( false, "Contact_Geometry::remap_nonconformal_rays - multiply_leader_dtangentprojectiondu function is not correct." );
            //                }
            //            }

            for ( uint idof = 0; idof < tNumDofs; idof++ )
            {
                tLeaderdTangentVelocitydu.get_column( idof ) =
                        tMultiplyLeaderTangentProjectionDu.get_column( idof )
                        + tLeaderTangentProjection * ( -tFollowerYgpdu.get_column( idof )             //
                                                       - tFollowerPreviousVgpdu.get_column( idof )    //
                                                       + tdGapdu( idof ) * tLeaderPreviousNormal );

                tLeaderdTangentVelocitydv.get_column( idof ) = tLeaderTangentProjection * (                           //
                                                                       -tFollowerYgpdv.get_column( idof )             //
                                                                       - tFollowerPreviousVgpdv.get_column( idof )    //
                                                                       + tdGapdv( idof ) * tLeaderPreviousNormal );

                tLeaderdTangentVelocityduprev.get_column( idof ) = tLeaderTangentProjection * (                        //
                                                                           tLeaderPreviousNUgp.get_column( idof ) +    //
                                                                           tGap * tLeaderPreviousdNormaldU.get_column( idof ) );

                tLeaderdTangentVelocitydvprev.get_column( idof ) = tLeaderTangentProjection * (    //
                                                                           -tFollowerPreviousNVgp.get_column( idof ) );
            }

            tFollowerPreviousdDisplacementsdu = trans( tFollowerPreviousUHat ) * tFollowerShapeFunctiondu;
            tFollowerPreviousdDisplacementsdv = trans( tFollowerPreviousUHat ) * tFollowerShapeFunctiondv;
        }

        // save results in GapData resorted such that nodes first
        aGapData->mGap                           = tGap;
        aGapData->mEta                           = tEta;
        aGapData->mLeaderNormal                  = tLeaderNormal;
        aGapData->mLeaderRefNormal               = tLeaderRefNormal;
        aGapData->mGapVec                        = tGapVector;
        aGapData->mLeaderTangentProjection       = tLeaderTangentProjection;
        aGapData->mLeaderTangentVelocity         = tLeaderTangentVelocity;
        aGapData->mFollowerDisplacements         = tFollowerVgp;
        aGapData->mFollowerPreviousDisplacements = tFollowerPreviousVgp;

        aGapData->set_first_order_derivatives(
                tSpaceDim,
                tNumNodes,
                tdGapdu,
                tdGapdv,
                tdEtadu,
                tdEtadv,
                tLeaderdNormaldU,
                tdGapvecdu,
                tdGapvecdv,
                tLeaderdTangentVelocitydu,
                tLeaderdTangentVelocitydv,
                tLeaderdTangentVelocityduprev,
                tLeaderdTangentVelocitydvprev,
                tFollowerdDisplacementsdu,
                tFollowerdDisplacementsdv,
                tFollowerPreviousdDisplacementsdu,
                tFollowerPreviousdDisplacementsdv,
                tFollowerShapeFunctiondu,
                tFollowerShapeFunctiondv );

        // gap data evaluation flag to false, i.e. does not need to be recomputed
        aGapData->mEval = false;

        // return parametric point on follower side
        return tFollowerParamPoint;
    }
}    // namespace moris::fem
