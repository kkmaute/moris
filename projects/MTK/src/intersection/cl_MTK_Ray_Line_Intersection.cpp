//
// Created by frank on 1/11/24.
//

#include "cl_MTK_Ray_Line_Intersection.hpp"
#include "armadillo"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_QuadraturePointMapper.hpp"
#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_Logger.hpp"
#include <chrono>
#include <cstdio>

namespace
{
    inline void assert_param_in_bounds_box( moris::Matrix< moris::DDRMat > const & aParam, const char* aContext )
    {
        if ( aParam.n_rows() >= 1 )
        {
            MORIS_ASSERT(
                    aParam( 0 ) >= -1.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                    "MTK Ray-Line parametric out of bounds (LINE domain, param[0]=%e, context=%s)",
                    aParam( 0 ),
                    aContext );
        }
        if ( aParam.n_rows() >= 2 )
        {
            MORIS_ASSERT(
                    aParam( 1 ) >= -1.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK Ray-Line parametric out of bounds (LINE domain, param[1]=%e, context=%s)",
                    aParam( 1 ),
                    aContext );
        }
    }

    inline void assert_param_in_bounds_simplex( moris::Matrix< moris::DDRMat > const & aParam, const char* aContext )
    {
        if ( aParam.n_rows() >= 2 )
        {
            MORIS_ASSERT(
                    aParam( 0 ) >= 0.0 - 1e-12 && aParam( 1 ) >= 0.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12 && aParam( 1 ) <= 1.0 + 1e-12 && aParam( 0 ) + aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK Ray-Line parametric out of bounds (TRI simplex, eta=%e zeta=%e, context=%s)",
                    aParam( 0 ),
                    aParam( 1 ),
                    aContext );
        }
    }
}    // namespace

namespace moris::mtk
{
    // NOT NEEDED ANYMORE
    /**
     * @brief Implementation according to https://stackoverflow.com/a/2932601
     */
    void Ray_Line_Intersection::perform_raytracing()
    {
        // Nomenclature:
        //     mOrigin: r
        //     mDirection: dr
        //        -> p(u) = r + u * dr
        //     mTargetOrigin: s
        //     mTargetSpan: ds = e - s (between line segment vertices s and e)
        //        -> q(v) = s + v * ds

        // distance between origins r and s
        // dOrigins = s - r

        Matrix< DDRMat > dOrigins = mTargetOrigin - mOrigin;

        // determinant of the matrix D = [dr, ds]
        MORIS_ASSERT( mSpatialDimension == 2,
                "Ray_Line_Intersection::perform: ray line intersection only implemented in 2D" );

        real const detD = mTargetSpan( 0 ) * mDirection( 1 ) - mTargetSpan( 1 ) * mDirection( 0 );

        // if detD is zero, the lines are parallel and no calculation is necessary
        if ( std::abs( detD ) > 1e-16 )
        {
            // scaling factor u and v for the ray and the line segment, respectively
            real u = NAN;
            real v = NAN;

            // u = (dOrigins_y * ds_x - dOrigins_x * ds_y) / detD
            u = ( dOrigins( 1 ) * mTargetSpan( 0 ) - dOrigins( 0 ) * mTargetSpan( 1 ) ) / detD;

            // v = (dOrigins_y * dr_x - dOrigins_x * dr_y) / detD
            // note: this is negative of matrix solution as positive of target span vector is used while is should be negative
            v = ( dOrigins( 1 ) * mDirection( 0 ) - dOrigins( 0 ) * mDirection( 1 ) ) / detD;

            // if v is between 0 and 1, the intersection point is on the line segment
            if ( v >= 0.0 && v <= 1.0 )
            {
                mHasIntersection = true;

                // calculate the (signed) ray length between the origin and the line segment
                // negative if the intersection point is in the negative direction of the ray
                mSignedRayLength = norm( u * mDirection ) * ( u > 0.0 ? 1.0 : -1.0 );

                // set the factors by which the span of the line segment has to be multiplied to get the intersection point
                // q(v) = s + v * ds
                this->set_intersection_factors( { v } );

                // the parametric coordinate goes from -1 to 1 and has the center in the middle of the line segment
                mParametricCoordinate = 2.0 * ( v - 0.5 );

                // xxxxxxxxxxxxxxxxxxxxx
                //                Matrix< DDRMat > tRefPoint = { { -2.457207927533254e-01 }, { -5.363658433292784e-04 } };
                //                if ( norm( mOrigin - tRefPoint ) < 1e-6 )
                //                {
                //                    Matrix< DDRMat > tJac =    //
                //                            { { mTargetSpan( 0 ), mDirection( 0 ) },
                //                                { mTargetSpan( 1 ), mDirection( 1 ) } };
                //                    print( mOrigin, "mOrigin" );
                //                    print( mTargetOrigin, "mTargetOrigin" );
                //                    print( mDirection, "mDirection" );
                //                    print( mTargetSpan, "mTargetSpan" );
                //                    print( tJac, "tJac" );
                //                    fprintf( stdout, "u = %e, v = %e, mSignedRayLength = %e\n", u, v, mSignedRayLength );
                //                }
                // xxxxxxxxxxxxxxxxxxxxx

                return;    // early exit
            }
        }
        mHasIntersection = false;
    }

    
    void Ray_Line_Intersection::perform_nonlinear_raytracing(
            Space_Interpolator&     tGeomSpaceInterpolator,
            Space_Interpolator&     tFieldSpaceInterpolator,
            Matrix< DDRMat >&       mOrigin,
            Matrix< DDRMat >&       mDirection,
            const Matrix< DDRMat >& tTargetLocalCoordinates )
    {
        // Newton solve for ray/line intersection in nonlinearly deformed configuration.
        int const        tProcRank = moris::Logger::logger_par_rank();
        auto const       tRaytraceStartTime = std::chrono::steady_clock::now();
        uint             tSpaceDim = mOrigin.n_rows();
        Matrix< DDRMat > tSolRay( tSpaceDim, 1, 0.0 );
        Matrix< DDRMat > tSolRayPrev( tSpaceDim, 1, 1.0 );
        real             tRayResRefNorm   = 0.0;
        real             tResRayNorm      = 0.0;
        bool             tNewtonConverged = false;
        uint             tNewtonIters     = 0;
        // Matrix< DDRMat > tTargetParamPoint( tSpaceDim - 1, 1, 0.0 );
        Matrix< DDRMat > tJacRayInv;
        const real       tRayResRefNormTol  = 1e-10;
        const real       tStagnationNormTol = 1e-12;
        const real       tStagnationAcceptTol = 1e-3;
        const uint       tMaxNewtonSteps    = 100;


        for ( uint inew = 0; inew < tMaxNewtonSteps; inew++ )
        {
            tNewtonIters = inew + 1;
            Matrix< DDRMat > tJacRay( tSpaceDim, tSpaceDim );
            if ( tSpaceDim == 2 )
            {
                // tSolRay( 0 ) = std::min( 1.0, std::max( -1.0, tSolRay( 0 ) ) );
                tSolRay( 1 ) = std::min( 1.0, std::max( -1.0, tSolRay( 1 ) ) );    // tSolRay = (s,eta)

                Matrix< DDRMat > tTargetEtagp = { { tSolRay( 1 ) } };
                assert_param_in_bounds_box( tTargetEtagp, "Ray_Line_Intersection::perform_nonlinear_2D_postclamp" );

                // Build a Lagrange interpolator along the target local coordinates on the side
                Interpolation_Rule tFacetInterpRule(
                        Geometry_Type::LINE,
                        Interpolation_Type::LAGRANGE,
                        Interpolation_Order::LINEAR,
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );
                Space_Interpolator tFacetInterpolator( tFacetInterpRule );

                // tTargetLocalCoordinates has rows = nodes, cols = parametric dims; transpose to (dim x nNodes)
                Matrix< DDRMat > tTargetLocalCoordsTrans = trans( tTargetLocalCoordinates );

                tFacetInterpolator.set_space_param_coeff( tTargetLocalCoordsTrans );
                tFacetInterpolator.set_space( tTargetEtagp );

                // Evaluate shape functions and interpolate the parametric coordinate along the side
                Matrix< DDRMat > tTargetNRsgp = trans( tFacetInterpolator.NXi() );
                Matrix< DDRMat > tTargetRsgp  = tTargetLocalCoordsTrans * tTargetNRsgp;
                tGeomSpaceInterpolator.set_space( tTargetRsgp );
                tFieldSpaceInterpolator.set_space( tTargetRsgp );

                // Get geometry and displacements for Target quadrature point
                const Matrix< DDRMat >& tTargetYgp      = tGeomSpaceInterpolator.valx();
                const Matrix< DDRMat >& tTargetdNYgpdrs = tGeomSpaceInterpolator.dNdXi();
                const Matrix< DDRMat >& tTargetVgp      = tFieldSpaceInterpolator.valx();

                // Chain rule for derivatives
                Matrix< DDRMat > tTargetdNrsdeta = trans( tFacetInterpolator.dNdXi() );
                Matrix< DDRMat > tTargetdrsdeta  = tTargetLocalCoordsTrans * tTargetdNrsdeta;

                // Derivatives with respect to r,s
                Matrix< DDRMat > tTargetdYgpdrs = tTargetdNYgpdrs * tGeomSpaceInterpolator.get_space_coeff();
                Matrix< DDRMat > tTargetdVgpdrs = tTargetdNYgpdrs * tFieldSpaceInterpolator.get_space_coeff();
                // Chain rule: dY/dxi = [dY/dr, dY/ds] * dxi/dxi_hat
                Matrix< DDRMat > tTargetdYgpdeta = tTargetdYgpdrs * tTargetdrsdeta;
                Matrix< DDRMat > tTargetdVgpdeta = tTargetdVgpdrs * tTargetdrsdeta;

                // Compute the residual
                Matrix< DDRMat > tResRay = mOrigin - trans( tTargetYgp ) + tSolRay( 0 ) * mDirection - trans( tTargetVgp );

                real tResRayNorm = norm( tResRay );
                if ( inew == 0 ) tRayResRefNorm = tResRayNorm;

                // Compute jacobian and its inverse
                tJacRay( { 0, tSpaceDim - 1 }, { 0, 0 } )             = mDirection.matrix_data();
                tJacRay( { 0, tSpaceDim - 1 }, { 1, tSpaceDim - 1 } ) = -tTargetdVgpdeta - tTargetdYgpdeta;
                if ( std::abs( det( tJacRay ) ) < MORIS_REAL_EPS ) break;
                tJacRayInv = inv( tJacRay );
    
                const real tConvThresh = std::max( tRayResRefNormTol * tRayResRefNorm, tRayResRefNormTol );
                if ( tResRayNorm < tConvThresh )
                {
                    tNewtonConverged = true;
                    break;
                }
               
                if ( norm( tSolRay - tSolRayPrev ) < tStagnationNormTol )
                {
                    if ( tResRayNorm < tStagnationAcceptTol )
                    {
                        tNewtonConverged = true;
                    }
                    break;
                }
                tSolRayPrev = tSolRay;
                tSolRay -= tJacRayInv * tResRay;
            }
            else if ( tSpaceDim == 3 )
            {

                real tEta  = tSolRay( 1 );
                real tZeta = tSolRay( 2 );

                // Project onto the triangular simplex:
                // eta >= 0, zeta >= 0, eta + zeta <= 1
                tEta  = std::max( 0.0, tEta );
                tZeta = std::max( 0.0, tZeta );

                if ( tEta + tZeta > 1.0 )
                {
                    real tSum = tEta + tZeta;
                    tEta  /= tSum;
                    tZeta /= tSum;
                }

                tSolRay( 1 ) = tEta;
                tSolRay( 2 ) = tZeta;
                Matrix< DDRMat > tTargetEtagp = { { tSolRay( 1 ) }, { tSolRay( 2 ) } };    // eta and zeta
                assert_param_in_bounds_simplex( tTargetEtagp, "Ray_Line_Intersection::perform_nonlinear_3D_postclamp" );


                // Build a Lagrange interpolator along the target local coordinates on the side
                Interpolation_Rule tTriangleInterpRule(
                        Geometry_Type::TRI,
                        Interpolation_Type::LAGRANGE,
                        Interpolation_Order::LINEAR,
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );
                Space_Interpolator tTriangleInterpolator( tTriangleInterpRule );

                // tTargetLocalCoordinates has rows = nodes, cols = parametric dims; transpose to (dim x nNodes)
                Matrix< DDRMat > tTargetLocalCoordsTrans = trans( tTargetLocalCoordinates );
                tTriangleInterpolator.set_space_param_coeff( tTargetLocalCoordsTrans );
                tTriangleInterpolator.set_space( tTargetEtagp );

                // Evaluate shape functions and interpolate the parametric coordinate along the side
                Matrix< DDRMat > tTargetNRsgp = trans( tTriangleInterpolator.NXi() );
                Matrix< DDRMat > tTargetRsgp  = tTargetLocalCoordsTrans * tTargetNRsgp;
                tGeomSpaceInterpolator.set_space( tTargetRsgp );
                tFieldSpaceInterpolator.set_space( tTargetRsgp );

                // Get geometry and displacements for Target quadrature point
                const Matrix< DDRMat >& tTargetYgp      = tGeomSpaceInterpolator.valx();
                const Matrix< DDRMat >& tTargetdNYgpdrs = tGeomSpaceInterpolator.dNdXi();
                const Matrix< DDRMat >& tTargetVgp      = tFieldSpaceInterpolator.valx();

                // Chain rule for derivatives
                Matrix< DDRMat > tTargetdNrsdeta = trans( tTriangleInterpolator.dNdXi() );
                Matrix< DDRMat > tTargetdrsdeta  = tTargetLocalCoordsTrans * tTargetdNrsdeta;

                // Derivatives with respect to r,s
                Matrix< DDRMat > tTargetdYgpdrs = tTargetdNYgpdrs * tGeomSpaceInterpolator.get_space_coeff();
                Matrix< DDRMat > tTargetdVgpdrs = tTargetdNYgpdrs * tFieldSpaceInterpolator.get_space_coeff();
                // Chain rule: dY/dxi = [dY/dr, dY/ds] * dxi/dxi_hat
    
                Matrix< DDRMat > tTargetdYgpdeta = trans( tTargetdYgpdrs )  * tTargetdrsdeta;
                Matrix< DDRMat > tTargetdVgpdeta = trans( tTargetdVgpdrs ) * tTargetdrsdeta;

                // Compute the residual
                Matrix< DDRMat > tResRay = mOrigin - trans( tTargetYgp ) + tSolRay( 0 ) * mDirection - trans( tTargetVgp );

                // DEBUG OUTPUT
                // PRINT( tTargetdVgpdeta );
                // PRINT( tTargetdYgpdeta );
                // PRINT( tTargetVgp );
                // PRINT( tTargetYgp );
                // PRINT( tResRay );
                real tResRayNorm = norm( tResRay );
                if ( inew == 0 ) tRayResRefNorm = tResRayNorm;
                // Compute jacobian and its inverse

                tJacRay( { 0, tSpaceDim - 1 }, { 0, 0 } )             = mDirection.matrix_data();
                tJacRay( { 0, tSpaceDim - 1 }, { 1, tSpaceDim - 1 } ) = -tTargetdVgpdeta - tTargetdYgpdeta;
                if ( std::abs( det( tJacRay ) ) < MORIS_REAL_EPS ) break;
                tJacRayInv = inv( tJacRay );
              
                const real tConvThresh = std::max( tRayResRefNormTol * tRayResRefNorm, tRayResRefNormTol );
                if ( tResRayNorm < tConvThresh )
                {
                    tNewtonConverged = true;
                    break;
                }
                
                if ( norm( tSolRay - tSolRayPrev ) < tStagnationNormTol )
                {
                    if ( tResRayNorm < tStagnationAcceptTol )
                    {
                        tNewtonConverged = true;
                    }
                    break;
                }
                tSolRayPrev = tSolRay;
                tSolRay -= tJacRayInv * tResRay;
            }
        }

        bool tPlotRays = true;

        auto const tRaytraceEndTime = std::chrono::steady_clock::now();
        auto const tRaytraceSeconds = std::chrono::duration< double >( tRaytraceEndTime - tRaytraceStartTime ).count();
        std::fprintf( stdout,
                "Rank %d MTK nonlinear raytrace timing converged=%d iters=%u elapsed=%.6f s\n",
                tProcRank,
                tNewtonConverged ? 1 : 0,
                tNewtonIters,
                tRaytraceSeconds );

        // Check if Newton converged
        if ( !tNewtonConverged )
        {
            // return a point with eta set to -2.0
            // mParametricCoordinate.fill( -2.0 );

            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            if ( tPlotRays )
            {
                sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
                Matrix< DDRMat > tRayCastPoint = mOrigin;
                Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * mDirection;

                if ( tSpaceDim == 2 )
                {
                    fprintf( stdout, "MTKNiter = %d Iter = %u NotConverged %e  %e  %e  %e res=%e\n",    //
                            tNiter,
                            tNewtonIters,
                            tRayCastPoint( 0 ),
                            tRayCastPoint( 1 ),
                            tTgtPoint( 0 ),
                            tTgtPoint( 1 ),
                            tResRayNorm );
                }
                else
                {
                    fprintf( stdout, "MTKNiter = %d Iter = %u NotConverged %e  %e  %e  %e  %e  %e res=%e\n",    //
                            tNiter,
                            tNewtonIters,
                            tRayCastPoint( 0 ),
                            tRayCastPoint( 1 ),
                            tRayCastPoint( 2 ),
                            tTgtPoint( 0 ),
                            tTgtPoint( 1 ),
                            tTgtPoint( 2 ),
                            tResRayNorm );
                }
            }
            mHasIntersection = false;
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        }
        else
        {
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            mHasIntersection   = true;
            mRayDirectionParam = tSolRay( 0 ) * mDirection;
            // Store signed scalar ray length (same convention as linear raytrace).
            mSignedRayLength      = norm( mRayDirectionParam ) * ( tSolRay( 0 ) > 0.0 ? 1.0 : -1.0 );
            mParametricCoordinate = tSolRay( { 1, tSpaceDim - 1 }, { 0, 0 } );

            if ( tSpaceDim == 3 )
            {
                assert_param_in_bounds_simplex( mParametricCoordinate, "Ray_Line_Intersection::perform_nonlinear_3D_final" );
            }
            else
            {
                assert_param_in_bounds_box( mParametricCoordinate, "Ray_Line_Intersection::perform_nonlinear_2D_final" );
            }

            if ( tSpaceDim == 3
                    && ( std::abs( mParametricCoordinate( 0 ) ) > 1.2 || std::abs( mParametricCoordinate( 1 ) ) > 1.2 ) )
            {
                std::cout << "RAYLINE_3D_PARAM_OOB"
                          << " s=" << tSolRay( 0 )
                          << " eta=" << mParametricCoordinate( 0 )
                          << " zeta=" << mParametricCoordinate( 1 )
                          << std::endl;
            }
        }
    }


    Matrix< DDRMat > Ray_Line_Intersection::get_intersection_parametric() const
    {
        return mParametricCoordinate;
    }

    real Ray_Line_Intersection::get_signed_ray_length() const
    {
        return mSignedRayLength;
    }
}    // namespace moris::mtk
