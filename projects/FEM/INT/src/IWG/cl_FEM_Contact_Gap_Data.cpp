/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include "cl_FEM_Contact_Gap_Data.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include <cmath>

namespace moris::fem
{
    void GapData::compute_outward_normal_at_gp_for_linear_deformed_geometry(
            Field_Interpolator*    aLeaderFieldInterpolator,
            Geometry_Interpolator* aLeaderIGGI,
            Matrix< DDRMat >&      aLeaderNormal,
            Matrix< DDRMat >&      aLeaderRefNormal,
            Matrix< DDRMat >&      aLeaderdNormaldU,
            Matrix< DDRMat >&      aLeaderdNormal2dU2,
            const bool             aEvaluateLinearization )
    {
        // get the space dimension
        const uint tSpaceDim = aLeaderFieldInterpolator->get_space_dim();

        // get coordinates (global and parametric) of IG nodes
        const Matrix< DDRMat >& tLeaderIGNodes      = aLeaderIGGI->get_space_coeff();
        const Matrix< DDRMat >& tLeaderIGNodesParam = aLeaderIGGI->get_space_param_coeff();

        // get displacements for leader
        const Matrix< DDRMat >& tLeaderUHat = aLeaderFieldInterpolator->get_coeff();

        // determine number of IP nodes and dofs (assume same for leader and follower)
        const uint tNumNodes = tLeaderUHat.n_rows();
        const uint tNumDofs  = tLeaderUHat.numel();

        const uint tNumIGNodes = tLeaderIGNodes.n_rows();

        MORIS_ASSERT( tNumDofs == tNumNodes * tSpaceDim,
                "IWG::remap_nonconformal_rays - Number of dofs does not match number of nodes." );

        // get displacements at IG Nodes of leader
        Matrix< DDRMat > tLeaderUIGNodes( tSpaceDim, tNumIGNodes );
        Matrix< DDRMat > tLeaderNUIGNodes( tSpaceDim, tNumDofs );

        const Matrix< DDRMat > tLeaderParamSpaceTime        = aLeaderFieldInterpolator->get_space_time();
        Matrix< DDRMat >       tLeaderIGNodesParamSpaceTime = tLeaderParamSpaceTime;

        for ( uint in = 0; in < tNumIGNodes; in++ )
        {
            tLeaderIGNodesParamSpaceTime( { 0, tSpaceDim - 1 }, { 0, 0 } ) =
                    trans( tLeaderIGNodesParam( { in, in }, { 0, tSpaceDim - 1 } ) );
            aLeaderFieldInterpolator->set_space_time( tLeaderIGNodesParamSpaceTime );

            tLeaderUIGNodes.get_column( in ) = aLeaderFieldInterpolator->val().matrix_data();

            const Matrix< DDRMat > tLeaderNUIGNodesTilde = aLeaderFieldInterpolator->N().matrix_data();

            tLeaderNUIGNodes( { in, in }, { 0, tNumDofs - 1 } ) = tLeaderNUIGNodesTilde( { 0, 0 }, { 0, tNumDofs - 1 } );
        }
        aLeaderFieldInterpolator->set_space_time( tLeaderParamSpaceTime );

        // compute difference of shape functions at IG nodes
        const Matrix< DDRMat > tLeaderDifNUIGNodes =
                tLeaderNUIGNodes( { 1, 1 }, { 0, tNumDofs - 1 } )    //
                - tLeaderNUIGNodes( { 0, 0 }, { 0, tNumDofs - 1 } );

        // compute position of IG nodes in the deformed configuration
        const Matrix< DDRMat > tLeaderDefIGNodes = trans( tLeaderIGNodes ) + tLeaderUIGNodes;

        //        const Matrix< DDRMat >& tLeaderUgp     = aLeaderFieldInterpolator->val();
        //        const Matrix< DDRMat >& tLeaderNUgp    = tLeaderIPGI->NXi();
        //        const Matrix< DDRMat >& tLeaderdNUgpdr = tLeaderIPGI->dNdXi();
        //
        //        // compute derivatives along IG element side
        //        const Matrix< DDRMat > tLeaderdXgpdxi = trans( tLeaderdNXgpdxi * tLeaderIGNodes );
        //        const Matrix< DDRMat > tLeaderRgp     = tLeaderNXgp * tLeaderIGNodesParam;
        //        const Matrix< DDRMat > tLeaderdRgpdxi = tLeaderdNXgpdxi * tLeaderIGNodesParam;
        //        const Matrix< DDRMat > tLeaderdUgpdr  = tLeaderdNUgpdr * tLeaderUHat;
        //
        //        const Matrix< DDRMat > tLeaderdUgpdxi  = trans( tLeaderdRgpdxi * tLeaderdUgpdr );
        //        const Matrix< DDRMat > tLeaderdNUgpdxi = trans( tLeaderdRgpdxi * tLeaderdNUgpdr );


        // compute outwards pointing normal vector at the leader quadrature point
        const Matrix< DDRMat > RotMat             = { { 0, 1 }, { -1, 0 } };    // FIXME: this is only valid for 2D line elements
        const Matrix< DDRMat > tLeaderNormalTilde = RotMat * ( tLeaderDefIGNodes.get_column( 1 ) - tLeaderDefIGNodes.get_column( 0 ) );

        aLeaderRefNormal = RotMat * trans( tLeaderIGNodes.get_row( 1 ) - tLeaderIGNodes.get_row( 0 ) );
        aLeaderRefNormal = 1.0 / norm( aLeaderRefNormal ) * aLeaderRefNormal;

        const real tNtildeSquared = dot( tLeaderNormalTilde, tLeaderNormalTilde );
        const real tNtildeOm12    = 1.0 / std::sqrt( tNtildeSquared );
        const real tNtildeOm32    = std::pow( tNtildeSquared, -1.5 );
        const real tNtilde3Om52   = 3.0 * std::pow( tNtildeSquared, -2.5 );

        aLeaderNormal = tNtildeOm12 * tLeaderNormalTilde;

        if ( aEvaluateLinearization )
        {
            // compute first derivative of the normal vector with respect to the leader dofs
            const Matrix< DDRMat > tIdentity = eye( tSpaceDim, tSpaceDim );

            const Matrix< DDRMat > tPreMultiply =
                    ( tNtildeOm12 * tIdentity - tNtildeOm32 * tLeaderNormalTilde * trans( tLeaderNormalTilde ) ) * RotMat;

            //        aLeaderdNormaldU =                                                           //
            //                tPreMultiply * ( tLeaderNUIGNodes( { tSpaceDim, tSpaceDim + 1 }, { 0, tNumDofs - 1 } ) -    //
            //                                 tLeaderNUIGNodes( { 0, tSpaceDim - 1 }, { 0, tNumDofs - 1 } ) );

            // Fixme: see for better solution above but requires reordering of the dof indices
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                aLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                        tLeaderDifNUIGNodes( in ) * tPreMultiply * tIdentity;
            }

            // compute second derivative of the normal vector with respect to the leader dofs
            const Matrix< DDRMat > tLeaderNormalAux = trans( RotMat ) * tLeaderNormalTilde;

            uint tCounter = 0;
            for ( uint idim = 0; idim < tSpaceDim; idim++ )
            {
                for ( uint in = 0; in < tNumNodes; in++ )
                {
                    for ( uint id = 0; id < tSpaceDim; id++ )
                    {
                        for ( uint jn = 0; jn < tNumNodes; jn++ )
                        {
                            for ( uint jd = 0; jd < tSpaceDim; jd++ )
                            {
                                aLeaderdNormal2dU2( idim, tCounter ) =
                                        ( tNtilde3Om52 * tLeaderNormalTilde( idim ) * tLeaderNormalAux( jd ) * tLeaderNormalAux( id )    //
                                                - tNtildeOm32 * ( RotMat( idim, id ) * tLeaderNormalAux( jd )                            //
                                                                  + RotMat( idim, jd ) * tLeaderNormalAux( id )                          //
                                                                  + tLeaderNormalTilde( idim ) * tIdentity( id, jd ) ) )                 //
                                        * tLeaderDifNUIGNodes( in ) * tLeaderDifNUIGNodes( jn );
                                tCounter++;
                            }
                        }
                    }
                }
                tCounter = 0;
            }
        }
    }

    //------------------------------------------------------------------------------

    void GapData::compute_outward_normal_at_gp_for_consistent_deformed_geometry(
            Field_Interpolator*    aLeaderFieldInterpolator,
            Geometry_Interpolator* aLeaderIGGI,
            Matrix< DDRMat >&      aLeaderNormal,
            Matrix< DDRMat >&      aLeaderRefNormal,
            Matrix< DDRMat >&      aLeaderdNormaldU,
            Matrix< DDRMat >&      aLeaderdNormal2dU2,
            const bool             aEvaluateLinearization )
    {
        // get geometry and displacements for leader quadrature point
        // const Matrix< DDRMat >& tLeaderXgp      = tLeaderIGGI->valx(); // global coordinates at quadrature points
        const Matrix< DDRMat >& tLeaderNXgp     = aLeaderIGGI->NXi();    // shape functions at quadrature points
        const Matrix< DDRMat >& tLeaderdNXgpdxi = aLeaderIGGI->dNdXi();

        // get IP geometry interpolator for leader
        Geometry_Interpolator* tLeaderIPGI = aLeaderFieldInterpolator->get_IG_geometry_interpolator();
        // const Matrix< DDRMat >& tLeaderUgp     = tLeaderFieldInterpolator->val();
        // const Matrix< DDRMat >& tLeaderNUgp    = tLeaderIPGI->NXi();
        const Matrix< DDRMat >& tLeaderdNUgpdr = tLeaderIPGI->dNdXi();

        // get coordinates (global and parametric) of IG nodes
        const Matrix< DDRMat >& tLeaderIGNodes      = aLeaderIGGI->get_space_coeff();
        const Matrix< DDRMat >& tLeaderIGNodesParam = aLeaderIGGI->get_space_param_coeff();

        // get displacements for leader
        const Matrix< DDRMat >& tLeaderUHat = aLeaderFieldInterpolator->get_coeff();

        // compute derivatives along IG element side
        const Matrix< DDRMat > tLeaderdXgpdxi = trans( tLeaderdNXgpdxi * tLeaderIGNodes );
        const Matrix< DDRMat > tLeaderRgp     = tLeaderNXgp * tLeaderIGNodesParam;
        const Matrix< DDRMat > tLeaderdRgpdxi = tLeaderdNXgpdxi * tLeaderIGNodesParam;
        const Matrix< DDRMat > tLeaderdUgpdr  = tLeaderdNUgpdr * tLeaderUHat;

        const Matrix< DDRMat > tLeaderdUgpdxi  = trans( tLeaderdRgpdxi * tLeaderdUgpdr );
        const Matrix< DDRMat > tLeaderdNUgpdxi = trans( tLeaderdRgpdxi * tLeaderdNUgpdr );

        // compute outwards pointing normal vector at the leader quadrature point
        const Matrix< DDRMat > RotMat             = { { 0, 1 }, { -1, 0 } };    // FIXME: this is only valid for 2D line elements
        const Matrix< DDRMat > tLeaderNormalTilde = RotMat * ( tLeaderdXgpdxi + tLeaderdUgpdxi );

        aLeaderRefNormal = RotMat * ( tLeaderdXgpdxi );
        aLeaderRefNormal = 1.0 / norm( aLeaderRefNormal ) * aLeaderRefNormal;

        const real tNtildeSquared = dot( tLeaderNormalTilde, tLeaderNormalTilde );
        const real tNtildeOm12    = 1.0 / std::sqrt( tNtildeSquared );
        const real tNtildeOm32    = std::pow( tNtildeSquared, -1.5 );
        const real tNtilde3Om52   = 3.0 * std::pow( tNtildeSquared, -2.5 );

        aLeaderNormal = tNtildeOm12 * tLeaderNormalTilde;

        if ( aEvaluateLinearization )
        {
            // get the space dimension
            const uint tSpaceDim = aLeaderFieldInterpolator->get_space_dim();

            // determine number of IP nodes
            const uint tNumNodes = tLeaderUHat.n_rows();

            // compute first derivative of the normal vector with respect to the leader dofs
            const Matrix< DDRMat > tIdentity = eye( tSpaceDim, tSpaceDim );

            const Matrix< DDRMat > tPreMultiply =
                    ( tNtildeOm12 * tIdentity - tNtildeOm32 * tLeaderNormalTilde * trans( tLeaderNormalTilde ) ) * RotMat;

            for ( uint in = 0; in < tNumNodes; in++ )
            {
                aLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                        tLeaderdNUgpdxi( in ) * tPreMultiply * tIdentity;
            }

            // compute second derivative of the normal vector with respect to the leader dofs
            const Matrix< DDRMat > tLeaderNormalAux = trans( RotMat ) * tLeaderNormalTilde;

            // std::cout << "aLeaderdNormaldU.n_cols(): " << aLeaderdNormaldU.n_cols() << std::endl;
            // std::cout << "aLeaderdNormal2dU2.n_cols(): " << aLeaderdNormal2dU2.n_cols() << std::endl;
            // if ( aLeaderdNormaldU.n_cols() == 0 or aLeaderdNormal2dU2.n_cols() == 0 )
            //   return;

            uint tCounter = 0;
            for ( uint idim = 0; idim < tSpaceDim; idim++ )
            {
                for ( uint in = 0; in < tNumNodes; in++ )
                {
                    for ( uint id = 0; id < tSpaceDim; id++ )
                    {
                        for ( uint jn = 0; jn < tNumNodes; jn++ )
                        {
                            for ( uint jd = 0; jd < tSpaceDim; jd++ )
                            {
                                aLeaderdNormal2dU2( idim, tCounter ) =
                                        ( tNtilde3Om52 * tLeaderNormalTilde( idim ) * tLeaderNormalAux( jd ) * tLeaderNormalAux( id )    //
                                                - tNtildeOm32 * ( RotMat( idim, id ) * tLeaderNormalAux( jd )                            //
                                                                  + RotMat( idim, jd ) * tLeaderNormalAux( id )                          //
                                                                  + tLeaderNormalTilde( idim ) * tIdentity( id, jd ) ) )                 //
                                        * tLeaderdNUgpdxi( in ) * tLeaderdNUgpdxi( jn );
                                tCounter++;
                            }
                        }
                    }
                }
                tCounter = 0;
            }
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > GapData::compute_tangential_plane_projector( const Matrix< DDRMat >& aNormal )
    {
        // compute tangential plane projector
        const uint tSpaceDim = aNormal.numel();

        // check if space dimension is two or three
        MORIS_ERROR( tSpaceDim == 2 || tSpaceDim == 3,
                "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_tangential_plane_projector - Only two and three dimensional cases are implemented." );

        // initialize tangential plane projector with identity matrix
        Matrix< DDRMat > tTangentialPlane = eye( tSpaceDim, tSpaceDim );

        for ( uint i = 0; i < tSpaceDim; i++ )
        {
            for ( uint j = 0; j < tSpaceDim; j++ )
            {
                tTangentialPlane( i, j ) = -aNormal( i ) * aNormal( j );
            }
        }

        return tTangentialPlane;
    }
}    // namespace moris::fem
