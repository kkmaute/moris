/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#ifndef MORIS_CL_FEM_CONTACT_GAP_DATA_HPP
#define MORIS_CL_FEM_CONTACT_GAP_DATA_HPP

#include <memory>
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "fn_dot.hpp"
#include "fn_trans.hpp"

namespace moris::fem
{
    class Field_Interpolator;
    class Geometry_Interpolator;

    /** Point-local contact geometry and its displacement derivatives.
     * The owning IWG controls evaluation invalidation and contact/stick locking
     * during finite differences. Derivatives use field/component-major ordering;
     * the setters convert the mapping calculation's node-major ordering.
     */
    struct GapData
    {
        bool mEval = true;    // The owning IWG invalidates all geometry orders with this flag.
        // A populated Hessian is valid only while mEval is false.
        bool mHasSecondOrderDerivatives = false;

        bool mIsLocked  = false;    // flag to indicate if the gap data is locked and should not be updated
        bool mIsContact = false;    // flag to indicate if the point is in contact
        bool mIsStick   = false;    // flag to indicate if the point is in stick

        real             mGap;
        Matrix< DDRMat > mdGapdu;
        Matrix< DDRMat > mdGap2du2;
        Matrix< DDRMat > mdGapdv;
        Matrix< DDRMat > mdGap2dv2;
        Matrix< DDRMat > mdGap2duv;

        Matrix< DDRMat > mEta;
        Matrix< DDRMat > mdEtadu;
        Matrix< DDRMat > mdEta2du2;
        Matrix< DDRMat > mdEtadv;
        Matrix< DDRMat > mdEta2dv2;
        Matrix< DDRMat > mdEta2duv;

        Matrix< DDRMat > mLeaderNormal;
        Matrix< DDRMat > mLeaderRefNormal;
        Matrix< DDRMat > mLeaderdNormaldu;
        Matrix< DDRMat > mLeaderdNormal2du2;

        Matrix< DDRMat > mGapVec;
        Matrix< DDRMat > mdGapvecdu;
        Matrix< DDRMat > mdGapvecdv;
        Matrix< DDRMat > mdGapvec2du2;
        Matrix< DDRMat > mdGapvec2dv2;
        Matrix< DDRMat > mdGapvec2duv;

        Matrix< DDRMat > mLeaderTangentProjection;
        Matrix< DDRMat > mLeaderTangentVelocity;
        Matrix< DDRMat > mLeaderdTangentVelocitydu;
        Matrix< DDRMat > mLeaderdTangentVelocitydv;
        Matrix< DDRMat > mLeaderdTangentVelocityduprev;
        Matrix< DDRMat > mLeaderdTangentVelocitydvprev;

        Matrix< DDRMat > mFollowerDisplacements;
        Matrix< DDRMat > mFollowerPreviousDisplacements;

        Matrix< DDRMat > mFollowerdDisplacementsdu;
        Matrix< DDRMat > mFollowerdDisplacementsdv;

        Matrix< DDRMat > mFollowerPreviousdDisplacementsdu;
        Matrix< DDRMat > mFollowerPreviousdDisplacementsdv;

        Matrix< DDRMat > mFollowerShapeFunctiondu;
        Matrix< DDRMat > mFollowerShapeFunctiondv;

        //----------------------------------------------------------------------------

        bool is_contact( const bool aIsContact )
        {
            if ( !mIsLocked )
            {
                mIsContact = aIsContact;
            }
            return mIsContact;
        }

        bool is_stick( const bool aIsStick )
        {
            if ( !mIsLocked )
            {
                mIsStick = aIsStick;
            }
            return mIsStick;
        }

        //----------------------------------------------------------------------------

        void copy( const std::unique_ptr< GapData >& tGapData )
        {
            mHasSecondOrderDerivatives = tGapData->mHasSecondOrderDerivatives;

            mGap      = tGapData->mGap;
            mdGapdu   = tGapData->mdGapdu;
            mdGap2du2 = tGapData->mdGap2du2;
            mdGapdv   = tGapData->mdGapdv;
            mdGap2dv2 = tGapData->mdGap2dv2;
            mdGap2duv = tGapData->mdGap2duv;

            mEta      = tGapData->mEta;
            mdEtadu   = tGapData->mdEtadu;
            mdEta2du2 = tGapData->mdEta2du2;
            mdEtadv   = tGapData->mdEtadv;
            mdEta2dv2 = tGapData->mdEta2dv2;
            mdEta2duv = tGapData->mdEta2duv;

            mLeaderNormal      = tGapData->mLeaderNormal;
            mLeaderRefNormal   = tGapData->mLeaderRefNormal;
            mLeaderdNormaldu   = tGapData->mLeaderdNormaldu;
            mLeaderdNormal2du2 = tGapData->mLeaderdNormal2du2;

            mGapVec      = tGapData->mGapVec;
            mdGapvecdu   = tGapData->mdGapvecdu;
            mdGapvecdv   = tGapData->mdGapvecdv;
            mdGapvec2du2 = tGapData->mdGapvec2du2;
            mdGapvec2dv2 = tGapData->mdGapvec2dv2;
            mdGapvec2duv = tGapData->mdGapvec2duv;

            mLeaderTangentProjection = tGapData->mLeaderTangentProjection;

            mLeaderTangentVelocity        = tGapData->mLeaderTangentVelocity;
            mLeaderdTangentVelocitydu     = tGapData->mLeaderdTangentVelocitydu;
            mLeaderdTangentVelocitydv     = tGapData->mLeaderdTangentVelocitydv;
            mLeaderdTangentVelocityduprev = tGapData->mLeaderdTangentVelocityduprev;
            mLeaderdTangentVelocitydvprev = tGapData->mLeaderdTangentVelocitydvprev;

            mFollowerDisplacements         = tGapData->mFollowerDisplacements;
            mFollowerPreviousDisplacements = tGapData->mFollowerPreviousDisplacements;

            mFollowerdDisplacementsdu = tGapData->mFollowerdDisplacementsdu;
            mFollowerdDisplacementsdv = tGapData->mFollowerdDisplacementsdv;

            mFollowerPreviousdDisplacementsdu = tGapData->mFollowerPreviousdDisplacementsdu;
            mFollowerPreviousdDisplacementsdv = tGapData->mFollowerPreviousdDisplacementsdv;

            mFollowerShapeFunctiondu = tGapData->mFollowerShapeFunctiondu;
            mFollowerShapeFunctiondv = tGapData->mFollowerShapeFunctiondv;
        }

        //----------------------------------------------------------------------------

        void set_matrix_sizes( const uint aSpaceDim, const uint aNumDofs, const uint aNumNodes, bool aEvaluateSecondDerivatives = true )
        {
            mdGapdu.set_size( 1, aNumDofs );
            mdGapdv.set_size( 1, aNumDofs );
            mdEtadu.set_size( 1, aNumDofs );
            mdEtadv.set_size( 1, aNumDofs );

            mLeaderdNormaldu.set_size( aSpaceDim, aNumDofs );
            mdGapvecdu.set_size( aSpaceDim, aNumDofs );
            mdGapvecdv.set_size( aSpaceDim, aNumDofs );

            if ( aEvaluateSecondDerivatives )
            {
                set_second_order_matrix_sizes( aSpaceDim, aNumDofs );
            }
            else
            {
                clear_second_order_derivatives();
            }

            mLeaderdTangentVelocitydu.set_size( aSpaceDim, aNumDofs );
            mLeaderdTangentVelocitydv.set_size( aSpaceDim, aNumDofs );
            mLeaderdTangentVelocityduprev.set_size( aSpaceDim, aNumDofs );
            mLeaderdTangentVelocitydvprev.set_size( aSpaceDim, aNumDofs );

            mFollowerdDisplacementsdu.set_size( aSpaceDim, aNumDofs );
            mFollowerdDisplacementsdv.set_size( aSpaceDim, aNumDofs );

            mFollowerPreviousdDisplacementsdu.set_size( aSpaceDim, aNumDofs );
            mFollowerPreviousdDisplacementsdv.set_size( aSpaceDim, aNumDofs );

            mFollowerShapeFunctiondu.set_size( aNumNodes, aNumDofs );
            mFollowerShapeFunctiondv.set_size( aNumNodes, aNumDofs );
        }

        //----------------------------------------------------------------------------

        void clear_second_order_derivatives()
        {
            mHasSecondOrderDerivatives = false;
            mdGap2du2.set_size( 0, 0 );
            mdGap2dv2.set_size( 0, 0 );
            mdGap2duv.set_size( 0, 0 );
            mdEta2du2.set_size( 0, 0 );
            mdEta2dv2.set_size( 0, 0 );
            mdEta2duv.set_size( 0, 0 );
            mLeaderdNormal2du2.set_size( 0, 0 );
            mdGapvec2du2.set_size( 0, 0 );
            mdGapvec2dv2.set_size( 0, 0 );
            mdGapvec2duv.set_size( 0, 0 );
        }

        //----------------------------------------------------------------------------

        void set_second_order_matrix_sizes( const uint aSpaceDim, const uint aNumDofs )
        {
            mHasSecondOrderDerivatives = false;
            mdGap2du2.set_size( aNumDofs, aNumDofs );
            mdGap2dv2.set_size( aNumDofs, aNumDofs );
            mdGap2duv.set_size( aNumDofs, aNumDofs );
            mdEta2du2.set_size( aNumDofs, aNumDofs );
            mdEta2dv2.set_size( aNumDofs, aNumDofs );
            mdEta2duv.set_size( aNumDofs, aNumDofs );

            mLeaderdNormal2du2.set_size( aSpaceDim, aNumDofs * aNumDofs );
            mdGapvec2du2.set_size( aSpaceDim, aNumDofs * aNumDofs );
            mdGapvec2dv2.set_size( aSpaceDim, aNumDofs * aNumDofs );
            mdGapvec2duv.set_size( aSpaceDim, aNumDofs * aNumDofs );
        }

        //----------------------------------------------------------------------------

        void set_first_order_derivatives(
                const uint              aSpaceDim,
                const uint              aNumNodes,
                const Matrix< DDRMat >& adGapdu,
                const Matrix< DDRMat >& adGapdv,
                const Matrix< DDRMat >& adEtadu,
                const Matrix< DDRMat >& adEtadv,
                const Matrix< DDRMat >& aLeaderdNormaldU,
                const Matrix< DDRMat >& adGapvecdu,
                const Matrix< DDRMat >& adGapvecdv,
                const Matrix< DDRMat >& aLeaderdTangentVelocitydu,
                const Matrix< DDRMat >& aLeaderdTangentVelocitydv,
                const Matrix< DDRMat >& aLeaderdTangentVelocityduprev,
                const Matrix< DDRMat >& aLeaderdTangentVelocitydvprev,
                const Matrix< DDRMat >& aFollowerdDisplacementsdu,
                const Matrix< DDRMat >& aFollowerdDisplacementsdv,
                const Matrix< DDRMat >& aFollowerPreviousdDisplacementsdu,
                const Matrix< DDRMat >& aFollowerPreviousdDisplacementsdv,
                const Matrix< DDRMat >& aFollowerShapeFunctiondu,
                const Matrix< DDRMat >& aFollowerShapeFunctiondv )
        {
            uint tIcounter = 0;

            for ( uint idim = 0; idim < aSpaceDim; idim++ )
            {
                for ( uint in = 0; in < aNumNodes; in++ )
                {
                    mdGapdu( tIcounter ) = adGapdu( idim + in * aSpaceDim );
                    mdGapdv( tIcounter ) = adGapdv( idim + in * aSpaceDim );
                    mdEtadu( tIcounter ) = adEtadu( idim + in * aSpaceDim );
                    mdEtadv( tIcounter ) = adEtadv( idim + in * aSpaceDim );

                    mLeaderdNormaldu.get_column( tIcounter ) = aLeaderdNormaldU.get_column( idim + in * aSpaceDim );

                    mdGapvecdu.get_column( tIcounter ) = adGapvecdu.get_column( idim + in * aSpaceDim );
                    mdGapvecdv.get_column( tIcounter ) = adGapvecdv.get_column( idim + in * aSpaceDim );

                    if ( aLeaderdTangentVelocitydu.numel() > 0 )
                    {
                        mLeaderdTangentVelocitydu.get_column( tIcounter )     = aLeaderdTangentVelocitydu.get_column( idim + in * aSpaceDim );
                        mLeaderdTangentVelocitydv.get_column( tIcounter )     = aLeaderdTangentVelocitydv.get_column( idim + in * aSpaceDim );
                        mLeaderdTangentVelocityduprev.get_column( tIcounter ) = aLeaderdTangentVelocityduprev.get_column( idim + in * aSpaceDim );
                        mLeaderdTangentVelocitydvprev.get_column( tIcounter ) = aLeaderdTangentVelocitydvprev.get_column( idim + in * aSpaceDim );
                    }

                    mFollowerdDisplacementsdu.get_column( tIcounter ) = aFollowerdDisplacementsdu.get_column( idim + in * aSpaceDim );
                    mFollowerdDisplacementsdv.get_column( tIcounter ) = aFollowerdDisplacementsdv.get_column( idim + in * aSpaceDim );

                    if ( aFollowerPreviousdDisplacementsdu.numel() > 0 )
                    {
                        mFollowerPreviousdDisplacementsdu.get_column( tIcounter ) = aFollowerPreviousdDisplacementsdu.get_column( idim + in * aSpaceDim );
                        mFollowerPreviousdDisplacementsdv.get_column( tIcounter ) = aFollowerPreviousdDisplacementsdv.get_column( idim + in * aSpaceDim );
                    }

                    mFollowerShapeFunctiondu.get_column( tIcounter ) = aFollowerShapeFunctiondu.get_column( idim + in * aSpaceDim );
                    mFollowerShapeFunctiondv.get_column( tIcounter ) = aFollowerShapeFunctiondv.get_column( idim + in * aSpaceDim );

                    tIcounter++;
                }
            }
        }

        //----------------------------------------------------------------------------

        void set_second_order_derivatives(
                const uint              aSpaceDim,
                const uint              tNumNodes,
                const uint              aNumDofs,
                const Matrix< DDRMat >& tdGap2du2,
                const Matrix< DDRMat >& tdGap2dv2,
                const Matrix< DDRMat >& tdGap2duv,
                const Matrix< DDRMat >& tdEta2du2,
                const Matrix< DDRMat >& tdEta2dv2,
                const Matrix< DDRMat >& tdEta2duv,
                const Matrix< DDRMat >& tLeaderdNormal2dU2,
                const Matrix< DDRMat >& tdGapvec2du2,
                const Matrix< DDRMat >& tdGapvec2dv2,
                const Matrix< DDRMat >& tdGapvec2duv )
        {
            uint tIcounter = 0;
            uint tJcounter = 0;

            for ( uint idim = 0; idim < aSpaceDim; idim++ )
            {
                for ( uint in = 0; in < tNumNodes; in++ )
                {
                    tJcounter = 0;
                    for ( uint jdim = 0; jdim < aSpaceDim; jdim++ )
                    {
                        for ( uint jn = 0; jn < tNumNodes; jn++ )
                        {
                            mdGap2du2( tIcounter, tJcounter ) = tdGap2du2( idim + in * aSpaceDim, jdim + jn * aSpaceDim );
                            mdGap2dv2( tIcounter, tJcounter ) = tdGap2dv2( idim + in * aSpaceDim, jdim + jn * aSpaceDim );
                            mdGap2duv( tIcounter, tJcounter ) = tdGap2duv( idim + in * aSpaceDim, jdim + jn * aSpaceDim );
                            mdEta2du2( tIcounter, tJcounter ) = tdEta2du2( idim + in * aSpaceDim, jdim + jn * aSpaceDim );
                            mdEta2dv2( tIcounter, tJcounter ) = tdEta2dv2( idim + in * aSpaceDim, jdim + jn * aSpaceDim );
                            mdEta2duv( tIcounter, tJcounter ) = tdEta2duv( idim + in * aSpaceDim, jdim + jn * aSpaceDim );

                            mLeaderdNormal2du2.get_column( tIcounter * aNumDofs + tJcounter ) =
                                    tLeaderdNormal2dU2.get_column( ( idim + in * aSpaceDim ) * aNumDofs + jdim + jn * aSpaceDim );
                            mdGapvec2du2.get_column( tIcounter * aNumDofs + tJcounter ) =
                                    tdGapvec2du2.get_column( ( idim + in * aSpaceDim ) * aNumDofs + jdim + jn * aSpaceDim );
                            mdGapvec2dv2.get_column( tIcounter * aNumDofs + tJcounter ) =
                                    tdGapvec2dv2.get_column( ( idim + in * aSpaceDim ) * aNumDofs + jdim + jn * aSpaceDim );
                            mdGapvec2duv.get_column( tIcounter * aNumDofs + tJcounter ) =
                                    tdGapvec2duv.get_column( ( idim + in * aSpaceDim ) * aNumDofs + jdim + jn * aSpaceDim );
                            tJcounter++;
                        }
                    }
                    tIcounter++;
                }
            }
            mHasSecondOrderDerivatives = true;
        }

        //----------------------------------------------------------------------------

        Matrix< DDRMat > multiply_leader_dnormal2du2( const Matrix< DDRMat >& aVector )
        {
            uint tSpaceDim = mLeaderdNormal2du2.n_rows();
            uint tNumDofs  = mdGapdu.numel();

            MORIS_ASSERT( aVector.n_rows() == tSpaceDim,
                    "GapData::multiply_leader_dnormal2du2 - spatial dimensions do not match." );

            Matrix< DDRMat > tResult( tNumDofs, tNumDofs );

            const Matrix< DDRMat > tVectorTrans = trans( aVector );

            for ( uint i = 0; i < tNumDofs; i++ )
            {
                tResult.get_row( i ) = tVectorTrans * mLeaderdNormal2du2( { 0, tSpaceDim - 1 }, { i * tNumDofs, ( i + 1 ) * tNumDofs - 1 } );
            }

            return tResult;
        }

        //----------------------------------------------------------------------------

        Matrix< DDRMat > multiply_leader_dtangentprojectiondu( const Matrix< DDRMat >& tVector )
        {
            MORIS_ASSERT( tVector.n_rows() == mLeaderNormal.n_rows(),
                    "GapData::multiply_leader_dtangentprojectiondu - dimension mismatch." );

            Matrix< DDRMat > tResult = -dot( tVector, mLeaderNormal ) * mLeaderdNormaldu    //
                                     - mLeaderNormal * trans( tVector ) * mLeaderdNormaldu;

            return tResult;
        }

        //----------------------------------------------------------------------------

        Matrix< DDRMat > multiply_follower_dshapefunctiondu( const Matrix< DDRMat >& tVector )
        {
            uint tSpaceDim = tVector.n_rows();
            uint tNumNodes = mFollowerShapeFunctiondu.n_rows();
            uint tNumDofs  = mFollowerShapeFunctiondu.n_cols();

            MORIS_ASSERT( tSpaceDim * tNumNodes == tNumDofs,
                    "GapData::multiply_follower_dshapefunctiondu - dimension mismatch." );

            Matrix< DDRMat > tResult( tNumDofs, tNumDofs );

            for ( uint idim = 0; idim < tSpaceDim; idim++ )
            {
                tResult( { idim * tNumNodes, ( idim + 1 ) * tNumNodes - 1 }, { 0, tNumDofs - 1 } ) =
                        tVector( idim ) * mFollowerShapeFunctiondu;
            }

            return tResult;
        }

        //----------------------------------------------------------------------------

        Matrix< DDRMat > multiply_follower_dshapefunctiondv( const Matrix< DDRMat >& tVector )
        {
            uint tSpaceDim = tVector.n_rows();
            uint tNumNodes = mFollowerShapeFunctiondv.n_rows();
            uint tNumDofs  = mFollowerShapeFunctiondv.n_cols();

            MORIS_ASSERT( tSpaceDim * tNumNodes == tNumDofs,
                    "GapData::multiply_follower_dshapefunctiondv - dimension mismatch." );

            Matrix< DDRMat > tResult( tNumDofs, tNumDofs );

            for ( uint idim = 0; idim < tSpaceDim; idim++ )
            {
                tResult( { idim * tNumNodes, ( idim + 1 ) * tNumNodes - 1 }, { 0, tNumDofs - 1 } ) =
                        tVector( idim ) * mFollowerShapeFunctiondv;
            }

            return tResult;
        }

        //----------------------------------------------------------------------------

        static void compute_outward_normal_at_gp_for_linear_deformed_geometry(
                Field_Interpolator*    aLeaderFieldInterpolatorManager,
                Geometry_Interpolator* aLeaderIGGI,
                Matrix< DDRMat >&      aLeaderNormal,
                Matrix< DDRMat >&      aLeaderRefNormal,
                Matrix< DDRMat >&      aLeaderdNormaldU,
                Matrix< DDRMat >&      aLeaderNormal2dU2,
                const bool             aEvaluateLinearization = true );

        //----------------------------------------------------------------------------

        static void compute_outward_normal_at_gp_for_consistent_deformed_geometry(
                Field_Interpolator*    aLeaderFieldInterpolatorManager,
                Geometry_Interpolator* aLeaderIGGI,
                Matrix< DDRMat >&      aLeaderNormal,
                Matrix< DDRMat >&      aLeaderRefNormal,
                Matrix< DDRMat >&      aLeaderdNormaldU,
                Matrix< DDRMat >&      aLeaderNormal2dU2,
                const bool             aEvaluateLinearization     = true,
                const bool             aEvaluateSecondDerivatives = true );

        //----------------------------------------------------------------------------

        static Matrix< DDRMat > compute_tangential_plane_projector( const Matrix< DDRMat >& aNormal );
    };
}    // namespace moris::fem

#endif
