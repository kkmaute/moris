/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_
// MRS/CNT/src
#include <utility>
#include <vector>

#include "cl_Vector.hpp"
// LNA/src
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "fn_vectorize.hpp"
#include "fn_isfinite.hpp"
// MRS/COR/src // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Enums.hpp"
#include "fn_FEM_FD_Scheme.hpp"
// FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
// GEN/src
#include "GEN_Data_Types.hpp"
#include <iomanip>

namespace moris::fem
{
    struct GapData
    {
        bool mEval = true;    // flag to indicate if the gap data needs to be evaluate

        real             mGap;
        Matrix< DDRMat > mdGapdu;
        Matrix< DDRMat > mdGap2du2;
        Matrix< DDRMat > mdGapdv;
        Matrix< DDRMat > mdGap2dv2;
        Matrix< DDRMat > mdGap2duv;

        Matrix< DDRMat >                mEta;
        Matrix< DDRMat >                mdEtadu;
        Matrix< DDRMat >                mdEta2du2;
        Matrix< DDRMat >                mdEtadv;
        Matrix< DDRMat >                mdEta2dv2;
        std::vector< Matrix< DDRMat > > mdEta2duv;        // [nEta], each is nDof x nDof 
        Matrix< DDRMat >                mdEta2duv_mat;    // flattened: (nEta * nDof) x nDof for backward compat

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

        //----------------------------------------------------------------------------

        void copy( const std::unique_ptr< GapData >& tGapData )
        {
            mEval     = tGapData->mEval;
            mGap      = tGapData->mGap;
            mdGapdu   = tGapData->mdGapdu;
            mdGap2du2 = tGapData->mdGap2du2;
            mdGapdv   = tGapData->mdGapdv;
            mdGap2dv2 = tGapData->mdGap2dv2;
            mdGap2duv = tGapData->mdGap2duv;

            mEta          = tGapData->mEta;
            mdEtadu       = tGapData->mdEtadu;
            mdEta2du2     = tGapData->mdEta2du2;
            mdEtadv       = tGapData->mdEtadv;
            mdEta2dv2     = tGapData->mdEta2dv2;
            mdEta2duv_mat = tGapData->mdEta2duv_mat;
            mdEta2duv.resize( tGapData->mdEta2duv.size() );
            for ( size_t i = 0; i < tGapData->mdEta2duv.size(); ++i )
            {
                mdEta2duv[ i ] = tGapData->mdEta2duv[ i ];
            }

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
        }

        //----------------------------------------------------------------------------

        void set_matrix_sizes( const uint aSpaceDim, const uint aLeaderNumDofs, const uint aFollowerNumDofs )
        {
            // Basic first-order sizes (gap and eta derivatives use leader DOFs by convention)
            mdGapdu.set_size( 1, aLeaderNumDofs );
            mdGapdv.set_size( 1, aFollowerNumDofs );
            mdEtadu.set_size( aSpaceDim - 1, aLeaderNumDofs );
            mdEtadv.set_size( aSpaceDim - 1, aFollowerNumDofs );

            mLeaderdNormaldu.set_size( aSpaceDim, aLeaderNumDofs );
            mdGapvecdu.set_size( aSpaceDim, aLeaderNumDofs );
            mdGapvecdv.set_size( aSpaceDim, aFollowerNumDofs );

            // For second-order matrices, size explicitly with leader and follower DOF counts
            // gap second derivatives
            mdGap2du2.set_size( aLeaderNumDofs, aLeaderNumDofs );
            mdGap2dv2.set_size( aFollowerNumDofs, aFollowerNumDofs );
            mdGap2duv.set_size( aLeaderNumDofs, aFollowerNumDofs );

            // eta second derivatives: du2 uses leader DOFs, dv2 uses follower DOFs
            mdEta2du2.set_size( ( aSpaceDim - 1 ) * aLeaderNumDofs, aLeaderNumDofs );
            mdEta2dv2.set_size( ( aSpaceDim - 1 ) * aFollowerNumDofs, aFollowerNumDofs );

            // flattened mdEta2duv matrix: (nEta * leaderDofs) x followerDofs
            mdEta2duv_mat.set_size( ( aSpaceDim - 1 ) * aLeaderNumDofs, aFollowerNumDofs );

            // per-eta storage: leader rows, follower cols
            mdEta2duv.resize( aSpaceDim - 1 );
            for ( uint i = 0; i < aSpaceDim - 1; ++i )
            {
                mdEta2duv[ i ].set_size( aLeaderNumDofs, aFollowerNumDofs );
            }

            // leader-normal and gap-vector second derivatives: sizes depend on leader/follower pairing
            mLeaderdNormal2du2.set_size( aSpaceDim, aLeaderNumDofs * aLeaderNumDofs );
            mdGapvec2du2.set_size( aSpaceDim, aLeaderNumDofs * aLeaderNumDofs );
            mdGapvec2dv2.set_size( aSpaceDim, aFollowerNumDofs * aFollowerNumDofs );
            mdGapvec2duv.set_size( aSpaceDim, aLeaderNumDofs * aFollowerNumDofs );
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
        const Matrix< DDRMat >& adGapvecdv )
        {
        uint tIcounter = 0;

        for ( uint idim = 0; idim < aSpaceDim; idim++ )
        {
                for ( uint in = 0; in < aNumNodes; in++ )
                {
                const uint tSrcDof = in * aSpaceDim + idim;

                mdGapdu( tIcounter ) = adGapdu( 0, tSrcDof );
                mdGapdv( tIcounter ) = adGapdv( 0, tSrcDof );
                for ( uint iEta = 0; iEta < aSpaceDim - 1; ++iEta )
                {
                        mdEtadu( iEta, tIcounter ) = adEtadu( iEta, tSrcDof );
                        mdEtadv( iEta, tIcounter ) = adEtadv( iEta, tSrcDof );
                }

                mLeaderdNormaldu( { 0, aSpaceDim - 1 }, { tIcounter, tIcounter } ) =
                        aLeaderdNormaldU( { 0, aSpaceDim - 1 }, { tSrcDof, tSrcDof } );
                mdGapvecdu( { 0, aSpaceDim - 1 }, { tIcounter, tIcounter } ) =
                        adGapvecdu( { 0, aSpaceDim - 1 }, { tSrcDof, tSrcDof } );
                mdGapvecdv( { 0, aSpaceDim - 1 }, { tIcounter, tIcounter } ) =
                        adGapvecdv( { 0, aSpaceDim - 1 }, { tSrcDof, tSrcDof } );

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
                const uint tNumEta = aSpaceDim - 1;

                // leader and follower DOF counts (aNumDofs is leader DOFs)
                const uint tLeaderDofs   = aNumDofs;
                const uint tFollowerDofs = mdEta2duv_mat.n_cols();

                // Sanity checks for expected matrix sizes using explicit leader/follower DOFs
                MORIS_ASSERT( mdEta2du2.n_rows() >= tNumEta * tLeaderDofs,
                        "GapData::set_second_order_derivatives - mdEta2du2 has unexpected number of rows." );
                MORIS_ASSERT( mdEta2du2.n_cols() >= tLeaderDofs,
                        "GapData::set_second_order_derivatives - mdEta2du2 has unexpected number of cols." );

                MORIS_ASSERT( mdEta2dv2.n_rows() >= tNumEta * tFollowerDofs,
                        "GapData::set_second_order_derivatives - mdEta2dv2 has unexpected number of rows." );
                MORIS_ASSERT( mdEta2dv2.n_cols() >= tFollowerDofs,
                        "GapData::set_second_order_derivatives - mdEta2dv2 has unexpected number of cols." );

                MORIS_ASSERT( mdEta2duv_mat.n_rows() >= tNumEta * tLeaderDofs,
                        "GapData::set_second_order_derivatives - mdEta2duv_mat has unexpected number of rows." );
                MORIS_ASSERT( mdEta2duv_mat.n_cols() >= tFollowerDofs,
                        "GapData::set_second_order_derivatives - mdEta2duv_mat has unexpected number of cols." );

                MORIS_ASSERT( tdEta2du2.n_rows() == tNumEta * tLeaderDofs && tdEta2du2.n_cols() == tLeaderDofs,
                        "GapData::set_second_order_derivatives - tdEta2du2 input has unexpected shape." );
                MORIS_ASSERT( tdEta2duv.n_rows() == tNumEta * tLeaderDofs && tdEta2duv.n_cols() == tFollowerDofs,
                        "GapData::set_second_order_derivatives - tdEta2duv input has unexpected shape." );

                MORIS_ASSERT( mdEta2duv.size() == tNumEta,
                        "GapData::set_second_order_derivatives - mdEta2duv has unexpected size." );
                for ( uint i = 0; i < mdEta2duv.size(); ++i )
                {
                        MORIS_ASSERT( mdEta2duv[ i ].n_rows() == tLeaderDofs && mdEta2duv[ i ].n_cols() == tFollowerDofs,
                                "GapData::set_second_order_derivatives - mdEta2duv[i] has unexpected shape." );
                }

                // Map second-order derivatives using explicit leader and follower DOF/node counts
                const uint tLeaderNumNodes   = tLeaderDofs / aSpaceDim;
                const uint tFollowerNumNodes = ( tFollowerDofs > 0 ) ? ( tFollowerDofs / aSpaceDim ) : 0;

                // 1) Leader-leader mappings (du2, leader-side blocks)
                uint tDstIDofL = 0;
                for ( uint idim = 0; idim < aSpaceDim; ++idim )
                {
                        for ( uint in = 0; in < tLeaderNumNodes; ++in, ++tDstIDofL )
                        {
                        const uint tSrcIDofL = in * aSpaceDim + idim;

                        uint tDstJDofL = 0;
                        for ( uint jdim = 0; jdim < aSpaceDim; ++jdim )
                        {
                                for ( uint jn = 0; jn < tLeaderNumNodes; ++jn, ++tDstJDofL )
                                {
                                const uint tSrcJDofL = jn * aSpaceDim + jdim;

                                mdGap2du2( tDstIDofL, tDstJDofL ) = tdGap2du2( tSrcIDofL, tSrcJDofL );

                                for ( uint iEta = 0; iEta < tNumEta; ++iEta )
                                {
                                        // leader-side eta du2 mapping
                                        if ( ( iEta * tLeaderDofs + tDstIDofL ) < mdEta2du2.n_rows() && tDstJDofL < mdEta2du2.n_cols() )
                                        {
                                        mdEta2du2( iEta * tLeaderDofs + tDstIDofL, tDstJDofL ) = tdEta2du2( iEta * tLeaderDofs + tSrcIDofL, tSrcJDofL );
                                        }
                                }

                                // leader-side flattened/vec mappings
                                mLeaderdNormal2du2( { 0, aSpaceDim - 1 }, { tDstIDofL * tLeaderDofs + tDstJDofL, tDstIDofL * tLeaderDofs + tDstJDofL } ) =
                                        tLeaderdNormal2dU2( { 0, aSpaceDim - 1 }, { tSrcIDofL * tLeaderDofs + tSrcJDofL, tSrcIDofL * tLeaderDofs + tSrcJDofL } );
                                mdGapvec2du2( { 0, aSpaceDim - 1 }, { tDstIDofL * tLeaderDofs + tDstJDofL, tDstIDofL * tLeaderDofs + tDstJDofL } ) =
                                        tdGapvec2du2( { 0, aSpaceDim - 1 }, { tSrcIDofL * tLeaderDofs + tSrcJDofL, tSrcIDofL * tLeaderDofs + tSrcJDofL } );
                                }
                        }
                        }
                }

                // 2) Follower-follower mappings (dv2, follower-side blocks)
                if ( tFollowerNumNodes > 0 )
                {
                        uint tDstIDofF = 0;
                        for ( uint idim = 0; idim < aSpaceDim; ++idim )
                        {
                        for ( uint in = 0; in < tFollowerNumNodes; ++in, ++tDstIDofF )
                        {
                                const uint tSrcIDofF = in * aSpaceDim + idim;

                                uint tDstJDofF = 0;
                                for ( uint jdim = 0; jdim < aSpaceDim; ++jdim )
                                {
                                for ( uint jn = 0; jn < tFollowerNumNodes; ++jn, ++tDstJDofF )
                                {
                                        const uint tSrcJDofF = jn * aSpaceDim + jdim;

                                        mdGap2dv2( tDstIDofF, tDstJDofF ) = tdGap2dv2( tSrcIDofF, tSrcJDofF );

                                        for ( uint iEta = 0; iEta < tNumEta; ++iEta )
                                        {
                                        if ( ( iEta * tFollowerDofs + tDstIDofF ) < mdEta2dv2.n_rows() && tDstJDofF < mdEta2dv2.n_cols() )
                                        {
                                                mdEta2dv2( iEta * tFollowerDofs + tDstIDofF, tDstJDofF ) = tdEta2dv2( iEta * tFollowerDofs + tSrcIDofF, tSrcJDofF );
                                        }
                                        }

                                        mdGapvec2dv2( { 0, aSpaceDim - 1 }, { tDstIDofF * tFollowerDofs + tDstJDofF, tDstIDofF * tFollowerDofs + tDstJDofF } ) =
                                                tdGapvec2dv2( { 0, aSpaceDim - 1 }, { tSrcIDofF * tFollowerDofs + tSrcJDofF, tSrcIDofF * tFollowerDofs + tSrcJDofF } );
                                }
                                }
                        }
                        }
                }

                // 3) Mixed leader-follower mapping for mdEta2duv (rows = leader DOFs, cols = follower DOFs)
                for ( uint idim = 0, tDstID = 0; idim < aSpaceDim; ++idim )
                {
                        for ( uint in = 0; in < tLeaderNumNodes; ++in, ++tDstID )
                        {
                        const uint tSrcID = in * aSpaceDim + idim;

                        // iterate over follower columns
                        for ( uint jdim = 0, tDstJ = 0; jdim < aSpaceDim; ++jdim )
                        {
                                for ( uint jn = 0; jn < tFollowerNumNodes; ++jn, ++tDstJ )
                                {
                                const uint tSrcJ = jn * aSpaceDim + jdim;

                                mdGap2duv( tDstID, tDstJ ) = tdGap2duv( tSrcID, tSrcJ );

                                for ( uint iEta = 0; iEta < tNumEta; ++iEta )
                                {
                                        uint tRowIdx = iEta * tLeaderDofs + tDstID;
                                        uint tColIdx = tDstJ;
                                        if ( tRowIdx < mdEta2duv_mat.n_rows() && tColIdx < mdEta2duv_mat.n_cols() )
                                        {
                                        mdEta2duv_mat( tRowIdx, tColIdx ) = tdEta2duv( iEta * tLeaderDofs + tSrcID, tSrcJ );
                                        }

                                        if ( mdEta2duv.size() > iEta )
                                        {
                                        if ( tDstID < mdEta2duv[ iEta ].n_rows() && tDstJ < mdEta2duv[ iEta ].n_cols() )
                                        {
                                                mdEta2duv[ iEta ]( tDstID, tDstJ ) = tdEta2duv( iEta * tLeaderDofs + tSrcID, tSrcJ );
                                        }
                                        }

                                        const uint tGapVecCol = tDstID * tFollowerDofs + tDstJ;
                                        if ( tGapVecCol < mdGapvec2duv.n_cols() )
                                        {
                                        mdGapvec2duv( { 0, aSpaceDim - 1 }, { tGapVecCol, tGapVecCol } ) =
                                                tdGapvec2duv( { 0, aSpaceDim - 1 }, { tSrcID * tFollowerDofs + tSrcJ, tSrcID * tFollowerDofs + tSrcJ } );
                                        }
                                }
                                }
                        }
                        }
                }
        }

        //----------------------------------------------------------------------------

        Matrix< DDRMat > multiply_leader_dnormal2du2( const Matrix< DDRMat >& tVector )
        {
            uint tSpaceDim = mLeaderdNormal2du2.n_rows();
            uint tNumDofs  = mdGapdu.numel();

            MORIS_ASSERT( tVector.n_rows() == tSpaceDim,
                    "GapData::multiply_leader_dnormal2du2 - spatial dimensions do not match." );

            Matrix< DDRMat > tResult( tNumDofs, tNumDofs );

            for ( uint i = 0; i < tNumDofs; i++ )
            {
                tResult.get_row( i ) = trans( tVector ) * mLeaderdNormal2du2( { 0, tSpaceDim - 1 }, { i * tNumDofs, ( i + 1 ) * tNumDofs - 1 } );
            }

            return tResult;
        }

        //----------------------------------------------------------------------------

        // static void compute_outward_normal_at_gp_for_linear_deformed_geometry(
        //         Field_Interpolator*    aLeaderFieldInterpolatorManager,
        //         Geometry_Interpolator* aLeaderIGGI,
        //         Matrix< DDRMat >&      aLeaderNormal,
        //         Matrix< DDRMat >&      aLeaderRefNormal,
        //         Matrix< DDRMat >&      aLeaderdNormaldU,
        //         Matrix< DDRMat >&      aLeaderNormal2dU2,
        //         const bool             aEvaluateLinearization = true );

        //----------------------------------------------------------------------------

        static void compute_outward_normal_at_gp_for_consistent_deformed_geometry(
                Field_Interpolator*    aLeaderFieldInterpolatorManager,
                Geometry_Interpolator* aLeaderIGGI,
                Matrix< DDRMat >&      aLeaderNormal,
                Matrix< DDRMat >&      aLeaderRefNormal,
                Matrix< DDRMat >&      aLeaderdNormaldU,
                Matrix< DDRMat >&      aLeaderNormal2dU2,
                const bool             aEvaluateLinearization = true );

        //----------------------------------------------------------------------------

        static Matrix< DDRMat > compute_tangential_plane_projector( const Matrix< DDRMat >& aNormal );
    };

    class Set;
    class Cluster;
    class Field_Interpolator_Manager;
    class FEM_Model;

    //------------------------------------------------------------------------------
    /**
     * Integrand of Weak Form of Governing Equations
     */
    class IWG
    {
      protected:
        // FEM set pointer
        fem::Set* mSet = nullptr;

        // cluster pointer
        fem::Cluster* mCluster = nullptr;

        // nodal weak BCs
        Matrix< DDRMat > mNodalWeakBCs;

        // list of parameters
        Vector< Matrix< DDRMat > > mParameters;

        // normal
        Matrix< DDRMat > mNormal;

        // residual dof type
        Vector< Vector< MSI::Dof_Type > > mResidualDofType;

        // leader and follower dof type lists
        Vector< Vector< MSI::Dof_Type > > mLeaderDofTypes;
        Vector< Vector< MSI::Dof_Type > > mFollowerDofTypes;

        // bool for building global dof type list and map
        bool mGlobalDofBuild   = true;
        bool mGlobalDvBuild    = true;
        bool mGlobalFieldBuild = true;

        // leader and follower global dof type lists
        Vector< Vector< MSI::Dof_Type > > mLeaderGlobalDofTypes;
        Vector< Vector< MSI::Dof_Type > > mFollowerGlobalDofTypes;

        // leader and follower requested global dof type lists
        Vector< Vector< MSI::Dof_Type > > mRequestedLeaderGlobalDofTypes;
        Vector< Vector< MSI::Dof_Type > > mRequestedFollowerGlobalDofTypes;

        // leader and follower field interpolator managers
        Field_Interpolator_Manager* mLeaderFIManager           = nullptr;
        Field_Interpolator_Manager* mFollowerFIManager         = nullptr;
        Field_Interpolator_Manager* mLeaderPreviousFIManager   = nullptr;
        Field_Interpolator_Manager* mFollowerPreviousFIManager = nullptr;

        // leader and follower dv type lists
        Vector< Vector< gen::PDV_Type > > mLeaderDvTypes;
        Vector< Vector< gen::PDV_Type > > mFollowerDvTypes;

        // leader and follower global dv type list
        Vector< Vector< gen::PDV_Type > > mLeaderGlobalDvTypes;
        Vector< Vector< gen::PDV_Type > > mFollowerGlobalDvTypes;

        // leader and follower field type lists
        Vector< Vector< mtk::Field_Type > > mLeaderFieldTypes;
        Vector< Vector< mtk::Field_Type > > mFollowerFieldTypes;

        // leader and follower global dv type list
        Vector< Vector< mtk::Field_Type > > mLeaderGlobalFieldTypes;
        Vector< Vector< mtk::Field_Type > > mFollowerGlobalFieldTypes;

        // leader and follower properties
        Vector< std::shared_ptr< Property > > mLeaderProp;
        Vector< std::shared_ptr< Property > > mFollowerProp;

        // local string to int map for properties
        std::map< std::string, uint > mPropertyMap;

        // leader and follower material models
        Vector< std::shared_ptr< fem::Material_Model > > mLeaderMM;
        Vector< std::shared_ptr< fem::Material_Model > > mFollowerMM;

        // Local string to int map for material models
        std::map< std::string, uint > mMaterialMap;

        // leader and follower constitutive models
        Vector< std::shared_ptr< fem::Constitutive_Model > > mLeaderCM;
        Vector< std::shared_ptr< fem::Constitutive_Model > > mFollowerCM;

        // Local string to int map for constitutive models
        std::map< std::string, uint > mConstitutiveMap;

        // stabilization parameters
        Vector< std::shared_ptr< fem::Stabilization_Parameter > > mStabilizationParam;

        // local string to int map for stabilizations
        std::map< std::string, uint > mStabilizationMap;

        // active cluster measure on IWG flag
        bool mActiveCMEAFlag = false;

        // interpolation order for IWG
        uint mOrder = MORIS_UINT_MAX;

        // tolerance for FD perturbation
        const real mToleranceFD = 1e-12;

        // bulk type
        fem::Element_Type mBulkType = fem::Element_Type::BULK;

        // strings for leader and follower phase name
        std::string mLeaderPhaseName;
        std::string mFollowerPhaseName;

        // bool for time continuity
        bool mTimeContinuity = false;

        // real for time final
        real mTimeFinal = -1.0;

        // bool for time boundary
        bool mTimeBoundary = false;

        // bool for ghost
        bool mIsGhost = false;

        // compute the jacobian using finite differencing (independent on the setting for the whole FEM set)
        bool mIsFDJacobian = false;

      protected:
        // string for IWG name
        std::string mName;

        //! string for IWG name
        enum moris::fem::IWG_Type mIWGType = moris::fem::IWG_Type::END_IWG_TYPE;

        // function pointers
        void ( IWG::*m_compute_jacobian_FD )(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType,
                bool               aUseAbsolutePerturbations ) = nullptr;
        void ( IWG::*m_compute_dRdp_FD_material )(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType ) = nullptr;
        void ( IWG::*m_compute_dRdp_FD_geometry )(
                real                          aWStar,
                real                          aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices ) = nullptr;

        // function pointer for building the perturbation size for FD
        real ( IWG::*m_build_perturbation_size )(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance ) = nullptr;

        // for non-conformal IWGs - gap data
        bool mUseDeformedGeometryForGap = false;
        // bool mUseConsistentDeformedGeometryForGap = false;

        std::unique_ptr< GapData > mGapData = nullptr;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * constructor
         */
        IWG() {};

        //------------------------------------------------------------------------------
        /**
         * virtual destructor
         */
        virtual ~IWG() {};

        //------------------------------------------------------------------------------
        /**
         * set name
         * param[ in ] aName a string for CM name
         */
        void
        set_name( std::string aName )
        {
            mName = std::move( aName );
        }

        //------------------------------------------------------------------------------
        /**
         * get name
         * param[ out ] mName a string for CM name
         */
        std::string
        get_name()
        {
            return mName;
        }

        //------------------------------------------------------------------------------
        /**
         * set time continuity flag
         * param[ in ] aTimeContinuity bool true if IWG for time continuity
         */
        void
        set_time_continuity( bool aTimeContinuity )
        {
            mTimeContinuity = aTimeContinuity;
        }

        //------------------------------------------------------------------------------
        /**
         * get time continuity flag
         * param[ out ] mTimeContinuity ool true if IWG for time continuity
         */
        bool
        get_time_continuity()
        {
            return mTimeContinuity;
        }

        //------------------------------------------------------------------------------
        /**
         * set time final value
         * param[ in ] aTimeFinal real for integration at final time
         */
        void
        set_time_final( real aTimeFinal )
        {
            mTimeFinal = aTimeFinal;
        }

        //------------------------------------------------------------------------------
        /**
         * get time final value
         * param[ in ] mTimeFinal real for integration at final time
         */
        real
        get_time_final()
        {
            return mTimeFinal;
        }

        //------------------------------------------------------------------------------
        /**
         * set time boundary flag
         * param[ in ] aTimeBoundary bool true if IWG for time boundary
         */
        void
        set_time_boundary( bool aTimeBoundary )
        {
            mTimeBoundary = aTimeBoundary;
        }

        //------------------------------------------------------------------------------
        /**
         * get time boundary flag
         * param[ out ] mTimeBoundary bool true if IWG for time boundary
         */
        bool
        get_time_boundary()
        {
            return mTimeBoundary;
        }

        //------------------------------------------------------------------------------
        /**
         * set ghost flag
         * param[ in ] aIsGhost bool true if IWG for ghost
         */
        void
        set_ghost_flag( bool aIsGhost )
        {
            mIsGhost = aIsGhost;
        }

        //------------------------------------------------------------------------------
        /**
         * get ghost flag
         * param[ out ] mIsGhost bool true if IWG for ghost
         */
        bool
        get_ghost_flag()
        {
            return mIsGhost;
        }

        //------------------------------------------------------------------------------
        /**
         * get IWG type
         * param[ out ] mIWGType an enum of the IWG type. type only implemented for TIME_CONTINUITY_DOF. All others return UNDEFINED.
         *              If needed you can implement the type for the others. Just folow the TIME_CONTINUITY_DOF
         */
        enum moris::fem::IWG_Type
        get_IWG_type()
        {
            return mIWGType;
        };

        //------------------------------------------------------------------------------
        /**
         * set phase name
         * param[ in ] aPhaseName a string for phase name
         * param[ in ] aIsLeader  an enum for leader or follower
         */
        void set_phase_name(
                const std::string&   aPhaseName,
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get phase name
         * param[ in ]  aIsLeader an enum for leader or follower
         * param[ out ] mName     a string for phase name
         */
        std::string get_phase_name( mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * print name
         */
        void print_names();

        //------------------------------------------------------------------------------

        bool is_fd_jacobian() const
        {
            return mIsFDJacobian;
        };

        //------------------------------------------------------------------------------

        void set_is_fd_jacobian( bool aIsFDJacobian )
        {
            mIsFDJacobian = aIsFDJacobian;
        };

        //------------------------------------------------------------------------------
        /*
         * set member set pointer
         * @param[ in ] aSetPointer a FEM set pointer
         */
        void
        set_set_pointer( Set* aSetPointer )
        {
            mSet = aSetPointer;

            // set function pointer
            this->set_function_pointers();
        }

        //------------------------------------------------------------------------------
        /*
         * set fem cluster pointer
         * @param[ in ] aClusterPointer a FEM cluster pointer
         */
        void
        set_cluster_pointer( fem::Cluster* aClusterPointer )
        {
            mCluster = aClusterPointer;
        }

        //------------------------------------------------------------------------------
        /*
         * set function pointers
         */
        void set_function_pointers();

        //------------------------------------------------------------------------------
        /*
         * get member set pointer
         * @param[ out ] aSetPointer a FEM set pointer
         */
        Set*
        get_set_pointer()
        {
            return mSet;
        }

        //------------------------------------------------------------------------------
        /*
         * set field interpolator manager
         * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
         * @param[ in ] aIsLeader                 an enum for leader or follower
         */
        void set_field_interpolator_manager(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * set field interpolator manager for previous time step
         * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
         * @param[ in ] aIsLeader                 an enum for leader or follower
         */
        void set_field_interpolator_manager_previous_time(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * get field interpolator manager
         * @param[ out ] aFieldInterpolatorManager a field interpolator manager pointer
         * @param[ in ]  aIsLeader                 an enum for leader or follower
         */
        Field_Interpolator_Manager* get_field_interpolator_manager(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * free memory
         */
        void
        free_memory()
        {
        }

        //------------------------------------------------------------------------------
        /**
         * set parameters
         * @param[ in ] aParameters a list of parameters
         */
        virtual void
        set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
        {
            // set a cluster
            mParameters = aParameters;
        }

        //------------------------------------------------------------------------------
        /**
         * set nodal weak BCs
         * @param[ in ] aNodalWeakBCs matrix with nodal values
         */
        void
        set_nodal_weak_bcs( const Matrix< DDRMat >& aNodalWeakBCs )
        {
            mNodalWeakBCs = aNodalWeakBCs;
        }

        //------------------------------------------------------------------------------
        /**
         * set normal
         * @param[ in ] aNormal normal vector
         */
        void set_normal( const Matrix< DDRMat >& aNormal );

        //------------------------------------------------------------------------------
        /**
         * set residual dof type
         * @param[ in ] aResidualdofType a cell of residual dof types
         */
        void
        set_residual_dof_type( const Vector< Vector< MSI::Dof_Type > >& aResidualDofType )
        {
            mResidualDofType = aResidualDofType;
        }

        //------------------------------------------------------------------------------
        /**
         * return a dof type for the residual
         * @param[ out ] aResidualdofType a cell of residual dof types
         */
        const Vector< Vector< MSI::Dof_Type > >&
        get_residual_dof_type() const
        {
            return mResidualDofType;
        }

        //------------------------------------------------------------------------------
        /**
         * set interpolation order for the residual dof type
         */
        void set_interpolation_order();

        void
        set_interpolation_order( uint aOrder )
        {
            // set order
            mOrder = aOrder;
        }

        //------------------------------------------------------------------------------
        /**
         * set bulk type
         * @param[ in ] aBulkType element type for the IWG
         */
        void
        set_bulk_type( fem::Element_Type aBulkType )
        {
            mBulkType = aBulkType;
        }

        //------------------------------------------------------------------------------
        /**
         * get bulk type
         * @param[ out ] mBulkType element type for the IWG
         */
        fem::Element_Type
        get_bulk_type()
        {
            return mBulkType;
        }

        //------------------------------------------------------------------------------
        /**
         * set IWG active dof types
         * @param[ in ] aDofTypes a list of group of dof types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                mtk::Leader_Follower                     aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * return a cell of dof types active for the IWG
         * @param[ in ] aIsLeader enum leader or follower
         * @param[ out ] aDofTypes a list of group of dof types
         */
        const Vector< Vector< MSI::Dof_Type > >& get_dof_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /**
         * set IWG active dv types
         * @param[ in ] aDvTypes a list of group of dv types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_dv_type_list(
                const Vector< Vector< gen::PDV_Type > >& aDvTypes,
                mtk::Leader_Follower                     aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * return a cell of dv types active for the IWG
         * @param[ in ]  aIsLeader enum leader or follower
         * @param[ out ] aDvTypes a list of group of dv types
         */
        const Vector< Vector< gen::PDV_Type > >& get_dv_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /**
         * return a cell of field types active for the IWG
         * @param[ in ]  aIsLeader enum leader or follower
         * @param[ out ] aFieldTypes a list of group of field types
         */
        const Vector< Vector< mtk::Field_Type > >& get_field_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /**
         * set IWG active field types
         * @param[ in ] aFieldTypes a list of group of field types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_field_type_list(
                const Vector< Vector< mtk::Field_Type > >& aDvTypes,
                mtk::Leader_Follower                       aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * check that field interpolators were assigned
         * @param[ in ]  aIsLeader enum leader or follower
         */
        void check_field_interpolators(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set property
         * @param[ in ] aProperty       a property pointer
         * @param[ in ] aPropertyString a string describing the property
         * @param[ in ] aIsLeader       enum leader or follower
         */
        void set_property(
                std::shared_ptr< Property > aProperty,
                const std::string&          aPropertyString,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get properties
         * @param[ in ]  aIsLeader   enum leader or follower
         * @param[ out ] aProperties cell of property pointers
         */
        Vector< std::shared_ptr< Property > >& get_properties(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set material model
         * @param[ in ] aMaterialModel       a material model pointer
         * @param[ in ] aMaterialModelString a string defining the material model
         * @param[ in ] aIsLeader            an enum for leader or follower
         */
        void set_material_model(
                std::shared_ptr< Material_Model > aMaterialModel,
                const std::string&                aMaterialModelString,
                mtk::Leader_Follower              aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get material models
         * @param[ in ]  aIsLeader           enum leader or follower
         * @param[ out ] aMaterialModels     cell of material model pointers
         */
        Vector< std::shared_ptr< Material_Model > >& get_material_models(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set constitutive model
         * @param[ in ] aConstitutiveModel  a constitutive model pointer
         * @param[ in ] aConstitutiveString a string defining the constitutive model
         * @param[ in ] aIsLeader           an enum for leader or follower
         */
        void set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                const std::string&                    aConstitutiveString,
                mtk::Leader_Follower                  aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get constitutive models
         * @param[ in ]  aIsLeader           enum leader or follower
         * @param[ out ] aConstitutiveModels cell of constitutive model pointers
         */
        Vector< std::shared_ptr< Constitutive_Model > >& get_constitutive_models(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set stabilization parameter
         * @param[ in ] aStabilizationParameter a stabilization parameter pointer
         * @param[ in ] aStabilizationString    a string defining the stabilization parameter
         */
        void set_stabilization_parameter(
                const std::shared_ptr< Stabilization_Parameter >& aStabilizationParameter,
                const std::string&                                aStabilizationString );

        //------------------------------------------------------------------------------
        /**
         * get stabilization parameters
         * @param[ out ] mStabilizationParam cell of stabilization parameter pointers
         */
        Vector< std::shared_ptr< Stabilization_Parameter > >&
        get_stabilization_parameters()
        {
            // return penalty parameter pointers
            return mStabilizationParam;
        }

        //------------------------------------------------------------------------------
        /**
         * create a global dof type list including
         * IWG, property, constitutive and stabilization dependencies
         */
        void build_global_dof_dv_and_field_type_list();

        //------------------------------------------------------------------------------
        /**
         * get a non unique list of dof type including
         * IWG, property, constitutive and stabilization dependencies
         * for both leader and follower
         */
        void get_non_unique_dof_dv_and_field_types(
                Vector< Vector< MSI::Dof_Type > >&   aDofTypes,
                Vector< Vector< gen::PDV_Type > >&   aDvTypes,
                Vector< Vector< mtk::Field_Type > >& aFieldTypes );

        //------------------------------------------------------------------------------
        /**
         * get global dof type list
         * @param[ in ]  aIsLeader       enum leader or follower
         * @param[ out ] mGlobalDofTypes global list of group of dof types
         */
        const Vector< Vector< MSI::Dof_Type > >& get_global_dof_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get global dv type list
         * @param[ in ]  aIsLeader       enum leader or follower
         * @param[ out ] mGlobalDvTypes global list of group of dv types
         */
        const Vector< Vector< gen::PDV_Type > >& get_global_dv_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get global field type list. TODO: Field types are only used by the IWG.
         * If a user wants to use them in a property or CM this cuntion has to be modified in the same way than get_global_dof_type_list()
         * @param[ in ]  aIsLeader    enum leader or follower
         * @param[ out ] mFieldTypes global list of group of dv types
         */
        const Vector< Vector< mtk::Field_Type > >& get_global_field_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * evaluate the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        virtual void compute_residual( real aWStar ) = 0;

        //------------------------------------------------------------------------------
        /**
         * evaluate the Jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        virtual void compute_jacobian( real aWStar ) = 0;

        virtual void
        compute_jacobian_previous( real aWStar )
        {
            MORIS_ERROR( false, "compute_jacobian_previous() not implemented" );
        };

        //------------------------------------------------------------------------------
        /**
         * evaluate the Jacobian by finite difference
         * @param[ in ] aPerturbation real to perturb for FD
         * @param[ in ] aWStar        weight associated to evaluation point
         * @param[ in ] aFDSchemeType enum for FD scheme
         */
        void
        compute_jacobian_FD(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType             = fem::FDScheme_Type::POINT_5,
                bool               aUseAbsolutePerturbations = false )
        {
            // compute jacobian by FD
            ( this->*m_compute_jacobian_FD )( aWStar, aPerturbation, aFDSchemeType, aUseAbsolutePerturbations );
        }

        void select_jacobian_FD(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType,
                bool               aUseAbsolutePerturbations );

        void select_jacobian_FD_double(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType,
                bool               aUseAbsolutePerturbations );

        //------------------------------------------------------------------------------
        /**
         * evaluate the residual and the Jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        virtual void compute_jacobian_and_residual( real aWStar ) = 0;

        //------------------------------------------------------------------------------
        /**
         * check the Jacobian with FD
         * @param[ in ] aPerturbation real to perturb for FD
         * @param[ in ] aEpsilon      real for check
         * @param[ in ] aWStar        real weight associated to evaluation point
         * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
         * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
         * @param[ in ] aErrorPrint   bool set to true to print non matching values in jacobian
         */
        bool check_jacobian(
                real              aPerturbation,
                real              aEpsilon,
                real              aWStar,
                Matrix< DDRMat >& aJacobians,
                Matrix< DDRMat >& aJacobiansFD,
                bool              aErrorPrint               = false,
                bool              aUseAbsolutePerturbations = false );

        //------------------------------------------------------------------------------
        /**
         * check Jacobian that uses multiple dof types with FD
         * @param[ in ] aPerturbation real to perturb for FD
         * @param[ in ] aEpsilon      real for check
         * @param[ in ] aWStar        real weight associated to evaluation point
         * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
         * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
         * @param[ in ] aErrorPrint   bool set to true to print non matching values in jacobian
         */
        bool check_jacobian_multi_residual(
                real              aPerturbation,
                real              aEpsilon,
                real              aWStar,
                Matrix< DDRMat >& aJacobians,
                Matrix< DDRMat >& aJacobiansFD,
                bool              aErrorPrint               = false,
                bool              aMaxErrorPrint            = false,
                moris::real       aFDtolerance              = -1.0,
                bool              aUseAbsolutePerturbations = false );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the residual wrt the design variables
         * @param[ in ] aWStar weight associated to evaluation point
         */
        virtual void compute_dRdp( real aWStar ) = 0;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the residual
         * wrt the material design variables by finite difference
         * @param[ in ] aWStar        weight associated to evaluation point
         * @param[ in ] aPerturbation real for dv perturbation
         * @param[ in ] aFDSchemeType enum for FD scheme
         */
        void
        compute_dRdp_FD_material(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL )
        {
            // compute jacobian by FD
            ( this->*m_compute_dRdp_FD_material )( aWStar, aPerturbation, aFDSchemeType );
        }

        void select_dRdp_FD_material(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType );

        void select_dRdp_FD_material_double(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the residual
         * wrt the geometry design variables by finite difference
         * @param[ in ] aWStar            weight associated to evaluation point
         * @param[ in ] aPerturbation     real for relative dv perturbation
         * @param[ in ] aGeoLocalAssembly matrix filled with pdv local assembly indices
         * @param[ in ] aFDSchemeType     enum for FD scheme
         */
        void
        compute_dRdp_FD_geometry(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices )
        {
            // compute jacobian by FD
            ( this->*m_compute_dRdp_FD_geometry )(
                    aWStar,
                    aPerturbation,
                    aFDSchemeType,
                    aGeoLocalAssembly,
                    aVertexIndices );
        }

        void select_dRdp_FD_geometry_bulk(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices );

        void select_dRdp_FD_geometry_sideset(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices );

        void select_dRdp_FD_geometry_time_sideset(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices );

        void select_dRdp_FD_geometry_double(
                moris::real                   aWStar,
                moris::real                   aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices );

        //------------------------------------------------------------------------------
        /**
         * \brief Applies the given perturbation amount on the coefficients of either follower or leader side and updates the geometry interpolators.
         * \param aPerturbationAmount
         * \param aCoefficients
         * \param aParametricCoefficients
         * \param aEvaluationPoint
         * \param aNodeIndex
         * \param aSpatialDirIndex
         * \param aLeaderFollowerType
         */
        void perturb_and_update_geometry_interpolators(
                real const                 aPerturbationAmount,
                Matrix< DDRMat > const &   aCoefficients,
                Matrix< DDRMat > const &   aParametricCoefficients,
                uint const                 aNodeIndex,
                uint const                 aSpatialDirIndex,
                mtk::Leader_Follower const aLeaderFollowerType ) const;

        /**
         * add the contribution of the cluster measure derivatives to the derivative of
         * the quantity of interest wrt to geometry dv by finite difference
         * @param[ in ] aWStar        weight associated to evaluation point
         * @param[ in ] aPerturbation pdv relative perturbation size
         * @param[ in ] aFDSchemeType enum for FD scheme
         */
        void add_cluster_measure_dRdp_FD_geometry(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType );

        /**
         * add the contribution of the cluster measure derivatives to the derivative of
         * the quantity of interest wrt to geometry dv by finite difference
         * for double sideset
         * @param[ in ] aWStar        weight associated to evaluation point
         * @param[ in ] aPerturbation pdv relative perturbation size
         * @param[ in ] aFDSchemeType enum for FD scheme
         */
        void add_cluster_measure_dRdp_FD_geometry_double(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType );

        //------------------------------------------------------------------------------
        /**
         * build perturbation size for finite difference
         * @param[ in ] aPerturbation         provided perturbation size from input
         * @param[ in ] aCoefficientToPerturb coefficient to perturb
         * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
         */
        real build_perturbation_size(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance );

        /**
         * build relative perturbation size for finite difference
         * @param[ in ] aPerturbation         provided perturbation size from input
         * @param[ in ] aCoefficientToPerturb coefficient to perturb
         * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
         */
        real build_perturbation_size_relative(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance );

        /**
         * build absolute perturbation size for finite difference
         * @param[ in ] aPerturbation         provided perturbation size from input
         * @param[ in ] aCoefficientToPerturb coefficient to perturb
         * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
         */
        real build_perturbation_size_absolute(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance );

        //------------------------------------------------------------------------------
        /**
         * check if ig node still inside ip element after perturbation in a specific
         * space direction, if not adapt the finite difference scheme used
         * @param[ in ] aPerturbation         provided perturbation size from input
         * @param[ in ] aCoefficientToPerturb coefficient to perturb
         * @param[ in ] aSpatialDirection     spatial direction in which we perturb
         * @param[ in ] aUsedFDScheme         FD scheme to be used, updated
         * @param[ in ] aIsLeader             if true leader IP element otherwise follower element is used (default is leader)
         * @param[ out ] aDeltaH              perturbation size built for finite difference
         */
        real check_ig_coordinates_inside_ip_element(
                const real&          aPerturbation,
                const real&          aCoefficientToPerturb,
                const uint&          aSpatialDirection,
                fem::FDScheme_Type&  aUsedFDScheme,
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * reset evaluation flags
         */
        void reset_eval_flags();

        //------------------------------------------------------------------------------
        /**
         * reset evaluation flags specific to child IWG
         */
        virtual void reset_spec_eval_flags() {};

        //------------------------------------------------------------------------------
        /**
         * build a list of dof types requested by the solver and owned by the IWG
         * @param[ in ] aIsResidual bool true if residual evaluation
         */
        void build_requested_dof_type_list( const bool aIsStaggered );

        //------------------------------------------------------------------------------

        Matrix< DDRMat > remap_nonconformal_rays(
                const bool aUseDeformedGeometryForGap,
                // const bool                     aUseConsistentDeformedGeometryForGap,
                const Vector< MSI::Dof_Type >& aDisplDofTypes,
                Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
                Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
                std::unique_ptr< GapData >&    aGapData ) const;

        static Matrix< DDRMat > compute_eta_pair_hessian_weights(
                const Matrix< DDRMat >& aDqgByEta,
                const Matrix< DDRMat >& aCompactHessian,
                const uint              aEtaIndexI,
                const uint              aEtaIndexJ );

        Matrix< DDRMat > remap_nonconformal_rays_consistent_deformed_geometry(
                const Vector< MSI::Dof_Type >& aDisplDofTypes,
                Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
                Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
                std::unique_ptr< GapData >&    aGapData ) const;

        // Matrix< DDRMat > remap_nonconformal_rays_linear_deformed_geometry(
        //         const Vector< MSI::Dof_Type >& aDisplDofTypes,
        //         Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
        //         Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
        //         std::unique_ptr< GapData >&    aGapData ) const;

        Matrix< DDRMat > remap_nonconformal_rays_undeformed_geometry(
                Field_Interpolator_Manager* aLeaderFieldInterpolatorManager,
                Field_Interpolator_Manager* aFollowerFieldInterpolatorManager ) const;

        const std::unique_ptr< GapData >& get_gap_data()
        {
            return mGapData;
        }
    };
    //------------------------------------------------------------------------------

}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
