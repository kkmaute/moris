/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <cmath>
#include "cl_FEM_IWG_Time_DG.hpp"
#include "cl_FEM_Time_DG_Operators.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "fn_trans.hpp"

namespace moris::fem
{
    IWG_Time_DG_Mass::IWG_Time_DG_Mass()
    {
        mIWGType        = IWG_Type::TIME_DG_MASS;
        mTimeContinuity = true;
        mLeaderProp.resize( COUNT, nullptr );
        mPropertyMap = { { "Density", DENSITY }, { "Thickness", THICKNESS }, { "InitialDisplacement", INITIAL_DISPLACEMENT }, { "InitialVelocity", INITIAL_VELOCITY }, { "InitialTime", INITIAL_TIME } };
    }

    std::string IWG_Time_DG_Mass::history_name() const
    {
        return "TimeDG/" + mName + "/Velocity";
    }

    real IWG_Time_DG_Mass::start_time() const
    {
        return mSet->get_equation_model()->get_time()( 0 );
    }

    real IWG_Time_DG_Mass::duration() const
    {
        const auto& tTime     = mSet->get_equation_model()->get_time();
        const real  tDuration = tTime( 1 ) - tTime( 0 );
        MORIS_ERROR( tDuration > 0.0, "Time DG requires a positive slab duration." );
        return tDuration;
    }

    bool IWG_Time_DG_Mass::initial_slab() const
    {
        const real tInitial = mLeaderProp( INITIAL_TIME ) ? mLeaderProp( INITIAL_TIME )->val()( 0 ) : 0.0;
        return start_time() == tInitial;
    }

    Matrix< DDRMat > IWG_Time_DG_Mass::initial_coefficients( uint aProperty )
    {
        auto*            tFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
        Matrix< DDRMat > tResult( tFI->get_number_of_space_bases(), tFI->get_number_of_fields(), 0.0 );
        if ( !mLeaderProp( aProperty ) )
        {
            return tResult;
        }

        // Prescribed data is sampled at spatial interpolation nodes. Callbacks
        // should use the IP geometry interpolator for spatially varying data.
        const Matrix< DDRMat > tSaved = tFI->get_space_time();
        const Matrix< DDRMat > tNodes = tFI->get_space_interpolation_points();
        for ( uint i = 0; i < tResult.n_rows(); ++i )
        {
            Matrix< DDRMat > tPoint = tSaved;
            for ( uint d = 0; d < tNodes.n_rows(); ++d )
            {
                tPoint( d ) = tNodes( d, i );
            }
            tPoint( tPoint.numel() - 1 ) = -1.0;
            mLeaderFIManager->set_space_time( tPoint );
            mLeaderProp( aProperty )->reset_eval_flags();
            const Matrix< DDRMat > tValue = mLeaderProp( aProperty )->val();
            MORIS_ERROR( tValue.numel() == tResult.n_cols(), "Time DG initial data requires one value per displacement component." );
            for ( uint c = 0; c < tResult.n_cols(); ++c )
            {
                tResult( i, c ) = tValue( c );
            }
        }
        mLeaderFIManager->set_space_time( tSaved );
        this->reset_eval_flags();
        return tResult;
    }

    Matrix< DDRMat > IWG_Time_DG_Mass::incoming_displacement()
    {
        if ( initial_slab() )
        {
            return initial_coefficients( INITIAL_DISPLACEMENT );
        }
        auto*      tPrevious = mLeaderPreviousFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
        const uint tSpace    = tPrevious->get_number_of_space_bases();
        MORIS_ERROR( tPrevious->get_number_of_time_bases() == 2, "Time DG requires two previous temporal bases." );
        return tPrevious->get_coeff()( { tSpace, 2 * tSpace - 1 }, { 0, tPrevious->get_number_of_fields() - 1 } );
    }

    Matrix< DDRMat > IWG_Time_DG_Mass::incoming_velocity()
    {
        if ( initial_slab() )
        {
            return initial_coefficients( INITIAL_VELOCITY );
        }
        return mCluster->get_cluster_ISV_data().mTimeHistory.values( history_name(), start_time() );
    }

    Matrix< DDRMat > IWG_Time_DG_Mass::spatial_shape()
    {
        auto* tFI = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
        MORIS_ERROR( tFI->get_number_of_time_bases() == 2, "Energy-corrected DG(1) requires two temporal displacement bases." );
        const uint tSpace = tFI->get_number_of_space_bases();
        // This IWG is evaluated at the current incoming endpoint, tau=-1.
        return tFI->NBuild()( { 0, 0 }, { 0, tSpace - 1 } );
    }

    void IWG_Time_DG_Mass::commit_time_history()
    {
        auto*      tFI     = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
        const uint tSpace  = tFI->get_number_of_space_bases();
        const uint tFields = tFI->get_number_of_fields();
        MORIS_ERROR( tFI->get_number_of_time_bases() == 2, "Time DG history requires DG(1) coefficients." );
        const auto&      tOp       = Time_DG_Operators::get();
        Matrix< DDRMat > tVelocity = -tOp.b( 1 ) * incoming_displacement();
        for ( uint b = 0; b < 2; ++b )
        {
            tVelocity += tOp.A( 1, b ) * tFI->get_coeff()( { b * tSpace, ( b + 1 ) * tSpace - 1 }, { 0, tFields - 1 } );
        }
        tVelocity = tVelocity / duration();

        auto&      tHistory = mCluster->get_cluster_ISV_data().mTimeHistory;
        const real tEnd     = mSet->get_equation_model()->get_time()( 1 );
        tHistory.stage( history_name(), start_time(), tEnd, tVelocity );
        tHistory.commit( history_name() );
    }

    void IWG_Time_DG_Mass::assemble( real aWeight, bool aResidual, bool aJacobian, bool aPrevious )
    {
        auto*       tFI     = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
        const uint  tSpace  = tFI->get_number_of_space_bases();
        const uint  tFields = tFI->get_number_of_fields();
        const uint  tIndex  = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ) );
        const uint  tStart  = mSet->get_res_dof_assembly_map()( tIndex )( 0 );
        const auto& tOp     = Time_DG_Operators::get();
        const real  tDt     = duration();
        MORIS_ERROR( mLeaderProp( DENSITY ) != nullptr, "Time DG requires Density." );
        // Velocity elimination assumes a common, state-independent mass in
        // both kinematic equations. Reject dependencies on any solution field,
        // including coupled temperature fields, rather than omit their tangent.
        // Spatial/design dependence is allowed; callbacks must be time independent.
        for ( const auto tProperty : { DENSITY, THICKNESS } )
        {
            if ( mLeaderProp( tProperty ) == nullptr )
            {
                continue;
            }
            const char* tName = tProperty == DENSITY ? "Density" : "Thickness";
            MORIS_ERROR( mLeaderProp( tProperty )->get_dof_type_list().size() == 0,
                    "Time DG requires state-independent %s.",
                    tName );
            const auto& tValue = mLeaderProp( tProperty )->val();
            MORIS_ERROR( tValue.numel() == 1 && std::isfinite( tValue( 0 ) ) && tValue( 0 ) > 0.0,
                    "Time DG requires a finite, strictly positive scalar %s.",
                    tName );
        }

        const real tDensity   = mLeaderProp( DENSITY )->val()( 0 );
        const real tThickness = mLeaderProp( THICKNESS ) ? mLeaderProp( THICKNESS )->val()( 0 ) : 1.0;

        const Matrix< DDRMat > tN      = spatial_shape();
        const real             tFactor = aWeight * tDensity * tThickness;

        if ( aResidual )
        {
            // R_a = M[(GA)_ab Q_b/dt - (Gb)_a q_minus/dt - e_a v_minus].
            // Local coefficient ordering is component, time basis, space basis.
            const Matrix< DDRMat > tQMinus = tN * incoming_displacement();
            const Matrix< DDRMat > tVMinus = tN * incoming_velocity();
            for ( uint a = 0; a < 2; ++a )
            {
                Matrix< DDRMat > tValue = -( tOp.Gb( a ) / tDt ) * tQMinus - tOp.e( a ) * tVMinus;
                for ( uint b = 0; b < 2; ++b )
                {
                    tValue += ( tOp.GA( a, b ) / tDt ) * tN * tFI->get_coeff()( { b * tSpace, ( b + 1 ) * tSpace - 1 }, { 0, tFields - 1 } );
                }
                for ( uint c = 0; c < tFields; ++c )
                {
                    for ( uint i = 0; i < tSpace; ++i )
                    {
                        mSet->get_residual()( 0 )( tStart + c * 2 * tSpace + a * tSpace + i ) += tFactor * tN( i ) * tValue( c );
                    }
                }
            }
        }

        if ( aJacobian && ( !aPrevious || !initial_slab() ) )
        {
            bool tRequested = false;
            for ( const auto& tGroup : mRequestedLeaderGlobalDofTypes )
            {
                tRequested = tRequested || tGroup( 0 ) == mResidualDofType( 0 )( 0 );
            }
            if ( !tRequested )
            {
                return;
            }
            const uint tColumnStart = mSet->get_jac_dof_assembly_map()( tIndex )( tIndex, 0 );
            const real tPreviousDt  = aPrevious ? start_time() - mCluster->get_cluster_ISV_data().mTimeHistory.start( history_name(), start_time() ) : tDt;
            for ( uint a = 0; a < 2; ++a )
            {
                for ( uint b = 0; b < 2; ++b )
                {
                    // Previous velocity depends on both previous slab coefficients.
                    // Its incoming-displacement dependence is propagated separately
                    // by add_adjoint_history (the second preceding slab).
                    const real tTemporal = aPrevious ? -( b == 1 ? tOp.Gb( a ) / tDt : 0.0 ) - tOp.e( a ) * tOp.A( 1, b ) / tPreviousDt : tOp.GA( a, b ) / tDt;
                    for ( uint c = 0; c < tFields; ++c )
                    {
                        for ( uint i = 0; i < tSpace; ++i )
                        {
                            for ( uint j = 0; j < tSpace; ++j )
                            {
                                mSet->get_jacobian()( tStart + c * 2 * tSpace + a * tSpace + i,
                                        tColumnStart + c * 2 * tSpace + b * tSpace + j ) += tFactor * tTemporal * tN( i ) * tN( j );
                            }
                        }
                    }
                }
            }
        }
    }

    void IWG_Time_DG_Mass::compute_residual( real aWeight ) { assemble( aWeight, true, false, false ); }
    void IWG_Time_DG_Mass::compute_jacobian( real aWeight ) { assemble( aWeight, false, true, false ); }
    void IWG_Time_DG_Mass::compute_jacobian_and_residual( real aWeight ) { assemble( aWeight, true, true, false ); }
    void IWG_Time_DG_Mass::compute_jacobian_previous( real aWeight ) { assemble( aWeight, false, true, true ); }

    void IWG_Time_DG_Mass::compute_dRdp( real aWeight )
    {
        this->compute_dRdp_FD_material( aWeight, mSet->get_finite_difference_perturbation_size(), mSet->get_finite_difference_scheme_for_sensitivity_analysis() );
    }

    void IWG_Time_DG_Mass::add_adjoint_history( real aWeight, uint aElement, uint aPoint, const Vector< Matrix< DDRMat > >& aAdjoint, Vector< Matrix< DDRMat > >& aRHS )
    {
        auto*                  tFI       = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
        const uint             tSpace    = tFI->get_number_of_space_bases();
        const uint             tFields   = tFI->get_number_of_fields();
        const uint             tIndex    = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ) );
        const uint             tStart    = mSet->get_res_dof_assembly_map()( tIndex )( 0 );
        const Matrix< DDRMat > tN        = spatial_shape();
        auto&                  tHistory  = mCluster->get_cluster_ISV_data().mTimeHistory;
        Matrix< DDRMat >       tIncoming = tHistory.take_adjoint( history_name(), start_time(), { aElement, aPoint } );
        Matrix< DDRMat >       tNext( tFields, aAdjoint.size(), 0.0 );
        const auto&            tOp            = Time_DG_Operators::get();
        const real             tDensity       = mLeaderProp( DENSITY )->val()( 0 );
        const real             tThickness     = mLeaderProp( THICKNESS ) ? mLeaderProp( THICKNESS )->val()( 0 ) : 1.0;
        const real             tPreviousStart = initial_slab() ? start_time() : tHistory.start( history_name(), start_time() );
        const real             tPreviousDt    = start_time() - tPreviousStart;

        for ( uint r = 0; r < aAdjoint.size(); ++r )
        {
            for ( uint c = 0; c < tFields; ++c )
            {
                for ( uint i = 0; i < tSpace; ++i )
                {
                    if ( tIncoming.numel() > 0 )
                    {
                        aRHS( r )( tStart + c * 2 * tSpace + tSpace + i ) += aWeight * tN( i ) * tIncoming( c, r );
                    }
                    if ( !initial_slab() )
                    {
                        for ( uint a = 0; a < 2; ++a )
                        {
                            // C_n = dR_n/dQ_(n-2). A weighted spatial covector
                            // is carried locally; quadrature weight is applied
                            // once, when assembling into the preceding topology.
                            tNext( c, r ) += tDensity * tThickness * tOp.e( a ) * tOp.b( 1 ) / tPreviousDt * tN( i ) * aAdjoint( r )( tStart + c * 2 * tSpace + a * tSpace + i );
                        }
                    }
                }
            }
        }
        if ( !initial_slab() )
        {
            tHistory.set_adjoint( history_name(), tPreviousStart, { aElement, aPoint }, tNext );
        }
    }

    IWG_Time_DG_Spatial::IWG_Time_DG_Spatial( std::shared_ptr< IWG > aSpatial )
            : mSpatial( std::move( aSpatial ) )
    {
        MORIS_ERROR( mSpatial && !mSpatial->has_time_history(), "A time-history IWG cannot be wrapped as a spatial force." );
        mPropertyMap        = mSpatial->mPropertyMap;
        mConstitutiveMap    = mSpatial->mConstitutiveMap;
        mMaterialMap        = mSpatial->mMaterialMap;
        mStabilizationMap   = mSpatial->mStabilizationMap;
        mLeaderProp         = mSpatial->mLeaderProp;
        mFollowerProp       = mSpatial->mFollowerProp;
        mLeaderCM           = mSpatial->mLeaderCM;
        mFollowerCM         = mSpatial->mFollowerCM;
        mLeaderMM           = mSpatial->mLeaderMM;
        mFollowerMM         = mSpatial->mFollowerMM;
        mStabilizationParam = mSpatial->mStabilizationParam;
        mActiveCMEAFlag     = mSpatial->mActiveCMEAFlag;
    }

    void IWG_Time_DG_Spatial::configure_spatial()
    {
        mSpatial->set_set_pointer( mSet );
        mSpatial->set_cluster_pointer( mCluster );
        mSpatial->mResidualDofType                 = mResidualDofType;
        mSpatial->mLeaderDofTypes                  = mLeaderDofTypes;
        mSpatial->mFollowerDofTypes                = mFollowerDofTypes;
        mSpatial->mLeaderGlobalDofTypes            = mLeaderGlobalDofTypes;
        mSpatial->mFollowerGlobalDofTypes          = mFollowerGlobalDofTypes;
        mSpatial->mRequestedLeaderGlobalDofTypes   = mRequestedLeaderGlobalDofTypes;
        mSpatial->mRequestedFollowerGlobalDofTypes = mRequestedFollowerGlobalDofTypes;
        mSpatial->mLeaderFIManager                 = mLeaderFIManager;
        mSpatial->mFollowerFIManager               = mFollowerFIManager;
        mSpatial->mLeaderProp                      = mLeaderProp;
        mSpatial->mFollowerProp                    = mFollowerProp;
        mSpatial->mLeaderCM                        = mLeaderCM;
        mSpatial->mFollowerCM                      = mFollowerCM;
        mSpatial->mLeaderMM                        = mLeaderMM;
        mSpatial->mFollowerMM                      = mFollowerMM;
        mSpatial->mStabilizationParam              = mStabilizationParam;
        mSpatial->mNormal                          = mNormal;
        mSpatial->mOrder                           = mOrder;
        mSpatial->mParameters                      = mParameters;
    }

    void IWG_Time_DG_Spatial::add_lifted_rows(
            const Matrix< DDRMat >& aMatrix,
            uint                    aTemporalRow,
            Matrix< DDRMat >&       aResult )
    {
        // Accumulate into caller-owned storage; each temporal evaluation only
        // contributes its active test rows, including leader/follower blocks.
        // Within each side, coefficient ordering is field/time/spatial basis.
        for ( auto tSide : { mtk::Leader_Follower::LEADER, mtk::Leader_Follower::FOLLOWER } )
        {
            const auto& tTypes = tSide == mtk::Leader_Follower::LEADER ? mLeaderDofTypes : mFollowerDofTypes;
            if ( tTypes.size() == 0 )
            {
                continue;
            }
            auto* tManager = tSide == mtk::Leader_Follower::LEADER ? mLeaderFIManager : mFollowerFIManager;
            for ( const auto& tGroup : mResidualDofType )
            {
                const sint tIndex = mSet->get_dof_index_for_type( tGroup( 0 ), tSide );
                if ( tIndex < 0 )
                {
                    continue;
                }
                auto* tFI = tManager->get_field_interpolators_for_type( tGroup( 0 ) );
                MORIS_ERROR( tFI->get_number_of_time_bases() == 2, "Time DG spatial terms require DG(1) trial functions." );
                const uint tSpace = tFI->get_number_of_space_bases();
                const uint tStart = mSet->get_res_dof_assembly_map()( tIndex )( 0 );
                for ( uint c = 0; c < tFI->get_number_of_fields(); ++c )
                {
                    for ( uint i = 0; i < tSpace; ++i )
                    {
                        const uint tBase = tStart + c * 2 * tSpace + i;
                        for ( uint j = 0; j < aMatrix.n_cols(); ++j )
                        {
                            // Sum over temporal test bases using partition of
                            // unity, extracting the ordinary spatial test row.
                            aResult( tBase + aTemporalRow * tSpace, j ) += aMatrix( tBase, j ) + aMatrix( tBase + tSpace, j );
                        }
                    }
                }
            }
        }
    }

    void IWG_Time_DG_Spatial::assemble( real aWeight, bool aResidual, bool aJacobian )
    {
        configure_spatial();
        auto& tResidual = mSet->get_residual()( 0 );
        auto& tJacobian = mSet->get_jacobian();

        // The wrapped IWG assembles in Set storage. Save and accumulate only
        // requested quantities, leaving the other assembly buffer untouched.
        // In particular, residual-only calls must avoid dense Jacobian work.
        Matrix< DDRMat > tSavedResidual;
        Matrix< DDRMat > tSavedJacobian;
        Matrix< DDRMat > tResultResidual;
        Matrix< DDRMat > tResultJacobian;

        if ( aResidual )
        {
            tSavedResidual = tResidual;
            tResultResidual.set_size( tResidual.n_rows(), tResidual.n_cols(), 0.0 );
        }

        if ( aJacobian )
        {
            tSavedJacobian = tJacobian;
            tResultJacobian.set_size( tJacobian.n_rows(), tJacobian.n_cols(), 0.0 );
        }

        const Matrix< DDRMat > tLeaderPoint   = mLeaderFIManager->get_IP_geometry_interpolator()->get_space_time();
        const Matrix< DDRMat > tLeaderIGPoint = mLeaderFIManager->get_IG_geometry_interpolator()->get_space_time();
        const bool             tFollower      = !mFollowerDofTypes.empty();
        Matrix< DDRMat >       tFollowerPoint;
        Matrix< DDRMat >       tFollowerIGPoint;

        if ( tFollower )
        {
            tFollowerPoint   = mFollowerFIManager->get_IP_geometry_interpolator()->get_space_time();
            tFollowerIGPoint = mFollowerFIManager->get_IG_geometry_interpolator()->get_space_time();
        }

        auto tRestore = [ & ]() {
            mLeaderFIManager->set_space_time( tLeaderPoint );
            mLeaderFIManager->get_IG_geometry_interpolator()->set_space_time( tLeaderIGPoint );
            if ( tFollower )
            {
                mFollowerFIManager->set_space_time( tFollowerPoint );
                mFollowerFIManager->get_IG_geometry_interpolator()->set_space_time( tFollowerIGPoint );
            }

            if ( aResidual )
            {
                tResidual = tSavedResidual;
            }

            if ( aJacobian )
            {
                tJacobian = tSavedJacobian;
            }

            this->reset_eval_flags();
        };

        try
        {
            for ( uint a = 0; a < 2; ++a )
            {
                if ( a == 1 )
                {
                    Matrix< DDRMat > tPoint      = tLeaderPoint;
                    tPoint( tPoint.numel() - 1 ) = -1.0;
                    mLeaderFIManager->set_space_time( tPoint );
                    mLeaderFIManager->get_IG_geometry_interpolator()->set_time( Matrix< DDRMat >{ { -1.0 } } );
                    if ( tFollower )
                    {
                        tPoint                       = tFollowerPoint;
                        tPoint( tPoint.numel() - 1 ) = -1.0;
                        mFollowerFIManager->set_space_time( tPoint );
                        mFollowerFIManager->get_IG_geometry_interpolator()->set_time( Matrix< DDRMat >{ { -1.0 } } );
                    }
                }
                mSpatial->reset_eval_flags();
                if ( aResidual )
                {
                    tResidual.fill( 0.0 );
                    mSpatial->compute_residual( aWeight );
                    add_lifted_rows( tResidual, a, tResultResidual );
                }
                if ( aJacobian )
                {
                    tJacobian.fill( 0.0 );
                    mSpatial->compute_jacobian( aWeight );
                    add_lifted_rows( tJacobian, a, tResultJacobian );
                }
            }
        } catch ( ... )
        {
            tRestore();
            throw;
        }
        tRestore();
        if ( aResidual )
        {
            tResidual += tResultResidual;
        }
        if ( aJacobian )
        {
            tJacobian += tResultJacobian;
        }
    }

    void IWG_Time_DG_Spatial::compute_residual( real aWeight ) { assemble( aWeight, true, false ); }
    void IWG_Time_DG_Spatial::compute_jacobian( real aWeight ) { assemble( aWeight, false, true ); }
    void IWG_Time_DG_Spatial::compute_jacobian_and_residual( real aWeight ) { assemble( aWeight, true, true ); }
    void IWG_Time_DG_Spatial::compute_jacobian_previous( real ) {}
    void IWG_Time_DG_Spatial::compute_dRdp( real aWeight )
    {
        this->compute_dRdp_FD_material( aWeight, mSet->get_finite_difference_perturbation_size(), mSet->get_finite_difference_scheme_for_sensitivity_analysis() );
    }
}    // namespace moris::fem
