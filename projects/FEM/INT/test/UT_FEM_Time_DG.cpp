/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <catch.hpp>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "cl_FEM_Time_DG_Operators.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"

#define protected public
#define private public
#include "cl_FEM_Cluster_Time_History.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Model.hpp"
#undef protected
#undef private
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_FEM_IWG_Time_DG.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"

using namespace moris;
using namespace moris::fem;

namespace
{
    using Pair  = std::array< real, 2 >;
    using Block = std::array< Pair, 2 >;

    Pair solve2( const Block& aJ, const Pair& aR )
    {
        const real tDet = aJ[ 0 ][ 0 ] * aJ[ 1 ][ 1 ] - aJ[ 0 ][ 1 ] * aJ[ 1 ][ 0 ];
        return { ( aJ[ 1 ][ 1 ] * aR[ 0 ] - aJ[ 0 ][ 1 ] * aR[ 1 ] ) / tDet,
            ( -aJ[ 1 ][ 0 ] * aR[ 0 ] + aJ[ 0 ][ 0 ] * aR[ 1 ] ) / tDet };
    }

    struct Slab
    {
        Pair  q;
        Pair  v;
        Pair  residual;
        Block jacobian;
    };

    Slab evaluate( Pair aQ, real aU, real aV, real aDt, real aDensity, real aHeight, bool aNonlinear )
    {
        Slab t;
        t.q               = aQ;
        const auto& tOp   = Time_DG_Operators::get();
        const real  tMass = aDensity * aHeight;
        const real  tK    = aHeight * aHeight * aHeight;
        const real  tBeta = aNonlinear ? 7.0 : 0.0;
        for ( uint a = 0; a < 2; ++a )
        {
            t.v[ a ] = ( tOp.A( a, 0 ) * aQ[ 0 ] + tOp.A( a, 1 ) * aQ[ 1 ] - tOp.b( a ) * aU ) / aDt;
        }
        for ( uint a = 0; a < 2; ++a )
        {
            // Mixed momentum evaluation is independent of condensed GA/Gb.
            t.residual[ a ] = tMass * ( tOp.G( a, 0 ) * t.v[ 0 ] + tOp.G( a, 1 ) * t.v[ 1 ] - tOp.e( a ) * aV );
            for ( uint b = 0; b < 2; ++b )
            {
                t.jacobian[ a ][ b ] = tMass * tOp.GA( a, b ) / aDt;
            }
        }
        for ( real x : { 0.5 - 0.5 / std::sqrt( 3.0 ), 0.5 + 0.5 / std::sqrt( 3.0 ) } )
        {
            const real q = ( 1.0 - x ) * aQ[ 0 ] + x * aQ[ 1 ];
            t.residual[ 0 ] += 0.5 * aDt * tK * ( q + tBeta * q * q * q );
            t.jacobian[ 0 ][ 0 ] += 0.5 * aDt * tK * ( 1.0 + 3.0 * tBeta * q * q ) * ( 1.0 - x );
            t.jacobian[ 0 ][ 1 ] += 0.5 * aDt * tK * ( 1.0 + 3.0 * tBeta * q * q ) * x;
        }
        t.residual[ 1 ] += aDt * tK * ( aQ[ 0 ] + tBeta * std::pow( aQ[ 0 ], 3 ) );
        t.jacobian[ 1 ][ 0 ] += aDt * tK * ( 1.0 + 3.0 * tBeta * aQ[ 0 ] * aQ[ 0 ] );
        return t;
    }

    std::vector< Slab > trajectory( real aDensity, real aHeight, bool aNonlinear )
    {
        std::vector< Slab > tSteps;
        real                u = 0.13;
        real                v = 0.7;
        for ( real dt : { 0.08, 0.11, 0.065, 0.09, 0.12 } )
        {
            Pair q = { u, u + dt * v };
            Slab t;
            for ( uint k = 0; k < 20; ++k )
            {
                t = evaluate( q, u, v, dt, aDensity, aHeight, aNonlinear );
                if ( std::hypot( t.residual[ 0 ], t.residual[ 1 ] ) < 1e-13 )
                {
                    break;
                }
                Pair d = solve2( t.jacobian, t.residual );
                q[ 0 ] -= d[ 0 ];
                q[ 1 ] -= d[ 1 ];
            }
            REQUIRE( std::hypot( t.residual[ 0 ], t.residual[ 1 ] ) < 1e-11 );
            tSteps.push_back( t );
            u = t.q[ 1 ];
            v = t.v[ 1 ];
        }
        return tSteps;
    }

    // All objects own their fixtures explicitly. Set and FI manager pointers
    // below are borrowed and detached before their stack owners are destroyed.
    struct Mass_Fixture
    {
        Set                     set;
        FEM_Model               model;
        Cluster                 cluster;
        mtk::Interpolation_Rule rule;
        Geometry_Interpolator   geometry;
        Geometry_Interpolator   previous_geometry;
        Field_Interpolator      current;
        Field_Interpolator      previous;
        // Managers borrow this list; it must outlive both managers.
        Vector< Vector< MSI::Dof_Type > > types = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
        Field_Interpolator_Manager        manager;
        Field_Interpolator_Manager        previous_manager;
        std::shared_ptr< IWG >            iwg;
        Matrix< DDRMat >                  time = { { 0.0 }, { 0.2 } };

        Mass_Fixture()
                : rule( mtk::Geometry_Type::QUAD, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR )
                , geometry( rule )
                , previous_geometry( rule )
                , current( 2, rule, &geometry, { MSI::Dof_Type::UX, MSI::Dof_Type::UY } )
                , previous( 2, rule, &previous_geometry, { MSI::Dof_Type::UX, MSI::Dof_Type::UY } )
                , manager( types, {}, {}, &set )
                , previous_manager( types, {}, {}, &set )
        {
            set.set_set_type( Element_Type::TIME_SIDESET );
            set.mEquationModel = &model;
            model.set_time( time );
            set.mUniqueDofTypeMap.set_size( uint( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            set.mLeaderDofTypeMap                              = set.mUniqueDofTypeMap;
            set.mUniqueDofTypeMap( uint( MSI::Dof_Type::UX ) ) = 0;
            set.mLeaderDofTypeMap( uint( MSI::Dof_Type::UX ) ) = 0;
            set.mUniqueDofTypeList                             = { MSI::Dof_Type::UX };
            manager.mDofTypeMap                                = set.mLeaderDofTypeMap;
            previous_manager.mDofTypeMap                       = set.mLeaderDofTypeMap;
            manager.mMaxNumDofFI                               = 1;
            previous_manager.mMaxNumDofFI                      = 1;
            set.mResDofAssemblyMap                             = { Matrix< DDSMat >{ { 0, 15 } } };
            set.mJacDofAssemblyMap                             = { Matrix< DDSMat >{ { 0, 15 } } };
            set.mResidual                                      = { Matrix< DDRMat >( 16, 1, 0.0 ) };
            set.mJacobian.set_size( 16, 16, 0.0 );
            geometry.set_space_coeff( Matrix< DDRMat >{ { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } } );
            geometry.set_time_coeff( time );
            geometry.set_param_coeff();
            previous_geometry.set_space_coeff( geometry.get_space_coeff() );
            previous_geometry.set_time_coeff( time );
            previous_geometry.set_param_coeff();
            Matrix< DDRMat > coeff( 8, 2 );
            for ( uint c = 0; c < 2; ++c )
            {
                for ( uint i = 0; i < 8; ++i )
                {
                    coeff( i, c ) = 0.01 * ( 1 + i + 2 * c );
                }
            }
            current.set_coeff( coeff );
            previous.set_coeff( coeff );
            manager.mFI                              = { &current };
            manager.mIPGeometryInterpolator          = &geometry;
            manager.mIGGeometryInterpolator          = &geometry;
            previous_manager.mFI                     = { &previous };
            previous_manager.mIPGeometryInterpolator = &previous_geometry;
            previous_manager.mIGGeometryInterpolator = &previous_geometry;
            iwg                                      = IWG_Factory().create_IWG( IWG_Type::TIME_DG_MASS );
            iwg->set_name( "mass" );
            iwg->set_set_pointer( &set );
            iwg->set_cluster_pointer( &cluster );
            iwg->set_residual_dof_type( types );
            iwg->set_dof_type_list( types );
            iwg->mRequestedLeaderGlobalDofTypes = types;
            auto density                        = std::make_shared< Property >();
            density->set_parameters( { Matrix< DDRMat >{ { 1.7 } } } );
            iwg->set_property( density, "Density" );
            iwg->set_field_interpolator_manager( &manager );
            iwg->set_field_interpolator_manager_previous_time( &previous_manager );
            manager.set_space_time( Matrix< DDRMat >{ { 0.23 }, { -0.31 }, { -1.0 } } );
        }

        ~Mass_Fixture()
        {
            manager.mFI.clear();
            previous_manager.mFI.clear();
            manager.mIPGeometryInterpolator          = nullptr;
            manager.mIGGeometryInterpolator          = nullptr;
            previous_manager.mIPGeometryInterpolator = nullptr;
            previous_manager.mIGGeometryInterpolator = nullptr;
            set.mEquationModel                       = nullptr;
        }
    };
}    // namespace

TEST_CASE( "Time DG temporal elimination and collocation", "[FEM][TimeDG]" )
{
    const auto& t = Time_DG_Operators::get();
    REQUIRE( norm( t.A - Matrix< DDRMat >{ { 3, 1 }, { -3, 1 } } ) < 1e-13 );
    REQUIRE( norm( t.b - Matrix< DDRMat >{ { 4 }, { -2 } } ) < 1e-13 );
    Matrix< DDRMat > q = { { 0.2, -0.3 }, { 0.4, 0.7 } };
    Matrix< DDRMat > u = { { 0.1, 0.8 } };
    Matrix< DDRMat > v = ( t.A * q - t.b * u ) / 0.17;
    REQUIRE( norm( t.G * q - 0.17 * t.H * v - t.e * u ) < 1e-13 );
    REQUIRE( std::abs( v( 1, 0 ) - ( q( 1, 0 ) - q( 0, 0 ) ) / 0.17 ) > 0.1 );
}

TEST_CASE( "Cluster time history commit rollback and adjoint isolation", "[FEM][TimeDG][ISV]" )
{
    Cluster_Time_History h;
    h.stage( "velocity", 0.0, 0.1, Matrix< DDRMat >{ { 2.0 }, { 3.0 } } );
    h.rollback( "velocity" );
    h.stage( "velocity", 0.0, 0.1, Matrix< DDRMat >{ { 4.0 }, { 5.0 } } );
    h.commit( "velocity" );
    h.stage( "velocity", 0.1, 0.3, Matrix< DDRMat >{ { 7.0 }, { 8.0 } } );
    REQUIRE( h.values( "velocity", 0.1 )( 0 ) == 4.0 );
    h.commit( "velocity" );
    REQUIRE( h.start( "velocity", 0.3 ) == 0.1 );
    h.set_adjoint( "velocity", 0.1, { 2, 3 }, Matrix< DDRMat >{ { 9.0, 10.0 } } );
    REQUIRE( h.take_adjoint( "velocity", 0.1, { 2, 4 } ).numel() == 0 );
    REQUIRE( h.take_adjoint( "velocity", 0.1, { 2, 3 } )( 1 ) == 10.0 );
    REQUIRE( h.take_adjoint( "velocity", 0.1, { 2, 3 } ).numel() == 0 );
    h.clear_adjoint();
    REQUIRE( h.values( "velocity", 0.3 )( 1 ) == 8.0 );
}

TEST_CASE( "Cluster adjoint history missing reads do not allocate entries", "[FEM][TimeDG][ISV][AdjointHistory]" )
{
    // Reverse assembly can request absent variables, times and quadrature points.
    // These reads must leave both empty and populated history banks unchanged.
    Cluster_Time_History tHistory;

    for ( uint i = 0; i < 128; ++i )
    {
        CHECK( tHistory.take_adjoint( "missing", i, { 2, 3 } ).numel() == 0 );
    }

    CHECK( tHistory.mAdjoint.empty() );
    const Matrix< DDRMat > tValue = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
    tHistory.set_adjoint( "velocity", 0.1, { 2, 3 }, tValue );
    CHECK( tHistory.take_adjoint( "other", 0.1, { 2, 3 } ).numel() == 0 );

    for ( uint i = 1; i <= 128; ++i )
    {
        CHECK( tHistory.take_adjoint( "velocity", i, { 2, 3 } ).numel() == 0 );
    }

    CHECK( tHistory.take_adjoint( "velocity", 0.1, { 2, 4 } ).numel() == 0 );
    CHECK( tHistory.take_adjoint( "velocity", 0.1, { 1, 3 } ).numel() == 0 );
    REQUIRE( tHistory.mAdjoint.size() == 1 );
    const auto& tTimes = tHistory.mAdjoint.at( "velocity" );
    REQUIRE( tTimes.size() == 1 );
    REQUIRE( tTimes.at( 0.1 ).size() == 1 );
    CHECK( norm( tTimes.at( 0.1 ).at( { 2, 3 } ) - tValue ) == 0.0 );
}

TEST_CASE( "Cluster adjoint history prunes consumed storage and preserves forward states", "[FEM][TimeDG][ISV][AdjointHistory]" )
{
    // Distinct element/point keys, times and names remain independent. Three RHS
    // columns must survive moving out the value and removing its map containers.
    Cluster_Time_History   tHistory;
    const Matrix< DDRMat > tValue  = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
    const Matrix< DDRMat > tSecond = 2.0 * tValue;
    const Matrix< DDRMat > tThird  = 3.0 * tValue;
    const Matrix< DDRMat > tLater  = 4.0 * tValue;
    const Matrix< DDRMat > tOther  = 5.0 * tValue;
    tHistory.stage( "velocity", 0.0, 0.1, tValue );
    tHistory.commit( "velocity" );
    tHistory.stage( "velocity", 0.1, 0.2, tSecond );
    tHistory.commit( "velocity" );
    tHistory.stage( "velocity", 0.2, 0.3, tThird );
    tHistory.set_adjoint( "velocity", 0.1, { 2, 3 }, tValue );
    tHistory.set_adjoint( "velocity", 0.1, { 2, 4 }, tSecond );
    tHistory.set_adjoint( "velocity", 0.1, { 1, 3 }, tThird );
    tHistory.set_adjoint( "velocity", 0.2, { 2, 3 }, tLater );
    tHistory.set_adjoint( "other", 0.1, { 2, 3 }, tOther );

    const Matrix< DDRMat > tTaken = tHistory.take_adjoint( "velocity", 0.1, { 2, 3 } );
    CHECK( norm( tTaken - tValue ) == 0.0 );
    REQUIRE( tHistory.mAdjoint.at( "velocity" ).at( 0.1 ).size() == 2 );
    CHECK( tHistory.take_adjoint( "velocity", 0.1, { 2, 3 } ).numel() == 0 );
    CHECK( norm( tHistory.take_adjoint( "velocity", 0.1, { 1, 3 } ) - tThird ) == 0.0 );
    CHECK( norm( tHistory.take_adjoint( "velocity", 0.1, { 2, 4 } ) - tSecond ) == 0.0 );
    REQUIRE( tHistory.mAdjoint.at( "velocity" ).size() == 1 );
    CHECK( tHistory.mAdjoint.at( "velocity" ).count( 0.1 ) == 0 );
    CHECK( norm( tHistory.take_adjoint( "velocity", 0.2, { 2, 3 } ) - tLater ) == 0.0 );
    CHECK( tHistory.mAdjoint.count( "velocity" ) == 0 );
    REQUIRE( tHistory.mAdjoint.size() == 1 );
    CHECK( tHistory.take_adjoint( "velocity", 0.2, { 2, 3 } ).numel() == 0 );
    CHECK( tHistory.mAdjoint.count( "velocity" ) == 0 );
    CHECK( norm( tHistory.take_adjoint( "other", 0.1, { 2, 3 } ) - tOther ) == 0.0 );
    CHECK( tHistory.mAdjoint.empty() );

    // Consuming adjoint carries must preserve both committed snapshots and a
    // separately staged forward state, which can still be committed afterwards.
    CHECK( norm( tTaken - tValue ) == 0.0 );
    CHECK( norm( tHistory.values( "velocity", 0.1 ) - tValue ) == 0.0 );
    CHECK( norm( tHistory.values( "velocity", 0.2 ) - tSecond ) == 0.0 );
    CHECK( tHistory.start( "velocity", 0.2 ) == 0.1 );
    tHistory.commit( "velocity" );
    CHECK( norm( tHistory.values( "velocity", 0.3 ) - tThird ) == 0.0 );
    CHECK( tHistory.start( "velocity", 0.3 ) == 0.2 );

    tHistory.set_adjoint( "velocity", 0.2, { 2, 3 }, tLater );
    tHistory.clear_adjoint();
    CHECK( tHistory.mAdjoint.empty() );
    CHECK( norm( tTaken - tValue ) == 0.0 );
    CHECK( norm( tHistory.values( "velocity", 0.3 ) - tThird ) == 0.0 );
}

TEST_CASE( "Time DG mass IWG current Jacobian and committed history", "[FEM][TimeDG][IWG]" )
{
    Mass_Fixture f;
    f.iwg->compute_jacobian( 0.63 );
    const Matrix< DDRMat > analytic = f.set.mJacobian;
    const Matrix< DDRMat > coeff    = f.current.get_coeff();
    for ( uint c = 0; c < 2; ++c )
    {
        for ( uint i = 0; i < 8; ++i )
        {
            Matrix< DDRMat > plus;
            for ( int sign : { 1, -1 } )
            {
                Matrix< DDRMat > perturbed = coeff;
                perturbed( i, c ) += sign * 1e-7;
                f.current.set_coeff( perturbed );
                f.set.mResidual( 0 ).fill( 0.0 );
                f.iwg->reset_eval_flags();
                f.iwg->compute_residual( 0.63 );
                if ( sign == 1 )
                {
                    plus = f.set.mResidual( 0 );
                }
                else
                {
                    Matrix< DDRMat > fd = ( plus - f.set.mResidual( 0 ) ) / 2e-7;
                    REQUIRE( norm( fd - analytic.get_column( c * 8 + i ) ) < 1e-7 );
                }
            }
        }
    }
    f.current.set_coeff( coeff );
    f.iwg->commit_time_history();
    f.cluster.get_cluster_ISV_data().reset();
    const auto& v = f.cluster.get_cluster_ISV_data().mTimeHistory.values( "TimeDG/mass/Velocity", 0.2 );
    for ( uint c = 0; c < 2; ++c )
    {
        for ( uint i = 0; i < 4; ++i )
        {
            REQUIRE( v( i, c ) == Approx( ( -3.0 * coeff( i, c ) + coeff( i + 4, c ) ) / 0.2 ) );
        }
    }
}

// The fixed-mass elimination must reject dependencies on other solution fields
// and invalid reference thickness before returning an incomplete state tangent.
TEST_CASE( "Time DG mass retains linear-time restriction", "[FEM][TimeDG][IWG]" )
{
    // General spatial lumping supports higher time orders; the energy-corrected
    // displacement/velocity elimination still requires exactly two time bases.
    Mass_Fixture tFixture;

    for ( auto tOrder : { mtk::Interpolation_Order::CONSTANT, mtk::Interpolation_Order::QUADRATIC } )
    {
        const mtk::Interpolation_Rule tRule( mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                tOrder == mtk::Interpolation_Order::CONSTANT ? mtk::Interpolation_Type::CONSTANT : mtk::Interpolation_Type::LAGRANGE,
                tOrder );
        Field_Interpolator            tField( 2, tRule, &tFixture.geometry, { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tField.set_coeff( Matrix< DDRMat >( 4 * tField.get_number_of_time_bases(), 2, 0.01 ) );
        tFixture.manager.mFI = { &tField };
        tFixture.manager.set_space_time( { { 0.23 }, { -0.31 }, { -1.0 } } );
        tFixture.iwg->reset_eval_flags();

        CHECK_THROWS_WITH( tFixture.iwg->compute_residual( 0.63 ),
                Catch::Contains( "Energy-corrected DG(1) requires two temporal displacement bases." ) );

        tFixture.manager.mFI = { &tFixture.current };
    }
}

TEST_CASE( "Time DG mass rejects state-dependent and invalid mass properties", "[FEM][TimeDG][IWG]" )
{
    Mass_Fixture f;
    auto         property = std::make_shared< Property >();
    property->set_parameters( { Matrix< DDRMat >{ { 1.0 } } } );
    std::string name = "Density";

    SECTION( "density depends on a coupled temperature field" )
    {
        property->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    }
    SECTION( "thickness depends on displacement" )
    {
        name = "Thickness";
        property->set_dof_type_list( f.types );
    }
    SECTION( "zero thickness makes mass singular" )
    {
        name = "Thickness";
        property->set_parameters( { Matrix< DDRMat >{ { 0.0 } } } );
    }
    SECTION( "vector thickness is not a scalar mass factor" )
    {
        name = "Thickness";
        property->set_parameters( { Matrix< DDRMat >{ { 1.0, 2.0 } } } );
    }

    f.iwg->set_property( property, name );
    REQUIRE_THROWS( f.iwg->compute_jacobian( 0.63 ) );
}

TEST_CASE( "Time DG nonlinear multislab adjoints include the second preceding slab", "[FEM][TimeDG][Adjoint]" )
{
    const std::array< real, 5 > dt = { 0.08, 0.11, 0.065, 0.09, 0.12 };
    const auto&                 op = Time_DG_Operators::get();
    for ( bool nonlinear : { false, true } )
    {
        const real            density = 1.3;
        const real            height  = 0.8;
        auto                  steps   = trajectory( density, height, nonlinear );
        std::array< Pair, 5 > lambda{};
        for ( int n = 4; n >= 0; --n )
        {
            Pair rhs = { 0.0, n == 4 ? 1.0 : 0.0 };
            for ( uint b = 0; b < 2; ++b )
            {
                for ( uint a = 0; a < 2; ++a )
                {
                    if ( n < 4 )
                    {
                        const real B = density * height * ( -( b == 1 ? op.Gb( a ) / dt[ n + 1 ] : 0.0 ) - op.e( a ) * op.A( 1, b ) / dt[ n ] );
                        rhs[ b ] -= B * lambda[ n + 1 ][ a ];
                    }
                    if ( n < 3 && b == 1 )
                    {
                        const real C = density * height * op.e( a ) * op.b( 1 ) / dt[ n + 1 ];
                        rhs[ b ] -= C * lambda[ n + 2 ][ a ];
                    }
                }
            }
            const auto& J = steps[ n ].jacobian;
            lambda[ n ]   = solve2( Block{ Pair{ J[ 0 ][ 0 ], J[ 1 ][ 0 ] }, Pair{ J[ 0 ][ 1 ], J[ 1 ][ 1 ] } }, rhs );
        }
        for ( uint parameter = 0; parameter < 2; ++parameter )
        {
            const real eps     = 1e-6;
            real       adjoint = 0.0;
            for ( uint n = 0; n < 5; ++n )
            {
                const real u     = n == 0 ? 0.13 : steps[ n - 1 ].q[ 1 ];
                const real v     = n == 0 ? 0.7 : steps[ n - 1 ].v[ 1 ];
                auto       plus  = evaluate( steps[ n ].q, u, v, dt[ n ], density + ( parameter == 0 ? eps : 0 ), height + ( parameter == 1 ? eps : 0 ), nonlinear );
                auto       minus = evaluate( steps[ n ].q, u, v, dt[ n ], density - ( parameter == 0 ? eps : 0 ), height - ( parameter == 1 ? eps : 0 ), nonlinear );
                for ( uint a = 0; a < 2; ++a )
                {
                    adjoint -= lambda[ n ][ a ] * ( plus.residual[ a ] - minus.residual[ a ] ) / ( 2 * eps );
                }
            }
            auto       plus  = trajectory( density + ( parameter == 0 ? eps : 0 ), height + ( parameter == 1 ? eps : 0 ), nonlinear );
            auto       minus = trajectory( density - ( parameter == 0 ? eps : 0 ), height - ( parameter == 1 ? eps : 0 ), nonlinear );
            const real fd    = ( plus.back().q[ 1 ] - minus.back().q[ 1 ] ) / ( 2 * eps );
            REQUIRE( adjoint == Approx( fd ).epsilon( 2e-5 ).margin( 1e-8 ) );
        }
    }
}

namespace
{
    class Spatial_Force_Probe : public IWG
    {
      public:
        void compute_residual( real aWeight ) override
        {
            auto*            fi    = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            Matrix< DDRMat > force = fi->val();
            for ( uint c = 0; c < force.numel(); ++c )
            {
                force( c ) += std::pow( force( c ), 3 );
            }
            mSet->get_residual()( 0 ) += aWeight * fi->N_trans() * force;
        }
        void compute_jacobian( real aWeight ) override
        {
            auto*            fi = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            Matrix< DDRMat > stiffness( 2, 2, 0.0 );
            for ( uint c = 0; c < 2; ++c )
            {
                stiffness( c, c ) = 1.0 + 3.0 * std::pow( fi->val()( c ), 2 );
            }
            mSet->get_jacobian() += aWeight * fi->N_trans() * stiffness * fi->N();
        }
        void compute_jacobian_and_residual( real aWeight ) override
        {
            compute_residual( aWeight );
            compute_jacobian( aWeight );
        }
        void compute_dRdp( real ) override {}
    };
}    // namespace

TEST_CASE( "Time DG spatial IWG reduced tests endpoint and tangent", "[FEM][TimeDG][IWG]" )
{
    Mass_Fixture f;
    auto         wrapped = std::make_shared< IWG_Time_DG_Spatial >( std::make_shared< Spatial_Force_Probe >() );
    wrapped->set_set_pointer( &f.set );
    wrapped->set_residual_dof_type( f.types );
    wrapped->set_dof_type_list( f.types );
    wrapped->set_field_interpolator_manager( &f.manager );
    f.manager.set_space_time( Matrix< DDRMat >{ { 0.13 }, { -0.27 }, { 0.41 } } );
    const Matrix< DDRMat > point = f.current.get_space_time();
    const Matrix< DDRMat > coeff = f.current.get_coeff();
    Matrix< DDRMat >       n( 1, 4 );
    for ( uint i = 0; i < 4; ++i )
    {
        n( i ) = f.current.NBuild()( i ) + f.current.NBuild()( i + 4 );
    }
    const Matrix< DDRMat > value = f.current.val();
    wrapped->compute_jacobian_and_residual( 0.73 );
    const Matrix< DDRMat > jac      = f.set.mJacobian;
    const Matrix< DDRMat > residual = f.set.mResidual( 0 );
    REQUIRE( norm( f.current.get_space_time() - point ) < 1e-14 );
    for ( uint c = 0; c < 2; ++c )
    {
        real incoming = 0.0;
        for ( uint i = 0; i < 4; ++i )
        {
            incoming += n( i ) * coeff( i, c );
        }
        for ( uint i = 0; i < 4; ++i )
        {
            CHECK( residual( c * 8 + i ) == Approx( 0.73 * n( i ) * ( value( c ) + std::pow( value( c ), 3 ) ) ) );
            CHECK( residual( c * 8 + 4 + i ) == Approx( 0.73 * n( i ) * ( incoming + std::pow( incoming, 3 ) ) ) );
        }
        for ( uint i = 0; i < 8; ++i )
        {
            Matrix< DDRMat > plus;
            for ( int sign : { 1, -1 } )
            {
                Matrix< DDRMat > perturbed = coeff;
                perturbed( i, c ) += sign * 1e-7;
                f.current.set_coeff( perturbed );
                f.set.mResidual( 0 ).fill( 0.0 );
                wrapped->compute_residual( 0.73 );
                if ( sign == 1 )
                {
                    plus = f.set.mResidual( 0 );
                }
                else
                {
                    Matrix< DDRMat > fd = ( plus - f.set.mResidual( 0 ) ) / 2e-7;
                    CHECK( norm( fd - jac.get_column( c * 8 + i ) ) < 1e-8 );
                }
            }
        }
    }
}

TEST_CASE( "Time DG history IWG previous tangent and delayed adjoint multiple RHS", "[FEM][TimeDG][Adjoint][IWG]" )
{
    Mass_Fixture           f;
    const Matrix< DDRMat > coeff       = f.previous.get_coeff();
    auto&                  history     = f.cluster.get_cluster_ISV_data().mTimeHistory;
    auto                   set_history = [ & ]( const Matrix< DDRMat >& q ) {
        history.clear();
        history.stage( "TimeDG/mass/Velocity", 0.0, 0.2, Matrix< DDRMat >( 4, 2, 0.0 ) );
        history.commit( "TimeDG/mass/Velocity" );
        Matrix< DDRMat > velocity( 4, 2 );
        for ( uint c = 0; c < 2; ++c )
        {
            for ( uint i = 0; i < 4; ++i )
            {
                velocity( i, c ) = ( -3.0 * q( i, c ) + q( i + 4, c ) + 0.12 ) / 0.3;
            }
        }
        history.stage( "TimeDG/mass/Velocity", 0.2, 0.5, velocity );
        history.commit( "TimeDG/mass/Velocity" );
    };
    f.time = { { 0.5 }, { 0.9 } };
    f.model.set_time( f.time );
    set_history( coeff );
    f.iwg->compute_jacobian_previous( 0.63 );
    const Matrix< DDRMat > jac = f.set.mJacobian;
    for ( uint c = 0; c < 2; ++c )
    {
        for ( uint i = 0; i < 8; ++i )
        {
            Matrix< DDRMat > plus;
            for ( int sign : { 1, -1 } )
            {
                Matrix< DDRMat > q = coeff;
                q( i, c ) += sign * 1e-7;
                f.previous.set_coeff( q );
                set_history( q );
                f.set.mResidual( 0 ).fill( 0.0 );
                f.iwg->compute_residual( 0.63 );
                if ( sign == 1 )
                {
                    plus = f.set.mResidual( 0 );
                }
                else
                {
                    Matrix< DDRMat > fd = ( plus - f.set.mResidual( 0 ) ) / 2e-7;
                    CHECK( norm( fd - jac.get_column( c * 8 + i ) ) < 1e-7 );
                }
            }
        }
    }
    set_history( coeff );
    Vector< Matrix< DDRMat > > adjoint = { Matrix< DDRMat >( 16, 1, 0.0 ), Matrix< DDRMat >( 16, 1, 0.0 ) };
    Vector< Matrix< DDRMat > > rhs     = adjoint;
    for ( uint r = 0; r < 2; ++r )
    {
        for ( uint i = 0; i < 16; ++i )
        {
            adjoint( r )( i ) = 0.1 * ( 1 + i + 3 * r );
        }
    }
    f.iwg->add_adjoint_history( 0.63, 2, 3, adjoint, rhs );
    REQUIRE( norm( rhs( 0 ) ) == 0.0 );
    f.time = { { 0.2 }, { 0.5 } };
    f.model.set_time( f.time );
    Vector< Matrix< DDRMat > > zero = { Matrix< DDRMat >( 16, 1, 0.0 ), Matrix< DDRMat >( 16, 1, 0.0 ) };
    f.iwg->add_adjoint_history( 0.63, 2, 3, zero, rhs );
    for ( uint r = 0; r < 2; ++r )
    {
        for ( uint c = 0; c < 2; ++c )
        {
            real contraction = 0.0;
            for ( uint i = 0; i < 4; ++i )
            {
                contraction += f.current.NBuild()( i ) * ( adjoint( r )( c * 8 + i ) + 4.0 * adjoint( r )( c * 8 + 4 + i ) );
            }
            for ( uint i = 0; i < 4; ++i )
            {
                CHECK( rhs( r )( c * 8 + i ) == 0.0 );
                CHECK( rhs( r )( c * 8 + 4 + i ) == Approx( -0.63 * 1.7 * 2.0 / 0.3 * f.current.NBuild()( i ) * contraction ) );
            }
        }
    }
}

TEST_CASE( "Time DG conservative energy identity on unequal slabs", "[FEM][TimeDG][Energy]" )
{
    const std::array< real, 5 > dt = { 0.08, 0.11, 0.065, 0.09, 0.12 };
    for ( bool nonlinear : { false, true } )
    {
        const real mass      = 1.3 * 0.8;
        const real stiffness = std::pow( 0.8, 3 );
        const real beta      = nonlinear ? 7.0 : 0.0;
        auto       potential = [ & ]( real q ) { return stiffness * ( 0.5 * q * q + 0.25 * beta * std::pow( q, 4 ) ); };
        auto       force     = [ & ]( real q ) { return stiffness * ( q + beta * std::pow( q, 3 ) ); };
        auto       steps     = trajectory( 1.3, 0.8, nonlinear );
        real       u         = 0.13;
        real       v         = 0.7;
        for ( uint n = 0; n < steps.size(); ++n )
        {
            const auto& s  = steps[ n ];
            const real  ju = s.q[ 0 ] - u;
            const real  jv = s.v[ 0 ] - v;
            const real  dE = 0.5 * mass * ( s.v[ 1 ] * s.v[ 1 ] - v * v ) + potential( s.q[ 1 ] ) - potential( u );
            const real  EJ = 0.5 * mass * jv * jv + force( s.q[ 0 ] ) * ju - ( potential( s.q[ 0 ] ) - potential( u ) );
            const real  ER = ( mass * ( s.v[ 1 ] - s.v[ 0 ] ) / dt[ n ] + force( s.q[ 0 ] ) ) * ju + mass * ( s.v[ 0 ] - ( s.q[ 1 ] - s.q[ 0 ] ) / dt[ n ] ) * jv;
            CHECK( std::abs( ER ) < 2e-10 );
            CHECK( std::abs( dE + EJ - ER ) < 2e-10 );
            CHECK( EJ >= -2e-12 );
            u = s.q[ 1 ];
            v = s.v[ 1 ];
        }
    }
}
