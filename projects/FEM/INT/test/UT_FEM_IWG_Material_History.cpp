/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */

#include <catch.hpp>

#define protected public
#define private public
#include "cl_FEM_Set.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#undef private
#undef protected

#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "fn_norm.hpp"

namespace moris::fem
{
    namespace
    {
        // A quadrature residual with independent current/previous material
        // contractions supplies exact derivative columns on either interface
        // side. It isolates the generic FD driver from a particular contact law.
        class Material_History_IWG : public IWG
        {
          public:
            bool mLeaderActive   = false;
            bool mFollowerActive = false;

            Material_History_IWG()
            {
                mTimeContinuity  = true;
                mResidualDofType = { { MSI::Dof_Type::TEMP } };
            }

            void compute_residual( real aWeight ) override
            {
                for ( uint s = 0; s < 2; ++s )
                {
                    if ( !( s == 0 ? mLeaderActive : mFollowerActive ) )
                    {
                        continue;
                    }
                    auto*       tManager  = s == 0 ? mLeaderFIManager : mFollowerFIManager;
                    auto*       tPrevious = s == 0 ? mLeaderPreviousFIManager : mFollowerPreviousFIManager;
                    const auto& tCurrent  = tManager->get_field_interpolators_for_type( gen::PDV_Type::ELASTIC_MODULUS )->get_coeff();
                    const auto& tHistory  = tPrevious->get_field_interpolators_for_type( gen::PDV_Type::ELASTIC_MODULUS )->get_coeff();
                    for ( uint i = 0; i < 4; ++i )
                    {
                        mSet->get_residual()( 0 )( 0 ) += aWeight * ( ( s == 0 ? 2.0 : 11.0 ) * ( i + 1 ) * tCurrent( i ) + ( s == 0 ? 3.0 : 13.0 ) * ( 4 - i ) * tHistory( i ) );
                        mSet->get_residual()( 0 )( 1 ) += aWeight * ( ( s == 0 ? 5.0 : 17.0 ) * ( i + 1 ) * tCurrent( i ) + ( s == 0 ? 7.0 : 19.0 ) * ( 4 - i ) * tHistory( i ) );
                    }
                }
            }

            void compute_jacobian( real ) override {}
            void compute_jacobian_and_residual( real ) override {}
            void compute_dRdp( real ) override {}
        };
    }    // namespace

    // Follower-only material dependencies used to return before differentiation.
    // History terms also need the same design perturbation on both slabs. Check
    // leader-only, follower-only and two-sided columns with all FD schemes, and
    // verify residual/coefficients are restored for later assembly operations.
    TEST_CASE( "Interface material partials include history and follower-only fields", "[FEM][MaterialHistory]" )
    {
        for ( uint tSides = 1; tSides <= 3; ++tSides )
        {
            Set                  tSet;
            Material_History_IWG tIWG;

            tSet.set_set_type( Element_Type::NONCONFORMAL_SIDESET );
            tIWG.mLeaderActive   = ( tSides & 1 ) != 0;
            tIWG.mFollowerActive = ( tSides & 2 ) != 0;
            tIWG.set_set_pointer( &tSet );

            const uint tDofMapSize = static_cast< uint >( MSI::Dof_Type::END_ENUM ) + 1;
            tSet.mLeaderDofTypeMap.set_size( tDofMapSize, 1, -1 );
            tSet.mFollowerDofTypeMap                                               = tSet.mLeaderDofTypeMap;
            tSet.mLeaderDofTypeMap( static_cast< uint >( MSI::Dof_Type::TEMP ) )   = 0;
            tSet.mFollowerDofTypeMap( static_cast< uint >( MSI::Dof_Type::TEMP ) ) = 0;
            tSet.mResDofAssemblyMap                                                = { Matrix< DDSMat >{ { 0, 0 } }, Matrix< DDSMat >{ { 1, 1 } } };
            tSet.mResidual                                                         = { Matrix< DDRMat >{ { 123.0 }, { -27.0 } } };

            const uint tDvMapSize = static_cast< uint >( gen::PDV_Type::UNDEFINED ) + 1;
            tSet.mLeaderDvTypeMap.set_size( tDvMapSize, 1, -1 );
            tSet.mFollowerDvTypeMap                        = tSet.mLeaderDvTypeMap;
            const Vector< Vector< gen::PDV_Type > > tTypes = { { gen::PDV_Type::ELASTIC_MODULUS } };
            if ( tIWG.mLeaderActive )
            {
                tSet.mLeaderDvTypes                                                            = tTypes;
                tSet.mLeaderDvTypeMap( static_cast< uint >( gen::PDV_Type::ELASTIC_MODULUS ) ) = 0;
            }
            if ( tIWG.mFollowerActive )
            {
                tSet.mFollowerDvTypes                                                            = tTypes;
                tSet.mFollowerDvTypeMap( static_cast< uint >( gen::PDV_Type::ELASTIC_MODULUS ) ) = 0;
            }
            const uint tBlocks = tIWG.mLeaderActive + tIWG.mFollowerActive;
            tSet.mPdvMatAssemblyMap.resize( tBlocks );
            for ( uint b = 0; b < tBlocks; ++b )
            {
                tSet.mPdvMatAssemblyMap( b ) = { { static_cast< sint >( 4 * b ), static_cast< sint >( 4 * b + 3 ) } };
            }
            tSet.mdRdp = { Matrix< DDRMat >( 2, 4 * tBlocks, 0.0 ) };

            mtk::Interpolation_Rule tGeometryRule( mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );
            Geometry_Interpolator   tGeometry( tGeometryRule );
            mtk::Interpolation_Rule tFieldRule( mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::CONSTANT,
                    mtk::Interpolation_Order::CONSTANT );
            const Matrix< DDRMat >  tCoeff = { { 2.0 }, { 3.0 }, { 4.0 }, { 5.0 } };

            Field_Interpolator_Manager tLeaderManager( {}, {}, {}, &tSet );
            Field_Interpolator_Manager tFollowerManager( {}, {}, {}, &tSet, mtk::Leader_Follower::FOLLOWER );
            Field_Interpolator_Manager tLeaderHistoryManager( {}, {}, {}, &tSet );
            Field_Interpolator_Manager tFollowerHistoryManager( {}, {}, {}, &tSet, mtk::Leader_Follower::FOLLOWER );
            // The managers own their heap-allocated interpolators. The Set is
            // deliberately not given these stack manager pointers.
            tLeaderManager.mDvFI          = { new Field_Interpolator( 1, tFieldRule, &tGeometry, tTypes( 0 ) ) };
            tFollowerManager.mDvFI        = { new Field_Interpolator( 1, tFieldRule, &tGeometry, tTypes( 0 ) ) };
            tLeaderHistoryManager.mDvFI   = { new Field_Interpolator( 1, tFieldRule, &tGeometry, tTypes( 0 ) ) };
            tFollowerHistoryManager.mDvFI = { new Field_Interpolator( 1, tFieldRule, &tGeometry, tTypes( 0 ) ) };
            for ( auto* tManager : { &tLeaderManager, &tFollowerManager, &tLeaderHistoryManager, &tFollowerHistoryManager } )
            {
                tManager->mDvFI( 0 )->set_coeff( tCoeff );
            }
            tIWG.set_field_interpolator_manager( &tLeaderManager );
            tIWG.set_field_interpolator_manager( &tFollowerManager, mtk::Leader_Follower::FOLLOWER );
            tIWG.set_field_interpolator_manager_previous_time( &tLeaderHistoryManager );
            tIWG.set_field_interpolator_manager_previous_time( &tFollowerHistoryManager, mtk::Leader_Follower::FOLLOWER );

            const real       tWeight = 1.7;
            Matrix< DDRMat > tExpected( 2, 4 * tBlocks, 0.0 );
            uint             tBlock = 0;
            for ( uint s = 0; s < 2; ++s )
            {
                if ( !( s == 0 ? tIWG.mLeaderActive : tIWG.mFollowerActive ) )
                {
                    continue;
                }
                for ( uint i = 0; i < 4; ++i )
                {
                    tExpected( 0, 4 * tBlock + i ) = tWeight * ( ( s == 0 ? 2.0 : 11.0 ) * ( i + 1 ) + ( s == 0 ? 3.0 : 13.0 ) * ( 4 - i ) );
                    tExpected( 1, 4 * tBlock + i ) = tWeight * ( ( s == 0 ? 5.0 : 17.0 ) * ( i + 1 ) + ( s == 0 ? 7.0 : 19.0 ) * ( 4 - i ) );
                }
                ++tBlock;
            }
            for ( auto tScheme : { FDScheme_Type::POINT_1_FORWARD, FDScheme_Type::POINT_1_BACKWARD, FDScheme_Type::POINT_3_CENTRAL } )
            {
                tSet.mdRdp( 0 ).fill( 0.0 );
                tIWG.select_dRdp_FD_material_double( tWeight, 1e-5, tScheme );
                CHECK( norm( tSet.mdRdp( 0 ) - tExpected ) / norm( tExpected ) < 1e-8 );
                CHECK( norm( tSet.mResidual( 0 ) - Matrix< DDRMat >{ { 123.0 }, { -27.0 } } ) == 0.0 );
                for ( auto* tManager : { &tLeaderManager, &tFollowerManager, &tLeaderHistoryManager, &tFollowerHistoryManager } )
                {
                    CHECK( norm( tManager->mDvFI( 0 )->get_coeff() - tCoeff ) == 0.0 );
                }
            }
        }
    }
}    // namespace moris::fem
