/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 */

#include <catch.hpp>
#include <memory>

#define protected public
#define private public
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_IQI_Volume.hpp"
#include "fn_FEM_FD_Scheme.hpp"

using namespace moris;
using namespace moris::fem;

TEST_CASE( "Three-point one-sided FD accuracy", "[FEM][FD]" )
{
    for ( FDScheme_Type tType : { FDScheme_Type::POINT_3_FORWARD, FDScheme_Type::POINT_3_BACKWARD } )
    {
        Vector< Vector< real > > tScheme;
        fd_scheme( tType, tScheme );
        const real tX = 0.7;
        const real tH = 0.01;
        real tDerivative = 0.0;
        for ( uint i = 0; i < tScheme( 0 ).size(); ++i )
        {
            const real tPoint = tX + tScheme( 0 )( i ) * tH;
            tDerivative += tScheme( 1 )( i ) * ( 2.0 + 3.0 * tPoint + tPoint * tPoint ) /
                           ( tScheme( 2 )( 0 ) * tH );
        }
        REQUIRE( tDerivative == Approx( 3.0 + 2.0 * tX ).margin( 1e-11 ) );
    }
}

TEST_CASE( "IQI geometry FD stays inside background element", "[FEM][FD]" )
{
    mtk::Interpolation_Rule tRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
    Geometry_Interpolator tGI( tRule );
    tGI.set_coeff( { { -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 } },
            { { 0.0 }, { 1.0 } } );

    Set tSet;
    tSet.set_set_type( Element_Type::BULK );
    tSet.set_perturbation_strategy( Perturbation_Type::ABSOLUTE );
    Vector< Vector< MSI::Dof_Type > > tDofs;
    Vector< Vector< gen::PDV_Type > > tDvs;
    Vector< Vector< mtk::Field_Type > > tFields;
    auto* tManager = new Field_Interpolator_Manager( tDofs, tDvs, tFields, &tSet );
    tManager->mIPGeometryInterpolator = &tGI;
    tSet.mLeaderFIManager = tManager;

    auto tIQI = std::make_shared< IQI_Volume >();
    tIQI->set_set_pointer( &tSet );

    for ( uint tDirection = 0; tDirection < 2; ++tDirection )
    {
        for ( real tX : { -1.0, -0.985, -0.5, 0.0, 0.5, 0.985, 1.0 } )
        {
            FDScheme_Type tScheme = FDScheme_Type::POINT_5;
            const real tH = tIQI->check_ig_coordinates_inside_ip_element( 0.01, tX, tDirection, tScheme );
            const FDScheme_Type tExpected = tX < -0.98 ? FDScheme_Type::POINT_3_FORWARD
                                                       : ( tX > 0.98 ? FDScheme_Type::POINT_3_BACKWARD : FDScheme_Type::POINT_5 );
            REQUIRE( tScheme == tExpected );

            Vector< Vector< real > > tStencil;
            fd_scheme( tScheme, tStencil );
            for ( real tOffset : tStencil( 0 ) )
            {
                REQUIRE( tX + tOffset * tH >= -1.0 );
                REQUIRE( tX + tOffset * tH <= 1.0 );
            }
        }
    }

    // The lower-order default retains its original boundary fallbacks.
    for ( real tX : { -1.0, 1.0 } )
    {
        FDScheme_Type tScheme = FDScheme_Type::POINT_3_CENTRAL;
        tIQI->check_ig_coordinates_inside_ip_element( 0.01, tX, 0, tScheme );
        REQUIRE( tScheme == ( tX < 0.0 ? FDScheme_Type::POINT_1_FORWARD : FDScheme_Type::POINT_1_BACKWARD ) );
    }

    // Large requested steps are capped so the two-step stencil still fits.
    FDScheme_Type tScheme = FDScheme_Type::POINT_5;
    const real tH = tIQI->check_ig_coordinates_inside_ip_element( 10.0, 0.0, 0, tScheme );
    REQUIRE( tScheme == FDScheme_Type::POINT_5 );
    REQUIRE( 2.0 * tH < 1.0 );
}
