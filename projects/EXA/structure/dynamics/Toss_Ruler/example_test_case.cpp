/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <memory>
#include <limits>
#include <string>
#include <vector>
#include "catch.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

using namespace moris;
extern "C" {
uint gTossElementsX        = 12;
uint gTossElementsY        = 2;
uint gTossTimeSteps        = 240;
real gTossFinalTime        = 0.012;
bool gPrintReferenceValues = false;
}

namespace
{
    constexpr real tLength        = 0.3;
    constexpr real tHeight        = 0.06;
    constexpr real tThickness     = 0.002;
    constexpr real tDensity       = 7800.0;
    constexpr real tYoungs        = 2.06e11;
    constexpr real tPulseDuration = 0.004;
    constexpr real tMass          = tDensity * tThickness * tLength * tHeight;

    struct Diagnostics
    {
        real px        = 0.0;
        real py        = 0.0;
        real angular   = 0.0;
        real kinetic   = 0.0;
        real strain    = 0.0;
        real cx        = 0.0;
        real cy        = 0.0;
        real min_det_f = std::numeric_limits< real >::infinity();
    };

    // Integrate the exported Q2 fields with MORIS interpolation and quadrature.
    // This is post-processing only; assembly, material tangents and time stepping
    // all use the normal MORIS FEM/SOL workflow in the input library.
    Diagnostics integrate( const Matrix< DDRMat >& aDisplacement, const Matrix< DDRMat >& aVelocity )
    {
        mtk::Interpolation_Function_Factory                 tFactory;
        std::unique_ptr< mtk::Interpolation_Function_Base > tInterpolation(
                tFactory.create_interpolation_function( mtk::Geometry_Type::QUAD,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::QUADRATIC ) );
        Matrix< DDRMat > tNodeXi;
        tInterpolation->get_param_coords( tNodeXi );
        mtk::Integration_Rule tRule( mtk::Geometry_Type::QUAD, mtk::Integration_Type::GAUSS, mtk::Integration_Order::QUAD_5x5, mtk::Geometry_Type::POINT, mtk::Integration_Type::GAUSS, mtk::Integration_Order::POINT );
        mtk::Integrator       tIntegrator( tRule );
        Matrix< DDRMat >      tPoints, tWeights;
        tIntegrator.get_points( tPoints );
        tIntegrator.get_weights( tWeights );

        const real  tDX        = tLength / gTossElementsX;
        const real  tDY        = tHeight / gTossElementsY;
        const uint  tGridWidth = 2 * gTossElementsX + 1;
        Diagnostics tResult;
        for ( uint ey = 0; ey < gTossElementsY; ++ey )
        {
            for ( uint ex = 0; ex < gTossElementsX; ++ex )
            {
                Matrix< DDRMat > tU( 9, 2 ), tV( 9, 2 );
                for ( uint n = 0; n < 9; ++n )
                {
                    const uint ix = 2 * ex + std::lround( tNodeXi( 0, n ) + 1.0 );
                    const uint iy = 2 * ey + std::lround( tNodeXi( 1, n ) + 1.0 );
                    for ( uint c = 0; c < 2; ++c )
                    {
                        tU( n, c ) = aDisplacement( ix + tGridWidth * iy, c );
                        tV( n, c ) = aVelocity( ix + tGridWidth * iy, c );
                    }
                }
                for ( uint q = 0; q < tWeights.numel(); ++q )
                {
                    Matrix< DDRMat > tXi = { { tPoints( 0, q ) }, { tPoints( 1, q ) } };
                    Matrix< DDRMat > tN, tDerivative;
                    tInterpolation->eval_N( tXi, tN );
                    tInterpolation->eval_dNdXi( tXi, tDerivative );
                    tDerivative.get_row( 0 ) *= 2.0 / tDX;
                    tDerivative.get_row( 1 ) *= 2.0 / tDY;
                    const Matrix< DDRMat > tGrad     = tDerivative * tU;
                    const Matrix< DDRMat > tDisp     = tN * tU;
                    const Matrix< DDRMat > tVelocity = tN * tV;
                    const real             x         = ( ex + 0.5 * ( tXi( 0 ) + 1.0 ) ) * tDX + tDisp( 0 );
                    const real             y         = -0.5 * tHeight + ( ey + 0.5 * ( tXi( 1 ) + 1.0 ) ) * tDY + tDisp( 1 );
                    const real             vx        = tVelocity( 0 );
                    const real             vy        = tVelocity( 1 );
                    const real             dv        = tWeights( q ) * tDX * tDY * tThickness / 4.0;
                    const real             dm        = tDensity * dv;
                    const real             f11       = 1.0 + tGrad( 0, 0 );
                    const real             f12       = tGrad( 1, 0 );
                    const real             f21       = tGrad( 0, 1 );
                    const real             f22       = 1.0 + tGrad( 1, 1 );
                    const real             e11       = 0.5 * ( f11 * f11 + f21 * f21 - 1.0 );
                    const real             e22       = 0.5 * ( f12 * f12 + f22 * f22 - 1.0 );
                    const real             e12       = 0.5 * ( f11 * f12 + f21 * f22 );

                    // Plane-stress StVK with nu=0, matching the active CM.
                    tResult.strain += dv * 0.5 * tYoungs * ( e11 * e11 + e22 * e22 + 2.0 * e12 * e12 );
                    tResult.kinetic += 0.5 * dm * ( vx * vx + vy * vy );
                    tResult.px += dm * vx;
                    tResult.py += dm * vy;
                    tResult.angular += dm * ( x * vy - y * vx );
                    tResult.cx += dm * x / tMass;
                    tResult.cy += dm * y / tMass;
                    tResult.min_det_f = std::min( tResult.min_det_f, f11 * f22 - f12 * f21 );
                }
            }
        }
        return tResult;
    }

    void check_results( bool aFullFlight )
    {
        const uint            tGridWidth = 2 * gTossElementsX + 1;
        const uint            tGridSize  = tGridWidth * ( 2 * gTossElementsY + 1 );
        const real            tDT        = gTossFinalTime / gTossTimeSteps;
        Matrix< DDRMat >      tPrevious( tGridSize, 2, 0.0 );
        Diagnostics           tFinal, tReleased;
        std::array< real, 2 > tTip = { { 0.0, 0.0 } };
        std::ofstream         tCSV;
        if ( par_rank() == 0 )
        {
            tCSV.open( "Toss_Ruler_history.csv", std::ios::trunc );
            REQUIRE( tCSV.good() );
            tCSV << std::setprecision( 16 )
                 << "time,cx,cy,px,py,angular_momentum,kinetic_energy,strain_energy,total_energy,min_det_f,tip_ux,tip_uy\n"
                 << "0,0.15,0,0,0,0,0,0,0,1,0,0\n";
        }

        // One reader per rank, closed before leaving this scope. VIS node IDs
        // need not be contiguous; map duplicate cluster vertices by coordinates.
        {
            mtk::Exodus_IO_Helper tExo( "Toss_Ruler.exo", 0 );
            REQUIRE( tExo.get_number_of_dimensions() == 2 );
            if ( !aFullFlight )
            {
                // Local VIS counts from the 12x2 Q2 regression, verified for
                // both MPI partitions. Both ranks own six columns of cells.
                REQUIRE( ( par_size() == 1 || par_size() == 2 ) );
                REQUIRE( tExo.get_number_of_nodes() == ( par_size() == 1 ? 125 : 65 ) );
                REQUIRE( tExo.get_number_of_elements() == ( par_size() == 1 ? 24 : 12 ) );
            }
            else
            {
                REQUIRE( tExo.get_number_of_nodes() > 0 );
                REQUIRE( tExo.get_number_of_elements() > 0 );
            }
            std::vector< uint > tIds( tExo.get_number_of_nodes() );
            std::vector< uint > tGridIndices( tIds.size() );
            for ( uint n = 0; n < tIds.size(); ++n )
            {
                tIds[ n ]                 = tExo.get_node_id( n );
                const Matrix< DDRMat > tX = tExo.get_nodal_coordinate( tIds[ n ] );
                const uint             ix = std::lround( tX( 0 ) * ( tGridWidth - 1 ) / tLength );
                const uint             iy = std::lround( ( tX( 1 ) + 0.5 * tHeight ) * ( 2 * gTossElementsY ) / tHeight );
                REQUIRE( ix < tGridWidth );
                REQUIRE( iy <= 2 * gTossElementsY );
                REQUIRE( tX( 0 ) == Approx( ix * tLength / ( tGridWidth - 1 ) ).margin( 1e-12 ) );
                REQUIRE( tX( 1 ) == Approx( -0.5 * tHeight + iy * tHeight / ( 2 * gTossElementsY ) ).margin( 1e-12 ) );
                tGridIndices[ n ] = ix + tGridWidth * iy;
            }

            for ( uint step = 0; step < gTossTimeSteps; ++step )
            {
                // Four fields plus a multiplicity count. Summation and averaging
                // make reconstruction independent of VIS duplication and MPI cuts.
                Matrix< DDRMat > tLocal( tGridSize, 5, 0.0 );
                for ( uint n = 0; n < tIds.size(); ++n )
                {
                    const uint i = tGridIndices[ n ];
                    for ( uint c = 0; c < 4; ++c )
                    {
                        tLocal( i, c ) += tExo.get_nodal_field_value( tIds[ n ], c + 2, step );
                    }
                    tLocal( i, 4 ) += 1.0;
                }
                Matrix< DDRMat > tGlobal = sum_all_matrix( tLocal );
                Matrix< DDRMat > tU( tGridSize, 2 ), tV( tGridSize, 2 );
                for ( uint i = 0; i < tGridSize; ++i )
                {
                    REQUIRE( tGlobal( i, 4 ) > 0.0 );
                    for ( uint c = 0; c < 2; ++c )
                    {
                        tU( i, c )           = tGlobal( i, c ) / tGlobal( i, 4 );
                        const real tIncoming = tGlobal( i, c + 2 ) / tGlobal( i, 4 );
                        // Report's alpha=4 DG(1) outgoing velocity. In particular,
                        // this is not the temporal derivative of displacement.
                        tV( i, c ) = ( -3.0 * tIncoming + tU( i, c ) + 2.0 * tPrevious( i, c ) ) / tDT;
                    }
                }
                tFinal           = integrate( tU, tV );
                const real tTime = ( step + 1 ) * tDT;
                CHECK( tExo.get_time_value() == Approx( tTime ).margin( 1e-12 ) );
                CHECK( std::isfinite( tFinal.kinetic + tFinal.strain ) );
                CHECK( tFinal.min_det_f > 0.0 );
                if ( std::abs( tTime - tPulseDuration ) < 0.5 * tDT )
                {
                    tReleased = tFinal;
                }
                const uint tTipIndex = 2 * gTossElementsX + tGridWidth * gTossElementsY;
                tTip                 = { { tU( tTipIndex, 0 ), tU( tTipIndex, 1 ) } };
                if ( par_rank() == 0 )
                {
                    tCSV << tTime << ',' << tFinal.cx << ',' << tFinal.cy << ','
                         << tFinal.px << ',' << tFinal.py << ',' << tFinal.angular << ','
                         << tFinal.kinetic << ',' << tFinal.strain << ',' << tFinal.kinetic + tFinal.strain << ','
                         << tFinal.min_det_f << ',' << tTip[ 0 ] << ',' << tTip[ 1 ] << '\n';
                }
                tPrevious = tU;
            }
            if ( gPrintReferenceValues )
            {
                MORIS_LOG_INFO( "Toss ruler rank %u: nodes=%d elements=%d, tip=(%.16e,%.16e), E=%.16e, P=(%.16e,%.16e), H=%.16e",
                        par_rank(),
                        tExo.get_number_of_nodes(),
                        tExo.get_number_of_elements(),
                        tTip[ 0 ],
                        tTip[ 1 ],
                        tFinal.kinetic + tFinal.strain,
                        tFinal.px,
                        tFinal.py,
                        tFinal.angular );
            }
        }

        // Use a final-step reader after closing the history reader: the helper
        // shares its current time index between nodal and global field caches.
        {
            mtk::Exodus_IO_Helper tFinalExo( "Toss_Ruler.exo", gTossTimeSteps - 1 );
            CHECK( tFinalExo.get_time_value() == Approx( gTossFinalTime ).margin( 1e-12 ) );
            CHECK( tFinalExo.get_global_variable( 0, gTossTimeSteps - 1 ) == Approx( tLength * tHeight ).epsilon( 1e-10 ) );
        }

        // Net force is 0.5*b*R(t) in +y. For a symmetric triangular pulse,
        // its impulse acts at Tpulse/2. These references are mesh-independent.
        const real tImpulse = 0.5 * tThickness * 0.5 * 4.0e4 * tPulseDuration;
        CHECK( tFinal.px == Approx( 0.0 ).margin( 1e-7 ) );
        CHECK( tFinal.py == Approx( tImpulse ).epsilon( 1e-5 ) );
        CHECK( tFinal.cx == Approx( 0.5 * tLength ).margin( 1e-7 ) );
        CHECK( tFinal.cy == Approx( tImpulse / tMass * ( gTossFinalTime - 0.5 * tPulseDuration ) ).epsilon( 1e-4 ).margin( 1e-8 ) );
        CHECK( tFinal.py == Approx( tReleased.py ).epsilon( 1e-7 ) );
        // Exact angular momentum/energy conservation is not assumed for this
        // dissipative DG scheme. Full histories report their drift and decay.
        CHECK( tFinal.kinetic > 0.0 );
        CHECK( tFinal.strain > 0.0 );
        CHECK( tFinal.kinetic + tFinal.strain <= ( tReleased.kinetic + tReleased.strain ) * ( 1.0 + 1e-5 ) );
        if ( !aFullFlight )
        {
            // References recorded from the verified one-rank 12x2 Q2 case,
            // dt=50 us and T=12 ms. Reuse them on two ranks; no reference solve
            // is run by the test. Tolerances allow roundoff in stiff free flight.
            CHECK( tTip[ 0 ] == Approx( -0.002021374514632963 ).epsilon( 1e-6 ).margin( 2e-9 ) );
            CHECK( tTip[ 1 ] == Approx( -0.02169343580232468 ).epsilon( 1e-6 ).margin( 2e-9 ) );
            CHECK( tFinal.kinetic + tFinal.strain == Approx( 0.3072041374185144 ).epsilon( 1e-6 ) );
            CHECK( tFinal.angular == Approx( -0.02399696451822027 ).epsilon( 1e-6 ) );
        }
        CHECK( std::abs( tTip[ 1 ] ) > ( aFullFlight ? 0.01 : 0.001 ) );
    }

    void run( bool aFullFlight )
    {
        gTossElementsX               = aFullFlight ? 30 : 12;
        gTossElementsY               = aFullFlight ? 4 : 2;
        gTossTimeSteps               = aFullFlight ? 800 : 240;
        gTossFinalTime               = aFullFlight ? 0.04 : 0.012;
        std::string         tLibrary = MORIS_TOSS_LIBRARY;
        std::vector< char > tArgument( tLibrary.begin(), tLibrary.end() );
        tArgument.push_back( '\0' );
        char  tName[]      = "";
        char* tArguments[] = { tName, tArgument.data() };
        if ( par_rank() == 0 )
        {
            std::remove( "Toss_Ruler_history.csv" );
        }
        barrier();
        REQUIRE( fn_WRK_Workflow_Main_Interface( 2, tArguments ) == 0 );
        check_results( aFullFlight );
    }
}    // namespace

// Free plane-stress StVK beam: asymmetric 4 ms pulse and 8 ms of free flight.
// Checks outgoing fields, reconstructed momentum, center-of-mass motion and IQI.
TEST_CASE( "Planar toss ruler: pulse and free flight", "[TossRuler]" )
{
    run( false );
}

// Finer spatial mesh, literature-sized time step and 40 ms duration. Explicitly
// selected demonstration case; writes every slab and the complete energy/momentum
// history for studying finite rotations and time-discretization dissipation.
TEST_CASE( "Planar toss ruler: full flight", "[.TossRulerFull]" )
{
    run( true );
}
