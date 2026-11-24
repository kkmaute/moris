/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <catch.hpp>
#include <cmath>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include <cstdio>
#include <string>
#include <vector>
#include "cl_Communication_Tools.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "cl_Logger.hpp"

using namespace moris;

extern "C" {
bool gDGNonlinear          = false;
uint gDGTimeSteps          = 6;
real gDGFinalTime          = 0.3;
real gDGBeamHeight         = 0.43;
real gDGDensity            = 1.0;
real gDGAmplitude          = 0.025;
bool gPrintReferenceValues = false;
}

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

namespace
{
    void run( bool aNonlinear )
    {
        gDGNonlinear                = aNonlinear;
        gDGTimeSteps                = 6;
        gDGFinalTime                = 0.3;
        gDGBeamHeight               = 0.43;
        gDGDensity                  = 1.0;
        gDGAmplitude                = 0.025;
        std::string         library = MORIS_DG_FORWARD_LIBRARY;
        std::vector< char > argument( library.begin(), library.end() );
        argument.push_back( '\0' );
        char  name[]      = "";
        char* arguments[] = { name, argument.data() };
        // Each rank removes its own Exodus output before the workflow so a
        // missing write cannot leave a stale field file for the checker.
        std::string output = "Time_DG_Cantilever_Forward.exo";
        if ( par_size() > 1 )
        {
            std::string rank = std::to_string( par_rank() );
            rank.insert( 0, std::to_string( par_size() ).size() - rank.size(), '0' );
            output += "." + std::to_string( par_size() ) + "." + rank;
        }
        std::remove( output.c_str() );
        barrier();
        REQUIRE( fn_WRK_Workflow_Main_Interface( 2, arguments ) == 0 );
    }

    void check_results()
    {
        // Analytical references: u=(0,-a*x^2*t), v=(0,-a*x^2), beam
        // [0,3]x[0,h]. Q2 space and DG(1) time reproduce these polynomials.
        // Each rank checks its local VIS file; node/element counts are local.
        uint tip_count = 0;
        {
            mtk::Exodus_IO_Helper exo( "Time_DG_Cantilever_Forward.exo", gDGTimeSteps - 1 );
            REQUIRE( exo.get_number_of_dimensions() == 2 );
            // Fixed 12x4 Q2 background, h=0.43: references from the verified
            // one- and two-rank runs, shared by both material models. Counts
            // are local to each Exodus file; VIS includes unused void nodes.
            REQUIRE( ( par_size() == 1 || par_size() == 2 ) );
            REQUIRE( exo.get_number_of_nodes() == ( par_size() == 1 ? 200 : 104 ) );
            REQUIRE( exo.get_number_of_elements() == ( par_size() == 1 ? 72 : 36 ) );
            REQUIRE( exo.get_time_value() == Approx( gDGFinalTime ).margin( 1e-12 ) );
            for ( uint i = 0; i < uint( exo.get_number_of_nodes() ); ++i )
            {
                const uint             id = uint( exo.get_node_id( i ) );
                const Matrix< DDRMat > x  = exo.get_nodal_coordinate( id );
                REQUIRE( x( 0 ) >= -1e-10 );
                REQUIRE( x( 0 ) <= 3.0 + 1e-10 );
                REQUIRE( x( 1 ) >= -1e-10 );
                // VIS retains the full cut-cluster vertex list, including
                // unused void vertices outside the active beam.
                if ( x( 1 ) > gDGBeamHeight + 1e-10 )
                {
                    continue;
                }
                // VIS fields 0 and 1 are Node_Id and Node_Index.
                CHECK( exo.get_nodal_field_value( id, 2, gDGTimeSteps - 1 ) == Approx( 0.0 ).margin( 2e-7 ) );
                CHECK( exo.get_nodal_field_value( id, 3, gDGTimeSteps - 1 ) == Approx( -gDGAmplitude * x( 0 ) * x( 0 ) * gDGFinalTime ).epsilon( 2e-5 ).margin( 2e-7 ) );
                if ( std::abs( x( 0 ) - 3.0 ) < 1e-10 && std::abs( x( 1 ) ) < 1e-10 )
                {
                    ++tip_count;
                }
            }
            if ( gPrintReferenceValues )
            {
                MORIS_LOG_INFO( "DG cantilever rank %u: local nodes %d, local elements %d, final time %.16e",
                        par_rank(),
                        exo.get_number_of_nodes(),
                        exo.get_number_of_elements(),
                        exo.get_time_value() );
            }
            // Global IQIs are already assembled over all ranks.
            CHECK( exo.get_global_variable( 0, gDGTimeSteps - 1 ) == Approx( -9.0 * gDGAmplitude * gDGFinalTime * gDGBeamHeight ).epsilon( 2e-5 ).margin( 2e-7 ) );
            CHECK( exo.get_global_variable( 1, gDGTimeSteps - 1 ) == Approx( 3.0 * gDGBeamHeight ).epsilon( 1e-8 ) );
        }
        // The Exodus reader is closed before collective checks or another reader.
        REQUIRE( sum_all( tip_count ) >= 1 );
    }

}    // namespace

// Linear manufactured bending on a level-set-cut cantilever; checks analytical
// nodal displacement, final integral displacement, volume, and mesh/time data.
TEST_CASE( "Time DG linear elastic cantilever", "[DGForwardLinear]" )
{
    run( false );
    check_results();
}

// Finite-deformation StVK cantilever with exact det(F)=1 manufactured shear;
// exercises nonlinear consistent tangents and temporal endpoint force assembly.
TEST_CASE( "Time DG nonlinear elastic cantilever", "[DGForwardNonlinear]" )
{
    run( true );
    check_results();
}
