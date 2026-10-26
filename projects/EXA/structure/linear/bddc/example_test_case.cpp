/* Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. */

#include <catch.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <petscconf.h>

#include "cl_Logger.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "fn_norm.hpp"

using namespace moris;

// Globals shared with the two input libraries, as in the neighboring EXA tests.
uint        gElementsPerSubdomain = 4;
bool        gUseMaterialJump      = false;
real        gYoungsModulusRatio   = 1000.0;
real        gLinearTolerance      = 1.e-9;
std::string gLinearSolver         = "bddc";
std::string gOutputFileName;

// Match PETSc's available factor package without configuring PETSc in the test.
#ifdef PETSC_HAVE_MUMPS
std::string gFactorSolverType = "mumps";
#else
std::string gFactorSolverType;
#endif

// Test selection and optional benchmark controls belong to the driver.
uint        gCaseIndex            = 0;
bool        gUseEdges             = true;
bool        gPrintReferenceValues = false;
bool        gRunBenchmark         = false;
std::string gPetscOptions;
std::string gCsvFileName;

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

// ----------------------------------------------------------------------------

void check_results()
{
    // Check the upper-right corner on the last rank, which contains that node
    // for HMR's Cartesian partition. No communication occurs inside this check.
    mtk::Exodus_IO_Helper tExoIO( gOutputFileName, 0, false, false );
    MORIS_LOG_INFO( "Checking BDDC results - test case %u on %i processors.", gCaseIndex, par_size() );

    // References were recorded from the Amesos solve with unit traction, Q1
    // elements (four per subdomain direction), and the same Nitsche clamp.
    // Cases: 2D serial, 2D parallel, 3D serial, 3D parallel; each has a
    // homogeneous case followed by the central inclusion with modulus 1000.
    const Vector< uint > tReferenceNumDims      = { 2, 2, 2, 2, 3, 3, 3, 3 };
    const Vector< uint > tReferenceNumNodes     = { 25, 25, 25, 25, 125, 125, 125, 125 };
    const Vector< uint > tReferenceNumElems     = { 16, 16, 16, 16, 64, 64, 64, 64 };
    const Vector< uint > tReferenceNodeId       = { 24, 24, 99, 99, 124, 124, 999, 999 };
    const Vector< real > tReferenceStrainEnergy = {
        1.6411834022095948e-08,
        2.9045171274205883e-01,
        1.6442332473594745e-08,
        3.0185790248479144e-01,
        1.6010246623920688e-08,
        3.1038450039872834e-01,
        1.6108614646464146e-08,
        3.3302715542548572e-01
    };

    const Vector< Matrix< DDRMat > > tReferenceDisplacement = {
        { { 3.2976653320264319e-08 }, { -5.2098618587134451e-09 } },
        { { 7.3418337223728269e-01 }, { -2.4778986955274088e-01 } },
        { { 3.3047002332297403e-08 }, { -5.1857799916392722e-09 } },
        { { 7.6844207393273234e-01 }, { -2.4911906855026938e-01 } },
        { { 3.2255949112611768e-08 }, { -5.1458372848119974e-09 }, { -5.1458372848120595e-09 } },
        { { 8.1227155646339577e-01 }, { -1.8652229389547650e-01 }, { -1.8652229389560621e-01 } },
        { { 3.2483577104741192e-08 }, { -5.1441576286936053e-09 }, { -5.1441576286935895e-09 } },
        { { 8.6564094037446915e-01 }, { -1.9898339297361192e-01 }, { -1.9898339297390155e-01 } }
    };

    uint             tNodeId     = tReferenceNodeId( gCaseIndex );
    Matrix< DDRMat > tCoordinate = tExoIO.get_nodal_coordinate( tNodeId );
    Matrix< DDRMat > tDisplacement( tReferenceNumDims( gCaseIndex ), 1 );

    // VIS reserves nodal fields 0 and 1 for Node_Id and Node_Index.
    for ( uint d = 0; d < tDisplacement.numel(); ++d )
    {
        tDisplacement( d ) = tExoIO.get_nodal_field_value( tNodeId, d + 2, 0 );
    }

    if ( gPrintReferenceValues )
    {
        std::cout << "Case index: " << gCaseIndex << '\n'
                  << "Number of dimensions: " << tExoIO.get_number_of_dimensions() << '\n'
                  << "Number of nodes: " << tExoIO.get_number_of_nodes() << '\n'
                  << "Number of elements: " << tExoIO.get_number_of_elements() << '\n'
                  << "Reference node ID: " << tNodeId << '\n';
        moris::print( tCoordinate, "Coordinates of reference point" );
        moris::print( tDisplacement, "Displacement at reference point" );
        std::cout << "Strain energy: " << std::scientific << std::setprecision( 16 )
                  << tExoIO.get_global_variable( 0, 0 ) << '\n';
        return;
    }

    REQUIRE( tExoIO.get_number_of_dimensions() == int( tReferenceNumDims( gCaseIndex ) ) );
    REQUIRE( tExoIO.get_number_of_nodes() == int( tReferenceNumNodes( gCaseIndex ) ) );
    REQUIRE( tExoIO.get_number_of_elements() == int( tReferenceNumElems( gCaseIndex ) ) );

    Matrix< DDRMat > tReferenceCoordinate( tReferenceNumDims( gCaseIndex ), 1, 1.0 );
    real             tCoordinateError = norm( tCoordinate - tReferenceCoordinate ) / norm( tReferenceCoordinate );
    REQUIRE( tCoordinateError < 1.e-12 );
    REQUIRE( tExoIO.get_time_value() == Approx( 1.0 ).epsilon( 1.e-12 ) );

    real tDisplacementError = norm( tDisplacement - tReferenceDisplacement( gCaseIndex ) )
                            / norm( tReferenceDisplacement( gCaseIndex ) );
    MORIS_LOG_INFO( "Check displacement: relative error %12.5e.", tDisplacementError );
    REQUIRE( tDisplacementError < 1.e-6 );

    real tEnergy      = tExoIO.get_global_variable( 0, 0 );
    real tEnergyError = std::abs( tEnergy / tReferenceStrainEnergy( gCaseIndex ) - 1.0 );
    MORIS_LOG_INFO( "Check strain energy: reference %12.5e, actual %12.5e, relative error %12.5e.",
            tReferenceStrainEnergy( gCaseIndex ),
            tEnergy,
            tEnergyError );
    REQUIRE( tEnergyError < 1.e-6 );
    REQUIRE( tExoIO.get_global_variable( 1, 0 ) == Approx( 1.0 ).epsilon( 1.e-10 ) );
    REQUIRE( tExoIO.get_global_variable( 2, 0 ) == Approx( 1.0 ).epsilon( 1.e-10 ) );
}

// ----------------------------------------------------------------------------

namespace
{
    void run_case(
            uint               aDimension,
            const std::string& aLibrary,
            bool               aJump     = false,
            bool               aUseEdges = true,
            const std::string& aSolver   = "bddc" )
    {
        // Each regression case supplies explicit globals before loading its input
        // library. Benchmark runs instead use globals populated by command-line options.
        if ( !gRunBenchmark )
        {
            REQUIRE( ( par_size() == 1 || par_size() == ( aDimension == 2 ? 4 : 8 ) ) );
            gElementsPerSubdomain = 4;
            gUseMaterialJump      = aJump;
            gYoungsModulusRatio   = 1000.0;
            gUseEdges             = aUseEdges;
            gLinearSolver         = aSolver;
            gLinearTolerance      = 1.e-9;
        }

        gCaseIndex      = ( aDimension == 2 ? 0 : 4 ) + ( par_size() == 1 ? 0 : 2 ) + uint( gUseMaterialJump );
        gOutputFileName = aLibrary + "_Case_" + std::to_string( gCaseIndex ) + "_" + gLinearSolver + ".exo";

        // Remove stale output on every rank before invoking the normal workflow.
        std::string tFile = gOutputFileName;
        if ( par_size() > 1 )
        {
            std::ostringstream tRank;
            tRank << std::setfill( '0' ) << std::setw( std::to_string( par_size() ).size() ) << par_rank();
            tFile += "." + std::to_string( par_size() ) + "." + tRank.str();
        }
        std::remove( tFile.c_str() );
        barrier();

        std::string tLibrary     = std::string( BDDC_INPUT_DIRECTORY ) + "/" + aLibrary + ".so";
        char        tProgram[]   = "";
        char*       tArguments[] = { tProgram, tLibrary.data() };

        MORIS_LOG_INFO( "Executing %s - test case %u on %i processors.", aLibrary.c_str(), gCaseIndex, par_size() );
        double tStart = MPI_Wtime();
        REQUIRE( fn_WRK_Workflow_Main_Interface( 2, tArguments ) == 0 );
        double tElapsed = MPI_Wtime() - tStart;
        double tWorkflowTime;
        MPI_Allreduce( &tElapsed, &tWorkflowTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

        // These integral/topology checks also apply to arbitrary benchmark meshes.
        {
            mtk::Exodus_IO_Helper tOutput( gOutputFileName, 0, false, false );
            uint                  tCells         = sum_all( uint( tOutput.get_number_of_elements() ) );
            uint                  tExpectedCells = par_size();
            for ( uint d = 0; d < aDimension; ++d )
            {
                tExpectedCells *= gElementsPerSubdomain;
                const auto& tField = tOutput.get_nodal_field_vector( d + 2, 0 );
                for ( uint i = 0; i < tField.numel(); ++i )
                {
                    REQUIRE( std::isfinite( tField( i ) ) );
                }
            }
            REQUIRE( tCells == tExpectedCells );
            REQUIRE( tOutput.get_number_of_dimensions() == int( aDimension ) );
            REQUIRE( tOutput.get_time_value() == Approx( 1.0 ) );

            real tEnergy = tOutput.get_global_variable( 0, 0 );
            real tVolume = tOutput.get_global_variable( 1, 0 );
            real tLoad   = tOutput.get_global_variable( 2, 0 );
            REQUIRE( std::isfinite( tEnergy ) );
            REQUIRE( tEnergy > 0.0 );
            REQUIRE( tVolume == Approx( 1.0 ).epsilon( 1.e-10 ) );
            REQUIRE( tLoad == Approx( 1.0 ).epsilon( 1.e-10 ) );

            if ( par_rank() == 0 && gRunBenchmark )
            {
                const std::string  tHeader = "dim,ranks,elements_per_subdomain_axis,cells,jump,sigma,solver,constraints,strain_energy,volume,load,workflow_seconds";
                std::ostringstream tRow;
                tRow << std::setprecision( 16 ) << aDimension << ',' << par_size() << ',' << gElementsPerSubdomain << ',' << tCells
                     << ',' << gUseMaterialJump << ',' << gYoungsModulusRatio << ',' << gLinearSolver << ',' << ( gUseEdges ? "ce" : "c" )
                     << ',' << tEnergy << ',' << tVolume << ',' << tLoad << ',' << tWorkflowTime;
                std::cout << "BDDC_RESULT " << tHeader << '\n'
                          << "BDDC_RESULT " << tRow.str() << '\n';

                if ( !gCsvFileName.empty() )
                {
                    std::ofstream tFile( gCsvFileName );
                    REQUIRE( tFile.good() );
                    tFile << tHeader << '\n'
                          << tRow.str() << '\n';
                    tFile.close();
                    REQUIRE( tFile.good() );
                }
            }
        }

        // Regression checks follow the neighboring EXA examples: a selected rank
        // reads Exodus and compares coordinates, displacement, and IQIs to constants.
        if ( !gRunBenchmark && par_rank() == par_size() - 1 )
        {
            gLogger.set_screen_output_rank( par_rank() );
            check_results();
            gLogger.set_screen_output_rank( 0 );
        }
    }
}    // namespace

// In regression mode, all cases use Q1 elements, a symmetric Nitsche clamp,
// and constant unit horizontal traction. They check mesh size, finite fields,
// time, unit volume, and unit total load. On the selected rank, they also check
// the reference-node coordinates and displacement and the global strain energy
// against stored Amesos results. Solver variants share the same physical references;
// iteration counts and timings are diagnostics, not regression assertions.
//
// PETSc cannot be reinitialized after finalization. CTest selects one hidden
// case per process, including distinct cases for C/CE and AIJ comparisons.
TEST_CASE( "BDDC_Plane_Stress", "[.][BDDC]" )
{
    // Homogeneous 2D plane stress with MATIS/BDDC and corner-plus-edge
    // constraints. Exercises the standard HMR/FEM workflow on one or four ranks.
    run_case( 2, "BDDC_Plane_Stress" );
}

TEST_CASE( "BDDC_Plane_Stress_Jump", "[.][BDDC]" )
{
    // Central square inclusion: E = 1000 inside and E = 1 outside.
    // Checks that spatially varying stiffness and BDDC reproduce the stored
    // inclusion solution on one or four ranks, retaining edge constraints.
    run_case( 2, "BDDC_Plane_Stress", true );
}

TEST_CASE( "BDDC_Plane_Stress_Corners", "[.][BDDC]" )
{
    // Homogeneous plane stress with edge constraints disabled. The four-rank
    // regression checks that corner-only BDDC reaches the same physical solution
    // as the corner-plus-edge and direct-solver variants.
    run_case( 2, "BDDC_Plane_Stress", false, false );
}

TEST_CASE( "BDDC_Plane_Stress_AIJ", "[.][BDDC]" )
{
    // Solve the homogeneous plane-stress model with PETSc AIJ and direct LU.
    // The four-rank regression protects the existing AIJ path while checking
    // the same displacement and energy references used by MATIS/BDDC.
    run_case( 2, "BDDC_Plane_Stress", false, true, "aij" );
}

TEST_CASE( "BDDC_Elasticity_3D", "[.][BDDC]" )
{
    // Homogeneous 3D elasticity with MATIS/BDDC and corner, edge, and face
    // constraints. Exercises HEX8 assembly and all displacement components on
    // one or eight ranks, including floating subdomains in the parallel case.
    run_case( 3, "BDDC_Elasticity_3D" );
}

TEST_CASE( "BDDC_Elasticity_3D_Jump", "[.][BDDC]" )
{
    // Central cubic inclusion: E = 1000 inside and E = 1 outside.
    // Checks the heterogeneous HEX8 model with corner, edge, and face constraints
    // against stored inclusion references on one or eight ranks.
    run_case( 3, "BDDC_Elasticity_3D", true );
}

TEST_CASE( "BDDC_Elasticity_3D_Corners", "[.][BDDC]" )
{
    // Homogeneous 3D elasticity with edge and face constraints disabled.
    // The eight-rank regression exercises corner-only BDDC on floating elasticity
    // subdomains and checks the same solution as the richer constraint variant.
    run_case( 3, "BDDC_Elasticity_3D", false, false );
}

TEST_CASE( "BDDC_Elasticity_3D_AIJ", "[.][BDDC]" )
{
    // Solve the homogeneous 3D model with PETSc AIJ and direct LU on eight ranks.
    // Checks that the existing assembled-matrix path reproduces the stored
    // HEX8 displacement and energy references independently of MATIS/BDDC.
    run_case( 3, "BDDC_Elasticity_3D", false, true, "aij" );
}
