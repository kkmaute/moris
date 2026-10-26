/* Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include "cl_Communication_Manager.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

extern moris::uint gElementsPerSubdomain;
extern bool        gUseMaterialJump;
extern moris::real gYoungsModulusRatio;
extern moris::real gLinearTolerance;
extern std::string gLinearSolver;
extern bool        gUseEdges;
extern bool        gPrintReferenceValues;
extern bool        gRunBenchmark;
extern std::string gPetscOptions;
extern std::string gCsvFileName;

// ----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    // Use the same MORIS/Catch initialization as the shared EXA main. Extra
    // command-line options populate globals for performance studies only.
    gMorisComm = moris::Comm_Manager( &argc, &argv );
    gLogger.initialize( 2 );
    moris::print_banner( argc, argv );

    int tResult;
    {
        Catch::Session tSession;
        auto           tCommandLine = tSession.cli()
                          | Catch::clara::Opt( gRunBenchmark )[ "--benchmark" ]( "Use benchmark globals instead of the fixed regression configuration" )
                          | Catch::clara::Opt( gElementsPerSubdomain, "count" )[ "--elements" ]( "Elements per subdomain direction" )
                          | Catch::clara::Opt( gUseMaterialJump )[ "--jump" ]( "Enable the central material inclusion" )
                          | Catch::clara::Opt( gYoungsModulusRatio, "modulus" )[ "--sigma" ]( "Young's modulus of the inclusion" )
                          | Catch::clara::Opt( gLinearTolerance, "tolerance" )[ "--rtol" ]( "Linear relative tolerance" )
                          | Catch::clara::Opt( gLinearSolver, "solver" )[ "--solver" ]( "bddc, aij, or amesos" )
                          | Catch::clara::Opt( gUseEdges, "bool" )[ "--use-edges" ]( "Enable BDDC edge and face constraints" )
                          | Catch::clara::Opt( gPrintReferenceValues )[ "--print-reference-values" ]( "Print regression reference values" )
                          | Catch::clara::Opt( gPetscOptions, "options" )[ "--petsc-options" ]( "Additional PETSc options" )
                          | Catch::clara::Opt( gCsvFileName, "filename" )[ "--csv" ]( "Benchmark measurement CSV" );
        tSession.cli( tCommandLine );
        tResult = tSession.run( argc, argv );
    }

    gMorisComm.finalize();
    return tResult;
}
