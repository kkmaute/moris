/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_WRK_Workflow_Main_Interface.cpp
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <fn_print.hpp>

#include <cstdio>    // nicer than streams in some respects
// C system files
#include <unistd.h>
// C++ system files
#include <stdio.h>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src
#include "typedefs.hpp"                    // COR/src
// other header files
// #include <catch.hpp>
// #include "fn_equal_to.hpp" //ALG
#include "cl_Tracer.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Stopwatch.hpp"    //CHR/src
#include "op_move.hpp"

#include "cl_Matrix.hpp"
#include "cl_Logger.hpp"    // MRS/IOS/src

#include "cl_WRK_Workflow_Factory.hpp"
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "cl_OPT_Manager.hpp"

#include "cl_Library_Factory.hpp"

#include "fn_stringify_matrix.hpp"

using namespace moris;

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] )
{
    // --------------------------------------------- //
    // check arguments provided

    // first, check if there are any arguments provided ...
    if ( argc < 2 )
    {
        // ... print an error if not
        std::cout << "\n Error: input file required\n"
                  << "\n";
        return -1;
    }

    // initialize file names and flags that could be provided
    std::string tSoFileName       = "";
    std::string tXmlFileName      = "";
    bool        tIsOnlyMeshing    = false;
    bool        tSoFileSpecified  = false;
    bool        tXmlFileSpecified = false;

    // go through user arguments and look for flags
    for ( int k = 1; k < argc; ++k )
    {
        // get the current argument
        std::string tArgString    = std::string( argv[ k ] );
        size_t      tArgStrLength = tArgString.length();

        if ( tArgString == "--help" || tArgString == "-h" )
        {
            tIsOnlyMeshing = true;
            MORIS_LOG( "Valid input arguments are:" );
            MORIS_LOG( "option       : --meshgen or -mg :  generate foreground and background meshes (for EXHUME project)" );
            MORIS_LOG( "option       : --pause   or -p  :  report process IDs and pause run upon user input (for debugging)" );
            MORIS_LOG( "last argument: <filename>.so    :  objective file with MORIS input" );
            MORIS_LOG( "last argument: <filename>.xml   :  xml file with MORIS input" );
            MORIS_LOG( " " );

            gMorisComm.finalize();
            exit( 1 );
        }

        // check if user requests to just generate foreground and background meshes (for EXHUME project)
        if ( tArgString == "--meshgen" || tArgString == "-mg" )
        {
            tIsOnlyMeshing = true;
            MORIS_LOG( "Only mesh generation and output requested." );
            continue;
        }

        // if an argument long enough to be a filename is provided
        if ( tArgStrLength > 4 )
        {
            // check if an .so file is provided
            if ( tArgString.substr( tArgStrLength - 3 ) == ".so" )
            {
                MORIS_ERROR( !tSoFileSpecified, "Multiple .so files specified. Specify no more than one!" );
                tSoFileSpecified = true;
                tSoFileName      = tArgString;
                MORIS_LOG( "Reading dynamically linked input file: %s", tArgString.c_str() );
                continue;
            }

            // check if an .xml file is provided
            if ( tArgString.substr( tArgStrLength - 4 ) == ".xml" || tArgString.substr( tArgStrLength - 4 ) == ".XML" )
            {
                MORIS_ERROR( !tXmlFileSpecified, "Multiple .xml files specified. Specify no more than one!" );
                tXmlFileSpecified = true;
                tXmlFileName      = tArgString;
                MORIS_LOG( "Reading parameters from static input file: %s", tArgString.c_str() );
                continue;
            }
        }

        // check whether the input argument is a flag
        if ( tArgString[ 0 ] == '-' )
        {
            MORIS_LOG( "Flag provided: %s", tArgString.c_str() );
            continue;
        }

        // TODO: exclude values that are provided after another flag and related to it, Should the Logger configuration be moved here?
        // if the argument is neither a known file type nor a flag, throw an error
        // MORIS_ERROR( false, "Argument provided is neither a known file-type nor a flag: %s", tArgString.c_str() );
    }

    // check that the arguments provided by user are complete and make sense
    if ( tIsOnlyMeshing )
    {
        // FIXME: uncomment once complete
        // MORIS_ERROR( tXmlFileSpecified, "An .xml input file must be specified for mesh generation." );
    }
    else
    {
        MORIS_ERROR( tSoFileSpecified || tXmlFileSpecified, "Neither an .so nor an .xml file has been provided. Provide at least one input file." );
    }

    // --------------------------------------------- //
    // create library according to inputs

    // prepare library to be generated
    moris::Library_Factory        tLibraryFactory;
    std::shared_ptr< Library_IO > tLibrary = nullptr;

    {
        // log & trace this set of operations
        Tracer tTracer( "WRK", "Main Interface", "Load Parameters" );

        // create the library based on the kind of workflow requested
        if ( tIsOnlyMeshing )
        {
            tLibrary = tLibraryFactory.create_Library( Library_Type::MESHGEN );
        }
        else    // no meshing workflow
        {
            tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
        }

        // load input parameters specified
        if ( tSoFileSpecified )
        {
            tLibrary->load_parameter_list( tSoFileName, File_Type::SO_FILE );
        }
        if ( tXmlFileSpecified )
        {
            tLibrary->load_parameter_list( tXmlFileName, File_Type::XML_FILE );
        }

        // finish initializing the library and lock it from modification
        tLibrary->finalize();
    }

    // --------------------------------------------- //
    // start workflow
    {
        // load the OPT parameter list
        ModuleParameterList tOPTParameterList = tLibrary->get_parameters_for_module( Parameter_List_Type::OPT );

        // Create performer manager
        wrk::Performer_Manager tPerformerManager( tLibrary );

        // FIXME: get this from parameter list "workflow"
        // get which workflow is to be used
        std::string tWRKFlowStr = tOPTParameterList( 0 )( 0 ).get< std::string >( "workflow" );

        // create and initialize (this includes running HMR) an instance of this workflow
        moris::Cell< std::shared_ptr< moris::opt::Criteria_Interface > > tWorkflows = { wrk::create_workflow( tWRKFlowStr, &tPerformerManager ) };

        if ( tOPTParameterList( 0 )( 0 ).get< bool >( "is_optimization_problem" ) )
        {
            moris::opt::Manager tManager( tOPTParameterList, tWorkflows );
            tManager.perform();
        }
        else
        {
            Matrix< DDRMat > tADVs( 0, 0 );
            Matrix< DDRMat > tDummyBounds;
            Matrix< IdMat >  tDummy1( 1, 1, 0.0 );
            tWorkflows( 0 )->initialize( tADVs, tDummyBounds, tDummyBounds, tDummy1 );
            Matrix< DDRMat > tIQIVal = tWorkflows( 0 )->get_criteria( tADVs );

            // print out matrix of IQI values
            MORIS_LOG_SPEC( "IQI values", ios::stringify_log( tIQIVal ) );
        }
    }

    // return success
    return 0;
}
