/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Query.cpp
 *
 */

#include "cl_Query.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>

// include moris assert functions
#include "fn_assert.hpp"

// Define Cells
#include "cl_Cell.hpp"

// Define uint, real, etc.
#include "typedefs.hpp"

// Define enums used
#include "cl_Tracer_Enums.hpp"

#include "Log_Constants.hpp"
#include "fn_stringify.hpp"

#include "fn_assert.hpp"

namespace moris
{
namespace ios
{

//-----------------------------------------------------------------------------------------------------------//
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------------------//
// run queries just with user input
void Query::run(int  & argc, char * argv[] )
{
    // initialize
    bool tReadFileIsDeclared = false;
    bool tWriteFileIsDeclared = false;
    bool tQueryTypeIsDeclared = false;
    bool tUseTreeQuery = false;
    bool tUseTableQuery = false;
    bool tSuppressText = false;

    std::string tFileNameRead = " ";
    std::string tFileNameWrite = " ";

    enum EntityBase   tEntityBase   =  EntityBase::Unknown;
    enum EntityType   tEntityType   =  EntityType::Unknown;
    enum EntityAction tEntityAction =  EntityAction ::Unknown;

    // go through user arguments and look for flags
    for( int k=0; k<argc; ++k )
    {
        // READ FILE
        if ( std::string( argv[ k ] ) == "--read" || std::string( argv[ k ] ) == "-rf" )
        {
            tReadFileIsDeclared = true;
            tFileNameRead =  std::string( argv[ k+1 ] );
            std::cout << "Reading from: " << tFileNameRead << "\n";
        }

        // WRITE FILE
        if ( std::string( argv[ k ] ) == "--write" || std::string( argv[ k ] ) == "-wf" )
        {
            tWriteFileIsDeclared = true;
            tFileNameWrite =  std::string( argv[ k+1 ] );
            std::cout << "Write to: " << tFileNameWrite << "\n";
        }

        // QUERY TYPE: TABLE
        if ( std::string( argv[ k ] ) == "--table" || std::string( argv[ k ] ) == "-tb" )
        {
            tQueryTypeIsDeclared = true;
            tUseTableQuery = true;
            std::cout << "Selected query type: " << "Table" << "\n";

            tEntityBase   =  get_entity_base_enum_from_str( std::string( argv[ k+1 ] ) );
            tEntityType   =  get_entity_type_enum_from_str( std::string( argv[ k+2 ] ) );
            tEntityAction =  get_entity_action_enum_from_str( std::string( argv[ k+3 ] ) );

            if (   (tEntityBase   == EntityBase::Unknown)
                || (tEntityType   == EntityType::Unknown)
                || (tEntityAction == EntityAction::Unknown))
            {
                std::cout << "Error: Entity unknown.";
                tQueryTypeIsDeclared = false;
            }
            else
            {
            std::cout << "Table query will look for instances of; " << get_enum_str(tEntityBase)
                                                           << " - " << get_enum_str(tEntityType)
                                                           << " - " << get_enum_str(tEntityAction) << "\n";
            }

        }

        // QUERY TYPE: TREE
        if ( std::string( argv[ k ] ) == "--tree" || std::string( argv[ k ] ) == "-tr" )
        {
            tQueryTypeIsDeclared = true;
            tUseTreeQuery = true;
            std::cout << "Selected query type: " << "Tree" << "\n";

            if ( std::string( argv[ k+1 ] ) == "SupText" )
            {
                tSuppressText = true;
                std::cout << "Text output will be removed from tree.\n";
            }
            else
            {
                tSuppressText = false;
                std::cout << "Text output will be written to tree. \n";
            }
        }

    } // end for each input argument

    // check user input
    if (!tReadFileIsDeclared)
    {
        std::cout << "Error: No .log file to read from provided. \n";
        std::cout << "Use flag -rf with file name, e.g.: -rf path/to/file/Log.log \n";
    }
    if (!tWriteFileIsDeclared)
    {
        std::cout << "Error: File name for output file not provided. \n";
        std::cout << "Use flag -wf with output file name, e.g.: -wf path/to/file/Out.log \n";
    }
    if (!tWriteFileIsDeclared)
    {
        std::cout << "Error: No query type or wrong arguments provided. \n";
        std::cout << "Examples: \n";
        std::cout << "-> To run a tree query with suppressed text output, enter:      --tree SupText \n";
        std::cout << "-> To run a table query looking for non-linear solvers, enter:  --table NonLinearSolver Arbitrary Solve \n";
    }

    // initialize query
    if (tReadFileIsDeclared)
        this->initialize(tFileNameRead);

    // run tree query
    if (tUseTreeQuery && tQueryTypeIsDeclared)
    {
        this->tree_query(tFileNameWrite, tSuppressText);
    }

    // run table query
    if (tUseTableQuery && tQueryTypeIsDeclared)
    {
        this->table_query(tFileNameWrite,
                          tEntityBase, tEntityType, tEntityAction);
    }

}

//-----------------------------------------------------------------------------------------------------------//
// function copies table from log file into cell arrays fed into function, returns number of lines in table
void Query::initialize(std::string aFileNameRead)
{
    std::cout << "Initializing query ... \n" << std::flush;

    // save file name
    mFileNameRead = aFileNameRead;

    std::cout << "Reading info from file: " << mFileNameRead << " ... " << std::flush;

    // read log file into memory
    mNumLines = this->extract_info_from_log_file();

    // get start and end lines for each function instance
    this->isolate_functions();

    // Prepare For Query ---------------------------------------------------------- //

    // open text file to read from
//    if (mLogFileRead.peek() == std::ifstream::traits_type::eof())
//            MORIS_ASSERT( false, "<MRS::IOS::cl_Query::initialize>: Read file is empty.");
    mLogFileRead.open(mFileNameRead);
    if( !mLogFileRead )
    {
        MORIS_ASSERT( false,
                "<MRS::IOS::cl_Query>: Read file could not be opened.");
    }

    std::cout << "Done. \n" << std::flush;

}

//-----------------------------------------------------------------------------------------------------------//
// function copies table from log file into cell arrays fed into function, returns number of lines in table
uint Query::extract_info_from_log_file()
{

    // Prepare reading file ---------------------------------------------------------- //

    // open text file to read from
//    if (mLogFileRead.peek() == std::ifstream::traits_type::eof())
//        MORIS_ASSERT( false, "<MRS::IOS::cl_Query::extract_info_from_text_file>: Read file is empty.");
    mLogFileRead.open(mFileNameRead);
    if( !mLogFileRead )
        MORIS_ASSERT( false, "<MRS::IOS::cl_Query.hpp::extract_info_from_text_file>: Read file could not be opened.");

    // Skip Header ------------------------------------------------------------------- //

    this->skip_header();

    // Read Table  ------------------------------------------------------------------- //

    // variables for saving the currently read data
    std::string tCurrentLine;

    uint tIndentLevel = 0;
    uint tFunctionID = 0;

    enum EntityBase tEntityBase = EntityBase::Unknown;
    enum EntityType tEntityType = EntityType::Unknown;
    enum EntityAction tEntityAction = EntityAction::Unknown;

    enum OutputSpecifier tOutputSpecifier = OutputSpecifier::Unknown;
    std::string tOutputValue = "";

    uint tLineCounter = 0;
    std::string tReadString =  "0";

    // get each line and take apart
    while (std::getline(mLogFileRead, tCurrentLine))
    {
        // count number of lines
        tLineCounter++;

        // copy current line to new stream
        std::istringstream tCurrentLineStream(tCurrentLine);

        // extract function identifiers
        std::getline(tCurrentLineStream, tReadString, ';');
        tIndentLevel = std::stoi(tReadString);
        std::getline(tCurrentLineStream, tReadString, ';');
        tFunctionID = std::stoi(tReadString);

        // extract entity identifiers
        std::getline(tCurrentLineStream, tReadString, ';');
        tEntityBase = get_entity_base_enum_from_str(tReadString);
        std::getline(tCurrentLineStream, tReadString, ';');
        tEntityType = get_entity_type_enum_from_str(tReadString);
        std::getline(tCurrentLineStream, tReadString, ';');
        tEntityAction = get_entity_action_enum_from_str(tReadString);

        // extract output information
        std::getline(tCurrentLineStream, tReadString, ';');
        tOutputSpecifier = get_output_spec_enum_from_str(tReadString);
        std::getline(tCurrentLineStream, tReadString, ';');
        tOutputValue = tReadString;

        // save extracted info to cell array
        if (tLineCounter == 0)
        {
            mIndentLevels.resize( 1 , tIndentLevel );
            mFunctionIDs.resize( 1 , tFunctionID );
            mEntityBases.resize( 1 , tEntityBase );
            mEntityTypes.resize( 1 , tEntityType );
            mEntityActions.resize( 1 , tEntityAction );
            mOutputSpecifiers.resize( 1 , tOutputSpecifier );
            mOutputValues.resize( 1 , tOutputValue );
        }
        else
        {
            mIndentLevels.push_back( tIndentLevel );
            mFunctionIDs.push_back( tFunctionID );
            mEntityBases.push_back( tEntityBase );
            mEntityTypes.push_back( tEntityType );
            mEntityActions.push_back( tEntityAction );
            mOutputSpecifiers.push_back( tOutputSpecifier );
            mOutputValues.push_back( tOutputValue );
        }
    }

    // close everything and end ------------------------------------------------------ //
    mLogFileRead.close();

    return tLineCounter;
}

//-----------------------------------------------------------------------------------------------------------//
// returns a two column list with the indices of the first and last line where information for this function
// is logged
void Query::isolate_functions()
{
    // get number of functions
    uint tNumFuncInstances = this->get_max_func_ID() + 1;

    // initialize output variable
    mInstanceStartEnd.resize(tNumFuncInstances);

    // initialize
//    uint tLastFuncID = 0;
//    uint tLastIndentLevel = 0;
    for(uint ii = 0; ii < tNumFuncInstances; ii++)
    {
        mInstanceStartEnd(ii).resize(2);
    }
    mInstanceStartEnd(0)(0) = 0;

    // go through all lines and collect start and end of each function instance
    for(uint ii = 1; ii < mNumLines; ii++)
    {
        // case: sign in: retrieve function id and save line index to first column
        if (mOutputSpecifiers(ii) == OutputSpecifier::SignIn)
            mInstanceStartEnd(mFunctionIDs(ii))(0) = ii;

        // case: sign out: retrieve function id and save line index to second column
        if (mOutputSpecifiers(ii) == OutputSpecifier::ElapsedTime)
            mInstanceStartEnd(mFunctionIDs(ii))(1) = ii;
    }

}

//-----------------------------------------------------------------------------------------------------------//
// Function to copy the header of an open read file to an open write file
void Query::copy_header(std::ofstream * aLogFileWrite)
{

    // check if files are open
    MORIS_ASSERT( mLogFileRead.is_open() && aLogFileWrite->is_open(),
                  "<MRS::IOS::cl_Query::copy_header>: Read or Write file not opened.");

    // initialize
    std::string tCurrentLine = "";

    // skip header and copy to file
    std::getline(mLogFileRead, tCurrentLine);
    while (tCurrentLine != LOGGER_HEADER_END )
    {
        *aLogFileWrite << tCurrentLine << "\n";
        std::getline(mLogFileRead, tCurrentLine);
    }

    // copy header ending
    *aLogFileWrite << tCurrentLine << "\n\n";

    // skip table header
    for (uint ii = 0; ii < 3; ii++)
    {
        std::getline(mLogFileRead, tCurrentLine);
    }
}

//-----------------------------------------------------------------------------------------------------------//
// Function to skip the header of an open read file
void Query::skip_header()
{
    // check if files are open
    MORIS_ASSERT( mLogFileRead.is_open(),
                  "<MRS::IOS::cl_Query::skip_header>: Read file not opened.");

    // initialize
    std::string tCurrentLine = "";

    // skip header and copy to file
    uint tLineCounter = 1;
    std::getline(mLogFileRead, tCurrentLine);
    while (tCurrentLine != LOGGER_HEADER_END)
    {
        std::getline(mLogFileRead, tCurrentLine);
        MORIS_ASSERT( tLineCounter < 1000, "query::skip_header: Header not found, check header markers." );
        tLineCounter++;
    }

    // skip table header
    for (uint ii = 0; ii < 3; ii++)
    {
        std::getline(mLogFileRead, tCurrentLine);
    }
}

//-----------------------------------------------------------------------------------------------------------//
// gets the maximum ID out of all function instance (= number of different function instances - 1)
uint Query::get_max_func_ID()
{
    // initialize
    uint tMaxFuncID = 0;

    // go through all lines and get max of cell array
    for(uint ii = 0; ii < mNumLines; ii++)
    {
        if (mFunctionIDs(ii) > tMaxFuncID)
            tMaxFuncID = mFunctionIDs(ii);
    }

    // return map that lists the first and last lines for logging for each instance by ID
    return tMaxFuncID;
}

//-----------------------------------------------------------------------------------------------------------//
// creates a line with vertical vertical markers according to indentation level
// (note: line break is NOT included)
std::string Query::create_empty_line(uint aIndentationLevel)
{
    // initialize
    std::string tEmptyLine = "|";

    // go through all lines and get max of cell array
    for(uint ii = 1; ii <= aIndentationLevel; ii++)
    {
        tEmptyLine.append("  |");
    }

    // return map that lists the first and last lines for logging for each instance by ID
    return tEmptyLine;
}

//-----------------------------------------------------------------------------------------------------------//
} // namespace std
} // namespace moris

