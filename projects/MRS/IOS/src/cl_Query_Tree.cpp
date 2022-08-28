/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Query_Tree.cpp
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

namespace moris
{
namespace ios
{

//-----------------------------------------------------------------------------------------------------------//
// TREE QUERY DEFINITION
//-----------------------------------------------------------------------------------------------------------//

// function to create tree log file
void Query::tree_query(std::string aFileNameWrite, bool aSuppressText)
{
    // Prepare reading file ---------------------------------------------------------- //

    std::cout << "Creating process tree ... " << std::flush;

    // create text file to write to
    std::ofstream tLogFileWrite;
    tLogFileWrite.open(aFileNameWrite);
    if( !tLogFileWrite )
    {
        MORIS_ASSERT( false,
                "<MRS::IOS::cl_Query::tree_query>: Write file could not be opened.");
    }

    // Read and Copy Header ---------------------------------------------------------- //
    copy_header(& tLogFileWrite);

    // Create Tree ------------------------------------------------------------------- //

    // initialize
    std::string tCurrentLine;
    uint tEndLine = 0;
    std::string tCurrentEntity = "";
    bool tPrevOutputLineIsEmpty = false;
    bool tPrevFuncInstaceIsShort = false;

    // logic loop for tree goes through every line of log information
    for (uint iLine = 0; iLine < mNumLines; iLine++)
    {
        // case: new function instance
        if (mOutputSpecifiers(iLine) == OutputSpecifier::SignIn)
        {
            // check if entity type is existent
            if ( mEntityTypes(iLine) == EntityType::NoType)
                tCurrentEntity = "";
            else
                tCurrentEntity = get_enum_str( mEntityTypes(iLine) ) + " ";

            // get line where the function instance is signed out
            tEndLine = mInstanceStartEnd( mFunctionIDs(iLine) )( 1 );

            // for the base level don't add vertical line
            if (mIndentLevels(iLine) > 0)
                tLogFileWrite << create_empty_line( mIndentLevels(iLine) - 1 ) << "__";

            // write function instance and execution time
            tLogFileWrite << tCurrentEntity
                          << get_enum_str( mEntityBases(iLine) ) << " - "
                          << get_enum_str( mEntityActions(iLine) ) << ": "
                          << mOutputValues( tEndLine ) << " seconds. \n";
            tPrevOutputLineIsEmpty = false;

            // add empty line if function instance is not terminated in next line
            if ( mInstanceStartEnd( mFunctionIDs(iLine) )( 1 ) > iLine + 1 )
            {
                tLogFileWrite << create_empty_line( mIndentLevels(iLine) ) << "\n";
                tPrevOutputLineIsEmpty = true;
                tPrevFuncInstaceIsShort = false;
            }
            else
                tPrevFuncInstaceIsShort = true;
        }

        // case: sign out line
        else if (mOutputSpecifiers(iLine) == OutputSpecifier::ElapsedTime)
        {
            // do nothing
        }

        // case: iteration or step
        else if ( (mOutputSpecifiers(iLine) == OutputSpecifier::Iteration)
                  || (mOutputSpecifiers(iLine) == OutputSpecifier::Step) )
        {
            // add empty line if previous line is not empty
            if ( !tPrevOutputLineIsEmpty )
                tLogFileWrite << create_empty_line( mIndentLevels(iLine) ) << "\n";

            // write output type and value to file
            tLogFileWrite << create_empty_line( mIndentLevels(iLine) ) << "_"
                          << get_enum_str( mOutputSpecifiers(iLine) ) << ": "
                          << mOutputValues( iLine ) << "\n";
            tPrevOutputLineIsEmpty = false;
        }

        // case: text output specifier and text is suppressed, don't do anything
        else if ( (    (mOutputSpecifiers(iLine) == OutputSpecifier::FreeText)
                    || (mOutputSpecifiers(iLine) == OutputSpecifier::InfoText)
                    || (mOutputSpecifiers(iLine) == OutputSpecifier::DebugText) )
                && aSuppressText )
        {
            // do nothing
        }

        // case: any other output specifier
        else
        {

            // write output type and value to file
            tLogFileWrite << create_empty_line( mIndentLevels(iLine) ) << "_"
                          << get_enum_str( mOutputSpecifiers(iLine) ) << " = "
                          <<  mOutputValues( iLine ) << "\n";

            tPrevOutputLineIsEmpty = false;
        }

        // step indentation level down in empty lines if previous output is of higher level
        if (iLine < mNumLines - 1)
        {
            if ( ( mIndentLevels(iLine) > mIndentLevels(iLine + 1) ) && ( !tPrevFuncInstaceIsShort )  )
            {
                tLogFileWrite << create_empty_line( mIndentLevels(iLine + 1) ) << "\n";
                tPrevOutputLineIsEmpty = true;
            }
            else if ( tPrevFuncInstaceIsShort )
                tPrevFuncInstaceIsShort = false;
        }

    } // end create tree logic loop

    // close everything and end ------------------------------------------------------ //
    tLogFileWrite.close();
    std::cout << "Success, Done. \n";
}

//-----------------------------------------------------------------------------------------------------------//
} // namespace std
} // namespace moris

