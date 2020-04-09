// contains methods for Query class
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

// needed?
#include "IO_Tools.hpp"

// Define enums used
#include "cl_Tracer_Enums.hpp"


namespace moris
{
namespace ios
{


//-----------------------------------------------------------------------------------------------------------//
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------//

// table and tree query in separate files
//-----------------------------------------------------------------------------------------------------------//
// function copies table from log file into cell arrays fed into function, returns number of lines in table
void Query::initialize()
{
    std::cout << "Reading info from file: " << mFileNameRead << " ... ";

    // read log file into memory
    mNumLines = this->extract_info_from_log_file();

    // get start and end lines for each function instance
    this->isolate_functions();

    // Prepare For Query ---------------------------------------------------------- //

    // open text file to read from
    if (mLogFileRead.peek() == std::ifstream::traits_type::eof())
            MORIS_ASSERT( false, "<MRS::IOS::cl_Query::initialize>: Read file is empty.");
    mLogFileRead.open(mFileNameRead);
    if( !mLogFileRead )
    {
        MORIS_ASSERT( false,
                "<MRS::IOS::cl_Query>: Read file could not be opened.");
    }

    std::cout << "Done. \n";

}

//-----------------------------------------------------------------------------------------------------------//
// function copies table from log file into cell arrays fed into function, returns number of lines in table
uint Query::extract_info_from_log_file()
{

    // Prepare reading file ---------------------------------------------------------- //

    // open text file to read from
    if (mLogFileRead.peek() == std::ifstream::traits_type::eof())
        MORIS_ASSERT( false, "<MRS::IOS::cl_Query::extract_info_from_text_file>: Read file is empty.");
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
    real tOutputValue = 0.0;

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
        tOutputValue = std::stod(tReadString);

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
    while (tCurrentLine != "<<<")
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
    std::getline(mLogFileRead, tCurrentLine);
    while (tCurrentLine != "<<<")
    {
        std::getline(mLogFileRead, tCurrentLine);
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
// QUERY DEFINITIONS
//-----------------------------------------------------------------------------------------------------------//

// function to create tree log file
void Query::tree_query(std::string aFileNameWrite)
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


        // case: any other output specifier
        else
        {
            // set outputs to floating point numbers

            // write output type and value to file
            tLogFileWrite << create_empty_line( mIndentLevels(iLine) ) << "_"
                          << get_enum_str( mOutputSpecifiers(iLine) ) << " = "
                          << std::scientific  << mOutputValues( iLine ) << "\n";

            // set outputs to whole numbers
            tLogFileWrite.unsetf(std::ios::fixed | std::ios::scientific);

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

// -------------------------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------------------------- //
// function to create tabled log file for instances of a given type
void Query::table_query(std::string aFileNameWrite,
                        enum EntityBase aEntityBase,
                        enum EntityType aEntityType,
                        enum EntityAction aEntityAction)
{


    // Prepare reading file ---------------------------------------------------------- //

    std::cout << "Creating process tables for: "
              << moris::get_enum_str(aEntityBase) << " - "
              << moris::get_enum_str(aEntityType) << " - "
              << moris::get_enum_str(aEntityAction) << " ... ";

    // create text file to write to
    std::ofstream tLogFileWrite;
    tLogFileWrite.open(aFileNameWrite);
    if( !tLogFileWrite )
        MORIS_ASSERT( false, "<MRS::IOS::fn_Query::table_query>: Write file could not be opened.");

    // Read and Copy Header ---------------------------------------------------------- //

    copy_header(& tLogFileWrite);


    // Starting Query ---------------------------------------------------------------- //

    // initialize
    std::string tCurrentLine;
    std::string tCurrentEntity = "";


    /////////////////////////////////
    // Look for function instances //
    /////////////////////////////////

    // initialize
    Cell<real> tInstanceIDs(0);
    Cell<real> tInstanceIndents(0);
    uint tNumInstances = 0;
    bool tTypesMatch = false;
    bool tActionsMatch = false;

    // get number of function instances
    //uint tMaxFuncID = get_max_func_ID(mNumLines, mFunctionIDs);
    uint tMaxFuncID = mInstanceStartEnd.size() - 1;

    // go through every function instance
    for (uint ID = 0; ID <= tMaxFuncID; ID++)
    {
        uint tFirstLineOfInstance = mInstanceStartEnd(ID)(0);

        // check if entities are the same
        if ( mEntityBases(tFirstLineOfInstance) == aEntityBase )
        {
            // check if entity types match
            if ( aEntityType == EntityType::Arbitrary )
                tTypesMatch = true;
            else
                tTypesMatch = ( aEntityType == mEntityTypes(tFirstLineOfInstance) );

            // check if entity actions match
            if ( aEntityType == EntityType::Arbitrary )
                tActionsMatch = true;
            else
                tActionsMatch = ( aEntityAction == mEntityActions(tFirstLineOfInstance) );

            // if everything matches, save instance
            if ( tTypesMatch && tActionsMatch )
            {
                tInstanceIDs.push_back(ID);
                tInstanceIndents.push_back(mIndentLevels(tFirstLineOfInstance));
                tNumInstances++;
            }
        } // end check
    } // end search


    /////////////////////////
    // Write Output Header //
    /////////////////////////

    // write out general header
    tLogFileWrite << "Looking for instances of function: " << get_enum_str(aEntityBase) << " - "
                                                           << get_enum_str(aEntityType) << " - "
                                                           << get_enum_str(aEntityAction) << "\n";
    tLogFileWrite << "Number of instances found: " << tNumInstances << ". \n";

    tLogFileWrite << "Indentation levels of respective instances: [" << tInstanceIndents(0);
    for (uint iInstance = 1; iInstance < tNumInstances; iInstance++)
    {
        tLogFileWrite << ", " << tInstanceIndents(iInstance);
    }
    tLogFileWrite << "]\n";

    tLogFileWrite << "IDs of respective instances: [" << tInstanceIDs(0);
    for (uint iInstance = 1; iInstance < tNumInstances; iInstance++)
    {
        tLogFileWrite << ", " << tInstanceIDs(iInstance);
    }
    tLogFileWrite << "]\n\n";


    ///////////////////////////
    // Record what is inside //
    ///////////////////////////

    // for each instance
    for (uint iInstance = 0; iInstance < tNumInstances; iInstance++)
    {


        /////////////////////////
        // Split up Iterations //
        /////////////////////////

        // get start end lines of function instance
        uint tCurrentInstanceID = tInstanceIDs(iInstance);
        uint tCurrentInstanceIndent = tInstanceIndents(iInstance);
        uint tCurrentInstanceStartLine = mInstanceStartEnd(tCurrentInstanceID)(0);
        uint tCurrentInstanceEndLine = mInstanceStartEnd(tCurrentInstanceID)(1);

        // initialize cell array storing start and end lines of each iteration
        uint tNumIterations = 0;
        Cell<Cell<uint>> tIterationStartEnd;
        tIterationStartEnd.resize(1);
        tIterationStartEnd(0).resize(2);
        tIterationStartEnd(0)(0) = tCurrentInstanceStartLine + 1;
        tIterationStartEnd(0)(1) = tCurrentInstanceEndLine - 1;

        // look through function instance and find iteration markers
        for (uint iLine = tCurrentInstanceStartLine; iLine < tCurrentInstanceEndLine; iLine++)
        {
            if ( ((mOutputSpecifiers(iLine) == OutputSpecifier::Iteration) || (mOutputSpecifiers(iLine) == OutputSpecifier::Step))
                    && (mIndentLevels(iLine) == tCurrentInstanceIndent) )
            {
                // append new iteration
                if (tNumIterations == 0)
                    tIterationStartEnd(0)(0) = iLine;
                else
                    tIterationStartEnd.append( {iLine, tCurrentInstanceEndLine - 1} );

                // save previous line as end line of previous iteration
                if (tNumIterations > 0)
                    tIterationStartEnd(tNumIterations-1)(1) = iLine - 1;

                tNumIterations++;
            }
        } // for each line in instance



        ///////////////////////////
        // Go through iterations //
        ///////////////////////////

        // initialize lists of outputs
        Cell<enum EntityBase> tListOfEntityBases;
        Cell<enum EntityType> tListOfEntityTypes;
        Cell<enum EntityAction> tListOfEntityActions;
        Cell<enum OutputSpecifier> tListOfOutputSpecs;
        Cell<Cell<real>> tListOfOutputValues(1);
        bool tFirstEntry = true;

        /////////////////////////////////////////////////////////
        // Go through first iteration and create reference set //
        /////////////////////////////////////////////////////////

        uint iLine = tIterationStartEnd(0)(0);

        // go through first iteration and collect initial information
        while (iLine <= tIterationStartEnd(0)(1))
        {
            // log information belongs to current function ID
            if ( mFunctionIDs(iLine) == tCurrentInstanceID )
            {
                // copy everything to list of what is happening inside of iteration
                if (tFirstEntry)
                {
                    tListOfEntityBases.resize( 1, mEntityBases(iLine) );
                    tListOfEntityTypes.resize( 1, mEntityTypes(iLine) );
                    tListOfEntityActions.resize( 1, mEntityActions(iLine) );
                    tListOfOutputSpecs.resize( 1, mOutputSpecifiers(iLine) );
                    tListOfOutputValues(0).resize( 1, mOutputValues(iLine) );
                    tFirstEntry = false;
                }
                else
                {
                    tListOfEntityBases.push_back( mEntityBases(iLine) );
                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
                    tListOfEntityActions.push_back( mEntityActions(iLine) );
                    tListOfOutputSpecs.push_back( mOutputSpecifiers(iLine) );
                    tListOfOutputValues(0).push_back( mOutputValues(iLine) );
                }
            }

            // check functions on indentation level one higher for important info
            if ( mIndentLevels(iLine) == tCurrentInstanceIndent + 1 )
            {
                // get end line
                uint tIndentedFunctionID = mFunctionIDs(iLine);
                uint tIndentedFunctionEnd = mInstanceStartEnd(tIndentedFunctionID)(1);

                // record function instance
                if (tFirstEntry)
                {
                    tListOfEntityBases.resize( 1, mEntityBases(iLine) );
                    tListOfEntityTypes.resize( 1, mEntityTypes(iLine) );
                    tListOfEntityActions.resize( 1, mEntityActions(iLine) );
                    tListOfOutputSpecs.resize( 1, mOutputSpecifiers(tIndentedFunctionEnd) );
                    tListOfOutputValues(0).resize( 1, mOutputValues(tIndentedFunctionEnd) );
                    tFirstEntry = false;
                }
                else
                {
                    tListOfEntityBases.push_back( mEntityBases(iLine) );
                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
                    tListOfEntityActions.push_back( mEntityActions(iLine) );
                    tListOfOutputSpecs.push_back( mOutputSpecifiers(tIndentedFunctionEnd) );
                    tListOfOutputValues(0).push_back( mOutputValues(tIndentedFunctionEnd) );
                }

                // initialize
                uint tNumInnerIterations = 0;
                uint tNumInnerSteps = 0;
                bool tInnerFunctionHasResidual = false;
                real tInitalInnerResidualNorm = - 1.0;
                real tFinalInnerResidualNorm = - 1.0;

                // increment line cursor to skip sign in line
                iLine++;

                // go through inner function
                while (iLine <= tIndentedFunctionEnd)
                {
                    // inner iteration is reported
                    if ( (mFunctionIDs(iLine) == tIndentedFunctionID) &&
                         (mOutputSpecifiers(iLine) == OutputSpecifier::Iteration) )
                    {
                        tNumInnerIterations = mOutputValues(iLine);
                    }

                    // inner step is reported
                    if ( (mFunctionIDs(iLine) == tIndentedFunctionID) &&
                         (mOutputSpecifiers(iLine) == OutputSpecifier::Step) )
                    {
                        tNumInnerSteps = mOutputValues(iLine);
                    }

                    // inner residual is reported
                    if ( (mFunctionIDs(iLine) == tIndentedFunctionID) && (mOutputSpecifiers(iLine) == OutputSpecifier::ResidualNorm) )
                    {
                        if (!tInnerFunctionHasResidual)
                        {
                            tInitalInnerResidualNorm = mOutputValues(iLine);
                            tInnerFunctionHasResidual = true;
                        }
                        else // if (tInnerFunctionHasResidual = true)
                            tFinalInnerResidualNorm = mOutputValues(iLine);
                    }

                    // increment line cursor
                    iLine++;

                } // end while: inner function

                // decrement line cursor to last line of indented function instance
                iLine--;

                // record iterations if present
                if (tNumInnerIterations != 0)
                {
                    // record function instance
                    tListOfEntityBases.push_back( mEntityBases(iLine) );
                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
                    tListOfEntityActions.push_back( mEntityActions(iLine) );
                    tListOfOutputSpecs.push_back( OutputSpecifier::Iteration );
                    tListOfOutputValues(0).push_back( tNumInnerIterations );
                }

                // record steps if present
                if (tNumInnerSteps != 0)
                {
                    // record function instance
                    tListOfEntityBases.push_back( mEntityBases(iLine) );
                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
                    tListOfEntityActions.push_back( mEntityActions(iLine) );
                    tListOfOutputSpecs.push_back( OutputSpecifier::Step );
                    tListOfOutputValues(0).push_back( tNumInnerIterations );
                }

                // record Residual Drop if present
                if (tInnerFunctionHasResidual)
                {
                    // record function instance
                    tListOfEntityBases.push_back( mEntityBases(iLine) );
                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
                    tListOfEntityActions.push_back( mEntityActions(iLine) );
                    tListOfOutputSpecs.push_back( OutputSpecifier::ResidualDrop );
                    tListOfOutputValues(0).push_back( tFinalInnerResidualNorm / tInitalInnerResidualNorm );
                }

            } // end if: indented function

            // increment line cursor
            iLine++;

        } // end while: go through iteration


        //////////////////
        // Write Tables //
        //////////////////

        uint tMaxColumnWidth = 18;
        uint tNumCols = tListOfEntityBases.size();


        // write general info
        tLogFileWrite << "Instance ID: " << tInstanceIDs(iInstance) << "\n";
        tLogFileWrite << "Execution Time: " << mOutputValues(mInstanceStartEnd(tInstanceIDs(iInstance))(1)) << " seconds.\n";

        // write header for a table
        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "|";
            for (uint iChar = 0; iChar < tMaxColumnWidth+2; iChar++ )
            tLogFileWrite << "-";
        }
        tLogFileWrite << "|\n";

        // write header for a table
        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfEntityBases(iCol)) << " ";
        }
        tLogFileWrite << "|\n";

        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfEntityTypes(iCol)) << " ";
        }
        tLogFileWrite << "|\n";

        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfEntityActions(iCol)) << " ";
        }
        tLogFileWrite << "|\n";

        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "|";
            for (uint iChar = 0; iChar < tMaxColumnWidth+2; iChar++ )
            tLogFileWrite << "-";
        }
        tLogFileWrite << "|\n";

        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfOutputSpecs(iCol)) << " ";
        }
        tLogFileWrite << "|\n";

        for (uint iCol = 0; iCol < tNumCols; iCol++ )
        {
            tLogFileWrite << "|";
            for (uint iChar = 0; iChar < tMaxColumnWidth+2; iChar++ )
            tLogFileWrite << "-";
        }
        tLogFileWrite << "|\n";

        // write output
        uint tNumRows = tListOfOutputValues.size();
        for (uint iRow = 0; iRow < tNumRows; iRow++ )
        {
            for (uint iCol = 0; iCol < tNumCols; iCol++ )
            {
                // if integer
                if (   (tListOfOutputSpecs(iCol) == OutputSpecifier::Iteration)
                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Step)
                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Count)
                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Error)
                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Restart) )
                {
                    tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << tListOfOutputValues(iRow)(iCol) << " ";
                }
                // else: float
                else
                {
                    tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << std::scientific << tListOfOutputValues(iRow)(iCol) << " ";
                    tLogFileWrite.unsetf(std::ios::fixed | std::ios::scientific);
                }
            }
            tLogFileWrite << "|\n";
        } // end: for each row

        // empty line between tables
        tLogFileWrite << "\n";

    } // end: for each instance




    // close everything and end ------------------------------------------------------ //
    tLogFileWrite.close();
    std::cout << "Success, Done. \n";
}


//-----------------------------------------------------------------------------------------------------------//
} // namespace std
} // namespace moris













