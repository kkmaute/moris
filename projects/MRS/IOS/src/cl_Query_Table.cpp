/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Query_Table.cpp
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

#include "fn_Logger.hpp"
#include "fn_stringify.hpp"
#include "fn_assert.hpp"

namespace moris
{
namespace ios
{

//-----------------------------------------------------------------------------------------------------------//
// TABLE QUERY DEFINITION
//-----------------------------------------------------------------------------------------------------------//

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
              << moris::get_enum_str(aEntityAction) << " ... " << std::flush;

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
    Cell<uint> tInstanceIDs(0);
    Cell<uint> tInstanceIndents(0);

    uint tNumInstances = this->find_instances( & tInstanceIDs,
                                               & tInstanceIndents,
                                                 aEntityBase,
                                                 aEntityType,
                                                 aEntityAction);

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

        // get current ID
        uint tCurrentInstanceID = tInstanceIDs(iInstance);

        // initialize cell array storing start and end lines of each iteration
        Cell<Cell<uint>> tIterationStartEnd = this->split_instances_into_iterations(tCurrentInstanceID);
        uint tNumIterations = tIterationStartEnd.size();

        ////////////////////////////////
        // Go through first iteration //
        ////////////////////////////////

        // initialize reference lists of outputs
        Cell<enum EntityBase> tRefListOfEntityBases;
        Cell<enum EntityType> tRefListOfEntityTypes;
        Cell<enum EntityAction> tRefListOfEntityActions;
        Cell<enum OutputSpecifier> tRefListOfOutputSpecs;
        Cell<std::string> tRefListOfOutputValues;

        // initialize full list of Output Values
        Cell<Cell<std::string>> tFullListOfOutputValues(1);

        // go through first iteration and create reference set of outputs
        this->extract_iteration( & tRefListOfEntityBases,
                                 & tRefListOfEntityTypes,
                                 & tRefListOfEntityActions,
                                 & tRefListOfOutputSpecs,
                                 & tRefListOfOutputValues,
                                   tIterationStartEnd,
                                   1,
                                   tCurrentInstanceID);

        // copy list of output values to table
        tFullListOfOutputValues(0) = tRefListOfOutputValues;

        /////////////////////////////////////
        // Go through subsequent iteration //
        /////////////////////////////////////

        // go through each subsequent iteration and modify reference set if necessary

        // initialize  lists of outputs
        Cell<enum EntityBase> tNewListOfEntityBases;
        Cell<enum EntityType> tNewListOfEntityTypes;
        Cell<enum EntityAction> tNewListOfEntityActions;
        Cell<enum OutputSpecifier> tNewListOfOutputSpecs;
        Cell<std::string> tNewListOfOutputValues;

        for (uint iIter = 2; iIter <= tNumIterations; iIter++)
        {

            // go through first iteration and get set of outputs
            this->extract_iteration( & tNewListOfEntityBases,
                                     & tNewListOfEntityTypes,
                                     & tNewListOfEntityActions,
                                     & tNewListOfOutputSpecs,
                                     & tNewListOfOutputValues,
                                       tIterationStartEnd,
                                       iIter,
                                       tCurrentInstanceID);

            // merge new list with old lists
            this->merge_iteration_outputs( & tRefListOfEntityBases,
                                           & tRefListOfEntityTypes,
                                           & tRefListOfEntityActions,
                                           & tRefListOfOutputSpecs,
                                           & tFullListOfOutputValues,
                                             tNewListOfEntityBases,
                                             tNewListOfEntityTypes,
                                             tNewListOfEntityActions,
                                             tNewListOfOutputSpecs,
                                             tNewListOfOutputValues);

        } // end for: each iteration

        //////////////////
        // Write Tables //
        //////////////////

        // write general info
        tLogFileWrite << "Instance ID: " << tInstanceIDs(iInstance) << "\n";
        tLogFileWrite << "Execution Time: " << mOutputValues(mInstanceStartEnd(tInstanceIDs(iInstance))(1)) << " seconds.\n";

        // write header for a table
        uint tNumCols = tRefListOfEntityBases.size();
        this->write_instance_table_header(& tLogFileWrite,
                                            tRefListOfEntityBases,
                                            tRefListOfEntityTypes,
                                            tRefListOfEntityActions,
                                            tRefListOfOutputSpecs);

        // write output
        uint tNumRows = tFullListOfOutputValues.size();
        for (uint iRow = 0; iRow < tNumRows; iRow++ )
        {
            for (uint iCol = 0; iCol < tNumCols; iCol++ )
            {
                tLogFileWrite << "| " << std::setw(QUERY_MAX_COLUMN_WIDTH) << tFullListOfOutputValues(iRow)(iCol) << " ";
            } // end: for each column

            tLogFileWrite << "|\n";
        } // end: for each row

        // empty line between tables
        tLogFileWrite << "\n";

    } // end: for each instance

    // close everything and end ------------------------------------------------------ //
    tLogFileWrite.close();
    std::cout << "Success, Done. \n" << std::flush;
}

//-----------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------//
// write header for a table of function instance
void Query::write_instance_table_header(std::ofstream * aLogFileWrite,
                                        Cell<enum EntityBase> aRefListOfEntityBases,
                                        Cell<enum EntityType> aRefListOfEntityTypes,
                                        Cell<enum EntityAction> aRefListOfEntityActions,
                                        Cell<enum OutputSpecifier> aRefListOfOutputSpecs)
{
    // initialize
    uint tNumCols = aRefListOfEntityBases.size();

    // initialize temporary modifiable strings
    std::stringstream tSeparationLine("");
    std::stringstream tEntityBases("");
    std::stringstream tEntityTypes("");
    std::stringstream tEntityActions("");
    std::stringstream tOutputSpecs("");

    for (uint iCol = 0; iCol < tNumCols; iCol++ )
    {
        tSeparationLine << "|";
        // create separation line
        for (uint iChar = 0; iChar < QUERY_MAX_COLUMN_WIDTH+2; iChar++ )
        {
            tSeparationLine << "-";
        }

        // write out each of the header lines to temporary stream
        tEntityBases   << "| " << std::setw(QUERY_MAX_COLUMN_WIDTH) << get_enum_str(aRefListOfEntityBases(iCol))   << " ";
        tEntityTypes   << "| " << std::setw(QUERY_MAX_COLUMN_WIDTH) << get_enum_str(aRefListOfEntityTypes(iCol))   << " ";
        tEntityActions << "| " << std::setw(QUERY_MAX_COLUMN_WIDTH) << get_enum_str(aRefListOfEntityActions(iCol)) << " ";
        tOutputSpecs   << "| " << std::setw(QUERY_MAX_COLUMN_WIDTH) << get_enum_str(aRefListOfOutputSpecs(iCol))   << " ";
    }
    tSeparationLine << "|";
    tEntityBases    << "|";
    tEntityTypes    << "|";
    tEntityActions  << "|";
    tOutputSpecs    << "|";

    // dump everything to file
    (*aLogFileWrite) << tSeparationLine.str() << "\n"
                     << tEntityBases.str() << "\n"
                     << tEntityTypes.str() << "\n"
                     << tEntityActions.str() << "\n"
                     << tSeparationLine.str() << "\n"
                     << tOutputSpecs.str() << "\n"
                     << tSeparationLine.str() << "\n";

}

//-----------------------------------------------------------------------------------------------------------//

// merges an old iteration table with a new one
void Query::merge_iteration_outputs( Cell<enum EntityBase>      * aRefListOfEntityBases,
                                     Cell<enum EntityType>      * aRefListOfEntityTypes,
                                     Cell<enum EntityAction>    * aRefListOfEntityActions,
                                     Cell<enum OutputSpecifier> * aRefListOfOutputSpecs,
                                     Cell<Cell<std::string>>    * aFullListOfOutputValues,

                                     Cell<enum EntityBase>      aNewListOfEntityBases,
                                     Cell<enum EntityType>      aNewListOfEntityTypes,
                                     Cell<enum EntityAction>    aNewListOfEntityActions,
                                     Cell<enum OutputSpecifier> aNewListOfOutputSpecs,
                                     Cell<std::string>          aNewListOfOutputValues)
{
    // initialize sorting lists
    Cell<uint> tRefSortingList(aRefListOfEntityBases->size());
    Cell<uint> tNewSortingList(aNewListOfEntityBases.size());

    // initialize counters
    uint tResultListCursor = 0;
    uint tNewListCursor = 0;
    uint tNewListIndex = 0;

    // --- Generate Sorting Lists ---------------------------------------------------- //
    // go through reference list of entities
    for (uint iRefListCursor = 0;  iRefListCursor < aRefListOfEntityBases->size();  iRefListCursor++)
    {
        // initialize
        bool tMatchFound = false;
        tNewListIndex = tNewListCursor;

        // for given reference list entry, go through new list and see if and where a corresponding entry in the new list can be found
        while ( ( tNewListIndex < aNewListOfEntityBases.size() ) && !tMatchFound )
        {
            // check for match
            if (    ((*aRefListOfEntityBases)(iRefListCursor)   == aNewListOfEntityBases(tNewListIndex))
                 && ((*aRefListOfEntityTypes)(iRefListCursor)   == aNewListOfEntityTypes(tNewListIndex))
                 && ((*aRefListOfEntityActions)(iRefListCursor) == aNewListOfEntityActions(tNewListIndex))
                 && ((*aRefListOfOutputSpecs)(iRefListCursor)   == aNewListOfOutputSpecs(tNewListIndex))   )
            {
                // indicate that match has been found
                tMatchFound = true;

                // save to both sorting lists
                tRefSortingList(iRefListCursor) = tResultListCursor;
                tNewSortingList(tNewListIndex) = tResultListCursor;
                tResultListCursor++;

                // check if there are new entries in new list that have been skipped
                if (tNewListIndex > tNewListCursor + 1)
                {
                    // save position
                    tNewSortingList(tNewListIndex) = tResultListCursor;
                    tResultListCursor++;
                }

                // indicate that the new list does not need to be searched before this instance
                tNewListCursor = tNewListIndex + 1;

            } // end if: match is found

            // increment
            tNewListIndex++;

        } // end while: go through new list

        // if match is not found in new list, create new entry
        if (!tMatchFound)
        {
            tRefSortingList(iRefListCursor) = tResultListCursor;
            tResultListCursor++;
        }

    } // end for: each entry in reference list

    // if end of reference list is reached, but end of new list is not
    if ( tNewListCursor < aNewListOfEntityBases.size() - 1)
    {
        // go through the rest and create new list entries
        while(tNewListCursor < aNewListOfEntityBases.size())
        {
            // save position
            tNewSortingList(tNewListCursor) = tResultListCursor;
            tResultListCursor++;

            // increment
            tNewListCursor++;
        }
    }

    // --- Create New Entity List ---------------------------------------------------- //

    // get number of entries in list
    uint tResultListSize = tResultListCursor;

    // initialize result entity lists
    Cell<enum EntityBase>      tResultListOfEntityBases(tResultListSize);
    Cell<enum EntityType>      tResultListOfEntityTypes(tResultListSize);
    Cell<enum EntityAction>    tResultListOfEntityActions(tResultListSize);
    Cell<enum OutputSpecifier> tResultListOfOutputSpecs(tResultListSize);

    //go reference list and copy to resulting list
    for (uint iRefListIndex = 0; iRefListIndex < tRefSortingList.size(); iRefListIndex++)
    {
        tResultListOfEntityBases(tRefSortingList(iRefListIndex)) = (*aRefListOfEntityBases)(iRefListIndex);
        tResultListOfEntityTypes(tRefSortingList(iRefListIndex)) = (*aRefListOfEntityTypes)(iRefListIndex);
        tResultListOfEntityActions(tRefSortingList(iRefListIndex)) = (*aRefListOfEntityActions)(iRefListIndex);
        tResultListOfOutputSpecs(tRefSortingList(iRefListIndex)) = (*aRefListOfOutputSpecs)(iRefListIndex);
    }

    //go new list and copy to resulting list
    for (uint iNewListIndex = 0; iNewListIndex < tNewSortingList.size(); iNewListIndex++)
    {
        tResultListOfEntityBases(tNewSortingList(iNewListIndex)) =   aNewListOfEntityBases(iNewListIndex);
        tResultListOfEntityTypes(tNewSortingList(iNewListIndex)) =   aNewListOfEntityTypes(iNewListIndex);
        tResultListOfEntityActions(tNewSortingList(iNewListIndex)) = aNewListOfEntityActions(iNewListIndex);
        tResultListOfOutputSpecs(tNewSortingList(iNewListIndex)) =   aNewListOfOutputSpecs(iNewListIndex);
    }

    // --- Merge Lists --------------------------------------------------------------- //

    // set cursor back to zero
    tResultListCursor = 0;
    // Go through old list and see where there are gaps, fill them
    for (uint tRefListIndex = 0; tRefListIndex < tResultListSize; tRefListIndex++ )
    {

        // if indices are missing
        if (tRefSortingList(tRefListIndex) > tResultListCursor + 1)
        {
            // insert missing columns
            for (uint ii = tResultListCursor; ii < tRefSortingList(tRefListIndex); ii++ )
            {
                this->insert_empty_column( aFullListOfOutputValues, ii );
            }
        }
        // update cursor
        tResultListCursor = tRefSortingList(tRefListIndex);
    }

    // insert missing columns into table of output values
    for (uint ii = tResultListCursor + 1; ii < tResultListSize; ii++ )
    {
        this->insert_empty_column( aFullListOfOutputValues, ii );
    }

    // go through new list and insert stuff where needed
    Cell<std::string> tResultListOfOutputValues(tResultListSize);
    for (uint iResultListIndex = 0; iResultListIndex < tResultListSize; iResultListIndex++ )
    {
        tResultListOfOutputValues(iResultListIndex) = "---";
    }

    for (uint iNewListIndex = 0; iNewListIndex < aNewListOfEntityBases.size(); iNewListIndex++ )
    {
        tResultListOfOutputValues(tNewSortingList(iNewListIndex)) = aNewListOfOutputValues(iNewListIndex);
    }

    // copy new lists to reference lists
    (*aRefListOfEntityBases)   = tResultListOfEntityBases;
    (*aRefListOfEntityTypes)   = tResultListOfEntityTypes;
    (*aRefListOfEntityActions) = tResultListOfEntityActions;
    (*aRefListOfOutputSpecs)   = tResultListOfOutputSpecs;

    // append new row to full list
    aFullListOfOutputValues->append({tResultListOfOutputValues});

}

//-----------------------------------------------------------------------------------------------------------//

// merges an old iteration table with a new one
void Query::insert_empty_column( Cell<Cell<std::string>> * aCellMatrix, uint aColumnIndex)
{
    uint tNumRows = aCellMatrix->size();
    uint tNumCols = (*aCellMatrix)(0).size();

    // insert column if index is in bounds
    if (aColumnIndex < tNumCols)
    {
        for (uint iRow = 0; iRow < tNumRows; iRow++ )
        {
            (*aCellMatrix)(iRow).insert(aColumnIndex, "---");
        }
    }

    // append column if index out of bounds
    else if (aColumnIndex == tNumCols)
    {
        for (uint iRow = 0; iRow < tNumRows; iRow++ )
        {
            (*aCellMatrix)(iRow).append({"---"});
        }
    }

    // error
    else
    {
        MORIS_ASSERT(false, "<MRS::IOS::cl_Query_Table::insert_empty_column>: invalid column index provided.");
    }
}

//-----------------------------------------------------------------------------------------------------------//

// goes through instance specified by ID, returns iteration start and end lines in cell matrix
Cell<Cell<uint>> Query::split_instances_into_iterations(uint aCurrentInstanceID)
{
    // get start end lines of function instance
    uint tCurrentInstanceStartLine = mInstanceStartEnd(aCurrentInstanceID)(0);
    uint tCurrentInstanceEndLine = mInstanceStartEnd(aCurrentInstanceID)(1);

    // get indentation level of current function instance
    uint tCurrentInstanceIndent = mIndentLevels(tCurrentInstanceStartLine);

    // initialize variables
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
                tIterationStartEnd.append( {{iLine, tCurrentInstanceEndLine - 1}} );

            // save previous line as end line of previous iteration
            if (tNumIterations > 0)
                tIterationStartEnd(tNumIterations-1)(1) = iLine - 1;

            // update number of iterations
            tNumIterations++;
        }
    } // for each line in instance

    return tIterationStartEnd;
}

//-----------------------------------------------------------------------------------------------------------//

// goes through logged info and finds instances of given type, returns number of instances
uint Query::find_instances(Cell<uint> * aInstanceIDs,
                           Cell<uint> * aInstanceIndents,
                           enum EntityBase aEntityBase,
                           enum EntityType aEntityType,
                           enum EntityAction aEntityAction)
{
    // initialize
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
                aInstanceIDs->push_back(ID);
                aInstanceIndents->push_back(mIndentLevels(tFirstLineOfInstance));
                tNumInstances++;
            }
        } // end check
    } // end search

    return tNumInstances;
}

//-----------------------------------------------------------------------------------------------------------//

void Query::extract_iteration(       Cell<enum EntityBase> * aListOfEntityBases,
                                     Cell<enum EntityType> * aListOfEntityTypes,
                                     Cell<enum EntityAction> * aListOfEntityActions,
                                     Cell<enum OutputSpecifier> * aListOfOutputSpecs,
                                     Cell<std::string> * aListOfOutputValues,
                               const Cell<Cell<uint>> aIterationStartEnd,
                               const uint aIteration,
                               const uint aCurrentInstanceID)
{

    //initialize
    bool tFirstEntry = true;
    uint iLine = aIterationStartEnd(aIteration-1)(0);

    // go through first iteration and collect initial information
    while (iLine <= aIterationStartEnd(aIteration-1)(1))
    {
        // log information belongs to current function ID, do not record if just text
        if ( (mFunctionIDs(iLine) == aCurrentInstanceID)
                &&  !(    (mOutputSpecifiers(iLine) == OutputSpecifier::FreeText)
                      || (mOutputSpecifiers(iLine) == OutputSpecifier::InfoText)
                      || (mOutputSpecifiers(iLine) == OutputSpecifier::DebugText)
                      || (mOutputSpecifiers(iLine) == OutputSpecifier::Warning) ) )
        {
            // copy everything to list of what is happening inside of iteration
            if (tFirstEntry)
            {
                aListOfEntityBases->resize( 1, mEntityBases(iLine) );
                aListOfEntityTypes->resize( 1, mEntityTypes(iLine) );
                aListOfEntityActions->resize( 1, mEntityActions(iLine) );
                aListOfOutputSpecs->resize( 1, mOutputSpecifiers(iLine) );
                aListOfOutputValues->resize( 1, mOutputValues(iLine) );
                tFirstEntry = false;
            }
            else
            {
                aListOfEntityBases->push_back( mEntityBases(iLine) );
                aListOfEntityTypes->push_back( mEntityTypes(iLine) );
                aListOfEntityActions->push_back( mEntityActions(iLine) );
                aListOfOutputSpecs->push_back( mOutputSpecifiers(iLine) );
                aListOfOutputValues->push_back( mOutputValues(iLine) );
            }
        }

        // check functions on indentation level one higher for important info
        if ( mIndentLevels(iLine) == mIndentLevels(mInstanceStartEnd(aCurrentInstanceID)(0)) + 1 )
        {
            // get end line
            uint tIndentedFunctionID = mFunctionIDs(iLine);
            uint tIndentedFunctionEnd = mInstanceStartEnd(tIndentedFunctionID)(1);

            // record function instance
            if (tFirstEntry)
            {
                aListOfEntityBases->resize( 1, mEntityBases(iLine) );
                aListOfEntityTypes->resize( 1, mEntityTypes(iLine) );
                aListOfEntityActions->resize( 1, mEntityActions(iLine) );
                aListOfOutputSpecs->resize( 1, mOutputSpecifiers(tIndentedFunctionEnd) );
                aListOfOutputValues->resize( 1, mOutputValues(tIndentedFunctionEnd) );
                tFirstEntry = false;
            }
            else
            {
                aListOfEntityBases->push_back( mEntityBases(iLine) );
                aListOfEntityTypes->push_back( mEntityTypes(iLine) );
                aListOfEntityActions->push_back( mEntityActions(iLine) );
                aListOfOutputSpecs->push_back( mOutputSpecifiers(tIndentedFunctionEnd) );
                aListOfOutputValues->push_back( mOutputValues(tIndentedFunctionEnd) );
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
                    tNumInnerIterations = std::stoi(mOutputValues(iLine));
                }

                // inner step is reported
                if ( (mFunctionIDs(iLine) == tIndentedFunctionID) &&
                        (mOutputSpecifiers(iLine) == OutputSpecifier::Step) )
                {
                    tNumInnerSteps = std::stoi(mOutputValues(iLine));
                }

                // inner residual is reported
                if ( (mFunctionIDs(iLine) == tIndentedFunctionID) && (mOutputSpecifiers(iLine) == OutputSpecifier::ResidualNorm) )
                {
                    if (!tInnerFunctionHasResidual)
                    {
                        tInitalInnerResidualNorm = std::stod(mOutputValues(iLine));
                        tInnerFunctionHasResidual = true;
                    }
                    else // if (tInnerFunctionHasResidual = true)
                        tFinalInnerResidualNorm = std::stod(mOutputValues(iLine));
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
                aListOfEntityBases->push_back( mEntityBases(iLine) );
                aListOfEntityTypes->push_back( mEntityTypes(iLine) );
                aListOfEntityActions->push_back( mEntityActions(iLine) );
                aListOfOutputSpecs->push_back( OutputSpecifier::Iteration );
                aListOfOutputValues->push_back( ios::stringify(tNumInnerIterations) );
            }

            // record steps if present
            if (tNumInnerSteps != 0)
            {
                // record function instance
                aListOfEntityBases->push_back( mEntityBases(iLine) );
                aListOfEntityTypes->push_back( mEntityTypes(iLine) );
                aListOfEntityActions->push_back( mEntityActions(iLine) );
                aListOfOutputSpecs->push_back( OutputSpecifier::Step );
                aListOfOutputValues->push_back( ios::stringify(tNumInnerIterations) );
            }

            // record Residual Drop if present
            if (tInnerFunctionHasResidual)
            {
                // record function instance
                aListOfEntityBases->push_back( mEntityBases(iLine) );
                aListOfEntityTypes->push_back( mEntityTypes(iLine) );
                aListOfEntityActions->push_back( mEntityActions(iLine) );
                aListOfOutputSpecs->push_back( OutputSpecifier::ResidualDrop );
                aListOfOutputValues->push_back( ios::stringify(tFinalInnerResidualNorm / tInitalInnerResidualNorm) );
            }

        } // end if: indented function

        // increment line cursor
        iLine++;

    } // end while: go through iteration
}

//-----------------------------------------------------------------------------------------------------------//
} // namespace std
} // namespace moris

