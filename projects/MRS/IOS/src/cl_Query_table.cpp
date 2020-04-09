//
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <cstdio>
//#include <string>
//#include <cstring>
//#include <sstream>
//
//// include moris assert functions
//#include "fn_assert.hpp"
//
//// Define Cells
//#include "cl_Cell.hpp"
//
//// Define uint, real, etc.
//#include "typedefs.hpp"
//
//// needed?
//#include "IO_Tools.hpp"
//
//// Define enums used
//#include "cl_Tracer_Enums.hpp"
//
//#include "cl_Query.hpp"
//
//
//namespace moris
//{
//namespace ios
//{
//
////-----------------------------------------------------------------------------------------------------------//
//// FUNCTION DEFINITIONS
////-----------------------------------------------------------------------------------------------------------//
//
//// function to create tabled log file for instances of a given type
//void Query::table_query(std::string aFileNameWrite,
//                        enum EntityBase aEntityBase,
//                        enum EntityType aEntityType,
//                        enum EntityAction aEntityAction)
//{
//
//
//
//    // Prepare reading file ---------------------------------------------------------- //
//
//    std::cout << "Starting query to create process tree ... ";
//
//    // create text file to write to
//    std::ofstream tLogFileWrite;
//    tLogFileWrite.open(aFileNameWrite);
//    if( !tLogFileWrite )
//        MORIS_ASSERT( false, "<MRS::IOS::fn_Query::table_query>: Write file could not be opened.");
//
//    // Read and Copy Header ---------------------------------------------------------- //
//
//    copy_header(& tLogFileWrite);
//
//
//    // Starting Query ---------------------------------------------------------------- //
//
//    // initialize
//    uint tColumnWidth = 20;
//    std::string tCurrentLine;
//    uint tEndLine = 0;
//    std::string tCurrentEntity = "";
//
//
//    /////////////////////////////////
//    // Look for function instances //
//    /////////////////////////////////
//
//    // initialize
//    Cell<real> tInstanceIDs(0);
//    Cell<real> tInstanceIndents(0);
//    uint tNumInstances = 0;
//    bool tTypesMatch = false;
//    bool tActionsMatch = false;
//
//    // get number of function instances
//    //uint tMaxFuncID = get_max_func_ID(mNumLines, mFunctionIDs);
//    uint tMaxFuncID = mInstanceStartEnd.size() - 1;
//
//    // go through every function instance
//    for (uint ID = 0; ID <= tMaxFuncID; ID++)
//    {
//        uint tFirstLineOfInstance = mInstanceStartEnd(ID)(0);
//
//        // check if entities are the same
//        if ( mEntityBases(tFirstLineOfInstance) == aEntityBase )
//        {
//            // check if entity types match
//            if ( aEntityType == EntityType::Arbitrary )
//                tTypesMatch = true;
//            else
//                tTypesMatch = ( aEntityType == mEntityTypes(tFirstLineOfInstance) );
//
//            // check if entity actions match
//            if ( aEntityType == EntityType::Arbitrary )
//                tActionsMatch = true;
//            else
//                tActionsMatch = ( aEntityAction == mEntityActions(tFirstLineOfInstance) );
//
//            // if everything matches, save instance
//            if ( tTypesMatch && tActionsMatch )
//            {
//                tInstanceIDs.push_back(ID);
//                tInstanceIndents.push_back(mIndentLevels(tFirstLineOfInstance));
//                tNumInstances++;
//            }
//        } // end check
//    } // end search
//
//
//    /////////////////////////
//    // Write Output Header //
//    /////////////////////////
//
//    // write out general header
//    tLogFileWrite << "Looking for instances of function: " << get_enum_str(aEntityBase) << " - "
//                                                           << get_enum_str(aEntityType) << " - "
//                                                           << get_enum_str(aEntityAction) << "\n";
//    tLogFileWrite << "Number of instances found: " << tNumInstances << ". \n";
//
//    tLogFileWrite << "Indentation levels of respective instances: [" << tInstanceIndents(0);
//    for (uint iInstance = 1; iInstance < tNumInstances; iInstance++)
//    {
//        tLogFileWrite << ", " << tInstanceIndents(iInstance);
//    }
//    tLogFileWrite << "]\n";
//
//    tLogFileWrite << "IDs of respective instances: [" << tInstanceIDs(0);
//    for (uint iInstance = 1; iInstance < tNumInstances; iInstance++)
//    {
//        tLogFileWrite << ", " << tInstanceIDs(iInstance);
//    }
//    tLogFileWrite << "]\n\n";
//
//
//    ///////////////////////////
//    // Record what is inside //
//    ///////////////////////////
//
//    // for each instance
//    for (uint iInstance = 0; iInstance < tNumInstances; iInstance++)
//    {
//
//
//        /////////////////////////
//        // Split up Iterations //
//        /////////////////////////
//
//        // get start end lines of function instance
//        uint tCurrentInstanceID = tInstanceIDs(iInstance);
//        uint tCurrentInstanceIndent = tInstanceIndents(iInstance);
//        uint tCurrentInstanceStartLine = mInstanceStartEnd(tCurrentInstanceID)(0);
//        uint tCurrentInstanceEndLine = mInstanceStartEnd(tCurrentInstanceID)(1);
//
//        // initialize cell array storing start and end lines of each iteration
//        uint tNumIterations = 0;
//        Cell<Cell<uint>> tIterationStartEnd;
//        tIterationStartEnd.resize(1);
//        tIterationStartEnd(0).resize(2);
//        tIterationStartEnd(0)(0) = tCurrentInstanceStartLine + 1;
//        tIterationStartEnd(0)(1) = tCurrentInstanceEndLine - 1;
//
//        // look through function instance and find iteration markers
//        for (uint iLine = tCurrentInstanceStartLine; iLine < tCurrentInstanceEndLine; iLine++)
//        {
//            if ( ((mOutputSpecifiers(iLine) == OutputSpecifier::Iteration) || (mOutputSpecifiers(iLine) == OutputSpecifier::Step))
//                    && (mIndentLevels(iLine) == tCurrentInstanceIndent) )
//            {
//                // append new iteration
//                if (tNumIterations == 0)
//                    tIterationStartEnd(0)(0) = iLine;
//                else
//                    tIterationStartEnd.append( {iLine, tCurrentInstanceEndLine - 1} );
//
//                // save previous line as end line of previous iteration
//                if (tNumIterations > 0)
//                    tIterationStartEnd(tNumIterations-1)(1) = iLine - 1;
//
//                tNumIterations++;
//            }
//        } // for each line in instance
//
//
//
//        ///////////////////////////
//        // Go through iterations //
//        ///////////////////////////
//
//        // initialize lists of outputs
//        Cell<enum EntityBase> tListOfEntityBases;
//        Cell<enum EntityType> tListOfEntityTypes;
//        Cell<enum EntityAction> tListOfEntityActions;
//        Cell<enum OutputSpecifier> tListOfOutputSpecs;
//        Cell<Cell<real>> tListOfOutputValues(1);
//        bool tFirstEntry = true;
//
//        /////////////////////////////////////////////////////////
//        // Go through first iteration and create reference set //
//        /////////////////////////////////////////////////////////
//
//        uint iLine = tIterationStartEnd(0)(0);
//
//        // go through first iteration and collect initial information
//        while (iLine <= tIterationStartEnd(0)(1))
//        {
//            // log information belongs to current function ID
//            if ( mFunctionIDs(iLine) == tCurrentInstanceID )
//            {
//                // copy everything to list of what is happening inside of iteration
//                if (tFirstEntry)
//                {
//                    tListOfEntityBases.resize( 1, mEntityBases(iLine) );
//                    tListOfEntityTypes.resize( 1, mEntityTypes(iLine) );
//                    tListOfEntityActions.resize( 1, mEntityActions(iLine) );
//                    tListOfOutputSpecs.resize( 1, mOutputSpecifiers(iLine) );
//                    tListOfOutputValues(0).resize( 1, mOutputValues(iLine) );
//                    tFirstEntry = false;
//                }
//                else
//                {
//                    tListOfEntityBases.push_back( mEntityBases(iLine) );
//                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
//                    tListOfEntityActions.push_back( mEntityActions(iLine) );
//                    tListOfOutputSpecs.push_back( mOutputSpecifiers(iLine) );
//                    tListOfOutputValues(0).push_back( mOutputValues(iLine) );
//                }
//            }
//
//            // check functions on indentation level one higher for important info
//            if ( mIndentLevels(iLine) == tCurrentInstanceIndent + 1 )
//            {
//                // get end line
//                uint tIndentedFunctionID = mFunctionIDs(iLine);
//                uint tIndentedFunctionEnd = mInstanceStartEnd(tIndentedFunctionID)(1);
//
//                // record function instance
//                if (tFirstEntry)
//                {
//                    tListOfEntityBases.resize( 1, mEntityBases(iLine) );
//                    tListOfEntityTypes.resize( 1, mEntityTypes(iLine) );
//                    tListOfEntityActions.resize( 1, mEntityActions(iLine) );
//                    tListOfOutputSpecs.resize( 1, mOutputSpecifiers(tIndentedFunctionEnd) );
//                    tListOfOutputValues(0).resize( 1, mOutputValues(tIndentedFunctionEnd) );
//                    tFirstEntry = false;
//                }
//                else
//                {
//                    tListOfEntityBases.push_back( mEntityBases(iLine) );
//                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
//                    tListOfEntityActions.push_back( mEntityActions(iLine) );
//                    tListOfOutputSpecs.push_back( mOutputSpecifiers(tIndentedFunctionEnd) );
//                    tListOfOutputValues(0).push_back( mOutputValues(tIndentedFunctionEnd) );
//                }
//
//                // initialize
//                uint tNumInnerIterations = 0;
//                uint tNumInnerSteps = 0;
//                bool tInnerFunctionHasResidual = false;
//                real tInitalInnerResidualNorm = - 1.0;
//                real tFinalInnerResidualNorm = - 1.0;
//
//                // increment line cursor to skip sign in line
//                iLine++;
//
//                // go through inner function
//                while (iLine <= tIndentedFunctionEnd)
//                {
//                    // inner iteration is reported
//                    if ( (mFunctionIDs(iLine) == tIndentedFunctionID) &&
//                         (mOutputSpecifiers(iLine) == OutputSpecifier::Iteration) )
//                    {
//                        tNumInnerIterations = mOutputValues(iLine);
//                    }
//
//                    // inner step is reported
//                    if ( (mFunctionIDs(iLine) == tIndentedFunctionID) &&
//                         (mOutputSpecifiers(iLine) == OutputSpecifier::Step) )
//                    {
//                        tNumInnerSteps = mOutputValues(iLine);
//                    }
//
//                    // inner residual is reported
//                    if ( (mFunctionIDs(iLine) == tIndentedFunctionID) && (mOutputSpecifiers(iLine) == OutputSpecifier::ResidualNorm) )
//                    {
//                        if (!tInnerFunctionHasResidual)
//                        {
//                            tInitalInnerResidualNorm = mOutputValues(iLine);
//                            tInnerFunctionHasResidual = true;
//                        }
//                        else // if (tInnerFunctionHasResidual = true)
//                            tFinalInnerResidualNorm = mOutputValues(iLine);
//                    }
//
//                    // increment line cursor
//                    iLine++;
//
//                } // end while: inner function
//
//                // decrement line cursor to last line of indented function instance
//                iLine--;
//
//                // record iterations if present
//                if (tNumInnerIterations != 0)
//                {
//                    // record function instance
//                    tListOfEntityBases.push_back( mEntityBases(iLine) );
//                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
//                    tListOfEntityActions.push_back( mEntityActions(iLine) );
//                    tListOfOutputSpecs.push_back( OutputSpecifier::Iteration );
//                    tListOfOutputValues(0).push_back( tNumInnerIterations );
//                }
//
//                // record steps if present
//                if (tNumInnerSteps != 0)
//                {
//                    // record function instance
//                    tListOfEntityBases.push_back( mEntityBases(iLine) );
//                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
//                    tListOfEntityActions.push_back( mEntityActions(iLine) );
//                    tListOfOutputSpecs.push_back( OutputSpecifier::Step );
//                    tListOfOutputValues(0).push_back( tNumInnerIterations );
//                }
//
//                // record Residual Drop if present
//                if (tInnerFunctionHasResidual)
//                {
//                    // record function instance
//                    tListOfEntityBases.push_back( mEntityBases(iLine) );
//                    tListOfEntityTypes.push_back( mEntityTypes(iLine) );
//                    tListOfEntityActions.push_back( mEntityActions(iLine) );
//                    tListOfOutputSpecs.push_back( OutputSpecifier::ResidualDrop );
//                    tListOfOutputValues(0).push_back( tFinalInnerResidualNorm / tInitalInnerResidualNorm );
//                }
//
//            } // end if: indented function
//
//            // increment line cursor
//            iLine++;
//
//        } // end while: go through iteration
//
//
//        //////////////////
//        // Write Tables //
//        //////////////////
//
//        uint tMaxColumnWidth = 18;
//        uint tNumCols = tListOfEntityBases.size();
//
//
//        // write general info
//        tLogFileWrite << "Instance ID: " << tInstanceIDs(iInstance) << "\n";
//        tLogFileWrite << "Execution Time: " << mOutputValues(mInstanceStartEnd(tInstanceIDs(iInstance))(1)) << " seconds.\n";
//
//        // write header for a table
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "|";
//            for (uint iChar = 0; iChar < tMaxColumnWidth+2; iChar++ )
//            tLogFileWrite << "-";
//        }
//        tLogFileWrite << "|\n";
//
//        // write header for a table
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfEntityBases(iCol)) << " ";
//        }
//        tLogFileWrite << "|\n";
//
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfEntityTypes(iCol)) << " ";
//        }
//        tLogFileWrite << "|\n";
//
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfEntityActions(iCol)) << " ";
//        }
//        tLogFileWrite << "|\n";
//
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "|";
//            for (uint iChar = 0; iChar < tMaxColumnWidth+2; iChar++ )
//            tLogFileWrite << "-";
//        }
//        tLogFileWrite << "|\n";
//
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << get_enum_str(tListOfOutputSpecs(iCol)) << " ";
//        }
//        tLogFileWrite << "|\n";
//
//        for (uint iCol = 0; iCol < tNumCols; iCol++ )
//        {
//            tLogFileWrite << "|";
//            for (uint iChar = 0; iChar < tMaxColumnWidth+2; iChar++ )
//            tLogFileWrite << "-";
//        }
//        tLogFileWrite << "|\n";
//
//        // write output
//        uint tNumRows = tListOfOutputValues.size();
//        for (uint iRow = 0; iRow < tNumRows; iRow++ )
//        {
//            for (uint iCol = 0; iCol < tNumCols; iCol++ )
//            {
//                // if integer
//                if (   (tListOfOutputSpecs(iCol) == OutputSpecifier::Iteration)
//                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Step)
//                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Count)
//                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Error)
//                    || (tListOfOutputSpecs(iCol) == OutputSpecifier::Restart) )
//                {
//                    tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << tListOfOutputValues(iRow)(iCol) << " ";
//                }
//                // else: float
//                else
//                {
//                    tLogFileWrite << "| " << std::setw(tMaxColumnWidth) << std::scientific << tListOfOutputValues(iRow)(iCol) << " ";
//                    tLogFileWrite.unsetf(std::ios::fixed | std::ios::scientific);
//                }
//            }
//            tLogFileWrite << "|\n";
//        } // end: for each row
//
//        // empty line between tables
//        tLogFileWrite << "\n";
//
//    } // end: for each instance
//
//
//
//
//    // close everything and end ------------------------------------------------------ //
//    tLogFileWrite.close();
//    std::cout << "Success, Done. \n";
//}
//
//
////-----------------------------------------------------------------------------------------------------------//
//} // namespace std
//} // namespace moris
//
//
//
//
//
//
//
//
//
//
//
//
//
