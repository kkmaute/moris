/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Query.hpp
 *
 */

#ifndef MORIS_IOS_CL_QUERY_HPP_
#define MORIS_IOS_CL_QUERY_HPP_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>

// Define Cells
#include "cl_Cell.hpp"

// Define uint, real, etc.
#include "moris_typedefs.hpp"

// Define enums used
#include "cl_Tracer_Enums.hpp"

namespace moris
{
namespace ios
{

//-----------------------------------------------------------------------------------------------------------//
class Query
    {
    private:

    //-----------------------------------------------------------------------------------------------------------//
    // PRIVATE MEMBERS
    //-----------------------------------------------------------------------------------------------------------//

    // files for reading and writing
    std::string mFileNameRead;
    std::ifstream mLogFileRead;

    // number of lines in read file
    uint mNumLines;

    // cell arrays containing info from log file in memory
    Cell<uint> mIndentLevels;
    Cell<uint> mFunctionIDs;
    Cell<enum EntityBase> mEntityBases;
    Cell<enum EntityType> mEntityTypes;
    Cell<enum EntityAction> mEntityActions;
    Cell<enum OutputSpecifier> mOutputSpecifiers;
    Cell<std::string> mOutputValues;

    Cell<Cell<uint>> mInstanceStartEnd;

    //-----------------------------------------------------------------------------------------------------------//
    // PRIVATE METHODS
    //-----------------------------------------------------------------------------------------------------------//

    uint get_max_func_ID();

    void isolate_functions();

    uint extract_info_from_log_file();

    std::string create_empty_line(uint aIndentationLevel);

    void copy_header(std::ofstream * aLogFileWrite);

    void skip_header();

    uint find_instances(Cell<uint> * aInstanceIDs,
                        Cell<uint> * aInstanceIndents,
                        enum EntityBase aEntityBase,
                        enum EntityType aEntityType,
                        enum EntityAction aEntityAction);

    Cell<Cell<uint>> split_instances_into_iterations(uint aCurrentInstanceID);

    void extract_iteration(      Cell<enum EntityBase> * aListOfEntityBases,
                                 Cell<enum EntityType> * aListOfEntityTypes,
                                 Cell<enum EntityAction> * aListOfEntityActions,
                                 Cell<enum OutputSpecifier> * aListOfOutputSpecs,
                                 Cell<std::string> * aListOfOutputValues,
                           const Cell<Cell<uint>> aIterationStartEnd,
                           const uint aIteration,
                           const uint aCurrentInstanceID);

    void merge_iteration_outputs( Cell<enum EntityBase>      * aRefListOfEntityBases,
                                  Cell<enum EntityType>      * aRefListOfEntityTypes,
                                  Cell<enum EntityAction>    * aRefListOfEntityActions,
                                  Cell<enum OutputSpecifier> * aRefListOfOutputSpecs,
                                  Cell<Cell<std::string>>    * aFullListOfOutputValues,

                                  Cell<enum EntityBase>      aNewListOfEntityBases,
                                  Cell<enum EntityType>      aNewListOfEntityTypes,
                                  Cell<enum EntityAction>    aNewListOfEntityActions,
                                  Cell<enum OutputSpecifier> aNewListOfOutputSpecs,
                                  Cell<std::string>          aNewListOfOutputValues);

    void insert_empty_column( Cell<Cell<std::string>> * aCellMatrix, uint aColumnIndex);

    void write_instance_table_header(std::ofstream * aLogFileWrite,
                                     Cell<enum EntityBase> aRefListOfEntityBases,
                                     Cell<enum EntityType> aRefListOfEntityTypes,
                                     Cell<enum EntityAction> aRefListOfEntityActions,
                                     Cell<enum OutputSpecifier> aRefListOfOutputSpecs);

    //-----------------------------------------------------------------------------------------------------------//
    // PUBLIC CONSTRUCTOR / DESTRUCTOR
    //-----------------------------------------------------------------------------------------------------------//
    public:

    Query()
    {
    }

    ~Query()
    {
        // close file
        mLogFileRead.close();

        std::cout << "Queried info dumped. \n";
    }

    //-----------------------------------------------------------------------------------------------------------//
    // PUBLIC METHODS
    //-----------------------------------------------------------------------------------------------------------//

    void initialize(std::string aFileNameRead);

    void run(int  & argc, char * argv[] );

    void tree_query(std::string aFileNameWrite, bool aSuppressText = true);

    void table_query(std::string aFileNameWrite,
                     enum EntityBase aEntityBase,
                     enum EntityType aEntityType,
                     enum EntityAction aEntityAction);

    //-----------------------------------------------------------------------------------------------------------//
    }; // end class

//-----------------------------------------------------------------------------------------------------------//
} // namespace std
} // namespace moris

#endif    /* MORIS_IOS_CL_QUERY_HPP_ */

