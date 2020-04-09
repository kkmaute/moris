#ifndef MORIS_IOS_CL_QUERY_HPP_
#define MORIS_IOS_CL_QUERY_HPP_

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
    Cell<real> mOutputValues;

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

    //-----------------------------------------------------------------------------------------------------------//
    // PUBLIC CONSTRUCTOR / DESTRUCTOR
    //-----------------------------------------------------------------------------------------------------------//
    public:

    Query(std::string aFileNameRead)
    {
        // save file name
        mFileNameRead = aFileNameRead;
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

    void initialize();

    void tree_query(std::string aFileNameWrite);

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












