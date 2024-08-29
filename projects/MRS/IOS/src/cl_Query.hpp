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
#include "cl_Vector.hpp"

// Define uint, real, etc.
#include "moris_typedefs.hpp"

// Define enums used
#include "cl_Tracer_Enums.hpp"

namespace moris::ios
{

    //-----------------------------------------------------------------------------------------------------------//
    class Query
    {
      private:
        //-----------------------------------------------------------------------------------------------------------//
        // PRIVATE MEMBERS
        //-----------------------------------------------------------------------------------------------------------//

        // files for reading and writing
        std::string   mFileNameRead;
        std::ifstream mLogFileRead;

        // number of lines in read file
        uint mNumLines;

        // cell arrays containing info from log file in memory
        Vector< uint >                 mIndentLevels;
        Vector< uint >                 mFunctionIDs;
        Vector< enum EntityBase >      mEntityBases;
        Vector< enum EntityType >      mEntityTypes;
        Vector< enum EntityAction >    mEntityActions;
        Vector< enum OutputSpecifier > mOutputSpecifiers;
        Vector< std::string >          mOutputValues;

        Vector< Vector< uint > > mInstanceStartEnd;

        //-----------------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //-----------------------------------------------------------------------------------------------------------//

        uint get_max_func_ID();

        void isolate_functions();

        uint extract_info_from_log_file();

        std::string create_empty_line( uint aIndentationLevel );

        void copy_header( std::ofstream *aLogFileWrite );

        void skip_header();

        uint find_instances( Vector< uint > *aInstanceIDs,
                Vector< uint >              *aInstanceIndents,
                enum EntityBase              aEntityBase,
                enum EntityType              aEntityType,
                enum EntityAction            aEntityAction );

        Vector< Vector< uint > > split_instances_into_iterations( uint aCurrentInstanceID );

        void extract_iteration( Vector< enum EntityBase > *aListOfEntityBases,
                Vector< enum EntityType >                 *aListOfEntityTypes,
                Vector< enum EntityAction >               *aListOfEntityActions,
                Vector< enum OutputSpecifier >            *aListOfOutputSpecs,
                Vector< std::string >                     *aListOfOutputValues,
                const Vector< Vector< uint > >            &aIterationStartEnd,
                const uint                                 aIteration,
                const uint                                 aCurrentInstanceID );

        void merge_iteration_outputs( Vector< enum EntityBase > *aRefListOfEntityBases,
                Vector< enum EntityType >                       *aRefListOfEntityTypes,
                Vector< enum EntityAction >                     *aRefListOfEntityActions,
                Vector< enum OutputSpecifier >                  *aRefListOfOutputSpecs,
                Vector< Vector< std::string > >                 *aFullListOfOutputValues,

                Vector< enum EntityBase >      aNewListOfEntityBases,
                Vector< enum EntityType >      aNewListOfEntityTypes,
                Vector< enum EntityAction >    aNewListOfEntityActions,
                Vector< enum OutputSpecifier > aNewListOfOutputSpecs,
                Vector< std::string >          aNewListOfOutputValues );

        void insert_empty_column( Vector< Vector< std::string > > *aCellMatrix, uint aColumnIndex );

        void write_instance_table_header( std::ofstream *aLogFileWrite,
                Vector< enum EntityBase >                aRefListOfEntityBases,
                Vector< enum EntityType >                aRefListOfEntityTypes,
                Vector< enum EntityAction >              aRefListOfEntityActions,
                Vector< enum OutputSpecifier >           aRefListOfOutputSpecs );

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

        void initialize( std::string aFileNameRead );

        void run( int &argc, char *argv[] );

        void tree_query( const std::string &aFileNameWrite, bool aSuppressText = true );

        void table_query( const std::string &aFileNameWrite,
                enum EntityBase              aEntityBase,
                enum EntityType              aEntityType,
                enum EntityAction            aEntityAction );

        //-----------------------------------------------------------------------------------------------------------//
    };    // end class

    //-----------------------------------------------------------------------------------------------------------//
}    // namespace moris::ios

#endif /* MORIS_IOS_CL_QUERY_HPP_ */
