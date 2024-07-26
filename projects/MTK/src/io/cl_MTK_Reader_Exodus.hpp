/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Reader_Exodus.hpp
 *
 */

#pragma once

#include <exodusII.h>
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

// TODO
#include "cl_MTK_Mesh.hpp"    // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Sets_Info.hpp"

namespace moris::mtk
{

    class Reader_Exodus
    {
      public:
        Mesh*                        mMesh;
        int                          mExoID;
        uint                         mNumSpatialDimensions;
        Matrix< DDRMat >             mNodeCoordinates;
        Matrix< IdMat >              mNodeMap;
        Matrix< IdMat >              mNodeOwner;
        MtkSetsInfo                  mMtkMeshSets;
        Vector< MtkBlockSetInfo >    mBlockSetInfo;
        Vector< std::string >        mBlockDescription;
        Vector< Matrix< IndexMat > > mElemConn;
        Vector< Matrix< IdMat > >    mLocaltoGlobalElemMap;
        Vector< MtkSideSetInfo >     mSideSetInfo;
        Vector< Matrix< IndexMat > > mSideSetElements;
        Vector< std::string >        mSideSetNames;

        /**
         * Constructor
         */
        Reader_Exodus();

        /** Destructor */
        ~Reader_Exodus();

        /**
         * Changes how Exodus handles errors
         * @param abort Causes fatal errors to force program exit.
         * @param debug Causes certain messages to print for debugging use.
         * @param verbose Causes all error messages to print when true, otherwise no error messages will print.
         */
        void set_error_options( bool abort, bool debug, bool verbose );

        /**
         * Reads an Exodus file and dumps the result into a specified mesh.
         * @param aFileName The name of the file to be read from
         */
        void read_file( std::string aFileName );

      private:
        /**
         *  Opens an Exodus file and stores the ID for future operations
         *  @param aExodusFileName Name of the Exodus file.
         *  @param aVersion Version of the database. Current version is 4.72 as of programming.
         */
        void open_file( std::string aExodusFileName, float aVersion = 4.72 );

        /**
         * Closes the open Exodus database *and* renames it to the permanent file name stored under mPermFileName. This
         * must be called in order for the Exodus file to be able to be read properly.
         */
        void close_file();

        /**
         * Gets the CellTopology from a character array describing the element type.
         * @param const char* the element type written to the Exodus file from Writer_Exodus
         * @return aCellTopology the type of element in MTK.
         */
        CellTopology get_cell_topology( char* aElementType );
    };

}    // namespace moris::mtk
