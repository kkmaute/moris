/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Exodus_IO_Helper.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_EXODUS_IO_HELPER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_EXODUS_IO_HELPER_HPP_

#include <exodusII.h>
#include <ne_nemesisI.h>

#include "cl_Matrix.hpp"

namespace moris
{
    namespace mtk
    {
        class Exodus_IO_Helper
        {
          private:
            int  mErrFlag     = true;
            bool mVerbose     = false;
            bool mBuildGlobal = false;

            // general mesh info
            int         mNumDim      = -1;
            int         mNumNodes    = -1;
            int         mNumElem     = -1;
            int         mNumElemBlk  = -1;
            int         mNumNodeSets = -1;
            int         mNumSideSets = -1;
            std::string mTitle;
            real        mCharLength;

            // Coordinates
            Matrix< DDRMat > mX;
            Matrix< DDRMat > mY;
            Matrix< DDRMat > mZ;

            // GLobal information
            int mNumNodesGlobal    = -1;
            int mNumElemsGlobal    = -1;
            int mNumElemBlksGlobal = -1;
            int mNumNodeSetsGlobal = -1;
            int mNumSideSetsGlobal = -1;

            std::vector< int > mGlobalElemBlkIds;
            std::vector< int > mGlobalElemBlkCnts;

            // Communication Map information
            int mExoFileId        = -1;
            int mNumInternalNodes = -1;
            int mNumBorderNodes   = -1;
            int mNumExternalNodes = -1;
            int mNumInternalElems = -1;
            int mNumBorderElems   = -1;
            int mNumNodeCmaps     = -1;
            int mNumElemCmaps     = -1;

            Matrix< IdMat > mNodeCmapIds;
            Matrix< IdMat > mNodeCmapNodeCnts;
            Matrix< IdMat > mElemCmapIds;
            Matrix< IdMat > mElemCmapElemCnts;

            Vector< Matrix< IdMat > > mNodeCmapNodeIds;
            Vector< Matrix< IdMat > > mNodeCmapProcIds;

            std::vector< int > mGlobalNodesetIds;
            std::vector< int > mNumGlobalNodeCounts;
            std::vector< int > mNumGlobalNodeDfCounts;

            Matrix< IdMat > mElemMapi;
            Matrix< IdMat > mElemMapb;
            Matrix< IdMat > mNodeMapi;
            Matrix< IdMat > mNodeMapb;
            Matrix< IdMat > mNodeMape;
            Matrix< IdMat > mNodeNumMap;
            Matrix< IdMat > mElemNumMap;

            // Node Sets
            std::vector< int >                mNodeSetIds;
            std::vector< int >                mNodeSetNEntries;
            std::vector< int >                mNodeSetDistFactors;
            std::vector< char >               mNodeSetNamesMemory;
            std::vector< char * >             mNodeSetNamePtrs;
            std::vector< Matrix< IndexMat > > mNodeSetNodeIds;

            // Side Sets
            std::vector< int >                mSideSetIds;
            std::vector< int >                mSideSetNEntries;
            std::vector< int >                mSideSetDistFactors;
            std::vector< char >               SideSetNamesMemory;
            std::vector< char * >             mSideSetNamePtrs;
            std::vector< Matrix< IndexMat > > mSideSetElemIds;
            std::vector< Matrix< IndexMat > > mSideSetSideOrd;

            // Block sets
            std::vector< int >                mBlockIds;
            std::vector< int >                mBlockSetNEntries;
            std::vector< int >                mBlockSetNNodesPerEntry;
            std::vector< int >                mBlockSetNedgesPerEntry;
            std::vector< int >                mBlockSetNfacesPerEntry;
            std::vector< int >                mBlockSetNattrPerEntry;
            std::vector< char >               mBlockNamesMemory;
            std::vector< char * >             mBlockNamesPtrs;
            std::vector< char >               mBlockElemTypeNamesMemory;
            std::vector< char * >             mBlockElemTypeNamesPtrs;
            std::vector< Matrix< IndexMat > > mBlockSetNodeConn;
            std::vector< Matrix< IndexMat > > mBlockSetEdgeConn;
            std::vector< Matrix< IndexMat > > mBlockSetFaceConn;

            // Timing information
            int  mNumTimeSteps  = -1;
            int  mTimeStepIndex = 0;
            real mTimeValue     = -1.0;

            // Nodal Fields
            int                             mNumNodalVars = -1;
            std::vector< char >             mNodeFieldNamesMemory;
            std::vector< char * >           mNodeFieldNamePtrs;
            std::vector< Matrix< DDRMat > > mFieldsNodalVars;

            // Global Variables
            int                   mNumGlobalVars = -1;
            std::vector< char >   mGlobalVariableNamesMemory;
            std::vector< char * > mGlobalVariableNamePtrs;
            Matrix< DDRMat >      mGlobalVariables;

            void
            get_init_mesh_data();

            void
            get_load_bal_parameters();

            void
            get_cmap_params();

            void
            get_node_cmap();

            void get_init_global();

            void
            get_eb_info_global();

            void
            get_node_coords();

            void
            get_ns_param_global();

            void
            get_node_map();

            void
            get_node_id_map();

            void
            get_elem_id_map();

            void
            get_set_information();

            void
            get_nodal_fields();

            void
            reload_nodal_fields();

            void
            get_global_variables();

            void
            reload_global_variables();

            void
            copy_coordinates( int tNewExoFileId );

            void
            copy_node_sets( int aNewExoFileId );

            void
            copy_side_sets( int aNewExoFileId );

            void
            copy_block_sets( int aNewExoFileId );

            void
            copy_nodal_fields(
                    int             aNewExoFileId,
                    ex_init_params &init_params );

            /*
             * @brief gets the length scale of the problem posed
             *
             */
            void
            get_characteristic_length();

            /*
             * from seacas test rd_wt_mesh.c
             *
             */
            std::string get_file_name(
                    const std::string &base,
                    const std::string &other = "" );

            static void setup_names(
                    int                    nnames,
                    std::vector< char >   &storage,
                    std::vector< char * > &ptrs );

          public:
            Exodus_IO_Helper(
                    const std::string &aExodusFile,
                    const int          aTimeStepIndex = 0,
                    const bool         aBuildGlobal   = false,
                    const bool         aVerbose       = false );

            ~Exodus_IO_Helper();

            //------------------------------------------------------------------------------
            /*
             * @brief returns number of dimension
             *
             */

            int get_number_of_dimensions()
            {
                return mNumDim;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns number of dimension
             *
             */

            int get_number_of_nodes()
            {
                return mNumNodes;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns number of elements
             *
             */

            int get_number_of_elements()
            {
                return mNumElem;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns number of elements
             *
             */

            int get_number_of_blocks()
            {
                return mNumElemBlk;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns number of node sets
             *
             */

            int get_number_of_node_sets()
            {
                return mNumNodeSets;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns number of node sets
             *
             */

            int get_number_of_side_sets()
            {
                return mNumSideSets;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns absolute time for currently loaded time step
             *
             */

            real get_time_value()
            {
                return mTimeValue;
            }

            //------------------------------------------------------------------------------
            /*
             * @brief returns vector with coordinates of node, size: number of dimensions x 1
             *
             * @param[ in ] aNodeId   id of node
             */

            moris::Matrix< DDRMat > get_nodal_coordinate( uint aNodeId );

            //------------------------------------------------------------------------------
            /*
             * @brief returns node index given a node Id
             *
             * @param[ in ] aNodeId   id of node
             */

            uint
            get_node_index_by_Id( uint aNodeId );

            //------------------------------------------------------------------------------
            /*
             * @brief returns nodal field value of node defined by its Id
             *
             * @param[ in ] aNodeId        id of node
             * @param[ in ] aFieldInd      index of field
             * @param[ in ] aTimeStepIndex index of time step; default 0
             */

            real
            get_nodal_field_value(
                    uint aNodeId,
                    uint aFieldIndex,
                    uint aTimeStepIndex = 0 );

            //------------------------------------------------------------------------------
            /*
             * @brief returns nodal field value based on coordinate location within an error threshold
             *        based on a characteristic length of the problem or a user defined error
             *
             * @param[ in ] aCoords         approximate nodal coordinate
             * @param[ in ] aFieldInd       index of field
             * @param[ in ] aTimeStepIndex  index of time step; default 0
             * @param[ in ] aThresholdScale Error threshold scaling for nodal coordinate;
             *                              default is 1e-6
             */

            real
            get_nodal_field_value_by_coords(
                    moris::Matrix< DDRMat > aCoords,
                    uint                    aFieldIndex,
                    uint                    aTimeStepIndex = 0,
                    real                    aThreshold     = 1.0e-6 );

            //------------------------------------------------------------------------------
            /*
             * @brief returns field value for all nodes
             *
             * @param[ in ] aFieldInd   index of field
             * @param[ in ] aTimeStepIndex index of time step; default 0
             */

            const Matrix< DDRMat > &
            get_nodal_field_vector(
                    uint aFieldIndex    = 0,
                    uint aTimeStepIndex = 0 );

            //------------------------------------------------------------------------------
            /*
             * @brief returns value of global variable for given a variable inex
             *
             * @param[ in ] aNodeId   id of node
             */

            real
            get_global_variable(
                    uint aGlobalVariableIndex,
                    uint aTimeStepIndex );

            //------------------------------------------------------------------------------
            /*
             * @brief returns name of global variable for given a variable inex
             *
             * @param[ in ] aNodeId   id of node
             */

            const char *
            get_global_variable_name( uint aGlobalVariableIndex );

            //------------------------------------------------------------------------------
            /*
             * Create a new exodus file and the information for an element communication map
             * Copy the exodus file in this Exodus_IO_Helper to a new one with the
             * communication map. This is because once the exodus file has been setup,
             * the element communication map cannot be appended
             *
             */
            void
            create_new_exo_with_elem_cmaps_from_existing_exo(
                    std::string     &aFileName,
                    Matrix< IdMat > &aElementIds,
                    Matrix< IdMat > &aElementSideOrds,
                    Matrix< IdMat > &aSharedProcIds );

            //------------------------------------------------------------------------------
            /*
             * @brief returns index of field given its name
             *
             * param[ in ] aFieldName   Name of field
             *
             */
            uint
            get_field_index_by_name( std::string aFileName );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the nodal field value based on the name of the field
             *
             * @param aNodeId
             * @param aFieldName
             * @param aTimeStepIndex
             * @return real
             */

            real
            get_nodal_field_value(
                    uint const        &aNodeId,
                    std::string const &aFieldName,
                    uint const        &aTimeStepIndex );
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_EXODUS_IO_HELPER_HPP_ */
