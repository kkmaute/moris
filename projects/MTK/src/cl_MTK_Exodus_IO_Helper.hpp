/*
 * cl_MTK_Exodus_IO_Helper.hpp
 *
 *  Created on: Nov 26, 2018
 *      Author: doble
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

                int  mErrFlag;
                bool mVerbose;

                // general mesh info
                int mNumDim;
                int mNumNodes;
                int mNumElem;
                int mNumElemBlk;
                int mNumNodeSets;
                int mNumSideSets;
                char * mTitle;

                // Coordinates
                Matrix<DDRMat> mX;
                Matrix<DDRMat> mY;
                Matrix<DDRMat> mZ;

                // GLobal information
                int mNumNodesGlobal;
                int mNumElemsGlobal;
                int mNumElemBlksGlobal;
                int mNumNodeSetsGlobal;
                int mNumSideSetsGlobal;

                std::vector<int> mGlobalElemBlkIds;
                std::vector<int> mGlobalElemBlkCnts;

                // Communication Map information
                int mExoFileId;
                int mNumInternalNodes;
                int mNumBorderNodes;
                int mNumExternalNodes;
                int mNumInternalElems;
                int mNumBorderElems;
                int mNumNodeCmaps;
                int mNumElemCmaps;
                Matrix<IdMat> mNodeCmapIds;
                Matrix<IdMat> mNodeCmapNodeCnts;
                Matrix<IdMat> mElemCmapIds;
                Matrix<IdMat> mElemCmapElemCnts;
                moris::Cell<Matrix<IdMat>> mNodeCmapNodeIds;
                moris::Cell<Matrix<IdMat>> mNodeCmapProcIds;

                std::vector<int> mGlobalNodesetIds;
                std::vector<int> mNumGlobalNodeCounts;
                std::vector<int> mNumGlobalNodeDfCounts;

                Matrix<IdMat> mElemMapi;
                Matrix<IdMat> mElemMapb;
                Matrix<IdMat> mNodeMapi;
                Matrix<IdMat> mNodeMapb;
                Matrix<IdMat> mNodeMape;
                Matrix<IdMat> mNodeNumMap;
                Matrix<IdMat> mElemNumMap;

                // Node Sets
                std::vector<int>              mNodeSetIds;
                std::vector<int>              mNodeSetNEntries;
                std::vector<int>              mNodeSetDistFactors;
                std::vector<char>             mNodeSetNamesMemory;
                std::vector<char*>            mNodeSetNamePtrs;
                std::vector<Matrix<IndexMat>> mNodeSetNodeIds;

                // Side Sets
                std::vector<int>              mSideSetIds;
                std::vector<int>              mSideSetNEntries;
                std::vector<int>              mSideSetDistFactors;
                std::vector<char>             SideSetNamesMemory;
                std::vector<char*>            mSideSetNamePtrs;
                std::vector<Matrix<IndexMat>> mSideSetElemIds;
                std::vector<Matrix<IndexMat>> mSideSetSideOrd;

                // Block sets
                std::vector<int>   mBlockIds;
                std::vector<int>   mBlockSetNEntries;
                std::vector<int>   mBlockSetNNodesPerEntry;
                std::vector<int>   mBlockSetNedgesPerEntry;
                std::vector<int>   mBlockSetNfacesPerEntry;
                std::vector<int>   mBlockSetNattrPerEntry;
                std::vector<char>  mBlockNamesMemory;
                std::vector<char*> mBlockNamesPtrs;
                std::vector<char>  mBlockElemTypeNamesMemory;
                std::vector<char*> mBlockElemTypeNamesPtrs;
                std::vector<Matrix<IndexMat>> mBlockSetNodeConn;
                std::vector<Matrix<IndexMat>> mBlockSetEdgeConn;
                std::vector<Matrix<IndexMat>> mBlockSetFaceConn;

                // Nodal Fields
                int                         mNumNodalVars;
                std::vector<char>           mNodeFieldNamesMemory;
                std::vector<char*>          mNodeFieldNamePtrs;
                std::vector<Matrix<DDRMat>> mFieldsNodalVars;

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
                copy_coordinates(int tNewExoFileId);

                void
                copy_node_sets(int aNewExoFileId);

                void
                copy_side_sets(int aNewExoFileId);


                void
                copy_block_sets(int aNewExoFileId);

                void
                copy_nodal_fields(
                        int              aNewExoFileId,
                        ex_init_params & init_params);

                /*
                 * from seacas test rd_wt_mesh.c
                 *
                 */
                void get_file_name(
                        const char *base,
                        const char *other,
                        char       *output);

                static void setup_names(
                        int nnames, std::vector<char>& storage,
                        std::vector<char*>           & ptrs);

            public:

                Exodus_IO_Helper(
                        const char * aExodusFile,
                        const bool   aVerbose = false);

                ~Exodus_IO_Helper();

                /*
                 * Create a new exodus file and the information for an element communication map
                 * Copy the exodus file in this Exodus_IO_Helper to a new one with the
                 * communication map. This is because once the exodus file has been setup,
                 * the element communication map cannot be appended
                 *
                 */
                void
                create_new_exo_with_elem_cmaps_from_existing_exo(
                        std::string    & aFileName,
                        Matrix<IdMat>  & aElementIds,
                        Matrix<IdMat>  & aElementSideOrds,
                        Matrix<IdMat>  & aSharedProcIds);
        };
    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_EXODUS_IO_HELPER_HPP_ */
