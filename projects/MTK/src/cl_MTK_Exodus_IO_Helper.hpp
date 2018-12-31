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
#include "fn_print.hpp"
#include "cl_Communication_Tools.hpp"


namespace moris
{
namespace mtk
{
  class Exodus_IO_Helper
  {
  public:

      Exodus_IO_Helper(const char * aExodusFile)
      {

          mTitle = new char[MAX_LINE_LENGTH+1];
          int   io_ws = 0;
          int   cpu_ws = 0;
          float version;

          get_file_name(aExodusFile, NULL, mTitle);

          mExoFileId = ex_open(mTitle,
                               EX_WRITE,
                               &cpu_ws,
                               &io_ws,
                               &version);

          MORIS_ASSERT(mExoFileId!=-1,"Exo open failure");


          mVerbose = false;

          get_init_mesh_data();
          get_init_global();
          get_eb_info_global();
          get_ns_param_global();
          get_load_bal_parameters();
          get_cmap_params();
          get_node_map();
          get_node_cmap();
          get_node_id_map();
          get_node_coords();
          get_elem_id_map();
          get_set_information();
          get_nodal_fields();
      }


      ~Exodus_IO_Helper()
      {
          delete mTitle;
      }

      /*
       * Create a new exodus file and the information for an element communication map
       * Copy the exodus file in this Exodus_IO_Helper to a new one with the
       * communication map. This is because once the exodus file has been setup,
       * the element communication map cannot be appended
       *
       */
      void
      create_new_exo_with_elem_cmaps_from_existing_exo(std::string    & aFileName,
                                                       Matrix<IdMat>  & aElementIds,
                                                       Matrix<IdMat>  & aElementSideOrds,
                                                       Matrix<IdMat>  & aSharedProcIds)
      {
          // retrieve initialization parameters from existing exodus file
          ex_init_params init_params;
          ex_get_init_ext(mExoFileId, &init_params);

          if (mVerbose)
          {
            std::cout << " init params:\n";
            std::cout << " Exodus ID: "     << mExoFileId << '\n';
            std::cout << " Title: "         << init_params.title << '\n';
            std::cout << " num_dim: "       << init_params.num_dim << '\n';
            std::cout << " num_nodes: "     << init_params.num_nodes << '\n';
            std::cout << " num_elem: "      << init_params.num_elem << '\n';
            std::cout << " num_elem_blk: "  << init_params.num_elem_blk << '\n';
            std::cout << " num_node_sets: " << init_params.num_node_sets << '\n';
            std::cout << " num_side_sets: " << init_params.num_side_sets << '\n';
          }

          // Word sizes
          int   io_ws  = 0;
          int   cpu_ws = 0;

          // Create a new exodus
          char * tNewTitle = new char[MAX_LINE_LENGTH+1];

          // Create the file name
          get_file_name(aFileName.c_str(), NULL, tNewTitle);

          // Create a new exodus file (clobber if already there)
          int tNewExoFileId = ex_create(tNewTitle,
                                        EX_CLOBBER,
                                        &cpu_ws,
                                        &io_ws);

          // Put file initialization information
          mErrFlag = ex_put_init_info(tNewExoFileId,
                                      par_size(),
                                      1,
                                      const_cast<char *>("p"));

          // Put global information on the file
          mErrFlag = ex_put_init_global( tNewExoFileId,
                              mNumNodesGlobal,
                              mNumElemsGlobal,
                              mNumElemBlksGlobal, /* I.   */
                              mNumNodeSetsGlobal, /* II.  */
                              mNumSideSetsGlobal );

          // Initialize information about the mesh (local information)
          mErrFlag = ex_put_init(tNewExoFileId,
                                 init_params.title,
                                 init_params.num_dim,
                                 init_params.num_nodes,
                                 init_params.num_elem,
                                 init_params.num_elem_blk,
                                 init_params.num_node_sets,
                                 init_params.num_side_sets);
          MORIS_ASSERT(!mErrFlag,"ex_put_init failed");



          // Put information about global blk ids and blk counts
          mErrFlag = ex_put_eb_info_global(tNewExoFileId,
                                           this->mGlobalElemBlkIds.data(),
                                           this->mGlobalElemBlkCnts.data());
          MORIS_ASSERT(!mErrFlag,"ex_put_eb_info_global failed");

          // Put the global node sets on the new file
          if (mGlobalNodesetIds.size())
          {
          mErrFlag = ex_put_ns_param_global(tNewExoFileId,
                                            this->mGlobalNodesetIds.data(),
                                            this->mNumGlobalNodeCounts.data(),
                                            this->mNumGlobalNodeDfCounts.data());
          MORIS_ASSERT(!mErrFlag,"ex_put_ns_param_global failed");
          }

          // TODO: SIDE SETS and BLOCKSET copying


          // Specify new load balance parameters but add the num_elem_cmaps of 1
          mNumElemCmaps = 1;
          mErrFlag = ex_put_loadbal_param( tNewExoFileId,
                                           this->mNumInternalNodes,
                                           this->mNumBorderNodes,
                                           this->mNumExternalNodes,
                                           this->mNumInternalElems,
                                           this->mNumBorderElems,
                                           this->mNumNodeCmaps,
                                           this->mNumElemCmaps,
                                           par_rank());

          MORIS_ASSERT(!mErrFlag,"ex_put_loadbal_param failed");

          mElemCmapIds.resize(1,1);
          mElemCmapIds(0,0) = 2;

          mElemCmapElemCnts.resize(1,1);
          mElemCmapElemCnts(0,0) = aElementIds.numel();

          // Put communication map parameters
          mErrFlag = ex_put_cmap_params(tNewExoFileId,
                                        mNodeCmapIds.data(),
                                        mNodeCmapNodeCnts.data(),
                                        mElemCmapIds.data(),
                                        mElemCmapElemCnts.data(),
                                        par_rank());
          // Put node communication map
          mErrFlag = ex_put_node_cmap(tNewExoFileId,
                                      mNodeCmapIds(0),
                                      mNodeCmapNodeIds(0).data(),
                                      mNodeCmapProcIds(0).data(),
                                      par_rank());

          // put node map
          mErrFlag =ne_put_node_map(tNewExoFileId,
                                    mNodeMapi.data(),
                                    mNodeMapb.data(),
                                    mNodeMape.data(),
                                      par_rank());

          // put node id map
          mErrFlag = ex_put_id_map(tNewExoFileId,EX_NODE_MAP, mNodeNumMap.data());

          // put the element communication maps
          mErrFlag = ex_put_elem_cmap(tNewExoFileId,
                                      2,
                                      aElementIds.data(),
                                      aElementSideOrds.data(),
                                      aSharedProcIds.data(),
                                      par_rank());

          // put the element id map
          mErrFlag = ex_put_id_map(tNewExoFileId,EX_ELEM_MAP, mElemNumMap.data());

          // Coordinates
          copy_coordinates(tNewExoFileId);

          // Copy sets
          copy_node_sets( tNewExoFileId );
          copy_side_sets( tNewExoFileId );
          copy_block_sets(tNewExoFileId );

          // Copy fields
          copy_nodal_fields(tNewExoFileId,  init_params);

          ex_close(tNewExoFileId);

      }



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

      void get_init_mesh_data()
      {

          mErrFlag = ex_get_init(mExoFileId,
                                 mTitle,
                                 &mNumDim,
                                 &mNumNodes,
                                 &mNumElem,
                                 &mNumElemBlk,
                                 &mNumNodeSets,
                                 &mNumSideSets);

//        MORIS_ERROR(mErrflag, "Error reading initial global data!");

        if (mVerbose)
          {
            std::cout  << "[" << par_rank() << "] " << "num_dim=" << mNumDim << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_nodes=" << mNumNodes << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_elem=" << mNumElem << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_elem_blk=" << mNumElemBlk << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_node_sets=" << mNumNodeSets << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_side_sets=" << mNumSideSets << std::endl;
          }
      }


      void
      get_load_bal_parameters()
      {
          mErrFlag =ex_get_loadbal_param(mExoFileId,
                                         &mNumInternalNodes,
                                         &mNumBorderNodes,
                                         &mNumExternalNodes,
                                         &mNumInternalElems,
                                         &mNumBorderElems,
                                         &mNumNodeCmaps,
                                         &mNumElemCmaps,
                                         par_rank() // The ID of the processor for which info is to be read
                                         );
          if (mErrFlag)
              printf("after ex_get_init, error = %d\n", mErrFlag);

         if (mVerbose)
           {
             std::cout << "[" << par_rank() << "] " << "num_internal_nodes=" << mNumInternalNodes << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_border_nodes=" << mNumBorderNodes << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_external_nodes=" << mNumExternalNodes << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_internal_elems=" << mNumInternalElems << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_border_elems=" << mNumBorderElems << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_node_cmaps=" << mNumNodeCmaps << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_elem_cmaps=" << mNumElemCmaps << std::endl;
           }
      }

      void
      get_cmap_params()
      {
          // Allocate space based on information from load balance parameter calls
           mNodeCmapIds.resize(1,mNumNodeCmaps);
           mNodeCmapNodeCnts.resize(1,mNumNodeCmaps);
           mElemCmapIds.resize(1,mNumElemCmaps);
           mElemCmapElemCnts.resize(1,mNumElemCmaps);

           // get the cmap parameters
           mErrFlag = ex_get_cmap_params( mExoFileId,
                                          mNodeCmapIds.data(),
                                          mNodeCmapNodeCnts.data(),
                                          mElemCmapIds.data(),
                                          mElemCmapElemCnts.data(),
                                          par_rank());
           MORIS_ERROR(!mErrFlag, "Error reading cmap parameters!");


           if (mVerbose)
             {
//               print(node_cmap_ids,"node_cmap_ids");
//               print(node_cmap_node_cnts,"node_cmap_node_cnts");
//               print(elem_cmap_ids,"elem_cmap_ids");
//               print(elem_cmap_elem_cnts,"elem_cmap_elem_cnts");
             }
      }


      void
      get_node_cmap()
      {
          mNodeCmapNodeIds.resize(mNumNodeCmaps);
          mNodeCmapProcIds.resize(mNumNodeCmaps);

           for (unsigned int i=0; i<mNodeCmapNodeIds.size(); ++i)
             {
               mNodeCmapNodeIds(i).resize(1,mNodeCmapNodeCnts(i));
               mNodeCmapProcIds(i).resize(1,mNodeCmapNodeCnts(i));

               mErrFlag = ex_get_node_cmap(mExoFileId,
                                           mNodeCmapIds(i),
                                           mNodeCmapNodeIds(i).data(),
                                           mNodeCmapProcIds(i).data(),
                                           par_rank());
               MORIS_ERROR(!mErrFlag, "Error reading node cmap node and processor ids!");

               if (mVerbose)
                 {
                   for (unsigned int j=0; j<mNodeCmapNodeIds.size(); ++j)
                       print(mNodeCmapNodeIds(j),"node_cmap_node_ids(j)");

                   // This is basically a vector, all entries of which are = node_cmap_ids[i]
                   // Not sure if it's always guaranteed to be that or what...
                   for (unsigned int j=0; j<mNodeCmapProcIds.size(); ++j)
                       print(mNodeCmapNodeIds(j),"node_cmap_node_ids(j)");
                 }
             }
      }

      void get_init_global()
      {
        mErrFlag =ex_get_init_global(mExoFileId,
                                      &mNumNodesGlobal,
                                      &mNumElemsGlobal,
                                      &mNumElemBlksGlobal,
                                      &mNumNodeSetsGlobal,
                                      &mNumSideSetsGlobal);
        MORIS_ERROR(!mErrFlag, "Error reading initial global data!");

        if (mVerbose)
          {
            std::cout << "[" << par_rank() << "] " << "num_nodes_global=" << mNumNodesGlobal << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_elems_global=" << mNumElemsGlobal << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_elem_blks_global=" << mNumElemBlksGlobal << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_node_sets_global=" << mNumNodeSetsGlobal << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_side_sets_global=" << mNumSideSetsGlobal << std::endl;
          }
      }


      void get_eb_info_global()
      {
        mGlobalElemBlkIds.resize(mNumElemBlksGlobal);
        mGlobalElemBlkCnts.resize(mNumElemBlksGlobal);

        if (mNumElemBlksGlobal > 0)
          {
            mErrFlag = ex_get_eb_info_global(mExoFileId,
                                             mGlobalElemBlkIds.data(),
                                             mGlobalElemBlkCnts.data());
            MORIS_ERROR(!mErrFlag, "Error reading global element block info!");
          }

        if (mVerbose)
          {
            std::cout << "[" << par_rank() << "] " << "Global Element Block IDs and Counts:" << std::endl;
            for (std::size_t bn=0; bn<mGlobalElemBlkIds.size(); ++bn)
              {
                std::cout << "  [" << par_rank() << "] "
                             << "global_elem_blk_ids["<<bn<<"]=" << mGlobalElemBlkIds[bn]
                             << ", global_elem_blk_cnts["<<bn<<"]=" << mGlobalElemBlkCnts[bn]
                             << std::endl;
              }
          }
      }

      void
      get_node_coords()
      {
          mX.resize(mNumNodes,1);
          mY.resize(mNumNodes,1);
          mZ.resize(mNumNodes,1);

          mErrFlag = ex_get_coord(mExoFileId,
                                      mX.data(),
                                      mY.data(),
                                      mZ.data());
      }

      void get_ns_param_global()
      {
        if (mNumNodeSetsGlobal > 0)
          {
            mGlobalNodesetIds.resize(mNumNodeSetsGlobal);
            mNumGlobalNodeCounts.resize(mNumNodeSetsGlobal);
            mNumGlobalNodeDfCounts.resize(mNumNodeSetsGlobal);

            mErrFlag =
              ex_get_ns_param_global(mExoFileId,
                                              mGlobalNodesetIds.data(),
                                              mNumGlobalNodeCounts.data(),
                                              mNumGlobalNodeDfCounts.data());
            MORIS_ERROR(!mErrFlag, "Error reading global nodeset parameters!");

            if (mVerbose)
              {
                std::cout << "[" << par_rank() << "] " << "Global Nodeset IDs, Node Counts, and DF counts:" << std::endl;
                for (std::size_t bn=0; bn<mGlobalNodesetIds.size(); ++bn)
                  {
                    std::cout << "  [" << par_rank() << "] "
                                 << "global_nodeset_ids["<<bn<<"]=" << mGlobalNodesetIds[bn]
                                 << ", num_global_node_counts["<<bn<<"]=" << mNumGlobalNodeCounts[bn]
                                 << ", num_global_node_df_counts["<<bn<<"]=" << mNumGlobalNodeDfCounts[bn]
                                 << std::endl;
                  }
              }
          }
      }

      void
      get_node_map()
      {
          mNodeMapi.resize(mNumInternalNodes,1);
          mNodeMapb.resize(mNumBorderNodes,1);
          mNodeMape.resize(mNumExternalNodes,1);
          mErrFlag =ne_get_node_map(mExoFileId,
                                     mNodeMapi.data(),
                                     mNodeMapb.data(),
                                     mNodeMape.data(),
                                      par_rank());
          MORIS_ERROR(!mErrFlag, "ne_get_node_map failed!");
      }


      void
      get_node_id_map()
      {
          mNodeNumMap.resize(1,mNumNodes);
          mErrFlag = ex_get_id_map(mExoFileId,EX_NODE_MAP, mNodeNumMap.data());
      }

      void
      get_elem_id_map()
      {
          mElemNumMap.resize(1,mNumElem);
          mErrFlag = ex_get_id_map(mExoFileId,EX_ELEM_MAP, mElemNumMap.data());
      }

      void
      get_set_information()
      {
          // Node set information
          mNodeSetIds.resize(mNumNodeSets);
          ex_get_ids(mExoFileId, EX_NODE_SET, mNodeSetIds.data());

          setup_names(int(mNumNodeSets), mNodeSetNamesMemory, mNodeSetNamePtrs);
          ex_get_names(mExoFileId, EX_NODE_SET, mNodeSetNamePtrs.data());

          mNodeSetNEntries.resize(mNumNodeSets);
          mNodeSetNodeIds.resize(mNumNodeSets);
          mNodeSetDistFactors.resize(mNumNodeSets);
          for (size_t i = 0; i < mNodeSetIds.size(); ++i)
          {
              ex_get_set_param(mExoFileId, EX_NODE_SET, mNodeSetIds[i], &mNodeSetNEntries[i], &mNodeSetDistFactors[i]);
              if (mVerbose) {
                  std::cout << "node set " << mNodeSetIds[i] << " has " << mNodeSetNEntries[i]<< " nodes\n";
              }

              mNodeSetNodeIds[i] = Matrix<IndexMat>(1,mNodeSetNEntries[i]);
              ex_get_set(mExoFileId, EX_NODE_SET, mNodeSetIds[i],
                         mNodeSetNodeIds[i].data(), nullptr);

          }

          // Copy over the side sets
          mSideSetIds.resize(mNumSideSets);
          ex_get_ids(mExoFileId, EX_SIDE_SET, mSideSetIds.data());

          // get/put side set names
          setup_names(int(mNumSideSets), SideSetNamesMemory, mSideSetNamePtrs);
          ex_get_names(mExoFileId, EX_SIDE_SET, mSideSetNamePtrs.data());

          // Put the side set data into the new file

          mSideSetNEntries.resize(mNumSideSets);
          mSideSetElemIds.resize(mNumSideSets);
          mSideSetSideOrd.resize(mNumSideSets);
          mSideSetDistFactors.resize(mNumSideSets);

          for (size_t i = 0; i < mSideSetIds.size(); ++i)
          {
              ex_get_set_param( mExoFileId, EX_SIDE_SET, mSideSetIds[i], &mSideSetNEntries[i], &mSideSetDistFactors[i] );
              if (mVerbose) {
                  std::cout << "side set #" << mSideSetIds[i] << " \"" << mSideSetNamePtrs[i]
                                                                                     << "\" has " << mSideSetNEntries[i] << " sides, will be surface "
                                                                                     << mSideSetIds[i] << "\n";
              }

              mSideSetElemIds[i] = Matrix<IndexMat> (1,mSideSetNEntries[i]);
              mSideSetSideOrd[i] = Matrix<IndexMat> (1,mSideSetNEntries[i]);
              ex_get_set(mExoFileId, EX_SIDE_SET, mSideSetIds[i],
                         mSideSetElemIds[i].data(), mSideSetSideOrd[i].data());
          }


          mBlockIds.resize(mNumElemBlk);
          ex_get_ids(mExoFileId, EX_ELEM_BLOCK, mBlockIds.data());

          setup_names(int(mNumElemBlk), mBlockNamesMemory, mBlockNamesPtrs);
          ex_get_names(mExoFileId, EX_ELEM_BLOCK, mBlockNamesPtrs.data());

          setup_names(int(mNumElemBlk), mBlockElemTypeNamesMemory, mBlockElemTypeNamesPtrs);

          mBlockSetNEntries.resize(mNumElemBlk);
          mBlockSetNNodesPerEntry.resize(mNumElemBlk);
          mBlockSetNedgesPerEntry.resize(mNumElemBlk);
          mBlockSetNfacesPerEntry.resize(mNumElemBlk);
          mBlockSetNattrPerEntry.resize(mNumElemBlk);
          mBlockSetNodeConn.resize(mNumElemBlk);
          mBlockSetEdgeConn.resize(mNumElemBlk);
          mBlockSetFaceConn.resize(mNumElemBlk);
          for (size_t i = 0; i < mBlockIds.size(); ++i)
          {
              ex_get_block(mExoFileId,
                           EX_ELEM_BLOCK,
                           mBlockIds[i],
                           mBlockElemTypeNamesPtrs[i],
                           &mBlockSetNEntries[i],
                           &mBlockSetNNodesPerEntry[i],
                           &mBlockSetNedgesPerEntry[i],
                           &mBlockSetNfacesPerEntry[i],
                           &mBlockSetNattrPerEntry[i]);

              if(mVerbose)
              {
              std::cout << "block " << mBlockIds[i] << " \"" << mBlockNamesPtrs[i] << "\""
                      << " has " << mBlockSetNEntries[i] << " elements of type " << mBlockElemTypeNamesPtrs[i]
                      << '\n';
              }
              mBlockSetNodeConn[i] = Matrix<IndexMat>(1,mBlockSetNEntries[i]*mBlockSetNNodesPerEntry[i]);
              mBlockSetEdgeConn[i] = Matrix<IndexMat>(1,mBlockSetNEntries[i]*mBlockSetNedgesPerEntry[i]);
              mBlockSetFaceConn[i] = Matrix<IndexMat>(1,mBlockSetNEntries[i]*mBlockSetNfacesPerEntry[i]);
              ex_get_conn(mExoFileId,
                          EX_ELEM_BLOCK,
                          mBlockIds[i],
                          mBlockSetNodeConn[i].data(),
                          mBlockSetEdgeConn[i].data(),
                          mBlockSetFaceConn[i].data());


              MORIS_ASSERT(!mErrFlag,"ex_put_conn failed");
          }


      }


      void
      get_nodal_fields()
      {
          int time_step = 0;
          ex_get_variable_param(mExoFileId, EX_NODAL, &mNumNodalVars);
          if(mVerbose) std::cout << mNumNodalVars << " nodal variables (in get_nodal_fields)\n";
          if(mNumNodalVars == 0) return ;

          setup_names(mNumNodalVars, mNodeFieldNamesMemory, mNodeFieldNamePtrs);
          ex_get_variable_names(mExoFileId, EX_NODAL, mNumNodalVars, mNodeFieldNamePtrs.data());

		mFieldsNodalVars.resize(mNumNodalVars);
          for (int i = 0; i < mNumNodalVars; ++i) {
            auto name = mNodeFieldNamePtrs[std::size_t(i)];
            if (mVerbose)
            {
              std::cout << "Loading nodal variable \"" << name << "\" at time step "
                        << time_step << '\n';
            }
            mFieldsNodalVars[i] = Matrix<DDRMat>(mNumNodes, 1);
            ex_get_var(mExoFileId, time_step + 1, EX_NODAL, i + 1, /*obj_id*/ 0,
                       mNumNodes, mFieldsNodalVars[i].data());

        }
      }


      void
      copy_coordinates(int tNewExoFileId)
      {
          const char *tmp[3] = {"x", "y", "z"};
          char *coord_names[3];
          memcpy(coord_names, tmp, sizeof coord_names);
          mErrFlag = ex_put_coord_names(tNewExoFileId, coord_names);

          mErrFlag = ex_put_coord(tNewExoFileId,
                                  mX.data(),
                                  mY.data(),
                                  mZ.data());
      }

      void
      copy_node_sets(int              aNewExoFileId)
      {

          // Put the names onto the new file
          ex_put_names(aNewExoFileId,EX_NODE_SET, mNodeSetNamePtrs.data());

          // Put the side set data into the new file
          for (size_t i = 0; i < mNodeSetIds.size(); ++i) {
              ex_put_set_param(aNewExoFileId, EX_NODE_SET, mNodeSetIds[i], mNodeSetNEntries[i], mNodeSetDistFactors[i]);
              if (par_rank() == 0) {
                  std::cout << "node set " << mNodeSetIds[i] << " has " << mNodeSetNEntries[i]
                          << " nodes\n";
              }

              ex_put_set(aNewExoFileId, EX_NODE_SET, mNodeSetIds[i],
                         mNodeSetNodeIds[i].data(), nullptr);
          }
      }

      void
      copy_side_sets(int              aNewExoFileId)
      {
          // Copy over the side sets
          ex_put_names(aNewExoFileId,EX_SIDE_SET, mSideSetNamePtrs.data());

          // Put the side set data into the new file
          for (size_t i = 0; i < mSideSetIds.size(); ++i) {
              ex_put_set_param(
                      aNewExoFileId, EX_SIDE_SET, mSideSetIds[i], mSideSetNEntries[i], mSideSetDistFactors[i]);
              if (mVerbose) {
                  std::cout << "side set #" << mSideSetIds[i] << " \"" << mSideSetNamePtrs[i]
                                                                                     << "\" has " << mSideSetNEntries[i] << " sides, will be surface "
                                                                                     << mSideSetIds[i] << "\n";
              }

              ex_put_set(aNewExoFileId, EX_SIDE_SET, mSideSetIds[i],
                         mSideSetElemIds[i].data(), mSideSetSideOrd[i].data());

          }
      }


      void
      copy_block_sets(int              aNewExoFileId)
      {
          ex_put_names(aNewExoFileId, EX_ELEM_BLOCK, mBlockNamesPtrs.data());

          for (size_t i = 0; i < mBlockIds.size(); ++i)
          {


              mErrFlag = ex_put_block(aNewExoFileId,
                                      EX_ELEM_BLOCK,
                                      mBlockIds[i],
                                      mBlockElemTypeNamesPtrs[i],
                                      mBlockSetNEntries[i],
                                      mBlockSetNNodesPerEntry[i],
                                      mBlockSetNedgesPerEntry[i],
                                      mBlockSetNfacesPerEntry[i],
                                      mBlockSetNattrPerEntry[i]);

              MORIS_ASSERT(!mErrFlag,"ex_put_block failed");
              mErrFlag = ex_put_conn(aNewExoFileId,
                                     EX_ELEM_BLOCK,
                                     mBlockIds[i],
                                     mBlockSetNodeConn[i].data(),
                                     mBlockSetEdgeConn[i].data(),
                                     mBlockSetFaceConn[i].data());

          }
      }

      void
      copy_nodal_fields(int              aNewExoFileId,
                        ex_init_params & init_params)
      {
          if(mVerbose) std::cout << mNumNodalVars << " nodal variables\n";
          if(mNumNodalVars == 0) return ;
          int time_step = 0;
          ex_put_variable_param(aNewExoFileId, EX_NODAL, mNumNodalVars);
          ex_put_variable_names(aNewExoFileId, EX_NODAL, mNumNodalVars, mNodeFieldNamePtrs.data());


          for (int i = 0; i < mNumNodalVars; ++i) {
            auto name = mNodeFieldNamePtrs[std::size_t(i)];
            if (mVerbose)
            {
              std::cout << "Loading nodal variable \"" << name << "\" at time step "
                        << time_step << '\n';
            }
            ex_put_var(aNewExoFileId, time_step + 1, EX_NODAL, i + 1, /*obj_id*/ 0,
                       mNumNodes, mFieldsNodalVars[i].data());

        }

      }

      /*
       * from seacas test rd_wt_mesh.c
       *
       */
      void get_file_name(const char *base,
                         const char *other,
                         char *output)
      {
        int nprocs = par_size();
        int rank = par_rank();
        int  i1, iTemp1;
        int  iMaxDigit = 0, iMyDigit = 0;
        char cTemp[128];

        output[0] = '\0';
        strcpy(output, base);
        if (other != NULL) {
          strcat(output, ".");
          strcat(output, other);
        }

        if (nprocs > 1) {
          /*
           * Find out the number of digits needed to specify the processor ID.
           * This allows numbers like 01-99, i.e., prepending zeros to the
           * name to preserve proper alphabetic sorting of the files.
           */

          iTemp1 = nprocs;
          do {
            iTemp1 /= 10;
            iMaxDigit++;
          } while (iTemp1 >= 1);

          iTemp1 = rank;
          do {
            iTemp1 /= 10;
            iMyDigit++;
          } while (iTemp1 >= 1);

          strcat(output, ".");
          sprintf(cTemp, "%d", nprocs);
          strcat(output, cTemp);
          strcat(output, ".");

          /*
           * Append the proper number of zeros to the filename.
           */
          for (i1 = 0; i1 < iMaxDigit - iMyDigit; i1++) {
            strcat(output, "0");
          }

          sprintf(cTemp, "%d", rank);
          strcat(output, cTemp);
      }
        std::cout<<"output = "<<output<<std::endl;
      }

      static void setup_names(
          int nnames, std::vector<char>& storage, std::vector<char*>& ptrs) {
        constexpr auto max_name_length = MAX_STR_LENGTH + 1;
        storage = std::vector<char>(std::size_t(nnames * max_name_length), '\0');
        ptrs = std::vector<char*>(std::size_t(nnames), nullptr);
        for (int i = 0; i < nnames; ++i) {
          ptrs[std::size_t(i)] = storage.data() + max_name_length * i;
        }
      }
  };
}
}




#endif /* PROJECTS_MTK_SRC_CL_MTK_EXODUS_IO_HELPER_HPP_ */
