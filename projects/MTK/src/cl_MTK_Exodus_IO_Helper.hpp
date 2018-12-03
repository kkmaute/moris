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

          title = new char[MAX_LINE_LENGTH+1];
          int   io_ws = 0;
          int   cpu_ws = 0;
          float version;

          get_file_name(aExodusFile, NULL, title);

          mExoFileId = ex_open(title,
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
      }


      ~Exodus_IO_Helper()
      {
          delete title;
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

          int   io_ws = 0;
          int   cpu_ws = 0;

          // Create a new exodus
          char * tNewTitle = new char[MAX_LINE_LENGTH+1];

          get_file_name(aFileName.c_str(), NULL, tNewTitle);

          // Create a new exodus file
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
                              num_nodes_global,
                              num_elems_global,
                              num_elem_blks_global, /* I.   */
                              num_node_sets_global, /* II.  */
                              num_side_sets_global );

          // Initialize information about the mesh (local information)
          mErrFlag = ex_put_init(tNewExoFileId,
                                 title,
                                 num_dim,
                                 num_nodes,
                                 num_elem,
                                 num_elem_blk,
                                 num_node_sets,
                                 num_side_sets);
          MORIS_ASSERT(!mErrFlag,"ex_put_init failed");



          // Put information about global blk ids and blk counts
          mErrFlag = ex_put_eb_info_global(tNewExoFileId,
                                           this->global_elem_blk_ids.data(),
                                           this->global_elem_blk_cnts.data());
          MORIS_ASSERT(!mErrFlag,"ex_put_eb_info_global failed");

          // Put the global node sets on the new file
          if (global_nodeset_ids.size())
          {
          mErrFlag = ex_put_ns_param_global(tNewExoFileId,
                                            this->global_nodeset_ids.data(),
                                            this->num_global_node_counts.data(),
                                            this->num_global_node_df_counts.data());
          MORIS_ASSERT(!mErrFlag,"ex_put_ns_param_global failed");
          }

          // TODO: SIDE SETS and BLOCKSET copying


          // Specify new load balance parameters but add the num_elem_cmaps of 1
          num_elem_cmaps = 1;
          mErrFlag = ex_put_loadbal_param( tNewExoFileId,
                                           this->num_internal_nodes,
                                           this->num_border_nodes,
                                           this->num_external_nodes,
                                           this->num_internal_elems,
                                           this->num_border_elems,
                                           this->num_node_cmaps,
                                           this->num_elem_cmaps,
                                           par_rank());

          MORIS_ASSERT(!mErrFlag,"ex_put_loadbal_param failed");

          elem_cmap_ids.resize(1,1);
          elem_cmap_ids(0,0) = 2;

          elem_cmap_elem_cnts.resize(1,1);
          elem_cmap_elem_cnts(0,0) = aElementIds.numel();


          mErrFlag = ex_put_cmap_params(tNewExoFileId,
                                        node_cmap_ids.data(),
                                        node_cmap_node_cnts.data(),
                                        elem_cmap_ids.data(),
                                        elem_cmap_elem_cnts.data(),
                                        par_rank());
          // Put node communication map
          mErrFlag = ex_put_node_cmap(tNewExoFileId,
                                      node_cmap_ids(0),
                                      node_cmap_node_ids(0).data(),
                                      node_cmap_proc_ids(0).data(),
                                      par_rank());
          node_mapi.resize(num_internal_nodes,1);
          node_mapb.resize(num_border_nodes,1);
          node_mape.resize(num_external_nodes,1);

          // put node map
          mErrFlag =ne_put_node_map(tNewExoFileId,
                                    node_mapi.data(),
                                    node_mapb.data(),
                                    node_mape.data(),
                                      par_rank());

          // put node id map
          mErrFlag = ex_put_id_map(tNewExoFileId,EX_NODE_MAP, node_num_map.data());

          // put the element communication maps
          mErrFlag = ex_put_elem_cmap(tNewExoFileId,
                                      2,
                                      aElementIds.data(),
                                      aElementSideOrds.data(),
                                      aSharedProcIds.data(),
                                      par_rank());

          // put the element id map
          mErrFlag = ex_put_id_map(tNewExoFileId,EX_ELEM_MAP, elem_num_map.data());

          //TODO: add coordinate names
          mErrFlag = ex_put_coord(tNewExoFileId,
                                      x.data(),
                                      y.data(),
                                      z.data());


          //FIXME: Add element connectivity block
          char * tStr = new char[MAX_LINE_LENGTH+1];
          for (unsigned int i=0; i<static_cast<unsigned>(this->num_elem_blks_global); ++i)
          {

              ex_block block{};
              mErrFlag = ex_get_block_param(mExoFileId, &block);
              std::cout<<"global_elem_blk_cnts[i]"<<global_elem_blk_cnts[i]<<std::endl;
              block.id   = 0;
              std::cout<<"Block id = "<<block.id;
              block.type = EX_ELEM_BLOCK;
              std::cout<<"block.type = "<<block.type;

              std::cout<<"Num entitiy = "<<block.num_entry <<std::endl;
              std::cout<<"num_nodes_per_elmt = "<<block.num_nodes_per_entry<<std::endl;

              Matrix<IndexMat> conn(1,block.num_entry*block.num_nodes_per_entry);
              mErrFlag = ex_get_conn(mExoFileId, EX_ELEM_BLOCK, block.id, conn.data(), nullptr,nullptr);
              block.id = 0;
              mErrFlag = ex_put_conn(tNewExoFileId, EX_ELEM_BLOCK,block.id, conn.data(), nullptr,nullptr);

          }

          delete tStr;
          ex_close(tNewExoFileId);


//          MORIS_ASSERT(mErrFlag,"ex_put_elem_cmap failed");

      }



  private:

      int  mErrFlag;
      bool mVerbose;

      // general mesh info
      int num_dim;
      int num_nodes;
      int num_elem;
      int num_elem_blk;
      int num_node_sets;
      int num_side_sets;
      char * title;

      // Coordinates
      Matrix<DDRMat> x;
      Matrix<DDRMat> y;
      Matrix<DDRMat> z;


      // GLobal information
      int num_nodes_global;
      int num_elems_global;
      int num_elem_blks_global;
      int num_node_sets_global;
      int num_side_sets_global;
      std::vector<int> global_elem_blk_ids;
      std::vector<int> global_elem_blk_cnts;

      // Communication Map information
      int mExoFileId;
      int num_internal_nodes;
      int num_border_nodes;
      int num_external_nodes;
      int num_internal_elems;
      int num_border_elems;
      int num_node_cmaps;
      int num_elem_cmaps;
      Matrix<IdMat> node_cmap_ids;
      Matrix<IdMat> node_cmap_node_cnts;
      Matrix<IdMat> elem_cmap_ids;
      Matrix<IdMat> elem_cmap_elem_cnts;
      moris::Cell<Matrix<IdMat>> node_cmap_node_ids;
      moris::Cell<Matrix<IdMat>> node_cmap_proc_ids;

      std::vector<int> global_nodeset_ids;
      std::vector<int> num_global_node_counts;
      std::vector<int> num_global_node_df_counts;


      Matrix<IdMat> elem_mapi;
      Matrix<IdMat> elem_mapb;
      Matrix<IdMat> node_mapi;
      Matrix<IdMat> node_mapb;
      Matrix<IdMat> node_mape;
      Matrix<IdMat> node_num_map;
      Matrix<IdMat> elem_num_map;



      void get_init_mesh_data()
      {

          mErrFlag = ex_get_init(mExoFileId,
                                 title,
                                 &num_dim,
                                 &num_nodes,
                                 &num_elem,
                                 &num_elem_blk,
                                 &num_node_sets,
                                 &num_side_sets);

//        MORIS_ERROR(mErrflag, "Error reading initial global data!");

        if (mVerbose)
          {
            std::cout  << "[" << par_rank() << "] " << "num_dim=" << num_dim << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_nodes=" << num_nodes << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_elem=" << num_elem << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_elem_blk=" << num_elem_blk << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_node_sets=" << num_node_sets << std::endl;
            std::cout  << "[" << par_rank() << "] " << "num_side_sets=" << num_side_sets << std::endl;
          }
      }


      void
      get_load_bal_parameters()
      {
          mErrFlag =ex_get_loadbal_param(mExoFileId,
                                         &num_internal_nodes,
                                         &num_border_nodes,
                                         &num_external_nodes,
                                         &num_internal_elems,
                                         &num_border_elems,
                                         &num_node_cmaps,
                                         &num_elem_cmaps,
                                         par_rank() // The ID of the processor for which info is to be read
                                         );
          if (mErrFlag)
              printf("after ex_get_init, error = %d\n", mErrFlag);

         if (mVerbose)
           {
             std::cout << "[" << par_rank() << "] " << "num_internal_nodes=" << num_internal_nodes << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_border_nodes=" << num_border_nodes << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_external_nodes=" << num_external_nodes << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_internal_elems=" << num_internal_elems << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_border_elems=" << num_border_elems << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_node_cmaps=" << num_node_cmaps << std::endl;
             std::cout << "[" << par_rank() << "] " << "num_elem_cmaps=" << num_elem_cmaps << std::endl;
           }
      }

      void
      get_cmap_params()
      {
          // Allocate space based on information from load balance parameter calls
           node_cmap_ids.resize(1,num_node_cmaps);
           node_cmap_node_cnts.resize(1,num_node_cmaps);
           elem_cmap_ids.resize(1,num_elem_cmaps);
           elem_cmap_elem_cnts.resize(1,num_elem_cmaps);

           // get the cmap parameters
           mErrFlag = ex_get_cmap_params( mExoFileId,
                                          node_cmap_ids.data(),
                                          node_cmap_node_cnts.data(),
                                          elem_cmap_ids.data(),
                                          elem_cmap_elem_cnts.data(),
                                          par_rank());
           MORIS_ERROR(!mErrFlag, "Error reading cmap parameters!");


           if (mVerbose && par_rank()==0)
             {
               print(node_cmap_ids,"node_cmap_ids");
               print(node_cmap_node_cnts,"node_cmap_node_cnts");
               print(elem_cmap_ids,"elem_cmap_ids");
               print(elem_cmap_elem_cnts,"elem_cmap_elem_cnts");
             }
      }


      void
      get_node_cmap()
      {
          node_cmap_node_ids.resize(num_node_cmaps);
          node_cmap_proc_ids.resize(num_node_cmaps);

           for (unsigned int i=0; i<node_cmap_node_ids.size(); ++i)
             {
               node_cmap_node_ids(i).resize(1,node_cmap_node_cnts(i));
               node_cmap_proc_ids(i).resize(1,node_cmap_node_cnts(i));

               mErrFlag = ex_get_node_cmap(mExoFileId,
                                           node_cmap_ids(i),
                                           node_cmap_node_ids(i).data(),
                                           node_cmap_proc_ids(i).data(),
                                           par_rank());
               MORIS_ERROR(!mErrFlag, "Error reading node cmap node and processor ids!");

               if (mVerbose && par_rank() == 0)
                 {
                   for (unsigned int j=0; j<node_cmap_node_ids.size(); ++j)
                       print(node_cmap_node_ids(j),"node_cmap_node_ids(j)");

                   // This is basically a vector, all entries of which are = node_cmap_ids[i]
                   // Not sure if it's always guaranteed to be that or what...
                   for (unsigned int j=0; j<node_cmap_proc_ids.size(); ++j)
                       print(node_cmap_node_ids(j),"node_cmap_node_ids(j)");
                 }
             }
      }

      void get_init_global()
      {
        mErrFlag =ex_get_init_global(mExoFileId,
                                      &num_nodes_global,
                                      &num_elems_global,
                                      &num_elem_blks_global,
                                      &num_node_sets_global,
                                      &num_side_sets_global);
        MORIS_ERROR(!mErrFlag, "Error reading initial global data!");

        if (mVerbose)
          {
            std::cout << "[" << par_rank() << "] " << "num_nodes_global=" << num_nodes_global << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_elems_global=" << num_elems_global << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_elem_blks_global=" << num_elem_blks_global << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_node_sets_global=" << num_node_sets_global << std::endl;
            std::cout << "[" << par_rank() << "] " << "num_side_sets_global=" << num_side_sets_global << std::endl;
          }
      }


      void get_eb_info_global()
      {
        global_elem_blk_ids.resize(num_elem_blks_global);
        global_elem_blk_cnts.resize(num_elem_blks_global);

        if (num_elem_blks_global > 0)
          {
            mErrFlag = ex_get_eb_info_global(mExoFileId,
                                             global_elem_blk_ids.data(),
                                             global_elem_blk_cnts.data());
            MORIS_ERROR(!mErrFlag, "Error reading global element block info!");
          }

        if (mVerbose)
          {
            std::cout << "[" << par_rank() << "] " << "Global Element Block IDs and Counts:" << std::endl;
            for (std::size_t bn=0; bn<global_elem_blk_ids.size(); ++bn)
              {
                std::cout << "  [" << par_rank() << "] "
                             << "global_elem_blk_ids["<<bn<<"]=" << global_elem_blk_ids[bn]
                             << ", global_elem_blk_cnts["<<bn<<"]=" << global_elem_blk_cnts[bn]
                             << std::endl;
              }
          }
      }

      void
      get_node_coords()
      {
          x.resize(num_nodes,1);
          y.resize(num_nodes,1);
          z.resize(num_nodes,1);

          mErrFlag = ex_get_coord(mExoFileId,
                                      x.data(),
                                      y.data(),
                                      z.data());
      }

      void get_ns_param_global()
      {
        if (num_node_sets_global > 0)
          {
            global_nodeset_ids.resize(num_node_sets_global);
            num_global_node_counts.resize(num_node_sets_global);
            num_global_node_df_counts.resize(num_node_sets_global);

            mErrFlag =
              ex_get_ns_param_global(mExoFileId,
                                              global_nodeset_ids.data(),
                                              num_global_node_counts.data(),
                                              num_global_node_df_counts.data());
            MORIS_ERROR(!mErrFlag, "Error reading global nodeset parameters!");

            if (mVerbose)
              {
                std::cout << "[" << par_rank() << "] " << "Global Nodeset IDs, Node Counts, and DF counts:" << std::endl;
                for (std::size_t bn=0; bn<global_nodeset_ids.size(); ++bn)
                  {
                    std::cout << "  [" << par_rank() << "] "
                                 << "global_nodeset_ids["<<bn<<"]=" << global_nodeset_ids[bn]
                                 << ", num_global_node_counts["<<bn<<"]=" << num_global_node_counts[bn]
                                 << ", num_global_node_df_counts["<<bn<<"]=" << num_global_node_df_counts[bn]
                                 << std::endl;
                  }
              }
          }
      }

      void
      get_node_map()
      {
          node_mapi.resize(num_internal_nodes,1);
          node_mapb.resize(num_border_nodes,1);
          node_mape.resize(num_external_nodes,1);
          mErrFlag =ne_get_node_map(mExoFileId,
                                     node_mapi.data(),
                                     node_mapb.data(),
                                     node_mape.data(),
                                      par_rank());
          MORIS_ERROR(!mErrFlag, "ne_get_node_map failed!");
      }


      void
      get_node_id_map()
      {
          node_num_map.resize(1,num_nodes);
          mErrFlag = ex_get_id_map(mExoFileId,EX_NODE_MAP, node_num_map.data());
      }

      void
      get_elem_id_map()
      {
          elem_num_map.resize(1,num_elem);
          mErrFlag = ex_get_id_map(mExoFileId,EX_ELEM_MAP, elem_num_map.data());
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
  };
}
}




#endif /* PROJECTS_MTK_SRC_CL_MTK_EXODUS_IO_HELPER_HPP_ */
