/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Mesh_Editor.hpp
 *
 */

#ifndef SRC_cl_MTK_Integration_Mesh_Editor_HPP_
#define SRC_cl_MTK_Integration_Mesh_Editor_HPP_

#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include <unordered_map>
#include "cl_TOL_Memory_Map.hpp"

namespace moris::mtk
{
    class Vertex;
    class Cluster;

    struct Integration_Mesh_Info
    {
        uint mSpatialDim;

        moris::Cell< mtk::Vertex const * > mVertices;

        // cell to vertex connectivity
        moris::Cell< moris_index > mCellToVertexIndicies;

        // cell cluster to primary ig cell connectivity
        moris::Cell< moris_index > mCellClusterToPrimaryIGCellIndices;

        // cell cluster to void ig cell connectivity
        moris::Cell< moris_index > mCellClusterToVoidIGCellIndices;

        // cell cluster to vertices connectivity
        moris::Cell< moris_index > mCellClusterToVertexIndices;

        // side cluster to the ig cells and ordinal connectivity
        moris::Cell< moris_index > mSideClusterToPrimaryIGCellIndices;

        moris::Cell< moris_index > mSideClusterToVoidIGCellIndices;

        moris::Cell< moris_index > mSideClusterToVertexOffSet;

        // side cluster to vertex connectivity
        moris::Cell< moris_index > mSideClusterToVertexIndices;

        // block - set info
        moris::Cell< moris_index > mBlockSetToCellClusterIndex;
        moris::Cell< moris_index > mBlockSetToCellClusterOffSet;

        // Side Set info
        moris::Cell< moris_index > mSideSetToSideClusterIndex;
        moris::Cell< moris_index > mSideSetToSideClusterOffset;

        // // Double Side Set info
        moris::Cell< moris_index > mDobleSideSetToDoubleSidedClusterIndex;
        moris::Cell< moris_index > mDobleSideSetToDoubleSidedClusterOffset;

        // the map used in constrcution of the double sided clusters
        std::unordered_map< Cluster const *, moris_index > mPreviousSideClusterToNewSideCluster;

        // Double sided cluster connectivity
        moris::Cell< moris_index >                           mDoubleSidedClusterToVertexOffSet;
        moris::Cell< std::pair< moris_index, moris_index > > mDoubleSidedClusterToNewSideClusterIndex;
        moris::Cell< moris_index >                           mDoubleSidedClusterToPairedVerticesIndex;

        moris::Cell< moris_index > mGhostToVertexPairIndices;
        moris::Cell< moris_index > mGhostToVertexOffset;

        // Temporary data needs to be fixed
        moris::Cell< moris_index > mGhostMasterToIPCellIndex;
        moris::Cell< moris_index > mGhostSlaveToIPCellIndex;

        moris::Cell< moris_index > mGhostMasterToIGCellIndex;
        moris::Cell< moris_index > mGhostSlaveToIGCellIndex;

        // Ghost side clusters and double sided clusters
        moris::Cell< moris_index > mDoubleSidedGhostToSideClusterGhostOffset;

        moris::Cell< moris_index > mGhostMasterSlaveVertexIndices;
    };
}    // namespace moris::mtk

namespace moris::mtk
{
    class Vertex;
    class Cluster;
    class Integration_Mesh;
    class Interpolation_Mesh_DataBase_IP;
    class Integration_Mesh_DataBase_IG;

    class Integration_Mesh_Editor
    {
      protected:
        moris::mtk::Integration_Mesh* mInputMesh = nullptr; /*!< the IG mesh we will use*/

        moris::mtk::Interpolation_Mesh_DataBase_IP* mIPMeshDataBase = nullptr; /*!< Detailed description after the member */

        moris::mtk::Integration_Mesh_DataBase_IG* mOutputMesh = nullptr;

        Integration_Mesh_Info* mIGMeshInfo = nullptr;

        bool mCheckMesh = true;

        // counters in order to add data later
        uint mNumPreviousVertices       = 0;
        uint mNumPreviousCells          = 0;
        uint mNumPreviousSideCluster    = 0;
        uint mNumPreviousDblSideCluster = 0;
        uint mNumPreviousDoubleSideSet  = 0;

      public:
        // ----------------------------------------------------------------------------

        /**
         * @brief Construct a new Integration_Mesh_Editor
         *
         * @param aMTKMesh
         * @param aIPMeshDataBase
         */

        Integration_Mesh_Editor() = default;

        // ----------------------------------------------------------------------------

        /**
         * @brief Construct a new Integration_Mesh_Editor
         *
         * @param aMTKMesh
         * @param aIPMeshDataBase
         */

        Integration_Mesh_Editor(
                moris::mtk::Integration_Mesh*               aMTKMesh,
                moris::mtk::Interpolation_Mesh_DataBase_IP* aIPMeshDataBase,
                bool                                        aCheckMesh = true );

        // ----------------------------------------------------------------------------

        /**
         * @brief Destroy the Integration_Mesh_Editor object
         *
         */

        virtual ~Integration_Mesh_Editor();

        // ----------------------------------------------------------------------------
        /**
         * @brief  generate vertex data
         *
         */
        void
        generate_vertex_data();

        // ----------------------------------------------------------------------------
        /**
         * @brief  generate cell data
         *
         */
        void
        generate_cell_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine sizes for veretx data
         *
         */

        void
        initialize_vertex_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine sizes for cell data
         *
         */
        void
        initialize_cell_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine sizes for cell cluster data
         *
         */
        void
        initialize_cell_cluster_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief generate cell cluster data
         *
         */

        void
        generate_cell_cluster_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine sizes for side cluster data
         *
         */

        void
        initialize_side_cluster_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief generate side cluster data
         *
         */
        void
        generate_side_cluster_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief perform function that generates vertex,cell,...data
         *
         */

        Integration_Mesh_DataBase_IG*
        perform();

        // ----------------------------------------------------------------------------

        /**
         * @brief Get the num side clusters
         *
         * @return uint
         */

        uint
        get_num_side_clusters();

        // ----------------------------------------------------------------------------

        /**
         * @brief generates the block set data
         *
         */

        void
        generate_block_set_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine the size of the block set data
         *
         */

        void
        initialize_block_set_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine the size of the double sided cluster data
         *
         */

        void
        initialize_double_sided_cluster_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief generate the double sided cluster data
         *
         */

        void
        generate_double_sided_cluster_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief determine the size of the ghost data
         *
         */

        void
        initialize_ghost_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief generate the ghost data
         *
         */


        void
        generate_ghost_data();

        // ----------------------------------------------------------------------------

        /**
         * @brief Get the memory usage
         *
         * @return moris::Memory_Map
         */

        moris::Memory_Map
        get_memory_usage();

        // ----------------------------------------------------------------------------

        /**
         * @brief delete unwanted memory at the end
         *
         */

        void
        free_memory();


        /**
         * @brief Create the vertices lagrange
         *
         */

        void create_vertices();

        // ----------------------------------------------------------------------------

        /**
         * @brief create lagrange ig cells
         *
         */

        void create_cells();

        // ----------------------------------------------------------------------------

        /**
         * @brief create the cell clusters
         *
         */

        void create_cell_clusters();

        // ----------------------------------------------------------------------------

        /**
         * @brief create side clusters
         *
         */

        void create_side_clusters();

        // ----------------------------------------------------------------------------

        /**
         * @brief create block sets
         *
         */

        void
        create_block_sets();


        // ----------------------------------------------------------------------------

        /**
         * @brief create side sets
         *
         */

        void
        create_side_sets();

        // ----------------------------------------------------------------------------

        /**
         * @brief create double side cluster
         *
         */

        void
        create_double_sided_clusters();

        // ----------------------------------------------------------------------------

        /**
         * @brief create double sided sets
         *
         */

        void
        create_double_sided_sets();

        // ----------------------------------------------------------------------------

        /**
         * @brief Create a ghost clusters objects
         *
         */

        void
        create_ghost_clusters();

        // ----------------------------------------------------------------------------

        /**
         * @brief Create a ghost sets object
         *
         */

        void
        create_ghost_sets();

        // ----------------------------------------------------------------------------

        /**
         * @brief Create the vertex glb id to loc vertex ind map
         *
         */
        void
        create_vertex_glb_id_to_loc_vertex_ind_map();

        // ----------------------------------------------------------------------------

        /**
         * @brief populate the map for the blockset name-cell topology
         *
         */

        void set_blockset_topology();

        // ----------------------------------------------------------------------------

        /**
         * @brief populate the map for the blockset name-cell topology
         *
         */

        void
        create_mesh();

        // ----------------------------------------------------------------------------

        /**
         * @brief populate the map for the blockset name-cell topology
         *
         */

        void
        add_vertices(
                moris::Cell< moris::Cell< moris_index > >& aSideClusterToVertexIndices,
                Matrix< DDRMat >                           aVerticesCoords );

        // ----------------------------------------------------------------------------

        /**
         * @brief populate the map for the blockset name-cell topology
         *
         */

        void
        add_cells( moris::Cell< moris::Cell< moris_index > >& aSideClusterToCells,
                moris::Cell< moris::Cell< moris_index > >&    aCellToVertexIndices );

        // ----------------------------------------------------------------------------

        virtual void
        construct_periodic_data_base(
                moris::Cell< moris::Cell< moris_index > >& aSideClusterToVertexIndices,
                Matrix< DDRMat >                           aVerticesCoords,
                moris::Cell< moris::Cell< moris_index > >& aSideClusterToCells,
                moris::Cell< moris::Cell< moris_index > >& aCellToVertexIndices,
                moris::Cell< moris_index >&                aSideClusterToIPCell,
                Matrix< DDRMat >&                          aVertexParametricCoords,
                moris::Cell< moris_index >&                aDoubleSidedClustersIndex,
                uint                                       mNumDblSideCluster,
                uint                                       aNumGeometry );

        // ----------------------------------------------------------------------------

        void
        create_parallel_consistent_new_vertex_ids( moris_index tNumPreviousVertices );

        // ----------------------------------------------------------------------------

        void
        add_side_clusters(
                moris::Cell< moris::Cell< moris_index > >& aSideClusterToCells,
                moris::Cell< moris_index >&                aSideClusterToIPCell,
                Matrix< DDRMat >&                          aVertexParametricCoords,
                moris::Cell< moris::Cell< moris_index > >& aSideClusterToVertexIndices );

        // ----------------------------------------------------------------------------

        void
        add_double_sided_clusters(
                uint                                       mNumDblSideCluster,
                moris::Cell< moris::Cell< moris_index > >& aSideClusterToVertexIndices );

        // ----------------------------------------------------------------------------

        void
        add_double_sided_set(
                moris::Cell< moris_index >& aDoubleSidedClustersIndex,
                uint                        aNumGeometry );

        // ----------------------------------------------------------------------------

        void
        reconstruct_connectivity();

        // ----------------------------------------------------------------------------

        void
        merge_meshes();

        // ----------------------------------------------------------------------------

        void
        recreate_side_sets();

        // ----------------------------------------------------------------------------

        void
        create_parallel_consistnet_cell_ids( moris_index aNumNewCells );

        // ----------------------------------------------------------------------------

        void
        generate_mesh_data();

        // ----------------------------------------------------------------------------

        friend class Integration_Mesh_Analysis;
        friend class Periodic2D_Analysis;

        /**
         * @name These function are for debug purposed
         */

        ////////----------------------------------------------------------------------------
        // Checking and debugging functions
        ////////----------------------------------------------------------------------------

        ///@{

        bool
        check_vertices();

        // ----------------------------------------------------------------------------

        bool
        check_maps();

        // ----------------------------------------------------------------------------

        bool
        check_cells();

        // ----------------------------------------------------------------------------

        bool
        check_cell_clusters();

        // ----------------------------------------------------------------------------

        bool
        check_side_clusters();

        // ----------------------------------------------------------------------------

        bool
        check_block_sets();

        // ----------------------------------------------------------------------------

        bool
        check_side_sets();

        // ----------------------------------------------------------------------------

        bool
        check_double_sided_clusters();

        // ----------------------------------------------------------------------------

        bool
        check_double_sided_sets();

        // ----------------------------------------------------------------------------

        bool
        check_ghost_clusters();

        // ----------------------------------------------------------------------------

        void
        check_input_output_mesh();

        ///@}
    };

}    // namespace moris::mtk


#endif /* cl_MTK_Integration_Mesh_Editor.hpp */
