/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_DataBase_IG.hpp
 *
 */

#ifndef SRC_cl_MTK_Mesh_DataBase_IG
#define SRC_cl_MTK_Mesh_DataBase_IG

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_TOL_Memory_Map.hpp"

namespace moris
{
    namespace mtk
    {
        // Forward declratios to
        class Integration_Mesh_DataBase;
        class Vertex_DataBase;
        class Vertex;
        class Cell_DataBase;
        // class Interpolation_Mesh;
        class Cell_Cluster_DataBase;
        class Side_Cluster_DataBase;
        class Interpolation_Mesh_Analysis;
        class Interpolation_Mesh_DataBase_IP;
        class Integration_Mesh_Info;

        class Integration_Mesh_DataBase_IG : public mtk::Integration_Mesh
        {
          private:
            Integration_Mesh_Info* mIGMeshInfo = nullptr;

            mtk::Interpolation_Mesh_DataBase_IP* mIPMesh = nullptr;

            // Vertex Information
            moris::Cell< Vertex_DataBase > mVertices;
            Matrix< moris::DDRMat >        mVertexCoordinates;

            // Cell to vertex connectivity
            moris::Cell< mtk::Vertex* > mCellToVertices;
            moris::Cell< moris_index >  mCellToVertexOffSet;

            // Cell cluster to primary IG cell connectivity
            moris::Cell< mtk::Cell* >   mCellClusterToPrimaryIGCell;
            moris::Cell< mtk::Cell* >   mCellClusterToVoidIGCell;
            moris::Cell< mtk::Vertex* > mCellClusterToVeretx;

            // offset data for cell cluster
            moris::Cell< moris_index > mCellClusterToPrimaryIGCellOffSet;
            moris::Cell< moris_index > mCellClusterToVoidIGCellOffset;
            moris::Cell< moris_index > mCellClusterToVertexOffset;

            // local coordinates of the cell cluster
            moris::Matrix< moris::DDRMat >* mCellClusterVertexCoords = nullptr;

            std::unordered_map< moris::moris_index, moris::moris_index > mCellClusterIndexToRowNumber;

            moris::Cell< mtk::Cell_Cluster_DataBase > mCellClusters;

            // number of complete side clusters
            uint mNumSideClusters;

            moris::Cell< mtk::Cell* >   mSideClusterToPrimaryIGCell;
            moris::Cell< moris_index >  mSideClusterToPrimaryIGCellSideOrd;
            moris::Cell< mtk::Cell* >   mSideClusterToVoidIGCell;
            moris::Cell< mtk::Vertex* > mSideClusterToVeretx;

            moris::Cell< moris_index > mSideClusterToPrimaryIGCellOffset;
            moris::Cell< moris_index > mSideClusterToVoidIGCellOffset;
            moris::Cell< moris_index > mSideClusterToVertexOffSet;

            // side cluster to IP cell
            moris::Cell< moris_index > mSideClusterToIPCell;

            moris::Cell< mtk::Side_Cluster_DataBase > mSideClusters;

            moris::Cell< mtk::Double_Side_Cluster > mDblSideClusters;

            moris::Cell< mtk::Side_Cluster_DataBase > mGhostLeader;
            moris::Cell< mtk::Side_Cluster_DataBase > mGhostFollower;
            moris::Cell< mtk::Double_Side_Cluster >   mGhostDblSidedSet;

            moris::Cell< moris_index > mGhostLeaderFollowerIPCellList;
            moris::Cell< mtk::Cell* >  mGhostLeaderFollowerIGCellList;

            moris::Cell< moris_index > mGhostLeaderFollowerOrd;

            moris::Cell< bool > mGhostLeaderFollowerIsTrivial;

            moris::Cell< moris_index > mGhostLeaderFollowerVertexOffSet;

            moris::Cell< mtk::Vertex* > mGhostLeaderFollowerToVertex;
            // Cell Information
            moris::Cell< Cell_DataBase > mCells;

            moris::map< std::string, enum CellTopology > mCellTopologyToNameMap;

            // vertex map ( used in GEN)
            std::unordered_map< moris_id, moris_index > mVertexGlobalIdToLocalIndex;

            moris::Cell< moris_id > mVertexIdList;
            moris::Cell< moris_id > mCellIdList;
            moris::Cell< moris_id > mVertexOwnerList;
            moris::Cell< moris_id > mCellOwnerList;
            moris::Cell< bool >     mCellClusterIsTrivial;
            moris::Cell< bool >     mSideClusterIsTrivial;

            moris::Cell< std::shared_ptr< moris::mtk::Cell_Info > > mCellInfoList;

            std::unordered_map< moris::moris_index, moris::moris_index > mSideClusterIndexToRowNumber;

            // coordinate of the ghost when they are not trivial
            moris::Matrix< moris::DDRMat >*                              mSecondaryClusterVertexCoords = nullptr;
            std::unordered_map< moris::moris_index, moris::moris_index > mSecondaryClusterIndexToRowNumber;

          public:
            // ----------------------------------------------------------------------------

            Integration_Mesh_DataBase_IG() = default;

            /**
             * @brief Construct a new Integration_Mesh_DataBase_IG object
             *
             * @param aIPDataBase an IP mesh data based containing the raw data
             * @param aIPMesh an already existing IP mesh
             */
            Integration_Mesh_DataBase_IG(
                    Integration_Mesh_DataBase*   aIGDataBase,
                    Integration_Mesh&            aIGMesh,
                    Interpolation_Mesh_Analysis& aIPMesh );

            // ----------------------------------------------------------------------------

            /**
             * @brief Destroy the Integration_Mesh_DataBase_IG object
             *
             */
            virtual ~Integration_Mesh_DataBase_IG();

            // ----------------------------------------------------------------------------

            /**
             * @name Get Functions for the FEM to work with
             *  These group of functions are implemented in the base class and some are unnecessary
             */

            ///@{

            /**
             * @brief Get the spatial dim
             *
             * @return uint the spatial dimension of the mesh
             */

            virtual uint get_spatial_dim() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the num entities object
             *
             * @param aEntityRank the entity rank (node, element)
             * @param aIndex the index of the mesh
             * @return uint
             */

            virtual uint get_num_entities(
                    enum EntityRank   aEntityRank,
                    const moris_index aIndex = 0 ) const override;
            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / non-writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex const&
             */

            virtual Vertex const & get_mtk_vertex( moris_index aVertexIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            virtual Vertex& get_mtk_vertex( moris_index aVertexIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            virtual mtk::Cell& get_mtk_cell( moris_index aCellIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            virtual mtk::Cell const & get_mtk_cell( moris_index aCellIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the communication table object
             *
             * @return Matrix< IdMat > communication table
             */

            virtual Matrix< IdMat >
            get_communication_table() const override;

            // ---------------------------------------------------------------------------

            /**
             * @brief Get the vertex glb id to loc vertex ind map object , used in gen
             *
             * @return std::unordered_map< moris_id, moris_index > glb id to loc vertex ind map
             */

            virtual std::unordered_map< moris_id, moris_index >
            get_vertex_glb_id_to_loc_vertex_ind_map() const override;

            // ----------------------------------------------------------------------------

            /**
             * Get a cell cluster related to an interpolation
             * cell
             */

            virtual Cell_Cluster const &
            get_cell_cluster( Cell const & aInterpCell ) const override;

            // ----------------------------------------------------------------------------

            /**
             * Get a cell cluster related to an interpolation
             * cell
             */

            virtual Cell_Cluster const &
            get_cell_cluster( moris_index aInterpCellIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the block set names object
             *
             * @return moris::Cell< std::string > block set names
             */

            virtual moris::Cell< std::string >
            get_block_set_names() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the cell clusters in set
             *
             * @param aBlockSetOrdinal
             * @return moris::Cell< Cluster const* > cell clusters
             */

            virtual moris::Cell< Cluster const * >
            get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the side set cluster
             *
             * @param aSideSetOrdinal
             * @return moris::Cell< Cluster const* > side clusters
             */

            virtual moris::Cell< Cluster const * >
            get_side_set_cluster( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * get number of side sets
             */

            virtual uint
            get_num_side_sets() const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the label
             */

            virtual std::string
            get_side_set_label( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the index given a label
             */

            virtual moris_index
            get_side_set_index( std::string aSideSetLabel ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the number of double sided side sets in the mesh
             */

            virtual uint
            get_num_double_sided_sets() const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the label
             */

            virtual std::string
            get_double_sided_set_label( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the double side clusters in the side set
             */

            virtual moris::Cell< Cluster const * >
            get_double_side_set_cluster( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------
            // Non Used Mesh Function (so far)
            // ----------------------------------------------------------------------------

            virtual MeshType get_mesh_type() const override;

            // ----------------------------------------------------------------------------

            virtual Matrix< IndexMat >
            get_entity_connected_to_entity_loc_inds(
                    moris_index       aEntityIndex,
                    enum EntityRank   aInputEntityRank,
                    enum EntityRank   aOutputEntityRank,
                    const moris_index aDiscretizationIndex = 0 ) const override;

            // ----------------------------------------------------------------------------

            Matrix< DDRMat >
            get_node_coordinate( moris_index aNodeIndex ) const override;

            // ----------------------------------------------------------------------------

            virtual moris_id
            get_glb_entity_id_from_entity_loc_index(
                    moris_index       aEntityIndex,
                    enum EntityRank   aEntityRank,
                    const moris_index aDiscretizationIndex = 0 ) const override;

            // ----------------------------------------------------------------------------

            virtual uint get_node_owner( moris_index aNodeIndex ) const override;

            // ----------------------------------------------------------------------------

            virtual Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex ) override;

            // ----------------------------------------------------------------------------

            virtual enum CellTopology
            get_blockset_topology( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------
            virtual enum CellShape
            get_IG_blockset_shape( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            virtual enum CellShape
            get_IP_blockset_shape( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            virtual uint get_element_owner( moris_index aElementIndex ) const override;

            // ----------------------------------------------------------------------------

            virtual Matrix< IndexMat >
            get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const override;

            ///@}

            // ----------------------------------------------------------------------------

            /**
             * @name Get functions for smaller entities
             *
             */

            ///@{

            ////////----------------------------------------------------------------------------
            // Checking and debugging functions
            ////////----------------------------------------------------------------------------

            /**
             * @brief Get the entity id
             *
             * @param aEntityRank
             * @param aEntityIndex
             * @return moris_id const
             */

            virtual moris_id
            get_entity_id(
                    enum EntityRank aEntityRank,
                    moris_index     aEntityIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the entity owner
             *
             * @param aEntityRank
             * @param aEntityIndex
             * @return moris_id const
             */

            virtual moris_id
            get_entity_owner(
                    enum EntityRank aEntityRank,
                    moris_index     aEntityIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the vertex coords ptr
             *
             * @param aVertexIndex
             * @return moris::real* const
             */

            virtual moris::real*
            get_vertex_coords_ptr( moris_index aVertexIndex );

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the cell vertices object
             *
             * @param aCellIndex
             * @return Vertex** the first pointer is a pointer to the location of the cell holding the vertex pointers
             */

            virtual Vertex**
            get_cell_vertices( moris_index aCellIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the ip cell in cluster
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return mtk::Cell* const
             */

            virtual mtk::Cell*
            get_ip_cell_in_cluster(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the ig cells in cluster
             *
             * @param aClusterType
             * @param aPrimaryOrVoid
             * @param aClusterIndex
             * @return mtk::Cell* const*  it is a pointer (location) to the ig cells
             */

            virtual mtk::Cell* const *
            get_ig_cells_in_cluster(
                    enum ClusterType       aClusterType,
                    enum mtk::Primary_Void aPrimaryOrVoid,
                    moris_index            aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the num cells in cluster
             *
             * @param aClusterType
             * @param aPrimaryOrVoid
             * @param aClusterIndex
             * @return uint const
             */

            virtual uint
            get_num_cells_in_cluster(
                    enum ClusterType       aClusterType,
                    enum mtk::Primary_Void aPrimaryOrVoid,
                    moris_index            aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the side ordinals in cluster object
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return moris_index* const a pointer to the moris_indices which are side ordinals
             */

            virtual moris_index*
            get_side_ordinals_in_cluster(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief if a cluster is trivial or nor
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return true cluster is trivial
             * @return false cluster is not trivial
             */

            virtual bool
            cluster_is_trivial(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the vertices in cluster
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return Vertex* const*
             */

            virtual Vertex* const *
            get_vertices_in_cluster(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the num vertices in cluster
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return uint const
             */

            virtual uint
            get_num_vertices_in_cluster(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the local coord matrix pointer
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return Matrix< DDRMat >* const
             */

            virtual Matrix< DDRMat >*
            get_local_coord_matrix_ptr(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the row number in the local coords matrix
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return uint const
             */

            virtual uint
            get_row_number_local_coords_matrix(
                    enum ClusterType aClusterType,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the associated cell cluster object
             *
             * @param aClusterIndex
             * @return mtk::Cell_Cluster const*
             */

            virtual mtk::Cell_Cluster const *
            get_associated_cell_cluster( moris_index aClusterIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the associated cell cluster object
             *
             * @param aClusterIndex
             * @return mtk::Cell_Cluster const*
             */

            virtual std::shared_ptr< mtk::Cell_Info >
            get_cell_info_sp( moris_index aEntityIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief this refres to the ghost side clusters and possibly other cluster( periodic o contact)
             *
             * @param aClusterIndex
             * @return true
             * @return false
             */

            virtual bool
            is_secondary_cluster( moris_index aClusterIndex ) const override;

            ///@}

            // ----------------------------------------------------------------------------

            /**
             * @brief deletes the unused data
             *
             */

            void
            free_memory();

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the memory usage
             *
             * @return moris::Memory_Map
             */

            moris::Memory_Map
            get_memory_usage();

            friend class Periodic2D_Analysis;
            friend class Integration_Mesh_Editor;
        };

    }    // namespace mtk
}    // namespace moris

#endif /* cl_MTK_Integration_Mesh_DataBase_IG.hpp */
