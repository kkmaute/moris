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
#include "cl_MTK_Nonconformal_Side_Cluster.hpp"

namespace moris::mig
{
    class Periodic_Mesh_Editor;
}

namespace moris
{
    class Memory_Map;
    namespace mtk
    {
        // Forward declarations to
        class Integration_Mesh_DataBase;
        class Vertex_DataBase;
        class Vertex;
        class Cell_DataBase;
        // class Interpolation_Mesh;
        class Cell_Cluster_DataBase;
        class Side_Cluster_DataBase;
        class Interpolation_Mesh_Analysis;
        class Interpolation_Mesh_DataBase_IP;

        struct Integration_Mesh_Info;

        class Integration_Mesh_DataBase_IG : public mtk::Integration_Mesh
        {
            friend class Periodic2D_Analysis;
            friend class Integration_Mesh_Editor;
            friend class mig::Periodic_Mesh_Editor;

          private:
            Integration_Mesh_Info* mIGMeshInfo = nullptr;

            mtk::Interpolation_Mesh_DataBase_IP* mIPMesh = nullptr;

            // Vertex Information
            Vector< Vertex_DataBase > mVertices;
            Matrix< moris::DDRMat >   mVertexCoordinates;

            // Cell to vertex connectivity
            Vector< mtk::Vertex* > mCellToVertices;
            Vector< moris_index >  mCellToVertexOffSet;

            // Cell cluster to primary IG cell connectivity
            Vector< mtk::Cell* >   mCellClusterToPrimaryIGCell;
            Vector< mtk::Cell* >   mCellClusterToVoidIGCell;
            Vector< mtk::Vertex* > mCellClusterToVeretx;

            // offset data for cell cluster
            Vector< moris_index > mCellClusterToPrimaryIGCellOffSet;
            Vector< moris_index > mCellClusterToVoidIGCellOffset;
            Vector< moris_index > mCellClusterToVertexOffset;

            // local coordinates of the cell cluster
            moris::Matrix< moris::DDRMat >* mCellClusterVertexCoords = nullptr;

            std::unordered_map< moris::moris_index, moris::moris_index > mCellClusterIndexToRowNumber;

            Vector< mtk::Cell_Cluster_DataBase > mCellClusters;

            // number of complete side clusters
            uint mNumSideClusters;

            Vector< mtk::Cell* >   mSideClusterToPrimaryIGCell;
            Vector< moris_index >  mSideClusterToPrimaryIGCellSideOrd;
            Vector< mtk::Cell* >   mSideClusterToVoidIGCell;
            Vector< mtk::Vertex* > mSideClusterToVeretx;

            Vector< moris_index > mSideClusterToPrimaryIGCellOffset;
            Vector< moris_index > mSideClusterToVoidIGCellOffset;
            Vector< moris_index > mSideClusterToVertexOffSet;

            // side cluster to IP cell
            Vector< moris_index > mSideClusterToIPCell;

            Vector< mtk::Side_Cluster_DataBase >     mSideClusters;
            Vector< mtk::Double_Side_Cluster >       mDblSideClusters;
            Vector< mtk::Nonconformal_Side_Cluster > mNonconformalSideClusters;

            Vector< mtk::Side_Cluster_DataBase > mGhostLeader;
            Vector< mtk::Side_Cluster_DataBase > mGhostFollower;
            Vector< mtk::Double_Side_Cluster >   mGhostDblSidedSet;

            Vector< moris_index > mGhostLeaderFollowerIPCellList;
            Vector< mtk::Cell* >  mGhostLeaderFollowerIGCellList;

            Vector< moris_index > mGhostLeaderFollowerOrd;

            Vector< bool > mGhostLeaderFollowerIsTrivial;

            Vector< moris_index > mGhostLeaderFollowerVertexOffSet;

            Vector< mtk::Vertex* > mGhostLeaderFollowerToVertex;
            // Cell Information
            Vector< Cell_DataBase > mCells;

            moris::map< std::string, enum CellTopology > mNameToCellTopologyMap;    // key: set name, value: cell topology

            // vertex map (used in GEN)
            std::unordered_map< moris_id, moris_index > mVertexGlobalIdToLocalIndex;

            Vector< moris_id > mVertexIdList;
            Vector< moris_id > mCellIdList;
            Vector< moris_id > mVertexOwnerList;
            Vector< moris_id > mCellOwnerList;
            Vector< bool >     mCellClusterIsTrivial;
            Vector< bool >     mSideClusterIsTrivial;

            Vector< std::shared_ptr< moris::mtk::Cell_Info > > mCellInfoList;

            std::unordered_map< moris::moris_index, moris::moris_index > mSideClusterIndexToRowNumber;

            // coordinate of the ghost when they are not trivial
            moris::Matrix< moris::DDRMat >*                              mSecondaryClusterVertexCoords = nullptr;
            std::unordered_map< moris::moris_index, moris::moris_index > mSecondaryClusterIndexToRowNumber;
            bool                                                         mMomentFittingFlag = false;

            // Moment fitting points
            Matrix< DDRMat >                                             mMomentFittingPoints;

          public:
            Integration_Mesh_DataBase_IG() = default;

            // ----------------------------------------------------------------------------
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
            ~Integration_Mesh_DataBase_IG() override;

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
            uint get_spatial_dim() const override;

            // ----------------------------------------------------------------------------

            void reserve_nonconformal_side_clusters( uint aNumNonconformalSideClusters );

            // ----------------------------------------------------------------------------

            void add_nonconformal_side_set( std::string const & aName,
                    Vector< Nonconformal_Side_Cluster > const & aClusters,
                    Matrix< IndexMat > const &                  aColors );

            // ----------------------------------------------------------------------------

            void reset_nonconformal_side_set();

            // ----------------------------------------------------------------------------

            bool get_moment_fitting_flag() const override;

            // ----------------------------------------------------------------------------

            void set_moment_fitting_flag( bool aMomentFittingFlag ) ;

            // ----------------------------------------------------------------------------

            const Matrix< DDRMat >& get_moment_fitting_points() const override;

            // ----------------------------------------------------------------------------

            void set_moment_fitting_points( Matrix< DDRMat > aMomentFittingPoints ) ;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the num entities object
             *
             * @param aEntityRank the entity rank (node, element)
             * @param aIndex the index of the mesh
             * @return uint
             */

            uint get_num_entities(
                    enum EntityRank   aEntityRank,
                    const moris_index aIndex = 0 ) const override;
            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / non-writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex const&
             */

            Vertex const & get_mtk_vertex( moris_index aVertexIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            Vertex& get_mtk_vertex( moris_index aVertexIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            mtk::Cell& get_mtk_cell( moris_index aCellIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            mtk::Cell const & get_mtk_cell( moris_index aCellIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the communication table object
             *
             * @return Matrix< IdMat > communication table
             */

            Matrix< IdMat >
            get_communication_table() const override;

            // ---------------------------------------------------------------------------

            /**
             * @brief Get the vertex glb id to loc vertex ind map object , used in gen
             *
             * @return std::unordered_map< moris_id, moris_index > glb id to loc vertex ind map
             */

            std::unordered_map< moris_id, moris_index >
            get_vertex_glb_id_to_loc_vertex_ind_map() const override;

            // ----------------------------------------------------------------------------

            /**
             * Get a cell cluster related to an interpolation
             * cell
             */

            Cell_Cluster const &
            get_cell_cluster( Cell const & aInterpCell ) const override;

            // ----------------------------------------------------------------------------

            /**
             * Get a cell cluster related to an interpolation
             * cell
             */

            Cell_Cluster const &
            get_cell_cluster( moris_index aInterpCellIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the block set names object
             *
             * @return Vector< std::string > block set names
             */

            Vector< std::string >
            get_block_set_names() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the cell clusters in set
             *
             * @param aBlockSetOrdinal
             * @return Vector< Cluster const* > cell clusters
             */

            Vector< Cluster const * >
            get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the side set cluster
             *
             * @param aSideSetOrdinal
             * @return Vector< Cluster const* > side clusters
             */

            Vector< Cluster const * >
            get_side_set_cluster( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * get number of side sets
             */

            uint
            get_num_side_sets() const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the label
             */

            std::string
            get_side_set_label( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the index given a label
             */

            moris_index
            get_side_set_index( std::string aSideSetLabel ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the number of double sided side sets in the mesh
             */

            uint
            get_num_double_sided_sets() const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the label
             */

            std::string
            get_double_sided_set_label( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------

            /*!
             * Returns the double side clusters in the side set
             */

            Vector< Cluster const * >
            get_double_side_set_cluster( moris_index aSideSetOrdinal ) const override;

            // ----------------------------------------------------------------------------
            // Non Used Mesh Function (so far)
            // ----------------------------------------------------------------------------

            MeshType get_mesh_type() const override;

            // ----------------------------------------------------------------------------

            Matrix< IndexMat >
            get_entity_connected_to_entity_loc_inds(
                    moris_index       aEntityIndex,
                    enum EntityRank   aInputEntityRank,
                    enum EntityRank   aOutputEntityRank,
                    const moris_index aDiscretizationIndex = 0 ) const override;

            // ----------------------------------------------------------------------------

            Matrix< DDRMat >
            get_node_coordinate( moris_index aNodeIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * Get a global entity ID from an entity rank and local index.
             *
             * @param aEntityIndex Local entity index
             * @param aEntityRank Entity rank
             * @param aDiscretizationIndex Discretization mesh index
             * @return Global entity ID
             */
            moris_id
            get_glb_entity_id_from_entity_loc_index(
                    moris_index aEntityIndex,
                    EntityRank  aEntityRank,
                    moris_index aDiscretizationIndex = 0 ) const override;

            /**
             * Get a local entity ID from an entity rank and global ID
             *
             * @param aEntityId Global entity ID
             * @param aEntityRank Entity rank
             * @param aDiscretizationIndex Discretization mesh index
             * @return Local entity index
             */
            moris_index
            get_loc_entity_ind_from_entity_glb_id(
                    moris_id    aEntityId,
                    EntityRank  aEntityRank,
                    moris_index aDiscretizationIndex = 0 ) const override;

            // ----------------------------------------------------------------------------

            uint get_node_owner( moris_index aNodeIndex ) const override;

            // ----------------------------------------------------------------------------

            Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex ) override;

            // ----------------------------------------------------------------------------

            enum CellTopology
            get_blockset_topology( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            enum CellShape
            get_IG_blockset_shape( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            enum CellShape
            get_IP_blockset_shape( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            uint get_element_owner( moris_index aElementIndex ) const override;

            // ----------------------------------------------------------------------------

            Matrix< IndexMat >
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
            moris_id
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
            moris_id
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
            moris::real*
            get_vertex_coords_ptr( moris_index aVertexIndex ) override;

            // ----------------------------------------------------------------------------
            /**
             * @brief Get the cell vertices object
             *
             * @param aCellIndex
             * @return Vertex** the first pointer is a pointer to the location of the cell holding the vertex pointers
             */
            Vertex**
            get_cell_vertices( moris_index aCellIndex ) override;

            // ----------------------------------------------------------------------------
            /**
             * @brief Get the ip cell in cluster
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return mtk::Cell* const
             */
            mtk::Cell*
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
            mtk::Cell* const *
            get_ig_cells_in_cluster(
                    enum ClusterType aClusterType,
                    Primary_Void     aPrimaryOrVoid,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------
            /**
             * @brief Get the num cells in cluster
             *
             * @param aClusterType
             * @param aPrimaryOrVoid
             * @param aClusterIndex
             * @return uint const
             */
            uint
            get_num_cells_in_cluster(
                    enum ClusterType aClusterType,
                    Primary_Void     aPrimaryOrVoid,
                    moris_index      aClusterIndex ) const override;

            // ----------------------------------------------------------------------------
            /**
             * @brief Get the side ordinals in cluster object
             *
             * @param aClusterType
             * @param aClusterIndex
             * @return moris_index* const a pointer to the moris_indices which are side ordinals
             */
            moris_index*
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
            bool
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
            Vertex* const *
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
            uint
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
            Matrix< DDRMat >*
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
            uint
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
            mtk::Cell_Cluster const *
            get_associated_cell_cluster( moris_index aClusterIndex ) const override;

            // ----------------------------------------------------------------------------
            /**
             * @brief Get the associated cell cluster object
             *
             * @param aClusterIndex
             * @return mtk::Cell_Cluster const*
             */
            std::shared_ptr< mtk::Cell_Info >
            get_cell_info_sp( moris_index aEntityIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief this refres to the ghost side clusters and possibly other cluster( periodic o contact)
             *
             * @param aClusterIndex
             * @return true
             * @return false
             */
            bool
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
        };

    }    // namespace mtk
}    // namespace moris

#endif /* cl_MTK_Integration_Mesh_DataBase_IG.hpp */
