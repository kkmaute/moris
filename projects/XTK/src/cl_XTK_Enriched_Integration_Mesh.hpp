/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_XTK_Enriched_Integration_Mesh.hpp  
 * 
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTEGRATION_MESH_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTEGRATION_MESH_HPP_

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Vertex.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_XTK_Field.hpp"
#include <unordered_map>
#include "cl_TOL_Memory_Map.hpp"
#include "cl_XTK_Subphase_Group.hpp"
#include "cl_MTK_Cluster_Group.hpp"

using namespace moris;

namespace xtk
{
class Model;
class Cell_Cluster;
class Side_Cluster;
class Interpolation_Cell_Unzipped;
class Ghost_Stabilization;
class Enrichment;
class Cut_Integration_Mesh;

class Enriched_Integration_Mesh : public mtk::Integration_Mesh
{
    //------------------------------------------------------------------------------

  protected:

    friend class Enrichment;
    friend class Ghost_Stabilization;

    Model*                   mModel;
    Cut_Integration_Mesh*    mCutIgMesh;

    // mesh indices
    Matrix< IndexMat > mBsplineMeshIndices;

    // Cell Clusters
    moris::Cell< std::shared_ptr< xtk::Cell_Cluster > > mCellClusters;

    // cluster groups, 
    // NOTE: first array index corresponds to index of the B-spline mesh in mBsplineMeshIndices 
    moris::Cell< moris::Cell< std::shared_ptr< mtk::Cluster_Group > > > mCellClusterGroups;
    moris::Cell< moris::Cell< std::shared_ptr< mtk::Cluster_Group > > > mSideClusterGroups;
    moris::Cell< moris::Cell< std::shared_ptr< mtk::Cluster_Group > > > mDblSideClusterGroups;

    // Vertex Set
    std::unordered_map< std::string, moris_index >     mVertexSetLabelToOrd;
    moris::Cell< std::string >                         mVertexSetNames;
    moris::Cell< moris::Cell< moris::mtk::Vertex * > > mVerticesInVertexSet;
    moris::Cell< moris::Matrix< IndexMat > >           mVertexSetColors;

    // Block sets containing Cell Clusters
    std::unordered_map< std::string, moris_index >          mBlockSetLabelToOrd;
    moris::Cell< std::string >                              mBlockSetNames;
    moris::Cell< enum CellTopology >                        mBlockSetTopology;
    moris::Cell< moris::Cell< xtk::Cell_Cluster const * > > mPrimaryBlockSetClusters;
    moris::Cell< moris::Matrix< IndexMat > >                mBlockSetColors; /*Bulk phases*/
    moris::Cell< moris::Cell< moris_index > >               mColorsBlockSets; /*transpose of mBlockSetColors*/

    // side sets
    std::unordered_map< std::string, moris_index >                     mSideSideSetLabelToOrd;
    moris::Cell< std::string >                                         mSideSetLabels;
    moris::Cell< moris::Cell< std::shared_ptr< xtk::Side_Cluster > > > mSideSets;
    moris::Cell< moris::Matrix< IndexMat > >                           mSideSetColors; /*Bulk phases of cells attached to side*/
    moris::Cell< moris::Cell< moris_index > >                          mColorsSideSets; /*transpose of mSideSetColors*/

    // double side sets
    std::unordered_map< std::string, moris_index >             mDoubleSideSetLabelToOrd;
    moris::Cell< std::string >                                 mDoubleSideSetLabels;
    moris::Cell< moris::Cell< std::shared_ptr< mtk::Double_Side_Cluster > > > mDoubleSideSets;
    moris::Cell< moris::Cell< moris_index > >                  mDoubleSideSetsMasterIndex;
    moris::Cell< moris::Cell< moris_index > >                  mDoubleSideSetsSlaveIndex;
    moris::Cell< std::shared_ptr< mtk::Double_Side_Cluster > > mDoubleSideClusters;
    moris::Cell< std::shared_ptr< xtk::Side_Cluster > >        mDoubleSideSingleSideClusters; /*lefts and rights of the double side sets*/
    moris::Matrix< moris::IndexMat >                           mBulkPhaseToDblSideIndex;
    moris::Cell< moris::Matrix< IndexMat > >                   mMasterDoubleSideSetColor;
    moris::Cell< moris::Matrix< IndexMat > >                   mSlaveDoubleSideSetColor;
    moris::Cell< moris::Cell< moris_index > >                  mColorMasterDoubleSideSet; /*transpose of mMasterDoubleSideSetColor*/
    moris::Cell< moris::Cell< moris_index > >                  mColorSlaveDoubleSideSet; /*transpose of mSlaveDoubleSideSetColor*/

    // Fields
    moris::Cell< xtk::Field >                                     mFields; /*Structure Node (0), Cell(1)*/
    moris::Cell< std::unordered_map< std::string, moris_index > > mFieldLabelToIndex;

    // Sub phase index to Cell Cluster Index (these only include the standard cluster i.e. non-ghost clusters.)
    moris::Matrix< moris::IndexMat > mSubphaseIndexToClusterIndex; // input: enr IP cell (= cluster) index || output: subphase index
    moris::Cell< moris::Cell< moris_index > > mClusterIndexToSubphaseIndices; // input: enr IP cell (= cluster) index || output: List of subphase indices in cluster

    // a connectivity pointer used by all transition cells
    moris::mtk::Cell_Info* mCellInfo;

    //------------------------------------------------------------------------------

  public:
    //------------------------------------------------------------------------------
    Enriched_Integration_Mesh( Model* aXTKModel );

    Enriched_Integration_Mesh( 
        Model*                   aXTKModel,
        const Matrix< IndexMat > aBsplineMeshIndices );
    //------------------------------------------------------------------------------
    ~Enriched_Integration_Mesh();
    //------------------------------------------------------------------------------
    // MTK Mesh Core Functionality (see base class mtk::Mesh for documentation)
    //------------------------------------------------------------------------------
    MeshType                                    get_mesh_type() const;
    moris::uint                                 get_spatial_dim() const;
    uint                                        get_num_entities( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
    Matrix< IndexMat >                          get_entity_connected_to_entity_loc_inds( moris_index aEntityIndex, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const;
    Matrix< IndexMat >                          get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const;
    moris_id                                    get_glb_entity_id_from_entity_loc_index( moris_index aEntityIndex, enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
    std::unordered_map< moris_id, moris_index > get_vertex_glb_id_to_loc_vertex_ind_map() const;
    moris_index                                 get_loc_entity_ind_from_entity_glb_id( moris_id aEntityId, enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
    Cell< mtk::Vertex const * >                 get_all_vertices() const;
    Matrix< IdMat >                             get_entity_connected_to_entity_glob_ids( moris_id aEntityId, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const;
    Matrix< DDRMat >                            get_node_coordinate( moris_index aNodeIndex ) const;
    mtk::Vertex &                               get_mtk_vertex( moris_index aVertexIndex );
    mtk::Vertex const &                         get_mtk_vertex( moris_index aVertexIndex ) const;
    mtk::Cell &                                 get_writable_mtk_cell( moris_index aElementIndex );
    mtk::Cell &                                 get_mtk_cell( moris_index aElementIndex );
    mtk::Cell const &                           get_mtk_cell( moris_index aElementIndex ) const;
    Matrix< IdMat >                             get_communication_table() const;
    moris::Cell< std::string >                  get_set_names( enum EntityRank aSetEntityRank ) const;
    enum CellTopology                           get_blockset_topology( const std::string &aSetName );
    enum CellShape                              get_IG_blockset_shape( const std::string &aSetName );
    enum CellShape                              get_IP_blockset_shape( const std::string &aSetName );
    Matrix< IndexMat >                          get_set_entity_loc_inds( enum EntityRank aSetEntityRank, std::string aSetName ) const;
    Matrix< IndexMat >                          get_element_indices_in_block_set( uint aSetIndex );
    void                                        get_sideset_elems_loc_inds_and_ords( const std::string &aSetName, Matrix< IndexMat > &aElemIndices, Matrix< IndexMat > &aSidesetOrdinals ) const;
    moris_id                                    get_max_entity_id( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
    uint                                        get_node_owner( moris_index aNodeIndex ) const;
    uint                                        get_element_owner( moris_index aElementIndex ) const;

    //------------------------------------------------------------------------------
    // end mesh core functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // MTK Integration Mesh Functions
    // see base class mtk::Integration_Mesh for documentation
    //------------------------------------------------------------------------------
    mtk::Cell_Cluster const &           get_cell_cluster( mtk::Cell const &aInterpCell ) const;
    mtk::Cell_Cluster const &           get_cell_cluster( moris_index aInterpCellIndex ) const;
    moris::Cell< std::string >          get_block_set_names() const;
    std::string                         get_block_set_label( moris_index aBlockSetOrdinal ) const;
    moris_index                         get_block_set_index( std::string aBlockSetLabel ) const;
    moris::Cell< mtk::Cluster const * > get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const;
    Matrix< IndexMat >                  get_block_set_colors( moris_index aBlockSetOrdinal ) const;
    moris::Cell< mtk::Cluster const * > get_side_set_cluster( moris_index aSideSetOrdinal ) const;
    Matrix< IndexMat >                  get_side_set_colors( moris_index aSideSetOrdinal ) const;
    uint                                get_num_side_sets() const;
    std::string                         get_side_set_label( moris_index aSideSetOrdinal ) const;
    moris_index                         get_side_set_index( std::string aSideSetLabel ) const;
    uint                                get_num_double_sided_sets() const;
    std::string                         get_double_sided_set_label( moris_index aSideSetOrdinal ) const;
    moris_index                         get_double_sided_set_index( std::string aDoubleSideSetLabel ) const;
    moris::Cell< mtk::Cluster const * > get_double_side_set_cluster( moris_index aSideSetOrdinal ) const;
    Matrix< IndexMat >                  get_double_side_set_colors( moris_index aSideSetOrdinal ) const;
    uint                                get_sidesets_num_faces( moris::Cell< moris_index > aSideSetIndex ) const;
    //------------------------------------------------------------------------------
    // end integration mesh functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Additional Set Functions
    //------------------------------------------------------------------------------
    moris::Cell< xtk::Cell_Cluster const * > const &
    get_xtk_cell_clusters_in_block_set( moris_index aBlockSetOrdinal ) const;

    /*!
     * @return Interface side set name
     */
    std::string
    get_interface_side_set_name( moris_index aGeomIndex,
        moris_index                          aBulkPhaseIndex0,
        moris_index                          aBulkPhaseIndex1 );

    //------------------------------------------------------------------------------
    /*!
     * @return Double sided interface side set name
     */
    std::string
    get_dbl_interface_side_set_name(
        moris_index aBulkPhaseIndex0,
        moris_index aBulkPhaseIndex1 );

    //------------------------------------------------------------------------------

    /*!
     * @return Primary cell local indices in a block set
     */
    Matrix< IndexMat >
    get_block_entity_loc_inds( std::string aSetName ) const;


    //------------------------------------------------------------------------------

    /*!
     * @brief Creates a double sided interface between the two bulk phases.
     * This function creates additional dbl sided interfaces. By default,
     * the enriched integrztion mesh creates only the low-master high-slave
     * dbl sided interfaces. This functions allows creation of  low-slave high-master
     * interfaces.
     * @param[in] aMasterBulkPhaseIndex Master bulk phase index
     * @param[in] aSlaveBulkPhaseIndex Slave bulk phase index
     */
    void
    create_dbl_sided_interface_set(
        moris_index aMasterBulkPhaseIndex,
        moris_index aSlaveBulkPhaseIndex );


    //------------------------------------------------------------------------------
    // Output/ Viz Functions
    //------------------------------------------------------------------------------
    /*!
     * @brief For cleanup when writing to an exodus file (note: in general this should not be used because
     * sets may not be always empty through an optimization run)
     */
    void
    deactivate_empty_sets();

    //------------------------------------------------------------------------------

    /*!
     * @brief Deactive empty side sets in mesh. 
     */
    void
    deactivate_empty_side_sets();


    //------------------------------------------------------------------------------

    /*!
     * @brief Deactive empty block sets in mesh. 
     */
    void
    deactivate_empty_block_sets();

    //------------------------------------------------------------------------------
    /*!
     * @brief Create basis support fields
     * @return Cell of field names for basis support (1xNumBasisFunc)
     */

    moris::Cell< std::string >
    create_basis_support_fields( moris::Matrix< moris::DDRMat > const &aProbeSpheres );

    void
    create_bg_cell_id_field();

    void create_subphase_fields();

    void
    write_subphase_neighborhood( std::string aFile );

    void
    create_cell_id_fields();

    void
    add_mesh_field_real_scalar_data_loc_inds( std::string const &aFieldName,
        enum EntityRank const &                                  aFieldEntityRank,
        Matrix< DDRMat > const &                                 aFieldData );
    //------------------------------------------------------------------------------

    /*!
     * @brief Write mesh
     */
    void
    write_mesh(
        moris::ParameterList *aParamList );


    void
    create_union_block(
        Cell< std::string > const &aBlocks,
        std::string                aNewBlock,
        Matrix< IndexMat > const & aNewBlockColor );

    void
    create_union_side_set(
        Cell< std::string > const &aSideSets,
        std::string                aNewSideSet,
        Matrix< IndexMat > const & aNewSideSetColor );


    void
    deactive_all_blocks_but_selected(
        Cell< std::string > const &aBlockSetsToKeep );

    void
    deactive_all_side_sets_but_selected(
        Cell< std::string > const &aSideSetsToKeep );


    //------------------------------------------------------------------------------
    // Memory Map
    //------------------------------------------------------------------------------
    /*!
     * @brief Memory map of the enriched integration mesh
     * @return Memory map
     */
    moris::Memory_Map
    get_memory_usage();


    //------------------------------------------------------------------------------
    // Additional Field Functions
    //------------------------------------------------------------------------------

    moris::Cell< std::string >
    get_field_names( enum moris::EntityRank aEntityRank );

    /*!
     * @brief Create a field
     * @param[in] aLabel Field label
     * @param[in] aEntityRank Field entity rank
     * @param[in]aBulkPhaseIndex Bulk phase field defined over
     * aBulkphaseIndex of MORIS_INDEX_MAX results in a field over all phases
     */
    moris::moris_index
    create_field( std::string  aLabel,
        enum moris::EntityRank aEntityRank,
        moris::moris_index     aBulkPhaseIndex = MORIS_INDEX_MAX );

    //------------------------------------------------------------------------------

    /*!
     * @brief Given a field name and rank, gets the field index(ordinal)
     * @return Field index
     */
    moris::moris_index
    get_field_index( std::string aLabel,
        enum moris::EntityRank   aEntityRank );

    //------------------------------------------------------------------------------

    /*!
     * @brief Add field data to created field
     * @param[in] aFieldIndex Field index (use fn:get_field_index(...) to get)
     * @param[in] aEntityRank Field entity rank
     * @param[in] aFieldData Field data
     */
    void
    add_field_data( moris::moris_index aFieldIndex,
        enum moris::EntityRank         aEntityRank,
        Matrix< DDRMat > const &       aFieldData );

    /*!
     * @brief Given a field index and field entity rank, get the field data
     * @param[in] aFieldIndex Field index (use fn:get_field_index(...) to get)
     * @param[in] aEntityRank Field entity rank
     * @return Field data
     */
    Matrix< DDRMat > const &
    get_field_data( moris::moris_index aFieldIndex,
        enum moris::EntityRank         aEntityRank ) const;
    //------------------------------------------------------------------------------
    /*!
     * @brief Convert entity indices to entity ids
     * @param[in] aIndices Enitity indices
     * @param[in] aEntityRank Entity rank
     * @return Entity ids
     */
    Matrix< IdMat > convert_indices_to_ids( Matrix< IndexMat > const &aIndices,
        enum EntityRank                                               aEntityRank ) const;
    //------------------------------------------------------------------------------
    /*!
     * @brief Convert entity ids to entity indices
     * @param[in] aIds Enitity ids
     * @param[in] aEntityRank Entity rank
     * @return Entity indices
     */
    Matrix< IndexMat > convert_ids_to_indices( Matrix< IdMat > const &aIds,
        enum EntityRank                                               aEntityRank ) const;
    //------------------------------------------------------------------------------

    /*!
     * @brief Get MTK cells from their indices
     * @param[in] aCellIndices Cell indices
     * @return MTK cells
     */
    moris::Cell< moris::mtk::Cell const * >
    get_mtk_cells_loc_inds( Matrix< IndexMat > const &aCellIndices );

    //------------------------------------------------------------------------------

    /*!
     * @brief Get MTK vertices from their indices
     * @param[in] aVertexIndices Vertex indices
     * @return MTK vertices
     */
    moris::Cell< moris::mtk::Vertex const * >
    get_mtk_vertices_loc_inds( Matrix< IndexMat > const &aVertexIndices );
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // Accessor functions of XTK specific data structures
    //------------------------------------------------------------------------------
    /*!
     * @brief Get the XTK cell cluster implementation
     * @param[in] aInterpCell Interpolation MTK cell
     * @return XTK cell cluster
     * NOTE: requires there is only one cell cluster associated with the interpolation cell
     */
    xtk::Cell_Cluster const &
    get_xtk_cell_cluster( mtk::Cell const &aInterpCell ) const;

    //------------------------------------------------------------------------------
    /*!
     * @brief Get the XTK cell cluster implementation
     * @param[in] aInterpCell Interpolation MTK cell
     * @return XTK cell cluster
     * NOTE: requires there is only one cell cluster associated with the interpolation cell
     */
    xtk::Cell_Cluster const &
    get_xtk_cell_cluster( moris_index aInterpCellIndex ) const;

    //------------------------------------------------------------------------------
    // Printing
    //------------------------------------------------------------------------------
    void print() const;
    void print_general() const;
    void print_cell_clusters( moris::uint aVerbosityLevel = 0 ) const;
    void print_block_sets( moris::uint aVerbosityLevel = 0 ) const;
    void print_side_sets( moris::uint aVerbosityLevel = 0 ) const;
    void print_double_side_sets( moris::uint aVerbosityLevel = 0 ) const;
    void print_double_side_clusters( moris::uint aVerbosityLevel = 0 ) const;

    //--------------------------------------------------------------------------------
    // Utilities for manipulating sets
    //--------------------------------------------------------------------------------
    /*!
     * @brief Create a single side set from a double sided side set
     * @param[in] aDblSideSetIndex Double side set index
     * @param[in] aSideSetName New side set name
     * @param[in] aCollectSets Perform collect operation upon completion of function. 
     * @return Side set index
     */
    moris_index
    create_side_set_from_dbl_side_set( moris_index const &aDblSideSetIndex,
        std::string const &                               aSideSetName,
        bool                                              aCollectSets = true );

    /*!
     * @brief Create a block set from a single sided side set
     * @param[in] aDblSideSetIndex Double side set index
     * @param[in] aBlockSetName New block set name
     * @param[in] aCellTopo New block set cell topology
     * @return Block set index
     */
    moris_index
    create_block_set_from_cells_of_side_set(
        moris_index const &      aSideSetIndex,
        std::string const &      aBlockSetName,
        enum CellTopology const &aCellTopo );

  protected:

    //------------------------------------------------------------------------------
    // Parallel functions
    //------------------------------------------------------------------------------
    moris_id allocate_entity_ids( moris::size_t aNumReqs, enum EntityRank aEntityRank );

    //------------------------------------------------------------------------------

    void
    commit_double_side_set( moris_index const &aDoubleSideSetIndex );

    //------------------------------------------------------------------------------

    void
    commit_double_side_set( const moris::Cell< moris_index > & aDoubleSideSetIndexList );

    //------------------------------------------------------------------------------

    void
    commit_side_set( moris_index const &aSideSetIndex );

    //------------------------------------------------------------------------------

    void
    commit_side_set( const moris::Cell< moris_index > & aSideSetIndexList );

    //------------------------------------------------------------------------------

    void
    commit_block_set( moris_index const &aBlockSetIndex );

  private:
    //------------------------------------------------------------------------------
    
    uint
    get_num_owned_cells() const;
    
    //------------------------------------------------------------------------------

    /**
     * @brief Fill enriched integration mesh with Cell clusters, setting their active and void IG cells
     * 
     */
    void
    setup_cell_clusters();

    void
    setup_cell_clusters_new();

    //------------------------------------------------------------------------------

    /**
     * @brief Group clusters with the same active bulk-phase into "block-" sets
     * 
     */
    void
    setup_blockset_with_cell_clusters();

    //------------------------------------------------------------------------------

    /**
     * @brief Set the up side set clusters object
     * 
     */
    void
    setup_side_set_clusters();

    //------------------------------------------------------------------------------

    void
    setup_double_side_set_clusters();

    //------------------------------------------------------------------------------

    void
    setup_color_to_set();

    //------------------------------------------------------------------------------

    void
    setup_cluster_groups();

    void
    setup_cell_cluster_groups();
    
    void
    setup_side_cluster_groups();
    
    void
    setup_dbl_side_cluster_groups();

    //------------------------------------------------------------------------------

    void
    visualize_cluster_measures(); 

    //------------------------------------------------------------------------------
    
    void 
    visualize_cluster_group_measures();

    //------------------------------------------------------------------------------

    void
    setup_double_sided_interface_sides();

    //------------------------------------------------------------------------------

    void
    declare_interface_double_side_sets();

    //------------------------------------------------------------------------------

    moris_index
    get_dbl_side_set_index( moris_index aPhase0,
        moris_index                     aPhase1 );

    //------------------------------------------------------------------------------

    void
    create_interface_double_side_sets_and_clusters();

    //------------------------------------------------------------------------------

    void
    add_side_to_cluster( std::shared_ptr< xtk::Side_Cluster > aSideCluster,
        moris_index                                           aCellIndex,
        moris_index                                           aSideOrdinal );

    //------------------------------------------------------------------------------

    moris::Cell< std::string >
    split_set_name_by_bulk_phase( std::string aBaseName );

    //------------------------------------------------------------------------------

    moris::Cell< std::string >
    split_set_name_by_child_no_child( std::string aBaseName );

    //------------------------------------------------------------------------------

    Cell< moris_index >
    register_vertex_set_names( moris::Cell< std::string > const &aVertexSetNames );

    //------------------------------------------------------------------------------

    Cell< moris_index >
    register_block_set_names_with_cell_topo( moris::Cell< std::string > const &aBlockSetNames,
        enum CellTopology                                                      aBlockTopology );

    //------------------------------------------------------------------------------

    void
    set_block_set_colors( moris_index const &aBlockSetIndex,
        Matrix< IndexMat > const &           aBlockSetColors );

    //------------------------------------------------------------------------------

    Cell< moris_index >
    register_side_set_names( moris::Cell< std::string > const &aSideSetNames );

    //------------------------------------------------------------------------------

    void
    set_side_set_colors( moris_index const &aSideSetIndex,
        Matrix< IndexMat > const &          aSideSetColors );

    //------------------------------------------------------------------------------

    Cell< moris_index >
    register_double_side_set_names( moris::Cell< std::string > const &aDblSideSetNames );

    //------------------------------------------------------------------------------

    void
    set_double_side_set_colors( moris_index const &aDblSideSetIndex,
        Matrix< IndexMat > const &                 aMasterSideColors,
        Matrix< IndexMat > const &                 aSlaveSideColors );

    //------------------------------------------------------------------------------

    void
    setup_interface_side_sets();

    //------------------------------------------------------------------------------

    void
    declare_interface_side_sets();

    //------------------------------------------------------------------------------

    void
    create_interface_side_sets_and_clusters();

    //------------------------------------------------------------------------------

    Cell< moris_index >
    declare_interface_vertex_sets();

    //------------------------------------------------------------------------------

    void
    set_vertex_set_color( moris_index const &aVertexSetIndex,
        Matrix< IndexMat > const &           aVertexSetColors );

    //------------------------------------------------------------------------------

    void
    construct_color_to_set_relationship( moris::Cell< moris::Matrix< IndexMat > > const &aSetColors,
        moris::Cell< moris::Cell< moris_index > > &                                      aColorToSetIndex );

    //------------------------------------------------------------------------------

    void
    create_interface_side_sets_from_interface_double_side_set( moris_index const &aBulkphase0,
        moris_index const &                                                       aBulkphase1 );

    //------------------------------------------------------------------------------
    // Internal Additional Field Functions
    //------------------------------------------------------------------------------
    /*
     * Returns an index in the data structure for a given entity rank (i.e. NODE = 0)
     */
    moris_index
    get_entity_rank_field_index( enum moris::EntityRank aEntityRank );
    //------------------------------------------------------------------------------
    /*!
     * @return  whether a field exists or not
     */
    bool
    field_exists( std::string  aLabel,
        enum moris::EntityRank aEntityRank );
};

}// namespace xtk


#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTEGRATION_MESH_HPP_ */
