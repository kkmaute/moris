/*
 * cl_XTK_Enriched_Integration_Mesh.cpp
 *
 *  Created on: Jul 22, 2019
 *      Author: doble
 */

#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "fn_isempty.hpp"
#include "cl_TOL_Memory_Map.hpp"
#include <memory>
#include "cl_Logger.hpp"

namespace xtk
{
//------------------------------------------------------------------------------
Enriched_Integration_Mesh::Enriched_Integration_Mesh( Model *aXTKModel,
    moris::moris_index                                       aInterpIndex ) :
    mModel( aXTKModel ),
    mCutIgMesh( mModel->get_cut_integration_mesh() ), mMeshIndexInModel( aInterpIndex ), mCellClusters( 0, nullptr ), mFields( 0 ), mFieldLabelToIndex( 2 ), mCellInfo( nullptr )
{
    Tracer tTracer( "XTK", "Enriched Integration Mesh", "Construction" );
    this->setup_cell_clusters();
    this->setup_blockset_with_cell_clusters();
    this->setup_side_set_clusters();
    this->setup_double_side_set_clusters();
    this->setup_interface_side_sets();
    // this->setup_interface_vertex_sets();
    this->setup_color_to_set();
    this->collect_all_sets();


    if ( this->get_spatial_dim() == 2 )
    {
        mCellInfo = new moris::mtk::Cell_Info_Quad4();
    }
    else if ( this->get_spatial_dim() == 3 )
    {
        mCellInfo = new moris::mtk::Cell_Info_Hex8();
    }
}

//------------------------------------------------------------------------------

Enriched_Integration_Mesh::~Enriched_Integration_Mesh()
{
    delete mCellInfo;

    for ( auto p : mListofBlocks )
    {
        delete p;
    }
    mListofBlocks.clear();

    for ( auto p : mListofSideSets )
    {
        delete p;
    }
    mListofSideSets.clear();

    for ( auto p : mListofDoubleSideSets )
    {
        delete p;
    }
    mListofDoubleSideSets.clear();

    for ( auto p : mDoubleSideClusters )
    {
        delete p;
    }
    mDoubleSideClusters.clear();

    mDoubleSideSingleSideClusters.clear();
}

//------------------------------------------------------------------------------

MeshType
Enriched_Integration_Mesh::get_mesh_type() const
{
    return MeshType::XTK;
}

//------------------------------------------------------------------------------

moris::uint
Enriched_Integration_Mesh::get_spatial_dim() const
{
    return mCutIgMesh->get_spatial_dim();
}

//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_num_entities( enum EntityRank aEntityRank, const moris_index aIndex ) const
{
    switch ( aEntityRank )
    {
    case EntityRank::NODE:
    {
        return mCutIgMesh->get_num_entities( EntityRank::NODE, 0 );
        break;
    }
    case EntityRank::ELEMENT:
    {
        return mCutIgMesh->get_num_entities( EntityRank::ELEMENT, 0 );
        break;
    }
    default:
    {
        MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
    }
        return 0;
    }
}

//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_num_owned_cells() const
{
    uint tNumEntities = this->get_num_entities( EntityRank::ELEMENT );

    uint tNumOwnedEntities = 0;

    moris::moris_id tParRank = moris::par_rank();

    // iterate and find out how many I own
    for ( moris::uint i = 0; i < tNumEntities; i++ )
    {
        mtk::Cell const &tCell = this->get_mtk_cell( (moris_index)i );
        if ( tCell.get_owner() == tParRank )
        {
            tNumOwnedEntities++;
        }
    }

    return tNumOwnedEntities;
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::get_entity_connected_to_entity_loc_inds(
    moris_index       aEntityIndex,
    enum EntityRank   aInputEntityRank,
    enum EntityRank   aOutputEntityRank,
    const moris_index aIndex ) const
{
    MORIS_ERROR( aInputEntityRank == EntityRank::ELEMENT && aOutputEntityRank == EntityRank::NODE,
        "Only support element to node connectivity" );

    return this->get_mtk_cell( aEntityIndex ).get_vertex_inds();
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
{
    MORIS_ERROR( 0,
        "XTK ENRICHED MESH ERROR: get_elements_connected_to_element_and_face_ind_loc_inds no implemented" );

    return Matrix< IndexMat >( 0, 0 );
}

//------------------------------------------------------------------------------

Cell< mtk::Vertex const * >
Enriched_Integration_Mesh::get_all_vertices() const
{
    moris::uint                 tNumNodes = this->get_num_entities( EntityRank::NODE );
    Cell< mtk::Vertex const * > tVertices( tNumNodes );

    for ( moris::uint i = 0; i < tNumNodes; i++ )
    {
        tVertices( i ) = &mCutIgMesh->get_mtk_vertex( i );
    }
    return tVertices;
}

//------------------------------------------------------------------------------

moris_id
Enriched_Integration_Mesh::get_glb_entity_id_from_entity_loc_index(
    moris_index       aEntityIndex,
    enum EntityRank   aEntityRank,
    const moris_index aIndex ) const
{
    return mCutIgMesh->get_glb_entity_id_from_entity_loc_index( aEntityIndex, aEntityRank );
}

//------------------------------------------------------------------------------

std::unordered_map< moris_id, moris_index >
Enriched_Integration_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
{
    return mCutIgMesh->get_vertex_glb_id_to_loc_vertex_ind_map();
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::get_loc_entity_ind_from_entity_glb_id(
    moris_id          aEntityId,
    enum EntityRank   aEntityRank,
    const moris_index aIndex ) const
{
    return mCutIgMesh->get_loc_entity_ind_from_entity_glb_id( aEntityId, aEntityRank );
}

//------------------------------------------------------------------------------

Matrix< IdMat >
Enriched_Integration_Mesh::get_entity_connected_to_entity_glob_ids(
    moris_id          aEntityId,
    enum EntityRank   aInputEntityRank,
    enum EntityRank   aOutputEntityRank,
    const moris_index aIndex ) const
{
    moris_index tEntityIndex = get_loc_entity_ind_from_entity_glb_id( aEntityId, aInputEntityRank );

    Matrix< IndexMat > tEntityToEntityLoc = this->get_entity_connected_to_entity_loc_inds( tEntityIndex, aInputEntityRank, aOutputEntityRank );

    return convert_indices_to_ids( tEntityToEntityLoc, aOutputEntityRank );
}

//------------------------------------------------------------------------------

Matrix< DDRMat >
Enriched_Integration_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
{
    mtk::Vertex const &tVertex = get_mtk_vertex( aNodeIndex );
    return tVertex.get_coords();
}

//------------------------------------------------------------------------------

mtk::Vertex &
Enriched_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex )
{
    return mCutIgMesh->get_mtk_vertex( aVertexIndex );
}

//------------------------------------------------------------------------------

mtk::Vertex const &
Enriched_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
{
    return mCutIgMesh->get_mtk_vertex( aVertexIndex );
}

//------------------------------------------------------------------------------

mtk::Cell &
Enriched_Integration_Mesh::get_writable_mtk_cell( moris_index aElementIndex )
{
    return mCutIgMesh->get_mtk_cell( aElementIndex );
}

//------------------------------------------------------------------------------

mtk::Cell &
Enriched_Integration_Mesh::get_mtk_cell( moris_index aElementIndex )
{
    return mCutIgMesh->get_mtk_cell( aElementIndex );
}

//------------------------------------------------------------------------------

mtk::Cell const &
Enriched_Integration_Mesh::get_mtk_cell( moris_index aElementIndex ) const
{
    return mCutIgMesh->get_mtk_cell( aElementIndex );
}

//------------------------------------------------------------------------------

Matrix< IdMat >
Enriched_Integration_Mesh::get_communication_table() const
{
    return mCutIgMesh->get_communication_table();
}

//------------------------------------------------------------------------------

moris::Cell< std::string >
Enriched_Integration_Mesh::get_set_names( enum EntityRank aSetEntityRank ) const
{
    switch ( aSetEntityRank )
    {
    case EntityRank::NODE:
    {
        return mVertexSetNames;
        break;
    }
    case EntityRank::EDGE:
    {
        MORIS_ASSERT( this->get_facet_rank() == EntityRank::EDGE, "side sets are defined on edges in 2d" );
        return mSideSetLabels;
        break;
    }
    case EntityRank::FACE:
    {
        MORIS_ASSERT( this->get_facet_rank() == EntityRank::FACE, "side sets are defined on faces in 3d" );
        return mSideSetLabels;
        break;
    }
    case EntityRank::ELEMENT:
    {
        return mBlockSetNames;
        break;
    }
    default:
    {
        MORIS_ERROR( 0, "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
    }
        return moris::Cell< std::string >( 0 );
        break;
    }
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::get_set_entity_loc_inds(
    enum EntityRank aSetEntityRank,
    std::string     aSetName ) const
{
    switch ( aSetEntityRank )
    {
    case EntityRank::NODE:
    {
        // get the vertex set index
        auto tSetIndex = mVertexSetLabelToOrd.find( aSetName );

        moris::Cell< moris::mtk::Vertex * > tVerticesInSet = mVerticesInVertexSet( tSetIndex->second );
        Matrix< IndexMat >                  tVerticesInSetMat( 1, tVerticesInSet.size() );
        for ( moris::uint i = 0; i < tVerticesInSet.size(); i++ )
        {
            tVerticesInSetMat( i ) = tVerticesInSet( i )->get_index();
        }

        return tVerticesInSetMat;
        break;
    }
    case EntityRank::EDGE:
    {
        MORIS_ASSERT( this->get_facet_rank() == EntityRank::EDGE, "side sets are defined on edges in 2d" );
        return Matrix< IndexMat >( 0, 0 );
        break;
    }
    case EntityRank::FACE:
    {
        MORIS_ASSERT( this->get_facet_rank() == EntityRank::FACE, "side sets are defined on faces in 3d" );
        return Matrix< IndexMat >( 0, 0 );
        break;
    }
    case EntityRank::ELEMENT:
    {
        return this->get_block_entity_loc_inds( aSetName );
        break;
    }
    default:
    {
        MORIS_ERROR( 0, "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
        return Matrix< IndexMat >( 0, 0 );
        break;
    }
    }
}

// ----------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::get_element_indices_in_block_set( uint aSetIndex )
{
    return this->get_block_entity_loc_inds( this->get_set_names( EntityRank::ELEMENT )( aSetIndex ) );
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::get_sideset_elems_loc_inds_and_ords(
    const std::string & aSetName,
    Matrix< IndexMat > &aElemIndices,
    Matrix< IndexMat > &aSidesetOrdinals ) const
{
    // get the index
    moris_index tSideSetIndex = this->get_side_set_index( aSetName );

    // get the cell clusters
    moris::Cell< mtk::Cluster const * > tSideClusters = this->get_side_set_cluster( tSideSetIndex );


    // iterate through side clusters and count number of sides in set
    moris::uint tNumSides = 0;
    for ( auto iCluster : tSideClusters )
    {
        tNumSides = tNumSides + iCluster->get_num_primary_cells();
    }

    // size outputs
    aElemIndices.resize( 1, tNumSides );
    aSidesetOrdinals.resize( 1, tNumSides );

    // reset count
    tNumSides = 0;
    for ( auto iCluster : tSideClusters )
    {
        Matrix< IndexMat > tSideOrdinals = iCluster->get_cell_side_ordinals();
        Matrix< IndexMat > tCellIndices  = iCluster->get_primary_cell_indices_in_cluster();

        aElemIndices( { 0, 0 }, { tNumSides, tNumSides + tCellIndices.numel() - 1 } )     = tCellIndices.matrix_data();
        aSidesetOrdinals( { 0, 0 }, { tNumSides, tNumSides + tCellIndices.numel() - 1 } ) = tSideOrdinals.matrix_data();

        tNumSides = tNumSides + tSideOrdinals.numel();
    }
}

//------------------------------------------------------------------------------

moris_id
Enriched_Integration_Mesh::get_max_entity_id( enum EntityRank aEntityRank, const moris_index aIndex ) const
{
    MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only Elements or Nodes have max id" );

    moris::uint tNumEntities = this->get_num_entities( aEntityRank );

    moris_id tLocalMaxId = 0;

    for ( moris::uint i = 0; i < tNumEntities; i++ )
    {
        moris_id tId = this->get_glb_entity_id_from_entity_loc_index( i, aEntityRank );

        if ( tId > tLocalMaxId )
        {
            tLocalMaxId = tId;
        }
    }

    moris_id tGlobalMaxId = moris::max_all( tLocalMaxId );
    return tGlobalMaxId;
}

//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_node_owner( moris_index aNodeIndex ) const
{
    return this->get_mtk_vertex( aNodeIndex ).get_owner();
}

//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_element_owner( moris_index aElementIndex ) const
{
    return this->get_mtk_cell( aElementIndex ).get_owner();
}

//------------------------------------------------------------------------------


Matrix< IndexMat >
Enriched_Integration_Mesh::get_block_entity_loc_inds( std::string aSetName ) const
{
    // ge tindex
    moris_index tBlockIndex = this->get_block_set_index( aSetName );

    // get clusters in block
    moris::Cell< mtk::Cluster const * > tCellClusters = this->get_cell_clusters_in_set( tBlockIndex );

    // iterate through clusters and count all primary integration cells
    moris::uint tCount = 0;
    for ( auto it : tCellClusters )
    {
        tCount = tCount + it->get_num_primary_cells();
    }

    // allocate output
    Matrix< IndexMat > tCellIndices( 1, tCount );

    // reset count to use it for something else
    tCount = 0;

    // iterate through clusters and collect all primary integration cell indices
    for ( auto it : tCellClusters )
    {
        // get primary cell clusters
        moris::Matrix< moris::IndexMat > tPrimaryCellIndices = it->get_primary_cell_indices_in_cluster();

        tCellIndices( { 0, 0 }, { tCount, tCount + tPrimaryCellIndices.numel() - 1 } ) = tPrimaryCellIndices.matrix_data();

        tCount = tCount + tPrimaryCellIndices.numel();
    }

    return tCellIndices;
}


void
Enriched_Integration_Mesh::create_dbl_sided_interface_set(
    moris_index aMasterBulkPhaseIndex,
    moris_index aSlaveBulkPhaseIndex )
{
    MORIS_ERROR( aMasterBulkPhaseIndex > aSlaveBulkPhaseIndex, "The master bulk phase needs to be lower than the slave bulk phase." );

    // get the name of this side set
    std::string tInterfaceDblSideName = this->get_dbl_interface_side_set_name( aMasterBulkPhaseIndex, aSlaveBulkPhaseIndex );

    // register it with the mesh
    Cell< moris_index > tDblSideSetOrds = this->register_double_side_set_names( { tInterfaceDblSideName } );

    // get the other
    moris_index tOtherInterfaceIndex = this->get_dbl_side_set_index( aSlaveBulkPhaseIndex, aMasterBulkPhaseIndex );

    // set the colors
    Matrix< IndexMat > tMasterColor = { { aMasterBulkPhaseIndex } };
    Matrix< IndexMat > tSlaveColor  = { { aSlaveBulkPhaseIndex } };
    this->set_double_side_set_colors( tDblSideSetOrds( 0 ), tMasterColor, tSlaveColor );

    // resize member data
    moris::uint tNumPairsInSet = mDoubleSideSetsMasterIndex( tOtherInterfaceIndex ).size();
    mDoubleSideSetsMasterIndex( tDblSideSetOrds( 0 ) ).resize( tNumPairsInSet );
    mDoubleSideSetsSlaveIndex( tDblSideSetOrds( 0 ) ).resize( tNumPairsInSet );

    for ( moris::uint i = 0; i < tNumPairsInSet; i++ )
    {
        // master is slave slave is master
        mDoubleSideSetsMasterIndex( tDblSideSetOrds( 0 ) )( i ) = mDoubleSideSetsSlaveIndex( tOtherInterfaceIndex )( i );
        mDoubleSideSetsSlaveIndex( tDblSideSetOrds( 0 ) )( i )  = mDoubleSideSetsMasterIndex( tOtherInterfaceIndex )( i );

        // get master ans slave clusters
        Side_Cluster *tMasterSideCluster = mDoubleSideSingleSideClusters( mDoubleSideSetsMasterIndex( tDblSideSetOrds( 0 ) )( i ) ).get();
        Side_Cluster *tSlaveSideCluster  = mDoubleSideSingleSideClusters( mDoubleSideSetsSlaveIndex( tDblSideSetOrds( 0 ) )( i ) ).get();

        // create double side set
        mtk::Double_Side_Cluster *tDblSideCluster = new mtk::Double_Side_Cluster( tMasterSideCluster, tSlaveSideCluster, tMasterSideCluster->get_vertices_in_cluster() );

        mDoubleSideClusters.push_back( tDblSideCluster );
        mDoubleSideSets( tDblSideSetOrds( 0 ) ).push_back( tDblSideCluster );
    }

    this->setup_color_to_set();
    this->commit_double_side_set( tDblSideSetOrds( 0 ) );
    this->collect_all_sets();
}


//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::deactivate_empty_sets()
{
    this->deactivate_empty_block_sets();
    this->deactivate_empty_side_sets();
    this->setup_color_to_set();
    this->collect_all_sets();
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::deactivate_empty_side_sets()
{
    // copy old data
    std::unordered_map< std::string, moris_index >                     tOldSetMap      = mSideSideSetLabelToOrd;
    moris::Cell< std::string >                                         tOldSetNames    = mSideSetLabels;
    moris::Cell< moris::Cell< std::shared_ptr< xtk::Side_Cluster > > > tOldSetClusters = mSideSets;

    // clear member data
    mSideSideSetLabelToOrd.clear();
    mSideSetLabels.clear();
    mSideSets.clear();

    for ( auto iB : mListofSideSets )
    {
        delete iB;
    }
    mListofSideSets.clear();

    // current index
    moris_index tSetIndex = 0;

    for ( moris::uint i = 0; i < tOldSetClusters.size(); i++ )
    {
        uint tMySize  = tOldSetClusters( i ).size();
        uint tAllSize = sum_all( tMySize );

        if ( tAllSize > 0 )
        {
            mSideSetLabels.push_back( tOldSetNames( i ) );
            mSideSets.push_back( tOldSetClusters( i ) );

            MORIS_ASSERT( mSideSideSetLabelToOrd.find( tOldSetNames( i ) ) == mSideSideSetLabelToOrd.end(), "Duplicate block set in mesh" );
            mSideSideSetLabelToOrd[tOldSetNames( i )] = tSetIndex;
            this->commit_side_set( tSetIndex );
            tSetIndex++;
        }
    }

    this->setup_color_to_set();
    this->collect_all_sets();
}

void
Enriched_Integration_Mesh::add_mesh_field_real_scalar_data_loc_inds(
    std::string const &     aFieldName,
    enum EntityRank const & aFieldEntityRank,
    Matrix< DDRMat > const &aFieldData )
{

    MORIS_ASSERT( aFieldEntityRank == EntityRank::ELEMENT, "Only tested for nodal and element scalar field" );

    moris_index tFieldIndex = this->create_field( aFieldName, EntityRank::ELEMENT, 0 );

    this->add_field_data( tFieldIndex, EntityRank::ELEMENT, aFieldData );
}

void
Enriched_Integration_Mesh::create_cell_id_fields()
{
    // Fields constructed here
    moris::Cell< std::string > tCellFields = { "cell_id" };

    moris_index tFieldIndex = this->create_field( tCellFields( 0 ), EntityRank::ELEMENT, 0 );

    moris::Matrix< moris::DDRMat > tCellIdField( 1, this->get_num_elems() );

    for ( moris::uint iCell = 0; iCell < this->get_num_elems(); iCell++ )
    {
        tCellIdField( iCell ) = (moris::real)this->get_mtk_cell( iCell ).get_id();
    }

    this->add_field_data( tFieldIndex, EntityRank::ELEMENT, tCellIdField );
}


//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::deactivate_empty_block_sets()
{
    std::unordered_map< std::string, moris_index >          tOldSetMap      = mBlockSetLabelToOrd;
    moris::Cell< std::string >                              tOldSetNames    = mBlockSetNames;
    moris::Cell< enum CellTopology >                        tOldSetTopo     = mBlockSetTopology;
    moris::Cell< moris::Cell< xtk::Cell_Cluster const * > > tOldSetClusters = mPrimaryBlockSetClusters;
    moris::Cell< moris::Matrix< IndexMat > >                tOldColors      = mBlockSetColors;

    // clear member data
    mBlockSetLabelToOrd.clear();
    mBlockSetNames.clear();
    mBlockSetTopology.clear();
    mPrimaryBlockSetClusters.clear();
    mBlockSetColors.clear();

    for ( auto iB : mListofBlocks )
    {
        delete iB;
    }
    mListofBlocks.clear();

    // current index
    moris_index tSetIndex = 0;

    for ( moris::uint i = 0; i < tOldSetClusters.size(); i++ )
    {
        uint tMySize  = tOldSetClusters( i ).size();
        uint tAllSize = sum_all( tMySize );

        if ( tAllSize > 0 )
        {
            mBlockSetNames.push_back( tOldSetNames( i ) );
            mPrimaryBlockSetClusters.push_back( tOldSetClusters( i ) );
            mBlockSetTopology.push_back( tOldSetTopo( i ) );
            mBlockSetColors.push_back( tOldColors( i ) );

            MORIS_ASSERT( mBlockSetLabelToOrd.find( tOldSetNames( i ) ) == mBlockSetLabelToOrd.end(), "Duplicate block set in mesh" );
            mBlockSetLabelToOrd[tOldSetNames( i )] = tSetIndex;
            this->commit_block_set( tSetIndex );
            tSetIndex++;
        }
    }

    this->setup_color_to_set();
    this->collect_all_sets();
}

//------------------------------------------------------------------------------

moris::Cell< std::string >
Enriched_Integration_Mesh::create_basis_support_fields()
{
    // get the enriched interpolation mesh
    Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( mMeshIndexInModel );

    // base string of field
    std::string tBaseStr = "Weight";

    // field names for output
    moris::Cell< std::string > tOutputFieldNames;

    // field information for internal use
    moris::Cell< moris::Cell< std::string > >      tFieldNames( tEnrInterpMesh->get_num_interpolation_types() );
    moris::Cell< moris::Cell< moris_index > >      tFieldIndices( tEnrInterpMesh->get_num_interpolation_types() );
    moris::Cell< moris::Cell< Matrix< DDRMat > > > tFieldData( tEnrInterpMesh->get_num_interpolation_types() );

    // iterate through interpolation types and for each basis declare the field in mesh
    for ( moris::uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
    {
        std::string tInterpTypeStr = "_mi_" + std::to_string( iBT );

        tFieldNames( iBT ).resize( tEnrInterpMesh->get_num_basis_functions() );
        tFieldIndices( iBT ).resize( tEnrInterpMesh->get_num_basis_functions() );
        tFieldData( iBT ).resize( tEnrInterpMesh->get_num_basis_functions(), Matrix< DDRMat >( this->get_num_entities( EntityRank::NODE ), 1 ) );

        // iterate through basis functions
        for ( moris::uint iB = 0; iB < tEnrInterpMesh->get_num_basis_functions(); iB++ )
        {
            // initialize the data to -1
            tFieldData( iBT )( iB ).fill( -1 );

            tFieldNames( iBT )( iB ) = tBaseStr + tInterpTypeStr + "_ind_" + std::to_string( iB );

            tOutputFieldNames.push_back( tFieldNames( iBT )( iB ) );

            // declare the field in this mesh
            tFieldIndices( iBT )( iB ) = this->create_field( tFieldNames( iBT )( iB ), EntityRank::NODE, 0 );
        }
    }

    // populate field data
    for ( moris::uint iCl = 0; iCl < this->mCellClusters.size(); iCl++ )
    {
        xtk::Cell_Cluster *tCluster = mCellClusters( iCl ).get();

        // get the interpolation cell
        moris::mtk::Cell const &tInterpCell = tCluster->get_interpolation_cell();

        // get the vertices attached to this cell
        moris::Cell< moris::mtk::Vertex * > tVertices = tInterpCell.get_vertex_pointers();

        // allocate place to put coefficients interpolating into these vertices
        moris::Cell< moris::Cell< moris_index > > tCoeffsIPIntoCluster( tEnrInterpMesh->get_num_interpolation_types() );

        // collect coefficients of this interpolation cell
        for ( moris::uint iV = 0; iV < tVertices.size(); iV++ )
        {
            for ( moris::uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
            {
                mtk::Vertex_Interpolation *tVertexIp = tVertices( iV )->get_interpolation( tEnrInterpMesh->get_interpolation_index( (moris_index)iBT ) );

                // get indices of coefficients
                Matrix< IndexMat > tCoeffInds = tVertexIp->get_indices();

                // resize data
                tCoeffsIPIntoCluster( iBT ).resize( tCoeffInds.numel() );

                for ( moris::uint iC = 0; iC < tCoeffInds.numel(); iC++ )
                {
                    tCoeffsIPIntoCluster( iBT )( iC ) = tCoeffInds( iC );
                }
            }
        }

        // get primary cells from the cluster
        moris::Cell< moris::mtk::Cell const * > const &tPrimaryCells = tCluster->get_primary_cells_in_cluster();

        // iterate through primary cells
        for ( moris::uint iCell = 0; iCell < tPrimaryCells.size(); iCell++ )
        {
            // get vertices attached to primary cells
            moris::Cell< moris::mtk::Vertex * > tVertices = tPrimaryCells( iCell )->get_vertex_pointers();

            // iterate through vertices and mark them as in support of all coefficients in tCoeffsIPIntoCluster
            for ( moris::uint iV = 0; iV < tVertices.size(); iV++ )
            {
                for ( moris::uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
                {
                    for ( moris::uint iC = 0; iC < tCoeffsIPIntoCluster( iBT ).size(); iC++ )
                    {
                        tFieldData( iBT )( tCoeffsIPIntoCluster( iBT )( iC ) )( tVertices( iV )->get_index() ) = 1;
                    }
                }
            }
        }
    }

    // add field data to mesh
    // iterate through interpolation
    for ( moris::uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
    {
        // iterate through basis functions
        for ( moris::uint iB = 0; iB < tEnrInterpMesh->get_num_basis_functions(); iB++ )
        {
            this->add_field_data( tFieldIndices( iBT )( iB ), EntityRank::NODE, tFieldData( iBT )( iB ) );
        }
    }

    return tOutputFieldNames;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::write_mesh( moris::ParameterList *aParamList )
{
    if ( aParamList->get< bool >( "deactivate_empty_sets" ) )
    {
        this->deactivate_empty_sets();
    }

    // get path to output XTK files to
    std::string tOutputPath = aParamList->get< std::string >( "output_path" );
    std::string tOutputFile = aParamList->get< std::string >( "output_file" );
    std::string tOutputBase = tOutputFile.substr( 0, tOutputFile.find( "." ) );
    std::string tOutputExt  = tOutputFile.substr( tOutputFile.find( "." ), tOutputFile.length() );

    MORIS_ASSERT( tOutputExt == ".exo" || tOutputExt == ".e", "Invalid file extension, needs to be .exo or .e" );

    // Write mesh
    moris::mtk::Writer_Exodus writer( this );

    // if user requests to keep XTK output for all iterations, add iteration count to output file name
    if ( aParamList->get< bool >( "keep_all_opt_iters" ) )
    {
        // get optimization iteration ( function returns zero if no optimization )
        uint tOptIter = gLogger.get_opt_iteration();

        writer.write_mesh(
            "", tOutputPath + tOutputBase + "." + std::to_string( tOptIter ) + tOutputExt, "", tOutputPath + "xtk_temp2." + std::to_string( tOptIter ) + tOutputExt );
    }
    // else: proceed as usual and overwrite xtk_temp.exo each iteration
    else
    {
        writer.write_mesh(
            "", tOutputPath + tOutputFile, "", tOutputPath + "xtk_temp2.exo" );
    }

    if ( aParamList->get< bool >( "write_enrichment_fields" ) )
    {
        // set up the nodal fields for basis support
        moris::Cell< std::string > tNodeFields = this->create_basis_support_fields();

        writer.set_nodal_fields( tNodeFields );

        for ( moris::uint iF = 0; iF < tNodeFields.size(); iF++ )
        {
            moris::moris_index tFieldIndex = this->get_field_index( tNodeFields( iF ), EntityRank::NODE );

            writer.write_nodal_field( tNodeFields( iF ), this->get_field_data( tFieldIndex, EntityRank::NODE ) );
        }
    }

    // Write the fields
    writer.set_time( 0.0 );
    writer.close_file();
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::create_union_block( Cell< std::string > const &aBlocks,
    std::string                                                           aNewBlock,
    Matrix< IndexMat > const &                                            aNewBlockColor )
{
    MORIS_ERROR( aBlocks.size() >= 2, "Union needs to happen between two blocks or more" );

    enum CellTopology tCellTopo = CellTopology::INVALID;

    moris::uint         tCount = 0;
    Cell< moris_index > tBlockIndices( aBlocks.size() );
    for ( moris::uint i = 0; i < aBlocks.size(); i++ )
    {
        tBlockIndices( i ) = this->get_block_set_index( aBlocks( i ) );
        tCount             = tCount + mPrimaryBlockSetClusters( tBlockIndices( i ) ).size();

        moris::mtk::Set *tSet = this->get_set_by_index( this->get_set_index_by_name( aBlocks( i ) ) );
        if ( i == 0 )
        {
            tCellTopo = tSet->get_cell_topology();
        }
        else
        {
            MORIS_ERROR( tCellTopo == tSet->get_cell_topology(), "Invalid merge detected, verify that all blocks have the same cell topology" );
        }
    }

    Cell< moris_index > tBlockSetIndex = this->register_block_set_names_with_cell_topo( { aNewBlock }, tCellTopo );
    mPrimaryBlockSetClusters( tBlockSetIndex( 0 ) ).reserve( tCount );

    for ( moris::uint i = 0; i < aBlocks.size(); i++ )
    {
        mPrimaryBlockSetClusters( tBlockSetIndex( 0 ) ).append( mPrimaryBlockSetClusters( tBlockIndices( i ) ) );
    }

    // Check compatibility of union

    this->commit_block_set( tBlockSetIndex( 0 ) );
    this->set_block_set_colors( tBlockSetIndex( 0 ), aNewBlockColor );
    this->setup_color_to_set();
    this->collect_all_sets( false );
}
//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::create_union_side_set( Cell< std::string > const &aSideSets,
    std::string                                                              aNewSideSet,
    Matrix< IndexMat > const &                                               aNewSideSetColor )
{
    MORIS_ERROR( aSideSets.size() >= 2, "Union needs to happen between two side sets or more" );

    moris::uint         tCount = 0;
    Cell< moris_index > tSideSetIndices( aSideSets.size() );
    for ( moris::uint i = 0; i < aSideSets.size(); i++ )
    {
        tSideSetIndices( i ) = this->get_side_set_index( aSideSets( i ) );
        tCount               = tCount + mSideSets( tSideSetIndices( i ) ).size();
    }

    Cell< moris_index > tSideSetIndex = this->register_side_set_names( { aNewSideSet } );
    mSideSets( tSideSetIndex( 0 ) ).reserve( tCount );

    for ( moris::uint i = 0; i < aSideSets.size(); i++ )
    {
        mSideSets( tSideSetIndex( 0 ) ).append( mSideSets( tSideSetIndices( i ) ) );
    }

    this->commit_side_set( tSideSetIndex( 0 ) );

    this->set_side_set_colors( tSideSetIndex( 0 ), aNewSideSetColor );

    this->setup_color_to_set();
    this->collect_all_sets( false );
}

void
Enriched_Integration_Mesh::deactive_all_blocks_but_selected( Cell< std::string > const &aBlockSetsToKeep )
{
    std::unordered_map< std::string, moris_index > tBlocksToKeepMap;
    for ( moris::uint i = 0; i < aBlockSetsToKeep.size(); i++ )
    {
        tBlocksToKeepMap[aBlockSetsToKeep( i )] = 1;
    }


    std::unordered_map< std::string, moris_index >          tOldSetMap      = mBlockSetLabelToOrd;
    moris::Cell< std::string >                              tOldSetNames    = mBlockSetNames;
    moris::Cell< enum CellTopology >                        tOldSetTopo     = mBlockSetTopology;
    moris::Cell< moris::Cell< xtk::Cell_Cluster const * > > tOldSetClusters = mPrimaryBlockSetClusters;
    moris::Cell< moris::Matrix< IndexMat > >                tOldColors      = mBlockSetColors;

    // clear member data
    mBlockSetLabelToOrd.clear();
    mBlockSetNames.clear();
    mBlockSetTopology.clear();
    mPrimaryBlockSetClusters.clear();
    mBlockSetColors.clear();

    for ( auto iB : mListofBlocks )
    {
        delete iB;
    }
    mListofBlocks.clear();

    // current index
    moris_index tSetIndex = 0;

    for ( moris::uint i = 0; i < tOldSetClusters.size(); i++ )
    {
        if ( tBlocksToKeepMap.find( tOldSetNames( i ) ) != tBlocksToKeepMap.end() )
        {
            mBlockSetNames.push_back( tOldSetNames( i ) );
            mPrimaryBlockSetClusters.push_back( tOldSetClusters( i ) );
            mBlockSetTopology.push_back( tOldSetTopo( i ) );
            mBlockSetColors.push_back( tOldColors( i ) );

            MORIS_ASSERT( mBlockSetLabelToOrd.find( tOldSetNames( i ) ) == mBlockSetLabelToOrd.end(), "Duplicate block set in mesh" );
            mBlockSetLabelToOrd[tOldSetNames( i )] = tSetIndex;
            this->commit_block_set( tSetIndex );
            tSetIndex++;
        }
    }
}
void
Enriched_Integration_Mesh::deactive_all_side_sets_but_selected( Cell< std::string > const &aSideSetsToKeep )
{
    std::unordered_map< std::string, moris_index > tSideSetsToKeepMap;
    for ( moris::uint i = 0; i < aSideSetsToKeep.size(); i++ )
    {
        tSideSetsToKeepMap[aSideSetsToKeep( i )] = 1;
    }


    // copy old data
    std::unordered_map< std::string, moris_index >                     tOldSetMap      = mSideSideSetLabelToOrd;
    moris::Cell< std::string >                                         tOldSetNames    = mSideSetLabels;
    moris::Cell< moris::Cell< std::shared_ptr< xtk::Side_Cluster > > > tOldSetClusters = mSideSets;

    // clear member data
    mSideSideSetLabelToOrd.clear();
    mSideSetLabels.clear();
    mSideSets.clear();

    for ( auto iB : mListofSideSets )
    {
        delete iB;
    }
    mListofSideSets.clear();

    // current index
    moris_index tSetIndex = 0;

    for ( moris::uint i = 0; i < tOldSetClusters.size(); i++ )
    {
        if ( tSideSetsToKeepMap.find( tOldSetNames( i ) ) != tSideSetsToKeepMap.end() )
        {
            mSideSetLabels.push_back( tOldSetNames( i ) );
            mSideSets.push_back( tOldSetClusters( i ) );

            MORIS_ASSERT( mSideSideSetLabelToOrd.find( tOldSetNames( i ) ) == mSideSideSetLabelToOrd.end(), "Duplicate block set in mesh" );
            mSideSideSetLabelToOrd[tOldSetNames( i )] = tSetIndex;
            this->commit_side_set( tSetIndex );
            tSetIndex++;
        }
    }
}

//------------------------------------------------------------------------------

moris::Memory_Map
Enriched_Integration_Mesh::get_memory_usage()
{
    // memory map of ig mesh
    moris::Memory_Map tMM;
    tMM.mMemoryMapData["mVertexSetNames"]            = moris::internal_capacity( mVertexSetNames );
    tMM.mMemoryMapData["mVerticesInVertexSet"]       = moris::internal_capacity( mVerticesInVertexSet );
    tMM.mMemoryMapData["mVertexSetColors"]           = moris::internal_capacity( mVertexSetColors );
    tMM.mMemoryMapData["mBlockSetNames"]             = moris::internal_capacity( mBlockSetNames );
    tMM.mMemoryMapData["mBlockSetTopology"]          = mBlockSetTopology.capacity();
    tMM.mMemoryMapData["mBlockSetNames"]             = moris::internal_capacity( mBlockSetNames );
    tMM.mMemoryMapData["mPrimaryBlockSetClusters"]   = moris::internal_capacity( mPrimaryBlockSetClusters );
    tMM.mMemoryMapData["mBlockSetColors"]            = moris::internal_capacity( mBlockSetColors );
    tMM.mMemoryMapData["mColorsBlockSets"]           = moris::internal_capacity( mColorsBlockSets );
    tMM.mMemoryMapData["mSideSetLabels"]             = moris::internal_capacity( mSideSetLabels );
    tMM.mMemoryMapData["mSideSets"]                  = moris::internal_capacity_nested_ptr( mSideSets );
    tMM.mMemoryMapData["mSideSetColors"]             = moris::internal_capacity( mSideSetColors );
    tMM.mMemoryMapData["mColorsSideSets"]            = moris::internal_capacity( mColorsSideSets );
    tMM.mMemoryMapData["mDoubleSideSetLabels"]       = moris::internal_capacity( mDoubleSideSetLabels );
    tMM.mMemoryMapData["mDoubleSideSets"]            = moris::internal_capacity( mDoubleSideSets );
    tMM.mMemoryMapData["mDoubleSideSetsMasterIndex"] = moris::internal_capacity( mDoubleSideSetsMasterIndex );
    tMM.mMemoryMapData["mDoubleSideSetsSlaveIndex"]  = moris::internal_capacity( mDoubleSideSetsSlaveIndex );
    //FIXME: Implement capacities down through MTK children
    //    tMM.mMemoryMapData["mDoubleSideClusters"] = moris::internal_capacity(mDoubleSideClusters);
    tMM.mMemoryMapData["mDoubleSideSingleSideClusters"] = moris::internal_capacity_ptr( mDoubleSideSingleSideClusters );
    tMM.mMemoryMapData["mBulkPhaseToDblSideIndex"]      = mBulkPhaseToDblSideIndex.capacity();
    tMM.mMemoryMapData["mMasterDoubleSideSetColor"]     = moris::internal_capacity( mMasterDoubleSideSetColor );
    tMM.mMemoryMapData["mSlaveDoubleSideSetColor"]      = moris::internal_capacity( mSlaveDoubleSideSetColor );
    tMM.mMemoryMapData["mColorMasterDoubleSideSet"]     = moris::internal_capacity( mColorMasterDoubleSideSet );
    tMM.mMemoryMapData["mColorSlaveDoubleSideSet"]      = moris::internal_capacity( mColorSlaveDoubleSideSet );
    return tMM;
}

//------------------------------------------------------------------------------

enum CellTopology
Enriched_Integration_Mesh::get_blockset_topology( const std::string &aSetName )
{
    moris_index tIndex = this->get_block_set_index( aSetName );
    return mBlockSetTopology( tIndex );
}

//------------------------------------------------------------------------------

enum CellShape
Enriched_Integration_Mesh::get_IG_blockset_shape( const std::string &aSetName )
{
    // get the clusters in the set
    moris::Cell< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

    // init cell shape
    CellShape tCellShape = CellShape::EMPTY;

    // if the set isn't empty exist
    if ( tSetClusters.size() > 0 )
    {
        // get the cells in the first cluster
        moris::Cell< moris::mtk::Cell const * > tClusterCells = tSetClusters( 0 )->get_primary_cells_in_cluster();

        // compute the cell shape of the first cell
        tCellShape = tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
    }

    // within debug, checking all cells to make sure that they are the same Cell Shape
    // if cells exist

    // looping through the clusters
    for ( uint iCluster = 0; iCluster < tSetClusters.size(); iCluster++ )
    {
        // get cell of cells in the cluster
        moris::Cell< moris::mtk::Cell const * > tClusterCellsCheck = tSetClusters( iCluster )->get_primary_cells_in_cluster();

        // looping through the cells in the cluster
        for ( uint iCheckCell = 0; iCheckCell < tClusterCellsCheck.size(); iCheckCell++ )
        {
            MORIS_ASSERT( tClusterCellsCheck( iCheckCell )->get_cell_info()->compute_cell_shape( tClusterCellsCheck( iCheckCell ) ) == tCellShape,
                "Mesh_Core_STK::get_IG_blockset_shape - cell shape is not consistent in the block" );
        }
    }

    return tCellShape;
}

//------------------------------------------------------------------------------

enum CellShape
Enriched_Integration_Mesh::get_IP_blockset_shape( const std::string &aSetName )
{
    // get the clusters in the set
    moris::Cell< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

    // init cell shape
    CellShape tCellShape = CellShape::EMPTY;

    // if the set isn't empty exist
    if ( tSetClusters.size() > 0 )
    {
        // get the cells in the first cluster
        mtk::Cell const &tClusterCell = tSetClusters( 0 )->get_interpolation_cell();

        // compute the cell shape of the first cell
        tCellShape = tClusterCell.get_cell_info()->compute_cell_shape( &tClusterCell );
    }

    // within debug, checking all cells to make sure that they are the same Cell Shape
    // if cells exist
    // looping through the clusters
    for ( uint iCluster = 1; iCluster < tSetClusters.size(); iCluster++ )
    {
        MORIS_ASSERT( tSetClusters( iCluster )->get_interpolation_cell().get_cell_info()->compute_cell_shape( &tSetClusters( iCluster )->get_interpolation_cell() )
                          == tCellShape,
            "Enriched_Integration_Mesh::get_IP_blockset_shape - cell shape is not consistent in the block" );
    }

    return tCellShape;
}

//------------------------------------------------------------------------------

Matrix< IdMat >
Enriched_Integration_Mesh::convert_indices_to_ids(
    Matrix< IndexMat > const &aIndices,
    enum EntityRank           aEntityRank ) const
{
    moris::uint     tNRow = aIndices.n_rows();
    moris::uint     tNCol = aIndices.n_cols();
    Matrix< IdMat > tIds( tNRow, tNCol );
    for ( moris::uint i = 0; i < tNRow; i++ )
    {
        for ( moris::uint j = 0; j < tNCol; j++ )
        {
            tIds( i, j ) = this->get_glb_entity_id_from_entity_loc_index( aIndices( i, j ), aEntityRank );
        }
    }
    return tIds;
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::convert_ids_to_indices(
    Matrix< IdMat > const &aIds,
    enum EntityRank        aEntityRank ) const
{
    moris::uint     tNRow = aIds.n_rows();
    moris::uint     tNCol = aIds.n_cols();
    Matrix< IdMat > tIndices( tNRow, tNCol );
    for ( moris::uint i = 0; i < tNRow; i++ )
    {
        for ( moris::uint j = 0; j < tNCol; j++ )
        {
            tIndices( i, j ) = this->get_loc_entity_ind_from_entity_glb_id( aIds( i, j ), aEntityRank );
        }
    }

    return tIndices;
}

//------------------------------------------------------------------------------

moris::Cell< moris::mtk::Cell const * >
Enriched_Integration_Mesh::get_mtk_cells_loc_inds( Matrix< IndexMat > const &aCellIndices )
{
    moris::uint                             tNumCells = aCellIndices.numel();
    moris::Cell< moris::mtk::Cell const * > tCells( tNumCells );

    for ( moris::uint i = 0; i < tNumCells; i++ )
    {
        tCells( i ) = &this->get_mtk_cell( aCellIndices( i ) );
    }

    return tCells;
}
//------------------------------------------------------------------------------

moris::Cell< moris::mtk::Vertex const * >
Enriched_Integration_Mesh::get_mtk_vertices_loc_inds( Matrix< IndexMat > const &aVertexIndices )
{
    moris::uint                               tNumVerts = aVertexIndices.numel();
    moris::Cell< moris::mtk::Vertex const * > tVertices( tNumVerts );

    for ( moris::uint i = 0; i < tNumVerts; i++ )
    {
        tVertices( i ) = &this->get_mtk_vertex( (moris_index)aVertexIndices( i ) );
    }

    return tVertices;
}

//------------------------------------------------------------------------------

xtk::Cell_Cluster const &
Enriched_Integration_Mesh::get_xtk_cell_cluster( mtk::Cell const &aInterpCell ) const
{
    return get_cell_cluster( aInterpCell.get_index() );
}

//------------------------------------------------------------------------------

mtk::Cell_Cluster const &
Enriched_Integration_Mesh::get_cell_cluster( mtk::Cell const &aInterpCell ) const
{
    return get_cell_cluster( aInterpCell.get_index() );
}
//------------------------------------------------------------------------------

Cell_Cluster const &
Enriched_Integration_Mesh::get_cell_cluster( moris_index aInterpCellIndex ) const
{
    MORIS_ASSERT( aInterpCellIndex < (moris_index)mCellClusters.size(), "Interpolation Cell index out of bounds" );
    return *mCellClusters( aInterpCellIndex );
}
//------------------------------------------------------------------------------

moris::Cell< std::string >
Enriched_Integration_Mesh::get_block_set_names() const
{
    return mBlockSetNames;
}

//------------------------------------------------------------------------------

std::string
Enriched_Integration_Mesh::get_block_set_label( moris_index aBlockSetOrdinal ) const
{
    MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mSideSetLabels.size(), "Block set ordinal out of bounds" );
    return mBlockSetNames( aBlockSetOrdinal );
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::get_block_set_index( std::string aBlockSetLabel ) const
{
    auto tIter = mBlockSetLabelToOrd.find( aBlockSetLabel );

    MORIS_ERROR( tIter != mBlockSetLabelToOrd.end(), "block set set label not found" );

    return tIter->second;
}

//------------------------------------------------------------------------------

moris::Cell< mtk::Cluster const * >
Enriched_Integration_Mesh::get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const
{
    MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mBlockSetNames.size(), "Requested block set ordinal out of bounds." );

    moris::Cell< xtk::Cell_Cluster const * > const &tXTKClustersInSet = mPrimaryBlockSetClusters( aBlockSetOrdinal );

    moris::Cell< mtk::Cluster const * > tClusterInSet( tXTKClustersInSet.size() );

    for ( moris::uint i = 0; i < tXTKClustersInSet.size(); i++ )
    {
        tClusterInSet( i ) = tXTKClustersInSet( i );
    }

    return tClusterInSet;
}
//------------------------------------------------------------------------------

moris::Cell< xtk::Cell_Cluster const * > const &
Enriched_Integration_Mesh::get_xtk_cell_clusters_in_block_set( moris_index aBlockSetOrdinal ) const
{
    MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mBlockSetNames.size(), "Requested block set ordinal out of bounds." );

    return mPrimaryBlockSetClusters( aBlockSetOrdinal );
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::get_block_set_colors( moris_index aBlockSetOrdinal ) const
{
    MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mBlockSetColors.size(), "Block set ordinal out of bounds" );
    return mBlockSetColors( aBlockSetOrdinal );
}

//------------------------------------------------------------------------------

moris::Cell< mtk::Cluster const * >
Enriched_Integration_Mesh::get_side_set_cluster( moris_index aSideSetOrdinal ) const
{
    MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSets.size(), "Side set ordinal out of bounds" );

    moris::uint tNumSideClustersInSet = mSideSets( aSideSetOrdinal ).size();

    moris::Cell< mtk::Cluster const * > tSideClustersInSet( tNumSideClustersInSet );

    for ( moris::uint i = 0; i < tNumSideClustersInSet; i++ )
    {
        tSideClustersInSet( i ) = mSideSets( aSideSetOrdinal )( i ).get();
    }

    return tSideClustersInSet;
}

//------------------------------------------------------------------------------

Matrix< IndexMat >
Enriched_Integration_Mesh::get_side_set_colors( moris_index aSideSetOrdinal ) const
{
    MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSetColors.size(), "Side set ordinal out of bounds" );
    return mSideSetColors( aSideSetOrdinal );
}
//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_num_side_sets() const
{
    return mSideSets.size();
}
//------------------------------------------------------------------------------

std::string
Enriched_Integration_Mesh::get_side_set_label( moris_index aSideSetOrdinal ) const
{
    MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSetLabels.size(), "Side set ordinal out of bounds" );
    return mSideSetLabels( aSideSetOrdinal );
}

//------------------------------------------------------------------------------
moris_index
Enriched_Integration_Mesh::get_side_set_index( std::string aSideSetLabel ) const
{
    auto tIter = mSideSideSetLabelToOrd.find( aSideSetLabel );

    MORIS_ERROR( tIter != mSideSideSetLabelToOrd.end(), "side side set label not found" );

    return tIter->second;
}

//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_num_double_sided_sets() const
{
    return mDoubleSideSetLabels.size();
}

//------------------------------------------------------------------------------

std::string
Enriched_Integration_Mesh::get_double_sided_set_label( moris_index aSideSetOrdinal ) const
{
    MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(), "Double side set ordinal out of bounds" );
    return mDoubleSideSetLabels( aSideSetOrdinal );
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::get_double_sided_set_index( std::string aDoubleSideSetLabel ) const
{
    auto tIter = mDoubleSideSetLabelToOrd.find( aDoubleSideSetLabel );

    MORIS_ERROR( tIter != mDoubleSideSetLabelToOrd.end(), "double side set label not found" );

    return tIter->second;
}

//------------------------------------------------------------------------------

moris::Cell< mtk::Cluster const * >
Enriched_Integration_Mesh::get_double_side_set_cluster( moris_index aSideSetOrdinal ) const
{
    MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(), "Double side set ordinal out of bounds" );
    moris::Cell< mtk::Cluster const * > tDblSideClusters;
    tDblSideClusters.data() = std::vector< mtk::Cluster const * >( mDoubleSideSets( aSideSetOrdinal ).cbegin(),
        mDoubleSideSets( aSideSetOrdinal ).cend() );

    return tDblSideClusters;
}

moris::Cell< std::string >
Enriched_Integration_Mesh::get_field_names( enum moris::EntityRank aEntityRank )
{
    moris::Cell< std::string > tOutputFieldNames;

    moris_index tRankFieldIndex = this->get_entity_rank_field_index( aEntityRank );

    for ( auto const &iter : mFieldLabelToIndex( tRankFieldIndex ) )
    {
        tOutputFieldNames.push_back( iter.first );
    }

    return tOutputFieldNames;
}


moris::moris_index
Enriched_Integration_Mesh::create_field(
    std::string            aLabel,
    enum moris::EntityRank aEntityRank,
    moris::moris_index     aBulkPhaseIndex )
{
    MORIS_ASSERT( !field_exists( aLabel, aEntityRank ), "Field already created" );

    moris::moris_index tFieldIndex                                                 = mFields.size();
    mFieldLabelToIndex( this->get_entity_rank_field_index( aEntityRank ) )[aLabel] = tFieldIndex;
    mFields.push_back( Field( aLabel, aBulkPhaseIndex ) );

    return tFieldIndex;
}

//------------------------------------------------------------------------------
Matrix< IndexMat >
Enriched_Integration_Mesh::get_double_side_set_colors( moris_index aSideSetOrdinal ) const
{
    MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(), "Double side set ordinal out of bounds" );
    return mMasterDoubleSideSetColor( aSideSetOrdinal );
}
//------------------------------------------------------------------------------

uint
Enriched_Integration_Mesh::get_sidesets_num_faces( moris::Cell< moris_index > aSideSetIndex ) const
{
    moris::uint tNumSideSetFaces = 0;

    for ( moris_index Ik = 0; Ik < (moris_index)aSideSetIndex.size(); ++Ik )
    {
        MORIS_ASSERT( aSideSetIndex( Ik ) < (moris_index)mSideSets.size(), "Side set index out of bounds" );

        // add up the sideset number of faces
        tNumSideSetFaces = tNumSideSetFaces + mSideSets( aSideSetIndex( Ik ) ).size();
    }

    return tNumSideSetFaces;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print() const
{
    this->print_general();
    this->print_block_sets();
    this->print_side_sets();
    this->print_double_side_sets();
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print_general() const
{
    moris::uint tNumCells = this->get_num_owned_cells();

    moris::uint tNumGlobalCells = sum_all( tNumCells );
    if ( par_rank() == 0 )
    {
        std::cout << "Num Cells: " << std::setw( 8 ) << tNumGlobalCells << std::endl;
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print_cell_clusters( moris::uint aVerbosityLevel ) const
{
    std::cout << "\nCell Clusters:" << std::endl;
    for ( moris::uint i = 0; i < mCellClusters.size(); i++ )
    {
        xtk::Cell_Cluster *tCluster    = mCellClusters( i ).get();
        std::string        tTrivialStr = "f";
        if ( tCluster->is_trivial() )
        {
            tTrivialStr = "t";
        }

        Interpolation_Cell_Unzipped const *tBaseInterpCell = tCluster->get_xtk_interpolation_cell();

        std::cout << "    Cluster Index: " << std::setw( 9 ) << i
                  << " | Interp Cell Id: " << std::setw( 9 ) << tCluster->get_interpolation_cell_id()
                  << " | Base Interp Cell Id: " << std::setw( 9 ) << tBaseInterpCell->get_base_cell()->get_id()
                  << " | Trivial: " << tTrivialStr
                  << " | Num Primary: " << std::setw( 9 ) << tCluster->get_num_primary_cells()
                  << " | Num Void: " << tCluster->get_num_void_cells()
                  << std::endl;

        moris::print( tCluster->get_vertices_local_coordinates_wrt_interp_cell(), "Local Coords" );

        if ( aVerbosityLevel > 0 )
        {
            moris::Cell< moris::mtk::Cell const * > const &tPrimaryCells = tCluster->get_primary_cells_in_cluster();
            std::cout << "        Primary Integration Cell Ids: ";
            for ( moris::uint i = 0; i < tCluster->get_num_primary_cells(); i++ )
            {
                std::cout << std::setw( 9 ) << tPrimaryCells( i )->get_id();
            }
            std::cout << "\n"
                      << std::endl;
        }
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print_block_sets( moris::uint aVerbosityLevel ) const
{
    std::cout << "\nBlock Sets:" << std::endl;
    std::cout << "    Num Block Sets: " << this->get_num_blocks();

    for ( moris::uint iBS = 0; iBS < this->get_num_blocks(); iBS++ )
    {
        std::cout << "\n    Block Name: " << std::setw( 20 ) << mBlockSetNames( iBS ) << " | Block Set Ord: " << std::setw( 9 ) << iBS << " | Num Cell Clusters: " << std::setw( 9 ) << mPrimaryBlockSetClusters( iBS ).size() << " | Bulk Phase: " << std::setw( 9 ) << mBlockSetColors( iBS )( 0 );

        if ( aVerbosityLevel > 0 )
        {
            moris::Cell< xtk::Cell_Cluster const * > tClusters = this->mPrimaryBlockSetClusters( iBS );
            std::cout << "\n            Cluster in set\n";
            for ( moris::uint i = 0; i < tClusters.size(); i++ )
            {
                xtk::Cell_Cluster const *tCluster    = tClusters( i );
                std::string              tTrivialStr = "f";
                if ( tCluster->is_trivial() )
                {
                    tTrivialStr = "t";
                }

                Interpolation_Cell_Unzipped const *tBaseInterpCell = tCluster->get_xtk_interpolation_cell();

                std::cout << "            Cluster Index: " << std::setw( 9 ) << i
                          << " | Interp Cell Id: " << std::setw( 9 ) << tCluster->get_interpolation_cell_id()
                          << " | Base Interp Cell Id: " << std::setw( 9 ) << tBaseInterpCell->get_base_cell()->get_id()
                          << " | Trivial: " << tTrivialStr
                          << " | Num Primary: " << std::setw( 9 ) << tCluster->get_num_primary_cells()
                          << " | Num Void: " << tCluster->get_num_void_cells()
                          << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print_side_sets( moris::uint aVerbosityLevel ) const
{
    std::cout << "\nSide Sets:" << std::endl;
    std::cout << "    Num Side Sets: " << this->get_num_side_sets() << std::endl;

    for ( moris::uint iSS = 0; iSS < this->get_num_side_sets(); iSS++ )
    {
        std::cout << "    Side Set Name: " << std::setw( 20 ) << mSideSetLabels( iSS ) << " | Side Set Ord: " << std::setw( 9 ) << iSS << " | Num Cell Clusters: " << std::setw( 9 ) << this->mSideSets( iSS ).size() << " | Bulk Phase: " << std::setw( 9 ) << mSideSetColors( iSS )( 0 ) << std::endl;
        if ( aVerbosityLevel > 0 )
        {
            for ( moris::uint iCL = 0; iCL < this->mSideSets( iSS ).size(); iCL++ )
            {
                const std::shared_ptr< xtk::Side_Cluster > tCluster = this->mSideSets( iSS )( iCL );

                moris::Cell< moris::mtk::Cell const * > const &tPrimaryCells = tCluster->get_primary_cells_in_cluster();

                for ( moris::uint i = 0; i < tCluster->get_num_primary_cells(); i++ )
                {
                    std::cout << "Integration Cell Id = " << std::setw( 9 ) << tPrimaryCells( i )->get_id() << std::endl;
                    moris::print( tCluster->get_cell_local_coords_on_side_wrt_interp_cell( i ), "Local Coords on Side" );
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print_double_side_sets( moris::uint aVerbosityLevel ) const
{
    std::cout << "\nDouble Side Sets:" << std::endl;
    std::cout << "    Num Side Sets: " << this->get_num_double_side_set() << std::endl;

    for ( moris::uint iSS = 0; iSS < this->get_num_double_side_set(); iSS++ )
    {
        std::cout << "    Dbl Side Set Name: " << std::setw( 20 ) << mDoubleSideSetLabels( iSS ) << " | Dbl Side Set Ord: " << std::setw( 9 ) << iSS << " | Num Cell Clusters: " << std::setw( 9 ) << this->mDoubleSideSets( iSS ).size() << " | Master Bulk Phase: " << std::setw( 9 ) << mMasterDoubleSideSetColor( iSS )( 0 ) << " | Slave Bulk Phase: " << std::setw( 9 ) << mSlaveDoubleSideSetColor( iSS )( 0 );

        if ( aVerbosityLevel > 0 )
        {
            for ( moris::uint i = 0; i < mDoubleSideSets( iSS ).size(); i++ )
            {
                std::cout << "\n      Master Interpolation Cell: " << std::setw( 9 ) << mDoubleSideSets( iSS )( i )->get_interpolation_cell( mtk::Master_Slave::MASTER ).get_id();
                std::cout << " | Slave Interpolation Cell: " << std::setw( 9 ) << mDoubleSideSets( iSS )( i )->get_interpolation_cell( mtk::Master_Slave::SLAVE ).get_id();
            }
        }

        std::cout << std::endl;
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::print_double_side_clusters( moris::uint aVerbosityLevel ) const
{
    std::cout << "\nDouble Side Clusters:" << std::endl;
    std::cout << "    Num Double Side Clusters: " << mDoubleSideClusters.size() << std::endl;

    for ( moris::uint i = 0; i < mDoubleSideClusters.size(); i++ )
    {
        std::cout << mDoubleSideClusters( i ) << std::endl;
    }
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::create_side_set_from_dbl_side_set( moris_index const &aDblSideSetIndex,
    std::string const &                                                          aSideSetName,
    bool                                                                         aCollectSets )
{
    Cell< moris_index > tSideSetIndex = this->register_side_set_names( { aSideSetName } );

    moris::Cell< mtk::Double_Side_Cluster * > &tDblSideClusters = mDoubleSideSets( aDblSideSetIndex );

    moris::uint tCount = 0;
    for ( moris::uint i = 0; i < tDblSideClusters.size(); i++ )
    {
        // get the index
        moris_index tMasterIndex = mDoubleSideSetsMasterIndex( aDblSideSetIndex )( i );
        moris_index tSlaveIndex  = mDoubleSideSetsSlaveIndex( aDblSideSetIndex )( i );

        mSideSets( tSideSetIndex( 0 ) ).push_back( mDoubleSideSingleSideClusters( tMasterIndex ) );
        mSideSets( tSideSetIndex( 0 ) ).push_back( mDoubleSideSingleSideClusters( tSlaveIndex ) );
        tCount++;
    }

    this->commit_side_set( tSideSetIndex( 0 ) );

    this->set_side_set_colors( tSideSetIndex( 0 ), this->get_double_side_set_colors( aDblSideSetIndex ) );

    if ( aCollectSets )
    {
        this->setup_color_to_set();
        this->collect_all_sets( false );
    }

    return tSideSetIndex( 0 );
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::create_block_set_from_cells_of_side_set(
    moris_index const &      aSideSetIndex,
    std::string const &      aBlockSetName,
    enum CellTopology const &aCellTopo )
{
    moris::Cell< std::shared_ptr< xtk::Side_Cluster > > &tSideClusters = mSideSets( aSideSetIndex );

    Cell< moris_index > tBlockSetIndex = this->register_block_set_names_with_cell_topo( { aBlockSetName }, aCellTopo );

    std::unordered_map< moris_index, moris_index > tIpCellInSet;

    for ( moris::uint i = 0; i < tSideClusters.size(); i++ )
    {
        // cast to xtk side cluster
        xtk::Side_Cluster *tSideCluster = tSideClusters( i ).get();

        if ( tIpCellInSet.find( tSideCluster->mIntegrationCells( 0 )->get_id() ) == tIpCellInSet.end() )
        {
            // create a new cell cluster
            std::shared_ptr< xtk::Cell_Cluster > tCellCluster = std::make_shared< Cell_Cluster >();

            // get the ith enriched interpolation cell
            tCellCluster->mInterpolationCell = tSideCluster->mInterpolationCell;

            // mark as trivial
            tCellCluster->mTrivial = true;

            // add interp cell as integration cell
            tCellCluster->mPrimaryIntegrationCells.append( tSideCluster->mIntegrationCells );

            mCellClusters.push_back( tCellCluster );

            mPrimaryBlockSetClusters( tBlockSetIndex( 0 ) ).push_back( tCellCluster.get() );

            tIpCellInSet[tSideCluster->mIntegrationCells( 0 )->get_id()] = i;
        }
    }

    this->commit_block_set( tBlockSetIndex( 0 ) );
    this->set_block_set_colors( tBlockSetIndex( 0 ), this->get_side_set_colors( aSideSetIndex ) );
    this->setup_color_to_set();
    this->collect_all_sets( false );

    return tBlockSetIndex( 0 );
}

//------------------------------------------------------------------------------

std::string
Enriched_Integration_Mesh::get_interface_side_set_name(
    moris_index aGeomIndex,
    moris_index aBulkPhaseIndex0,
    moris_index aBulkPhaseIndex1 )
{
    MORIS_ASSERT( aGeomIndex < (moris_index)mModel->get_geom_engine()->get_num_geometries(), "Geometry index out of bounds" );
    MORIS_ASSERT( aBulkPhaseIndex0 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 0 out of bounds" );
    MORIS_ASSERT( aBulkPhaseIndex1 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 1 out of bounds" );

    return "iside_b0_" + std::to_string( aBulkPhaseIndex0 ) + "_b1_" + std::to_string( aBulkPhaseIndex1 );
}

//------------------------------------------------------------------------------

std::string
Enriched_Integration_Mesh::get_dbl_interface_side_set_name(
    moris_index aBulkPhaseIndex0,
    moris_index aBulkPhaseIndex1 )
{
    MORIS_ASSERT( aBulkPhaseIndex0 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 0 out of bounds" );
    MORIS_ASSERT( aBulkPhaseIndex1 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 1 out of bounds" );

    return "dbl_iside_p0_" + std::to_string( aBulkPhaseIndex0 ) + "_p1_" + std::to_string( aBulkPhaseIndex1 );
}

//------------------------------------------------------------------------------

moris::moris_index
Enriched_Integration_Mesh::get_field_index( std::string aLabel,
    enum moris::EntityRank                              aEntityRank )
{
    MORIS_ASSERT( field_exists( aLabel, aEntityRank ), "Field does not exist in mesh" );

    moris_index tIndex = get_entity_rank_field_index( aEntityRank );
    auto        tIter  = mFieldLabelToIndex( tIndex ).find( aLabel );
    return tIter->second;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::add_field_data(
    moris::moris_index      aFieldIndex,
    enum moris::EntityRank  aEntityRank,
    Matrix< DDRMat > const &aFieldData )
{
    mFields( aFieldIndex ).mFieldData = aFieldData.copy();
}

//------------------------------------------------------------------------------

Matrix< DDRMat > const &
Enriched_Integration_Mesh::get_field_data(
    moris::moris_index     aFieldIndex,
    enum moris::EntityRank aEntityRank ) const
{
    return mFields( aFieldIndex ).mFieldData;
}

//------------------------------------------------------------------------------
moris_id
Enriched_Integration_Mesh::allocate_entity_ids(
    moris::size_t   aNumReqs,
    enum EntityRank aEntityRank )
{
    MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only Elements or Nodes have ids" );

    moris_id tGlobalMax = this->get_max_entity_id( aEntityRank );

    int tProcRank = par_rank();
    int tProcSize = par_size();

    moris::Cell< moris::moris_id > aGatheredInfo;
    moris::Cell< moris::moris_id > tFirstId( 1 );
    moris::Cell< moris::moris_id > tNumIdsRequested( 1 );

    tNumIdsRequested( 0 ) = (moris::moris_id)aNumReqs;

    moris::gather( tNumIdsRequested, aGatheredInfo );

    moris::Cell< moris::moris_id > tProcFirstID( tProcSize );

    if ( tProcRank == 0 )
    {
        // Loop over entities print the number of entities requested by each processor
        for ( int iProc = 0; iProc < tProcSize; ++iProc )
        {
            // Give each processor their desired amount of IDs
            tProcFirstID( iProc ) = tGlobalMax;

            // Increment the first available node ID
            tGlobalMax = tGlobalMax + aGatheredInfo( iProc );
        }
    }

    moris::scatter( tProcFirstID, tFirstId );

    return tFirstId( 0 );
}


//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::commit_double_side_set( moris_index const &aDoubleSideSetIndex )
{

    MORIS_ASSERT( mListofDoubleSideSets.size() == (uint)aDoubleSideSetIndex,
        "Committing double side set failed. aDoubleSideSetIndex needs to be equivalent to the size of the list of double side sets" );

    mListofDoubleSideSets.resize( mListofDoubleSideSets.size() + 1, nullptr );

    mListofDoubleSideSets( aDoubleSideSetIndex ) = new moris::mtk::Double_Side_Set( mDoubleSideSetLabels( aDoubleSideSetIndex ),
        this->get_double_side_set_cluster( aDoubleSideSetIndex ),
        this->get_double_side_set_colors( aDoubleSideSetIndex ),
        this->get_spatial_dim() );
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::commit_side_set( moris_index const &aSideSetIndex )
{
    MORIS_ASSERT( mListofSideSets.size() == (uint)aSideSetIndex,
        "Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of single side sets" );

    mListofSideSets.resize( mListofSideSets.size() + 1, nullptr );

    mListofSideSets( aSideSetIndex ) = new moris::mtk::Side_Set(
        mSideSetLabels( aSideSetIndex ),
        this->get_side_set_cluster( aSideSetIndex ),
        this->get_side_set_colors( aSideSetIndex ),
        this->get_spatial_dim() );
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::commit_block_set( moris_index const &aBlockSetIndex )
{
    MORIS_ASSERT( mListofBlocks.size() == (uint)aBlockSetIndex,
        "Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of double side sets" );

    mListofBlocks.resize( mListofBlocks.size() + 1, nullptr );

    mListofBlocks( aBlockSetIndex ) = new moris::mtk::Block( mBlockSetNames( aBlockSetIndex ),
        this->get_cell_clusters_in_set( aBlockSetIndex ),
        this->get_block_set_colors( aBlockSetIndex ),
        this->get_spatial_dim() );
}
//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_cell_clusters()
{
    Tracer tTracer( "XTK", "Enriched Integration Mesh", "setup_cell_clusters" );

    Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( mMeshIndexInModel );

    // Number of interpolation cells
    moris::uint tNumInterpCells = tEnrInterpMesh->get_num_entities( EntityRank::ELEMENT );

    // allocate cell cluster member data
    mCellClusters.resize( tNumInterpCells, nullptr );

    // Allocate subphase index to cluster index
    mSubphaseIndexToClusterIndex.resize( 1, tNumInterpCells );

    // reference the enriched cells
    Cell< Interpolation_Cell_Unzipped * > const &tEnrichedInterpCells = tEnrInterpMesh->get_enriched_interpolation_cells();

    // iterate through interpolation cells to create cell clusters
    for ( moris::uint i = 0; i < tNumInterpCells; i++ )
    {
        // index
        moris_index tInterpCellIndex = tEnrichedInterpCells( i )->get_index();

        // create a new cell cluster
        mCellClusters( tInterpCellIndex ) = std::make_shared< Cell_Cluster >();

        // get the ith enriched interpolation cell
        mCellClusters( tInterpCellIndex )->mInterpolationCell = tEnrichedInterpCells( i );

        // ad subphase index to cluster index to subphase index data
        mSubphaseIndexToClusterIndex( i ) = mCellClusters( tInterpCellIndex )->get_xtk_interpolation_cell()->get_subphase_index();

        // base cell
        moris::mtk::Cell const *tBaseInterpCell = mCellClusters( tInterpCellIndex )->mInterpolationCell->get_base_cell();

        // ask background mesh if the base cell has children (the opposite answer to this question is the trivial flag)
        mCellClusters( tInterpCellIndex )->mTrivial = !mCutIgMesh->parent_cell_has_children( tBaseInterpCell->get_index() );

        // if it has children get a pointer to the child mesh
        if ( !mCellClusters( tInterpCellIndex )->mTrivial )
        {
            // subphase index
            moris_index tProcSubphaseIndex = mCellClusters( tInterpCellIndex )->mInterpolationCell->get_subphase_index();

            std::shared_ptr< IG_Cell_Group > tSubphaseCells = mCutIgMesh->get_subphase_ig_cells( tProcSubphaseIndex );

            // number of subphase associated with the same background cell
            moris::Cell< moris_index > const &tParentCellSubphases = mCutIgMesh->get_parent_cell_subphases( tBaseInterpCell->get_index() );

            // vertex group
            std::shared_ptr< IG_Vertex_Group > tVertexGroupForCluster = mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tBaseInterpCell->get_index() ) );


            moris::Cell< std::shared_ptr< IG_Cell_Group > > tVoidSubphases;
            tVoidSubphases.reserve( tParentCellSubphases.size() - 1 );

            // get the subphases in the void region
            for ( moris::uint iVoid = 0; iVoid < tParentCellSubphases.size(); iVoid++ )
            {
                if ( tParentCellSubphases( iVoid ) != tProcSubphaseIndex )
                {
                    tVoidSubphases.push_back( mCutIgMesh->get_subphase_ig_cells( tParentCellSubphases( iVoid ) ) );
                }
            }

            mCellClusters( tInterpCellIndex )->set_primary_integration_cell_group( tSubphaseCells );
            mCellClusters( tInterpCellIndex )->set_void_integration_cell_groups( tVoidSubphases );
            mCellClusters( tInterpCellIndex )->set_ig_vertex_group( tVertexGroupForCluster );
        }

        // trivial case, the base of the enriched interpolation cell becomes the primary cell
        else
        {
            mCellClusters( tInterpCellIndex )->mPrimaryIntegrationCells.push_back( mCellClusters( tInterpCellIndex )->mInterpolationCell->get_base_cell() );
            Matrix< IndexMat > tVertexIndices                     = mCellClusters( tInterpCellIndex )->mPrimaryIntegrationCells( 0 )->get_vertex_inds();
            mCellClusters( tInterpCellIndex )->mVerticesInCluster = this->get_mtk_vertices_loc_inds( tVertexIndices );
            tBaseInterpCell->get_cell_info()->get_loc_coords_of_cell( mCellClusters( tInterpCellIndex )->mLocalCoords );
        }
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_blockset_with_cell_clusters()
{
    // get background mesh
    Background_Mesh &tBackgroundMesh = mModel->get_background_mesh();

    // enriched interpolation mesh
    Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( mMeshIndexInModel );

    // my proc rank
    moris_index tProcRank = par_rank();

    // get block sets (in background mesh data)
    Cell< std::string > tBlockSetsNames = tBackgroundMesh.get_mesh_data().get_set_names( EntityRank::ELEMENT );

    // for each block set construct
    for ( moris::uint iBS = 0; iBS < tBlockSetsNames.size(); iBS++ )
    {
        // split set into child and no child as we need to have the same type of integration cell in each set
        moris::Cell< std::string > tChildNoChildSetNames = this->split_set_name_by_child_no_child( tBlockSetsNames( iBS ) );

        // split child and no child sets by phases
        moris::Cell< std::string > tPhaseChildBlockSetNames   = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 0 ) );
        moris::Cell< std::string > tPhaseNoChildBlockSetNames = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 1 ) );

        // topology enums
        enum CellTopology tChildTopo  = mModel->get_cut_mesh().get_child_element_topology();
        enum CellTopology tParentTopo = mModel->get_background_mesh().get_parent_cell_topology();

        // add block set names to member data
        Cell< moris_index > tChildBlockSetOrds   = this->register_block_set_names_with_cell_topo( tPhaseChildBlockSetNames, tChildTopo );
        Cell< moris_index > tNoChildBlockSetOrds = this->register_block_set_names_with_cell_topo( tPhaseNoChildBlockSetNames, tParentTopo );

        // set block set colors
        for ( moris_index i = 0; i < (moris_index)tChildBlockSetOrds.size(); i++ )
        {
            this->set_block_set_colors( tChildBlockSetOrds( i ), { { i } } );
            this->set_block_set_colors( tNoChildBlockSetOrds( i ), { { i } } );
        }

        // get the cells in this block
        moris::Cell< moris::mtk::Cell const * > tCellsInBlock = tBackgroundMesh.get_mesh_data().get_block_set_cells( tBlockSetsNames( iBS ) );

        // get the enriched interpolation cells in this block
        moris::Cell< xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsInBlock = tEnrInterpMesh->get_enriched_cells_from_base_cells( tCellsInBlock );

        // iterate through and add cluster associated with enriched cell to block set
        for ( moris::uint iC = 0; iC < tEnrichedCellsInBlock.size(); iC++ )
        {
            // get the bulk phase
            moris_index tBulkPhaseIndex = tEnrichedCellsInBlock( iC )->get_bulkphase_index();

            // get cluster associated with enriched cell
            xtk::Cell_Cluster const &tCluster = this->get_cell_cluster( tEnrichedCellsInBlock( iC )->get_index() );

            if ( tEnrichedCellsInBlock( iC )->get_owner() == tProcRank )
            {
                // set ord
                moris_index tSetOrd = MORIS_INDEX_MAX;

                if ( tCluster.is_trivial() )
                {
                    tSetOrd = tNoChildBlockSetOrds( tBulkPhaseIndex );
                }

                else
                {
                    tSetOrd = tChildBlockSetOrds( tBulkPhaseIndex );
                }

                // add to member data
                mPrimaryBlockSetClusters( tSetOrd ).push_back( &tCluster );
            }
        }
    }

    for ( moris::uint Ik = mListofBlocks.size(); Ik < mPrimaryBlockSetClusters.size(); Ik++ )
    {
        this->commit_block_set( Ik );
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_side_set_clusters()
{
    // get data for easy access
    Enriched_Interpolation_Mesh *tEnrInterpMesh  = mModel->mEnrichedInterpMesh( mMeshIndexInModel );
    Background_Mesh &            tBackgroundMesh = mModel->mBackgroundMesh;
    Integration_Mesh_Generator   tIGMeshGen;

    // rank enum for facets
    enum EntityRank tFacetRank = mModel->mBackgroundMesh.get_mesh_data().get_facet_rank();

    // get side sets (in background mesh data)
    Cell< std::string > tSideSetNames = tBackgroundMesh.get_mesh_data().get_set_names( tFacetRank );

    tSideSetNames = mModel->check_for_and_remove_internal_seacas_side_sets( tSideSetNames );

    // my proc rank
    moris_index tParRank = par_rank();

    // access background facet to child facet connectivity
    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > const &tBGFacetToChildFacet = mCutIgMesh->get_background_facet_to_child_facet_connectivity();

    // for each side set construct
    for ( moris::uint iSS = 0; iSS < tSideSetNames.size(); iSS++ )
    {
        // split set into child and no child as we need to have the same type of integration cell in each set
        moris::Cell< std::string > tChildNoChildSetNames = this->split_set_name_by_child_no_child( tSideSetNames( iSS ) );

        // split child and no child sets by phases
        moris::Cell< std::string > tPhaseChildSideSetNames   = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 0 ) );
        moris::Cell< std::string > tPhaseNoChildSideSetNames = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 1 ) );

        // add side set names to member data
        Cell< moris_index > tChildSideSetOrds   = this->register_side_set_names( tPhaseChildSideSetNames );
        Cell< moris_index > tNoChildSideSetOrds = this->register_side_set_names( tPhaseNoChildSideSetNames );

        // set side set colors
        for ( moris_index i = 0; i < (moris_index)tChildSideSetOrds.size(); i++ )
        {
            this->set_side_set_colors( tChildSideSetOrds( i ), { { i } } );
            this->set_side_set_colors( tNoChildSideSetOrds( i ), { { i } } );
        }

        // get the cells in this side set and their side ordinals
        moris::Cell< mtk::Cell const * > tCellsInSideSet;
        Matrix< IndexMat >               tCellOrdsInSideSet;

        tBackgroundMesh.get_mesh_data().get_sideset_cells_and_ords(
            tSideSetNames( iSS ),
            tCellsInSideSet,
            tCellOrdsInSideSet );

        // iterate through cells in side set
        for ( moris::uint iC = 0; iC < tCellsInSideSet.size(); iC++ )
        {
            mtk::Cell const *tBaseCell  = tCellsInSideSet( iC );
            moris_index      tSideOrd   = tCellOrdsInSideSet( iC );
            moris_index      tSideIndex = tBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds( tBaseCell->get_index(), EntityRank::ELEMENT, tBackgroundMesh.get_mesh_data().get_facet_rank() )( tSideOrd );


            // only place cluster's related to the background cells owned by current proc in sets
            if ( tBaseCell->get_owner() == tParRank )
            {
                // get the enriched interpolation cells associated with base cell
                moris::Cell< xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsOfBaseCell = tEnrInterpMesh->get_enriched_cells_from_base_cell( tBaseCell );

                // get the subphase indices associated with the current parent cell
                moris::Cell< moris_index > const &tSubphasesWrtParentCell = mCutIgMesh->get_parent_cell_subphases( tBaseCell->get_index() );

                MORIS_ASSERT( tSubphasesWrtParentCell.size() == tEnrichedCellsOfBaseCell.size(), "Discrepency in number of interpolation cells and subphases related to base cell." );

                // check some things out
                for ( moris::uint iSP = 0; iSP < tSubphasesWrtParentCell.size(); iSP++ )
                {
                    MORIS_ASSERT( tSubphasesWrtParentCell( iSP ) == tEnrichedCellsOfBaseCell( iSP )->get_subphase_index(), "Discrepency in subphase index." );
                }

                // if this is a null ptr, the bg facet does not have any child facets associated with it, this indicates that the side set is on a facet not intersected by the geometry
                if ( tBGFacetToChildFacet( tSideIndex ) != nullptr )
                {
                    // for the base cell, get the associated vertices
                    std::shared_ptr< IG_Vertex_Group > tVertexGroupForCluster = mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tBaseCell->get_index() ) );

                    // collect all the integration cells on the current facet and the side ordinal
                    moris::Cell< moris::mtk::Cell * > tIgCellsOnBgFacet;
                    moris::Cell< moris_index >        tIgCellsSideOrdsOnBgFacet;
                    tIGMeshGen.collect_ig_cells_and_side_ords_on_bg_facet( mCutIgMesh, tSideIndex, tIgCellsOnBgFacet, tIgCellsSideOrdsOnBgFacet );

                    // iterate through the subphases
                    for ( moris::uint iSP = 0; iSP < tSubphasesWrtParentCell.size(); iSP++ )
                    {

                        moris::Cell< moris::mtk::Cell const * > tIgCellsInSubphase;
                        moris::Cell< moris_index >              tIgCellsSideOrdsInSubphase;
                        tIgCellsInSubphase.reserve( tIgCellsOnBgFacet.size() );
                        tIgCellsInSubphase.reserve( tIgCellsSideOrdsInSubphase.size() );

                        // bulk phase of this subphase
                        moris_index tBulkPhase = mCutIgMesh->get_subphase_bulk_phase( tSubphasesWrtParentCell( iSP ) );

                        // iterate through ig cells on bg facet
                        for ( moris::uint iIGCell = 0; iIGCell < tIgCellsOnBgFacet.size(); iIGCell++ )
                        {
                            if ( mCutIgMesh->get_ig_cell_subphase_index( tIgCellsOnBgFacet( iIGCell )->get_index() ) == tSubphasesWrtParentCell( iSP ) )
                            {
                                tIgCellsInSubphase.push_back( tIgCellsOnBgFacet( iIGCell ) );
                                tIgCellsSideOrdsInSubphase.push_back( tIgCellsSideOrdsOnBgFacet( iIGCell ) );
                            }
                        }

                        // convert the side ordinal to a matrix
                        Matrix< IndexMat > tSideOrdsMat( 1, tIgCellsSideOrdsInSubphase.size() );
                        for ( moris::uint iCopy = 0; iCopy < tIgCellsSideOrdsInSubphase.size(); iCopy++ )
                        {
                            tSideOrdsMat( iCopy ) = tIgCellsSideOrdsInSubphase( iCopy );
                        }

                        if ( tIgCellsInSubphase.size() > 0 )
                        {
                            // make a side cluster
                            std::shared_ptr< Side_Cluster > tSideCluster = std::make_shared< Side_Cluster >();
                            tSideCluster->mInterpolationCell             = tEnrichedCellsOfBaseCell( iSP );
                            tSideCluster->mTrivial                       = false;
                            tSideCluster->mIntegrationCellSideOrdinals   = tSideOrdsMat;
                            tSideCluster->mAssociatedCellCluster         = &this->get_xtk_cell_cluster( *tEnrichedCellsOfBaseCell( iSP ) );
                            tSideCluster->mIntegrationCells              = tIgCellsInSubphase;

                            tSideCluster->set_ig_vertex_group( tVertexGroupForCluster );

                            // add to the integration mesh
                            mSideSets( tChildSideSetOrds( tBulkPhase ) ).push_back( tSideCluster );
                        }
                    }
                }
                else
                {
                    MORIS_ASSERT( tEnrichedCellsOfBaseCell.size() == 1,
                        "For the trivial case, a base cell should have 1 enriched interpolation cell associated with it" );

                    // phase of cell
                    moris_index tBulkPhase = tEnrichedCellsOfBaseCell( 0 )->get_bulkphase_index();

                    // side set ordinal in mesh
                    moris_index tSideSetOrd = tNoChildSideSetOrds( tBulkPhase );

                    // create a new side cluster in the side set assoicate with this bulk phase
                    moris_index tIndex = mSideSets( tSideSetOrd ).size();
                    mSideSets( tSideSetOrd ).push_back( std::make_shared< Side_Cluster >() );
                    std::shared_ptr< Side_Cluster > tSideCluster = mSideSets( tSideSetOrd )( tIndex );

                    // set trivial flag
                    tSideCluster->mTrivial = true;

                    // get the set enriched interpolation cell
                    tSideCluster->mInterpolationCell = tEnrichedCellsOfBaseCell( 0 );

                    // mark child mesh as nullptr
                    tSideCluster->mChildMesh = nullptr;

                    // integration cell is the same as the interpolation cell in this case
                    tSideCluster->mIntegrationCells = { tSideCluster->mInterpolationCell->get_base_cell() };

                    // side ordinal
                    tSideCluster->mIntegrationCellSideOrdinals = Matrix< IndexMat >( { { tSideOrd } } );

                    // add vertices
                    tSideCluster->mVerticesInCluster.append( tSideCluster->mIntegrationCells( 0 )->get_vertices_on_side_ordinal( tSideOrd ) );

                    tSideCluster->mAssociatedCellCluster = &this->get_xtk_cell_cluster( *tEnrichedCellsOfBaseCell( 0 ) );
                }
            }
        }
    }

    for ( moris::uint Ik = mListofSideSets.size(); Ik < mSideSets.size(); Ik++ )
    {
        this->commit_side_set( Ik );
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_double_side_set_clusters()
{
    this->setup_double_sided_interface_sides();
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_color_to_set()
{
    this->construct_color_to_set_relationship( mBlockSetColors, mColorsBlockSets );
    this->construct_color_to_set_relationship( mSideSetColors, mColorsSideSets );
    this->construct_color_to_set_relationship( mMasterDoubleSideSetColor, mColorMasterDoubleSideSet );
    this->construct_color_to_set_relationship( mSlaveDoubleSideSetColor, mColorSlaveDoubleSideSet );
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_double_sided_interface_sides()
{
    this->declare_interface_double_side_sets();

    this->create_interface_double_side_sets_and_clusters();
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::declare_interface_double_side_sets()
{
    uint tNumBulkPhases = mModel->get_geom_engine()->get_num_bulk_phase();

    Cell< std::string > tDoubleInterfaceSideNames;

    mBulkPhaseToDblSideIndex.resize( tNumBulkPhases, tNumBulkPhases );
    mBulkPhaseToDblSideIndex.fill( MORIS_INDEX_MAX );

    moris_index tCount = 0;

    Cell< Matrix< IndexMat > > tInterfaceMasterSideColors;
    Cell< Matrix< IndexMat > > tInterfaceSlaveSideColors;

    for ( moris::moris_index iP0 = 0; iP0 < (moris_index)tNumBulkPhases; iP0++ )
    {
        for ( moris::moris_index iP1 = iP0 + 1; iP1 < (moris_index)tNumBulkPhases; iP1++ )
        {

            std::string tInterfaceSideSetName = this->get_dbl_interface_side_set_name( iP0, iP1 );

            tDoubleInterfaceSideNames.push_back( tInterfaceSideSetName );
            tInterfaceMasterSideColors.push_back( { { iP0 } } );
            tInterfaceSlaveSideColors.push_back( { { iP1 } } );

            mBulkPhaseToDblSideIndex( iP0, iP1 ) = tCount;
            mBulkPhaseToDblSideIndex( iP1, iP0 ) = tCount;
            tCount++;
        }
    }

    Cell< moris_index > tDblSideSetOrds = this->register_double_side_set_names( tDoubleInterfaceSideNames );

    // set interface side set colors
    for ( moris_index iSS = 0; iSS < (moris_index)tDblSideSetOrds.size(); iSS++ )
    {
        this->set_double_side_set_colors( tDblSideSetOrds( iSS ), tInterfaceMasterSideColors( iSS ), tInterfaceSlaveSideColors( iSS ) );
    }
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::get_dbl_side_set_index(
    moris_index aPhase0,
    moris_index aPhase1 )
{
    MORIS_ASSERT( aPhase0 < aPhase1, "Double side sets are defined from low phase index to high" );

    MORIS_ASSERT( mDoubleSideSetLabels( mBulkPhaseToDblSideIndex( aPhase0, aPhase1 ) ) == this->get_dbl_interface_side_set_name( aPhase0, aPhase1 ),
        "Interface double side set not showing up in correct index" );

    return mBulkPhaseToDblSideIndex( aPhase0, aPhase1 );
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::create_interface_double_side_sets_and_clusters()
{
    Tracer tTracer( "XTK", "Enriched Integration Mesh", "create_interface_double_side_sets_and_clusters" );

    // tool for generating double sided interface
    Integration_Mesh_Generator tIGMeshGen;

    // access interpolation mesh
    Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( mMeshIndexInModel );

    // get the enriched interpolation cell
    Cell< Interpolation_Cell_Unzipped * > &tEnrIpCells = tEnrInterpMesh->get_enriched_interpolation_cells();

    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > const &tDoubleSidedInterface = mCutIgMesh->get_bulk_phase_to_bulk_phase_dbl_side_interface();

    // for a subphase to subphase side cluster - value in map is the location in tSideClusters
    moris::Cell< std::unordered_map< moris_index, moris_index > > tSubphaseToSubphaseSideClusterIndex( mCutIgMesh->get_num_subphases() );

    moris::Cell< std::shared_ptr< xtk::Side_Cluster > > tSideClusters;
    moris::Cell< moris::Cell< moris_index > >           tSideClusterSideOrdinals;

    for ( moris::uint iBP0 = 0; iBP0 < tDoubleSidedInterface.size(); iBP0++ )
    {
        for ( moris::uint iBP1 = 0; iBP1 < tDoubleSidedInterface.size(); iBP1++ )
        {
            if ( tDoubleSidedInterface( iBP0 )( iBP1 ) != nullptr )
            {
                // tDoubleSidedInterface(iBP0)(iBP1)->print();

                std::unordered_map< moris_index, moris_index > tSubphaseIndexToClusterInterfaceOrd;

                // acess pointer
                std::shared_ptr< IG_Cell_Double_Side_Group > tDblSideGroup = tDoubleSidedInterface( iBP0 )( iBP1 );

                // iterate through the facet pairs
                for ( moris::uint iDblFacet = 0; iDblFacet < tDoubleSidedInterface( iBP0 )( iBP1 )->mLeaderIgCells.size(); iDblFacet++ )
                {
                    moris_index tLeaderSubphaseIndex   = mCutIgMesh->get_ig_cell_subphase_index( tDblSideGroup->mLeaderIgCells( iDblFacet )->get_index() );
                    moris_index tFollowerSubphaseIndex = mCutIgMesh->get_ig_cell_subphase_index( tDblSideGroup->mFollowerIgCells( iDblFacet )->get_index() );

                    moris_index tLeaderCellSideOrdinal   = tDblSideGroup->mLeaderIgCellSideOrdinals( iDblFacet );
                    moris_index tFollowerCellSideOrdinal = tDblSideGroup->mFollowerIgCellSideOrdinals( iDblFacet );

                    moris_index tLeaderBulkPhase   = mCutIgMesh->get_subphase_bulk_phase( tLeaderSubphaseIndex );
                    moris_index tFollowerBulkPhase = mCutIgMesh->get_subphase_bulk_phase( tFollowerSubphaseIndex );

                    moris_index tLeaderOwner = tEnrIpCells( tLeaderSubphaseIndex )->get_owner();

                    if ( tLeaderBulkPhase < tFollowerBulkPhase && tLeaderOwner == moris::par_rank() )
                    {

                        auto tLeaderClusterIter = tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).find( tFollowerSubphaseIndex );

                        // do we have a double side cluster between these subphases yet?

                        //if not we construct the two side clusters
                        if ( tLeaderClusterIter == tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).end() )
                        {
                            // index of double side set
                            moris_index tDoubleSideSetIndex = this->get_dbl_side_set_index( tLeaderBulkPhase, tFollowerBulkPhase );

                            // if not construct one from leader to follower
                            moris_index tNewClusterIndex = tSideClusters.size();
                            tSideClusters.push_back( std::make_shared< xtk::Side_Cluster >() );
                            tSideClusterSideOrdinals.push_back( moris::Cell< moris_index >() );

                            // add master side cluster
                            mDoubleSideSetsMasterIndex( tDoubleSideSetIndex ).push_back( mDoubleSideSingleSideClusters.size() );
                            mDoubleSideSingleSideClusters.push_back( tSideClusters( tNewClusterIndex ) );

                            tSideClusters( tNewClusterIndex )->mInterpolationCell     = tEnrIpCells( tLeaderSubphaseIndex );
                            tSideClusters( tNewClusterIndex )->mTrivial               = false;
                            tSideClusters( tNewClusterIndex )->mAssociatedCellCluster = &this->get_xtk_cell_cluster( *tSideClusters( tNewClusterIndex )->mInterpolationCell );

                            // leader vertex group
                            std::shared_ptr< IG_Vertex_Group > tLeaderVertexGroup =
                                mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tEnrIpCells( tLeaderSubphaseIndex )->get_base_cell()->get_index() ) );

                            tSideClusters( tNewClusterIndex )->set_ig_vertex_group( tLeaderVertexGroup );

                            // add to the map
                            tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex )[tFollowerSubphaseIndex] = tNewClusterIndex;

                            // construct one from follower to leader
                            tNewClusterIndex = tSideClusters.size();
                            tSideClusters.push_back( std::make_shared< xtk::Side_Cluster >() );
                            tSideClusterSideOrdinals.push_back( moris::Cell< moris_index >() );

                            // add master side cluster to double side set
                            mDoubleSideSetsSlaveIndex( tDoubleSideSetIndex ).push_back( mDoubleSideSingleSideClusters.size() );
                            mDoubleSideSingleSideClusters.push_back( tSideClusters( tNewClusterIndex ) );


                            tSideClusters( tNewClusterIndex )->mInterpolationCell     = tEnrIpCells( tFollowerSubphaseIndex );
                            tSideClusters( tNewClusterIndex )->mTrivial               = false;
                            tSideClusters( tNewClusterIndex )->mAssociatedCellCluster = &this->get_xtk_cell_cluster( *tSideClusters( tNewClusterIndex )->mInterpolationCell );

                            // follower vertex group
                            std::shared_ptr< IG_Vertex_Group > tFollowerVertexGroup =
                                mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tEnrIpCells( tFollowerSubphaseIndex )->get_base_cell()->get_index() ) );

                            tSideClusters( tNewClusterIndex )->set_ig_vertex_group( tFollowerVertexGroup );

                            MORIS_ASSERT( tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).find( tLeaderSubphaseIndex ) == tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).end(),
                                "Nonconcurrent leader follower interface construction occured." );

                            // add to the map
                            tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex )[tLeaderSubphaseIndex] = tNewClusterIndex;

                            // create double side set
                            mtk::Double_Side_Cluster *tDblSideCluster = new mtk::Double_Side_Cluster(
                                tSideClusters( tNewClusterIndex - 1 ).get(),
                                tSideClusters( tNewClusterIndex ).get(),
                                tSideClusters( tNewClusterIndex - 1 )->get_vertices_in_cluster() );

                            mDoubleSideClusters.push_back( tDblSideCluster );
                            mDoubleSideSets( tDoubleSideSetIndex ).push_back( tDblSideCluster );
                        }

                        // get the relevant side cluster indices
                        moris_index tLeaderToFollowerSideClusterIndex = tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).find( tFollowerSubphaseIndex )->second;
                        moris_index tFollowerToLeaderSideClusterIndex = tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).find( tLeaderSubphaseIndex )->second;

                        MORIS_ASSERT( tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).find( tFollowerSubphaseIndex ) != tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).end(), "Side cluster not constructed" );
                        MORIS_ASSERT( tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).find( tLeaderSubphaseIndex ) != tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).end(), "Side cluster not constructed" );

                        // place these in a coincident struct which we will convert to matrix and store in the cluster later
                        tSideClusterSideOrdinals( tLeaderToFollowerSideClusterIndex ).push_back( tLeaderCellSideOrdinal );
                        tSideClusterSideOrdinals( tFollowerToLeaderSideClusterIndex ).push_back( tFollowerCellSideOrdinal );

                        // add integration cell pointers
                        tSideClusters( tLeaderToFollowerSideClusterIndex )->mIntegrationCells.push_back( tDblSideGroup->mLeaderIgCells( iDblFacet ) );
                        tSideClusters( tFollowerToLeaderSideClusterIndex )->mIntegrationCells.push_back( tDblSideGroup->mFollowerIgCells( iDblFacet ) );
                    }
                }
            }
        }
    }

    // convert the cells of side ords to matrix and add to clusters
    for ( moris::uint iSC = 0; iSC < tSideClusterSideOrdinals.size(); iSC++ )
    {
        Matrix< IndexMat > tSideOrds( 1, tSideClusterSideOrdinals( iSC ).size() );

        for ( moris::uint iSide = 0; iSide < tSideClusterSideOrdinals( iSC ).size(); iSide++ )
        {
            tSideOrds( iSide ) = tSideClusterSideOrdinals( iSC )( iSide );
        }

        tSideClusters( iSC )->mIntegrationCellSideOrdinals = tSideOrds;
    }

    // construct double side set interfaces
    for ( moris::uint Ik = mListofDoubleSideSets.size(); Ik < mDoubleSideSets.size(); Ik++ )
    {
        this->commit_double_side_set( Ik );
    }
}
//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::add_side_to_cluster(
    std::shared_ptr< xtk::Side_Cluster > aSideCluster,
    moris_index                          aCellIndex,
    moris_index                          aSideOrdinal )
{
    // number of current sides in cluster
    uint tNumCurrentSides = aSideCluster->get_num_sides_in_cluster();

    aSideCluster->mIntegrationCells.push_back( &this->get_mtk_cell( aCellIndex ) );

    // add sides
    aSideCluster->mIntegrationCellSideOrdinals.resize( 1, tNumCurrentSides + 1 );
    aSideCluster->mIntegrationCellSideOrdinals( tNumCurrentSides ) = aSideOrdinal;
}

//------------------------------------------------------------------------------
void
Enriched_Integration_Mesh::setup_side_cluster_vertices( std::shared_ptr< xtk::Side_Cluster > aMasterSideCluster,
    std::shared_ptr< xtk::Side_Cluster >                                                     aSlaveSideCluster )
{
    moris::Cell< mtk::Cell const * > const &tMasterCellsInCluster = aMasterSideCluster->get_cells_in_side_cluster();
    moris::Cell< mtk::Cell const * > const &tSlaveCellsInCluster  = aSlaveSideCluster->get_cells_in_side_cluster();

    moris::Matrix< moris::IndexMat > tMasterSideOrds = aMasterSideCluster->get_cell_side_ordinals();
    moris::Matrix< moris::IndexMat > tSlaveSideOrds  = aSlaveSideCluster->get_cell_side_ordinals();
    MORIS_ASSERT( aMasterSideCluster->get_num_sides_in_cluster() == aSlaveSideCluster->get_num_sides_in_cluster(), "Number of sides in side cluster mismatch" );

    // vector of vertices on side ordinals
    moris::Cell< moris::mtk::Vertex const * > tMasterVerticesInCluster;
    moris::Cell< moris::mtk::Vertex const * > tSlaveVerticesInCluster;

    // map to ensure vertices are added only one time
    std::unordered_map< moris_index, moris_index > tMasterUniqueVertexMap;
    std::unordered_map< moris_index, moris_index > tSlaveUniqueVertexMap;


    // add integration cells to cluster
    for ( moris::uint iF = 0; iF < aMasterSideCluster->get_num_sides_in_cluster(); iF++ )
    {
        // iterate through vertices and keep track of the unique ones
        moris::Cell< moris::mtk::Vertex const * > tMasterVerticesOnSide = tMasterCellsInCluster( iF )->get_vertices_on_side_ordinal( tMasterSideOrds( iF ) );
        moris::Cell< moris::mtk::Vertex const * > tSlaveVerticesOnSide  = tSlaveCellsInCluster( iF )->get_vertices_on_side_ordinal( tSlaveSideOrds( iF ) );

        for ( moris::uint iVoS = 0; iVoS < tMasterVerticesOnSide.size(); iVoS++ )
        {
            if ( tMasterUniqueVertexMap.find( tMasterVerticesOnSide( iVoS )->get_id() ) == tMasterUniqueVertexMap.end() )
            {
                tMasterUniqueVertexMap[tMasterVerticesOnSide( iVoS )->get_id()] = 1;

                tMasterVerticesInCluster.push_back( tMasterVerticesOnSide( iVoS ) );
            }

            if ( tSlaveUniqueVertexMap.find( tSlaveVerticesOnSide( iVoS )->get_id() ) == tSlaveUniqueVertexMap.end() )
            {
                tSlaveUniqueVertexMap[tSlaveVerticesOnSide( iVoS )->get_id()] = 1;

                tSlaveVerticesInCluster.push_back( tSlaveVerticesOnSide( iVoS ) );
            }
        }
    }

    // add the vertices to the cluster
    aMasterSideCluster->mVerticesInCluster = tMasterVerticesInCluster;
    aSlaveSideCluster->mVerticesInCluster  = tMasterVerticesInCluster;// intentionally using left here so ordering is consistent
}
//------------------------------------------------------------------------------

moris::Cell< std::string >
Enriched_Integration_Mesh::split_set_name_by_bulk_phase( std::string aBaseName )
{
    moris::uint                tNumPhases = mModel->mGeometryEngine->get_num_bulk_phase();
    moris::Cell< std::string > tSetNames( tNumPhases );
    for ( moris::uint i = 0; i < tNumPhases; i++ )
    {
        tSetNames( i ) = aBaseName + "_p" + std::to_string( i );
    }

    return tSetNames;
}

//------------------------------------------------------------------------------

moris::Cell< std::string >
Enriched_Integration_Mesh::split_set_name_by_child_no_child( std::string aBaseName )
{
    moris::Cell< std::string > tSetNames( 2 );
    tSetNames( 0 ) = aBaseName + "_c";
    tSetNames( 1 ) = aBaseName + "_n";
    return tSetNames;
}

//------------------------------------------------------------------------------
Cell< moris_index >
Enriched_Integration_Mesh::register_vertex_set_names( moris::Cell< std::string > const &aVertexSetNames )
{
    uint tNumSetsToRegister = aVertexSetNames.size();


    // block set ords
    Cell< moris_index > tVertexSetOrds( tNumSetsToRegister );

    // iterate and add sets
    for ( moris::uint i = 0; i < tNumSetsToRegister; i++ )
    {
        tVertexSetOrds( i ) = mVertexSetNames.size();

        mVertexSetNames.push_back( aVertexSetNames( i ) );
        MORIS_ASSERT( mVertexSetLabelToOrd.find( aVertexSetNames( i ) ) == mVertexSetLabelToOrd.end(),
            "Duplicate vertex set in mesh" );

        mVertexSetLabelToOrd[aVertexSetNames( i )] = tVertexSetOrds( i );
    }

    mVerticesInVertexSet.resize( mVerticesInVertexSet.size() + tNumSetsToRegister );
    mVertexSetColors.resize( mVertexSetColors.size() + tNumSetsToRegister );

    return tVertexSetOrds;
}

//------------------------------------------------------------------------------

Cell< moris_index >
Enriched_Integration_Mesh::register_block_set_names_with_cell_topo(
    moris::Cell< std::string > const &aBlockSetNames,
    enum CellTopology                 aBlockTopology )
{
    uint tNumSetsToRegister = aBlockSetNames.size();

    // block set ords
    Cell< moris_index > tBlockSetOrds( tNumSetsToRegister );

    // iterate and add sets
    for ( moris::uint i = 0; i < tNumSetsToRegister; i++ )
    {
        tBlockSetOrds( i ) = mBlockSetNames.size();

        mBlockSetNames.push_back( aBlockSetNames( i ) );
        mBlockSetTopology.push_back( aBlockTopology );
        MORIS_ASSERT( mBlockSetLabelToOrd.find( aBlockSetNames( i ) ) == mBlockSetLabelToOrd.end(), "Duplicate block set in mesh" );
        mBlockSetLabelToOrd[aBlockSetNames( i )] = tBlockSetOrds( i );
    }

    mPrimaryBlockSetClusters.resize( mPrimaryBlockSetClusters.size() + tNumSetsToRegister );
    mBlockSetColors.resize( mBlockSetColors.size() + tNumSetsToRegister );
    return tBlockSetOrds;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::set_block_set_colors(
    moris_index const &       aBlockSetIndex,
    Matrix< IndexMat > const &aBlockSetColors )
{
    MORIS_ASSERT( moris::isempty( mBlockSetColors( aBlockSetIndex ) ), "Attempting to overwrite colors of a block set" );

    mBlockSetColors( aBlockSetIndex ) = aBlockSetColors;
}

//------------------------------------------------------------------------------

Cell< moris_index >
Enriched_Integration_Mesh::register_side_set_names( moris::Cell< std::string > const &aSideSetNames )
{
    uint tNumSetsToRegister = aSideSetNames.size();

    // block set ords
    Cell< moris_index > tSideSetOrds( tNumSetsToRegister );

    // iterate and add sets
    for ( moris::uint i = 0; i < tNumSetsToRegister; i++ )
    {
        tSideSetOrds( i ) = mSideSetLabels.size();

        mSideSetLabels.push_back( aSideSetNames( i ) );
        MORIS_ASSERT( mSideSideSetLabelToOrd.find( aSideSetNames( i ) ) == mSideSideSetLabelToOrd.end(),
            "Duplicate block set in mesh" );

        mSideSideSetLabelToOrd[aSideSetNames( i )] = tSideSetOrds( i );
    }

    mSideSets.resize( mSideSets.size() + tNumSetsToRegister );
    mSideSetColors.resize( mSideSetColors.size() + tNumSetsToRegister );

    return tSideSetOrds;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::set_side_set_colors(
    moris_index const &       aSideSetIndex,
    Matrix< IndexMat > const &aSideSetColors )
{

    MORIS_ASSERT( moris::isempty( mSideSetColors( aSideSetIndex ) ),
        "Attempting to overwrite colors of a side set" );

    mSideSetColors( aSideSetIndex ) = aSideSetColors;
}
//------------------------------------------------------------------------------
Cell< moris_index >
Enriched_Integration_Mesh::register_double_side_set_names( moris::Cell< std::string > const &aDblSideSetNames )
{
    uint tNumSetsToRegister = aDblSideSetNames.size();

    // block set ords
    Cell< moris_index > tDblSideSetOrds( tNumSetsToRegister );

    // iterate and add sets
    for ( moris::uint i = 0; i < tNumSetsToRegister; i++ )
    {
        tDblSideSetOrds( i ) = mDoubleSideSets.size() + i;

        mDoubleSideSetLabels.push_back( aDblSideSetNames( i ) );
        MORIS_ASSERT( mDoubleSideSetLabelToOrd.find( aDblSideSetNames( i ) ) == mDoubleSideSetLabelToOrd.end(), "Duplicate double side set in mesh" );
        mDoubleSideSetLabelToOrd[aDblSideSetNames( i )] = tDblSideSetOrds( i );
    }

    mDoubleSideSets.resize( mDoubleSideSets.size() + tNumSetsToRegister );
    mDoubleSideSetsMasterIndex.resize( mDoubleSideSetsMasterIndex.size() + tNumSetsToRegister );
    mDoubleSideSetsSlaveIndex.resize( mDoubleSideSetsSlaveIndex.size() + tNumSetsToRegister );
    mMasterDoubleSideSetColor.resize( mMasterDoubleSideSetColor.size() + tNumSetsToRegister );
    mSlaveDoubleSideSetColor.resize( mSlaveDoubleSideSetColor.size() + tNumSetsToRegister );

    return tDblSideSetOrds;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::set_double_side_set_colors(
    moris_index const &       aDblSideSetIndex,
    Matrix< IndexMat > const &aMasterSideColors,
    Matrix< IndexMat > const &aSlaveSideColors )
{
    MORIS_ASSERT( moris::isempty( mMasterDoubleSideSetColor( aDblSideSetIndex ) ),
        "Attempting to overwrite colors of a master side of double side set" );

    MORIS_ASSERT( moris::isempty( mSlaveDoubleSideSetColor( aDblSideSetIndex ) ),
        "Attempting to overwrite colors of a slave side of double side set" );

    mMasterDoubleSideSetColor( aDblSideSetIndex ) = aMasterSideColors;
    mSlaveDoubleSideSetColor( aDblSideSetIndex )  = aSlaveSideColors;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::setup_interface_side_sets()
{
    this->declare_interface_side_sets();

    this->create_interface_side_sets_and_clusters();
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::declare_interface_side_sets()
{
    uint tNumBulkPhases = mModel->get_geom_engine()->get_num_bulk_phase();

    Cell< std::string >        tInterfaceSideNames;
    Cell< Matrix< IndexMat > > tInterfaceSideColors;
    for ( moris::moris_index iP0 = 0; iP0 < (moris_index)tNumBulkPhases; iP0++ )
    {
        for ( moris::moris_index iP1 = 0; iP1 < (moris_index)tNumBulkPhases; iP1++ )
        {
            if ( iP1 != iP0 )
            {
                std::string tInterfaceSideSetName = get_interface_side_set_name( 0, iP0, iP1 );

                tInterfaceSideNames.push_back( tInterfaceSideSetName );
                tInterfaceSideColors.push_back( { { iP0 } } );
            }
        }
    }

    Cell< moris_index > tSideSetOrds = this->register_side_set_names( tInterfaceSideNames );

    // set interface side set colors
    for ( moris_index iSS = 0; iSS < (moris_index)tSideSetOrds.size(); iSS++ )
    {
        this->set_side_set_colors( tSideSetOrds( iSS ), tInterfaceSideColors( iSS ) );
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::create_interface_side_sets_and_clusters()
{
    uint tNumBulkPhases = mModel->get_geom_engine()->get_num_bulk_phase();

    // iterate through bulk phases
    for ( moris::moris_index iBP0 = 0; iBP0 < (moris_index)tNumBulkPhases; iBP0++ )
    {
        // iterate through bulk phase +1 (high to low
        for ( moris::moris_index iBP1 = iBP0 + 1; iBP1 < (moris_index)tNumBulkPhases; iBP1++ )
        {
            // create the interface side sets
            this->create_interface_side_sets_from_interface_double_side_set( iBP0, iBP1 );
        }
    }

    for ( moris::uint i = mListofSideSets.size(); i < mSideSets.size(); i++ )
    {
        this->commit_side_set( i );
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::construct_color_to_set_relationship(
    moris::Cell< moris::Matrix< IndexMat > > const &aSetColors,
    moris::Cell< moris::Cell< moris_index > > &     aColorToSetIndex )
{
    moris_index tMaxColor = 0;
    for ( moris::uint i = 0; i < aSetColors.size(); i++ )
    {
        moris_index tLocMax = aSetColors( i ).max();
        if ( tLocMax > tMaxColor )
        {
            tMaxColor = tLocMax;
        }
    }

    // size
    aColorToSetIndex.clear();
    aColorToSetIndex.resize( tMaxColor + 1 );

    for ( moris::uint i = 0; i < aSetColors.size(); i++ )
    {
        for ( moris::uint iC = 0; iC < aSetColors( i ).numel(); iC++ )
        {
            aColorToSetIndex( aSetColors( i )( iC ) ).push_back( (moris_index)i );
        }
    }
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::create_interface_side_sets_from_interface_double_side_set(
    moris_index const &aBulkphase0,
    moris_index const &aBulkphase1 )
{
    // get the side set names one going from 0->1 and another 1->0
    std::string tISideNameMasterToSlave = this->get_interface_side_set_name( 0, aBulkphase0, aBulkphase1 );
    std::string tISideNameSlaveToMaster = this->get_interface_side_set_name( 0, aBulkphase1, aBulkphase0 );

    // get the corresponding indices
    moris_index tISideIndexMasterToSlave = this->get_side_set_index( tISideNameMasterToSlave );
    moris_index tISideIndexSlaveToMaster = this->get_side_set_index( tISideNameSlaveToMaster );

    // get the double side set index
    moris_index tDblSideSetIndex = this->get_dbl_side_set_index( aBulkphase0, aBulkphase1 );

    moris::Cell< mtk::Double_Side_Cluster * > &tDblSideClusters = mDoubleSideSets( tDblSideSetIndex );

    // place the clusters in the two side sets
    for ( moris::uint i = 0; i < tDblSideClusters.size(); i++ )
    {
        // get the index
        moris_index tMasterIndex = mDoubleSideSetsMasterIndex( tDblSideSetIndex )( i );
        moris_index tSlaveIndex  = mDoubleSideSetsSlaveIndex( tDblSideSetIndex )( i );

        mSideSets( tISideIndexMasterToSlave ).push_back( mDoubleSideSingleSideClusters( tMasterIndex ) );
        mSideSets( tISideIndexSlaveToMaster ).push_back( mDoubleSideSingleSideClusters( tSlaveIndex ) );
    }
}

//------------------------------------------------------------------------------

Cell< moris_index >
Enriched_Integration_Mesh::declare_interface_vertex_sets()
{
    // number of geometries in the mesh
    moris::uint tNumGeometries = mModel->get_geom_engine()->get_num_geometries();

    // allocate a cell of strings
    moris::Cell< std::string > tInterfaceVertexSetNames( tNumGeometries );

    // base set name (interface vertex geometry #)
    std::string tSetNameBase = "iv_g_";

    for ( moris::uint i = 0; i < tNumGeometries; i++ )
    {
        // add the vertex set to the cell
        tInterfaceVertexSetNames( i ) = std::string( tSetNameBase + std::to_string( i ) );
    }

    // register vertex sets
    Cell< moris_index > tVertexSetOrds = this->register_vertex_set_names( tInterfaceVertexSetNames );

    // make the geometric index the color
    for ( moris::uint i = 0; i < tNumGeometries; i++ )
    {
        // the color of the interface ndoe sets is the geometric index
        this->set_vertex_set_color( tVertexSetOrds( i ), Matrix< IndexMat >( { { (moris_index)i } } ) );
    }

    return tVertexSetOrds;
}

//------------------------------------------------------------------------------

void
Enriched_Integration_Mesh::set_vertex_set_color(
    moris_index const &       aVertexSetIndex,
    Matrix< IndexMat > const &aVertexSetColors )
{
    MORIS_ASSERT( moris::isempty( mVertexSetColors( aVertexSetIndex ) ),
        "Attempting to overwrite colors of a side set" );

    mVertexSetColors( aVertexSetIndex ) = aVertexSetColors;
}

//------------------------------------------------------------------------------

bool
Enriched_Integration_Mesh::field_exists(
    std::string            aLabel,
    enum moris::EntityRank aEntityRank )
{
    moris::moris_index tIndex = this->get_entity_rank_field_index( aEntityRank );
    return mFieldLabelToIndex( tIndex ).find( aLabel ) != mFieldLabelToIndex( tIndex ).end();
}

//------------------------------------------------------------------------------

moris_index
Enriched_Integration_Mesh::get_entity_rank_field_index( enum moris::EntityRank aEntityRank )
{
    MORIS_ERROR( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT, "Only node and cell fields are supported" );

    moris_index tIndex = MORIS_INDEX_MAX;
    if ( aEntityRank == EntityRank::NODE )
    {
        tIndex = 0;
    }

    else if ( aEntityRank == EntityRank::ELEMENT )
    {
        tIndex = 1;
    }

    return tIndex;
}
}// namespace xtk
