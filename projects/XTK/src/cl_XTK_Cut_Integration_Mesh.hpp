#ifndef MORIS_cl_XTK_Cut_Integration_Mesh_HPP_
#define MORIS_cl_XTK_Cut_Integration_Mesh_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Vertex_XTK_Impl.hpp"
#include "cl_Tracer.hpp"

#include "cl_Communication_Tools.hpp"
#include <stdio.h>
#include <iostream>
#include <iomanip>
using namespace moris;

namespace xtk
{
struct IG_Cell_Group
{
    IG_Cell_Group( moris_index aNumCellsInGroup );

    IG_Cell_Group();

    moris::Cell< moris::mtk::Cell* > mIgCellGroup;
};

struct IG_Cell_Side_Group
{
    IG_Cell_Side_Group( moris_index aEstimatedNumCells );

    moris::Cell< moris::mtk::Cell* > mIgCells;                  // over allocated
    moris::Cell< moris_index >       mIgCellSideOrdinals;       // over allocated
};

struct IG_Cell_Double_Side_Group
{
    IG_Cell_Double_Side_Group( moris_index aEstimatedNumCells );

    moris::Cell< moris::mtk::Cell* > mLeaderIgCells;               // over allocated
    moris::Cell< moris_index >       mLeaderIgCellSideOrdinals;    // over allocated

    moris::Cell< moris::mtk::Cell* > mFollowerIgCells;             // over allocated
    moris::Cell< moris_index >       mFollowerIgCellSideOrdinals;  // over allocated

    void
    print()
    {
        std::cout << "Number of Leaders:   " << mLeaderIgCells.size() << std::endl;
        std::cout << "Number of Followers: " << mFollowerIgCells.size() << std::endl;

        int tStrLen = std::string( "Lead Cell Id   | " ).size();

        std::cout << "Lead Cell Id   | ";
        std::cout << "Side Ord       | ";
        std::cout << "Follow Cell Id | ";
        std::cout << "Side Ord       " << std::endl;
        // iterate through pairs
        for ( moris::uint i = 0; i < mLeaderIgCells.size(); i++ )
        {
            std::cout << std::setw( tStrLen ) << mLeaderIgCells( i )->get_id();
            std::cout << std::setw( tStrLen ) << mLeaderIgCellSideOrdinals( i );
            std::cout << std::setw( tStrLen ) << mFollowerIgCells( i )->get_id();
            std::cout << std::setw( tStrLen ) << mFollowerIgCellSideOrdinals( i ) << std::endl;
        }
    }
};

struct IG_Vertex_Group
{
    IG_Vertex_Group( moris_index aNumVerticesInGroup );

    std::size_t
    size();

    void
    reserve( std::size_t aReserveSize );

    void
    add_vertex(
        moris::mtk::Vertex const*                 aVertex,
        std::shared_ptr< Matrix< DDRMat > > aVertexLocalCoord );

    void
    add_vertex_local_coord_pointers();

    moris::mtk::Vertex const *
    get_vertex( moris_index aGroupVertexOrdinal );

    moris_index
    get_vertex_group_ordinal( moris_index aVertex );

    std::shared_ptr< Matrix< DDRMat > >
    get_vertex_local_coords( moris_index aVertex );

    moris::uint
    get_vertex_local_coords_dim() const;

    bool
    vertex_is_in_group(moris_index aVertex);

    void
    print();

  private:
    moris::Cell< moris::mtk::Vertex const * >          mIgVertexGroup;
    std::unordered_map< moris_index, moris_index >     mIgVertexIndexToVertexOrdinal;
    moris::Cell< std::shared_ptr< Matrix< DDRMat > > > mIgVertexLocalCoords;
};

struct Edge_Based_Connectivity
{
    moris::Cell< moris::Cell< moris::mtk::Vertex* > > mEdgeVertices;          // input: edge || output: list of vertices on edge
    moris::Cell< moris::Cell< moris::mtk::Cell* > >   mEdgeToCell;            // input: edge || output: list of cells attached to edge
    moris::Cell< moris::Cell< moris::moris_index > >  mEdgeToCellEdgeOrdinal; // input: edge || output: ?
    moris::Cell< moris::Cell< moris_index > >         mCellToEdge;            // input: cell || output: list of edge indices on cell
};

struct Edge_Based_Ancestry
{
    moris::Cell< moris::moris_index > mEdgeParentEntityIndex;
    moris::Cell< moris::moris_index > mEdgeParentEntityRank;
    moris::Cell< moris::moris_index > mEdgeParentEntityOrdinalWrtBackgroundCell;
};

struct Vertex_Ancestry
{
    Vertex_Ancestry() = default;
    Vertex_Ancestry(
        moris::Cell< moris::moris_index > const& aVertexParentEntityIndices,
        moris::Cell< enum EntityRank > const&    aVertexParentEntityRank ) :
        mVertexParentEntityIndex( aVertexParentEntityIndices ),
        mVertexParentEntityRank( aVertexParentEntityRank )
    {
    }

    moris::Cell< moris::moris_index > mVertexParentEntityIndex;
    moris::Cell< enum EntityRank >    mVertexParentEntityRank;

    moris_index
    get_vertex_parent_index( moris_index aVertexIndex ) const
    {
        return mVertexParentEntityIndex( aVertexIndex );
    }

    enum EntityRank
    get_vertex_parent_rank( moris_index aVertexIndex ) const
    {
        return mVertexParentEntityRank( aVertexIndex );
    }
};

struct Facet_Based_Connectivity
{
    moris::Cell< moris::Cell< moris::mtk::Vertex* > > mFacetVertices;            // over allocated
    moris::Cell< moris::Cell< moris::mtk::Cell* > >   mFacetToCell;              // over allocated
    moris::Cell< moris::Cell< moris::moris_index > >  mFacetToCellEdgeOrdinal;   // over allocated
    moris::Cell< moris::Cell< moris_index > >         mCellToFacet;              // over allocated
    std::unordered_map< moris_index, moris_index >    mCellIndexToCellOrdinal;   // over allocated

    moris_index
    get_cell_ordinal( const moris_index& aCellIndex )
    {
        auto tIter = mCellIndexToCellOrdinal.find( aCellIndex );
        MORIS_ASSERT( tIter != mCellIndexToCellOrdinal.end(), "Cell not in facet connectivity" );
        return tIter->second;
    }
};

struct Facet_Based_Ancestry
{
    moris::Cell< moris::moris_index > mFacetParentEntityIndex;
    moris::Cell< moris::moris_index > mFacetParentEntityRank;
    moris::Cell< moris::moris_index > mFacetParentEntityOrdinalWrtBackgroundCell;
};

struct Cell_Neighborhood_Connectivity
{
    moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Cell* > > > mNeighborCells;
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > >       mMySideOrdinal;
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > >       mNeighborSideOrdinal;
};

struct Cell_Connectivity
{
    Cell_Connectivity(){};
    Cell_Connectivity(
        Matrix< IndexMat > const& aCellVertexInds,
        Matrix< IndexMat > const& aCellEdgesInds,
        Matrix< IndexMat > const& aCellFacesInds ) :
        mCellVertexInds( aCellVertexInds ),
        mCellEdgesInds( aCellEdgesInds ), mCellFacesInds( aCellFacesInds )
    {
    }

    moris_index
    get_entity_index(
        moris_index     aEntityOrdinal,
        enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::FACE:
            return mCellFacesInds( aEntityOrdinal );
            break;
        case EntityRank::EDGE:
            return mCellEdgesInds( aEntityOrdinal );
            break;
        case EntityRank::NODE:
            return mCellVertexInds( aEntityOrdinal );
            break;
        default:
            MORIS_ERROR( 0, "INVALID RANK" );
            return 0;
            break;
        }
    }

    moris_index
    get_entity_ordinal(
        moris_index const& aEntityIndex,
        enum EntityRank    aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::FACE:
            for ( moris::uint iEnt = 0; iEnt < mCellFacesInds.numel(); iEnt++ )
            {
                if ( mCellFacesInds( iEnt ) == aEntityIndex )
                {
                    return iEnt;
                }
            }
            MORIS_ERROR( 0, "Face not found" );
            return MORIS_INDEX_MAX;
            break;
        case EntityRank::EDGE:
            for ( moris::uint iEnt = 0; iEnt < mCellEdgesInds.numel(); iEnt++ )
            {
                if ( mCellEdgesInds( iEnt ) == aEntityIndex )
                {
                    return iEnt;
                }
            }
            MORIS_ERROR( 0, "Edge not found" );
            return MORIS_INDEX_MAX;
            break;
        case EntityRank::NODE:
            for ( moris::uint iEnt = 0; iEnt < mCellVertexInds.numel(); iEnt++ )
            {
                if ( mCellVertexInds( iEnt ) == aEntityIndex )
                {
                    return iEnt;
                }
            }
            MORIS_ERROR( 0, "Node not found" );
            return MORIS_INDEX_MAX;
            break;
        default:
            MORIS_ERROR( 0, "INVALID RANK" );
            return 0;
            break;
        }
    }

    const Matrix< IndexMat > mCellVertexInds;
    const Matrix< IndexMat > mCellEdgesInds;
    const Matrix< IndexMat > mCellFacesInds;
};

struct Subphase_Neighborhood_Connectivity
{
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mSubphaseToSubPhase;
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mSubphaseToSubPhaseMySideOrds;
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mSubphaseToSubPhaseNeighborSideOrds;
    moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mTransitionNeighborCellLocation;

    void
    print_subphase_neighborhood()
    {

        std::cout << "Subphases" << std::endl;
        for ( moris::uint iC = 0; iC < mSubphaseToSubPhase.size(); iC++ )
        {
            std::cout << std::setw( 6 ) << iC << " | ";

            for ( moris::uint iN = 0; iN < mSubphaseToSubPhase( iC )->size(); iN++ )
            {
                std::cout << std::setw( 6 ) << (*mSubphaseToSubPhase( iC ))( iN );
            }
            std::cout << std::endl;
        }

        std::cout << "Subphases My Side Ordinals" << std::endl;
        for ( moris::uint iC = 0; iC < mSubphaseToSubPhaseMySideOrds.size(); iC++ )
        {
            std::cout << std::setw( 6 ) << iC << " | ";

            for ( moris::uint iN = 0; iN < mSubphaseToSubPhaseMySideOrds( iC )->size(); iN++ )
            {
                std::cout << std::setw( 6 ) << (*mSubphaseToSubPhaseMySideOrds( iC ))( iN );
            }
            std::cout << std::endl;
        }

        std::cout << "Subphases Neighbor Side Ordinals" << std::endl;
        for ( moris::uint iC = 0; iC < mSubphaseToSubPhaseNeighborSideOrds.size(); iC++ )
        {
            std::cout << std::setw( 6 ) << iC << " | ";

            for ( moris::uint iN = 0; iN < mSubphaseToSubPhaseNeighborSideOrds( iC )->size(); iN++ )
            {
                std::cout << std::setw( 6 ) << (*mSubphaseToSubPhaseNeighborSideOrds( iC ))( iN );
            }
            std::cout << std::endl;
        }

        std::cout << "Transition Neighbor Locations" << std::endl;
        for ( moris::uint iC = 0; iC < mTransitionNeighborCellLocation.size(); iC++ )
        {
            std::cout << std::setw( 6 ) << iC << " | ";

            for ( moris::uint iN = 0; iN < mTransitionNeighborCellLocation( iC )->size(); iN++ )
            {
                std::cout << std::setw( 12 ) << (*mTransitionNeighborCellLocation( iC ))( iN );
            }
            std::cout << std::endl;
        }
    }
};

class Child_Mesh_Experimental;
class Model;
class Cell_XTK_No_CM;
class Cut_Integration_Mesh : public moris::mtk::Mesh
{
    friend class Integration_Mesh_Generator;
    friend class Model;

  protected:
    moris::uint mSpatialDimension;

    bool mSameLevelChildMeshes = true;

    // integration cells
    moris_index                                           mFirstControlledCellIndex;
    moris::Cell< moris::mtk::Cell* >                      mIntegrationCells;
    moris::Cell< std::shared_ptr< xtk::Cell_XTK_No_CM > > mControlledIgCells;

    // quantities related to integration cells
    moris::Cell< moris::Cell< moris_index > > mIntegrationCellToCellGroupIndex;
    moris::Cell< moris_index >                mIntegrationCellToSubphaseIndex;
    moris::Cell< moris::moris_index >         mIntegrationCellBulkPhase;

    // integration vertices
    moris_index                                              mFirstControlledVertexIndex;
    moris::Cell< moris::mtk::Vertex* >                       mIntegrationVertices;
    moris::Cell< std::shared_ptr< moris::mtk::Vertex_XTK > > mControlledIgVerts;

    // vertex ancestry
    moris::Cell< moris::moris_index > mIgVertexParentEntityIndex;
    moris::Cell< moris::moris_index > mIgVertexParentEntityRank;
    moris::Cell< moris::moris_index > mIgVertexConnectedCell;

    // vertex quantities
    moris::Cell< std::shared_ptr< Matrix< DDRMat > > > mVertexCoordinates;

    // all data is stored in the current mesh. pointers are in the child mesh
    // as well as accessor functions are provided there
    moris::Cell< std::shared_ptr< Child_Mesh_Experimental > > mChildMeshes;
    moris::Cell< std::shared_ptr< Child_Mesh_Experimental > > mOwnedChildMeshes;
    moris::Cell< std::shared_ptr< Child_Mesh_Experimental > > mNotOwnedChildMeshes;

    // communication map
    moris::Matrix< IdMat > mCommunicationMap;

    // group of all integration cells in a single parent cell
    moris::Cell< std::shared_ptr< IG_Cell_Group > >   mIntegrationCellGroups;
    moris::Cell< std::shared_ptr< IG_Vertex_Group > > mIntegrationVertexGroups;
    moris::Cell< moris::mtk::Cell* >                  mIntegrationCellGroupsParentCell;
    moris::Cell< moris_index >                        mParentCellCellGroupIndex;

    moris::Cell< moris_index > mOwnedIntegrationCellGroupsInds;
    moris::Cell< moris_index > mNotOwnedIntegrationCellGroups;

    // subphase groupings
    moris::Cell< moris_index >                      mSubPhaseIds;
    moris::Cell< std::shared_ptr< IG_Cell_Group > > mSubPhaseCellGroups;
    moris::Cell< moris::moris_index >               mSubPhaseBulkPhase;
    moris::Cell< moris::mtk::Cell* >                mSubPhaseParentCell;
    moris::Cell< moris::Cell< moris_index > >       mParentCellToSubphase;
    moris::Cell< moris_index >                      mParentCellHasChildren;

    moris::Cell< moris_index >                                mOwnedSubphaseGroupsInds;
    moris::Cell< moris_index >                                mNotOwnedSubphaseGroupsInds;
    std::unordered_map< moris::moris_id, moris::moris_index > mGlobalToLocalSubphaseMap;

    std::shared_ptr< Subphase_Neighborhood_Connectivity > mSubphaseNeighborhood;

    // face connectivity
    std::shared_ptr< Facet_Based_Connectivity > mIgCellFaceConnectivity;

    // face ancestry
    std::shared_ptr< Facet_Based_Ancestry > mIgCellFaceAncestry;

    // interface facets - indexed based on mIgCellFaceConnectivity facet indices
    moris::Cell< moris_index > mInterfaceFacets;

    // double side interface groups
    // outer cell - bulk phase 0
    // inner cell - bulk phase 1
    // IG_Cell_Double_Side_Group pairings between integration cells
    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > mBptoBpDblSideInterfaces;

    // background facet to child facet connectivity
    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > mBGFacetToChildFacet;

    // block set data
    std::unordered_map< std::string, moris_index >  mBlockSetLabelToOrd;
    moris::Cell< std::string >                      mBlockSetNames;
    moris::Cell< std::shared_ptr< IG_Cell_Group > > mBlockSetCellGroup;
    moris::Cell< enum CellTopology >                mBlockCellTopo;

    // Side Set Data
    std::unordered_map< std::string, moris_index >       mSideSideSetLabelToOrd;
    Cell< std::string >                                  mSideSetLabels;
    moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > mSideSetCellSides;

    // connectivity from vertex to child mesh
    // outer cell vertex
    // inner cell child meshes associated with the vertex
    moris::Cell< moris::Cell< moris_index > > mVertexToChildMeshIndex;

    // connectivity from vertex to child mesh
    // outer cell integration cell index
    // inner cell child meshes associated with the integration cell
    moris::Cell< moris::Cell< moris_index > > mCellToChildMeshIndex;

    // outer cell geomtry index
    // if vertex index is in the map then it is a member of the geometric interface
    moris::Cell< std::unordered_map< moris_index, moris_index > > mGeometryInterfaceVertexIndices;

    moris::moris_index mGlobalMaxVertexId;
    moris::moris_index mGlobalMaxCellId;

    std::unordered_map< moris_id, moris_index > mIntegrationCellIdToIndexMap;
    std::unordered_map< moris_id, moris_index > mIntegrationVertexIdToIndexMap;

    moris::Cell< moris::moris_index > mIntegrationCellIndexToId;
    moris::Cell< moris::moris_index > mIntegrationVertexIndexToId;

    moris::mtk::Mesh*                        mBackgroundMesh;
    Model*                                   mXTKModel;

  public:

    // ----------------------------------------------------------------------------------

    Cut_Integration_Mesh(
        moris::mtk::Mesh* aBackgroundMesh,
        Model*            aXTKModel );

    // ----------------------------------------------------------------------------------

    ~Cut_Integration_Mesh();

    // ----------------------------------------------------------------------------------

    // Core Mesh Functions
    uint
    get_spatial_dim() const;

    // ----------------------------------------------------------------------------------

    MeshType
    get_mesh_type() const;

    // ----------------------------------------------------------------------------------

    uint
    get_num_entities(
        enum EntityRank   aEntityRank,
        const moris_index aIndex ) const;

    // ----------------------------------------------------------------------------------

    moris::uint get_num_sets() const;

    // ----------------------------------------------------------------------------------

    Matrix< DDRMat >
    get_node_coordinate( moris_index aNodeIndex ) const;

    // ----------------------------------------------------------------------------------

    uint
    get_node_owner( moris_index aNodeIndex ) const;

    // ----------------------------------------------------------------------------------

    uint
    get_element_owner( moris_index aElementIndex ) const;

    // ----------------------------------------------------------------------------------

    Matrix< IdMat >
    get_communication_table() const;

    // ----------------------------------------------------------------------------------

    void
    add_proc_to_comm_table(moris_index aProcRank);

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex );

    // ----------------------------------------------------------------------------------

    enum CellTopology
    get_blockset_topology( const std::string& aSetName );

    // ----------------------------------------------------------------------------------

    enum CellShape
    get_IG_blockset_shape( const std::string& aSetName );

    // ----------------------------------------------------------------------------------

    enum CellShape
    get_IP_blockset_shape( const std::string& aSetName );

    // ----------------------------------------------------------------------------------

    moris_id
    get_glb_entity_id_from_entity_loc_index(
        moris_index       aEntityIndex,
        enum EntityRank   aEntityRank,
        const moris_index aDiscretizationIndex = 0 ) const;

    // ----------------------------------------------------------------------------------

    moris_index
    get_loc_entity_ind_from_entity_glb_id(
        moris_id        aEntityId,
        enum EntityRank aEntityRank ) const;

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    get_entity_connected_to_entity_loc_inds(
        moris_index       aEntityIndex,
        enum EntityRank   aInputEntityRank,
        enum EntityRank   aOutputEntityRank,
        const moris_index aDiscretizationIndex = 0 ) const;

    // ----------------------------------------------------------------------------------

    moris::Cell< std::string >
    get_set_names( enum EntityRank aSetEntityRank ) const;

    // ----------------------------------------------------------------------------------

    moris_index
    get_block_set_index( std::string aBlockSetLabel ) const;

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    get_block_entity_loc_inds( std::string aSetName ) const;

    // ----------------------------------------------------------------------------------

    moris_index
    get_side_set_index( std::string aSideSetLabel ) const;

    // ----------------------------------------------------------------------------------

    void
    get_sideset_elems_loc_inds_and_ords(
        const std::string&  aSetName,
        Matrix< IndexMat >& aElemIndices,
        Matrix< IndexMat >& aSidesetOrdinals ) const;

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    get_set_entity_loc_inds(
        enum EntityRank aSetEntityRank,
        std::string     aSetName ) const;

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< Matrix< DDRMat > > >*
    get_all_vertex_coordinates_loc_inds();

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< IG_Cell_Group > >&
    get_all_cell_groups();

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Child_Mesh_Experimental >
    get_child_mesh( moris_index aChildMeshIndex );

    // ----------------------------------------------------------------------------------

    uint
    get_num_child_meshes() const;

    // ----------------------------------------------------------------------------------

    std::unordered_map< moris_id, moris_index >
    get_vertex_glb_id_to_loc_vertex_ind_map() const;

    // ----------------------------------------------------------------------------------

    bool
    vertex_exists( moris_index tId ) const;

    // ----------------------------------------------------------------------------------

    mtk::Cell const&
    get_mtk_cell( moris_index aElementIndex ) const;

    // ----------------------------------------------------------------------------------

    mtk::Cell&
    get_mtk_cell( moris_index aElementIndex );

    // ----------------------------------------------------------------------------------

    mtk::Vertex&
    get_mtk_vertex( moris_index aVertexIndex );

    // ----------------------------------------------------------------------------------

    mtk::Vertex const&
    get_mtk_vertex( moris_index aVertexIndex ) const;

    // ----------------------------------------------------------------------------------

    moris::mtk::Vertex*
    get_mtk_vertex_pointer( moris_index aVertexIndex );

    // ----------------------------------------------------------------------------------

    std::shared_ptr< IG_Vertex_Group >
    get_vertex_group( moris_index aVertexGroupIndex );

    // ----------------------------------------------------------------------------------

    moris_index
    get_parent_cell_group_index( moris_index aParentCellIndex );

    // ----------------------------------------------------------------------------------

    void
    replace_controlled_ig_cell(
        moris_index                              aCellIndex,
        moris_id                                 aCellId,
        std::shared_ptr< moris::mtk::Cell_Info > aCellInfo,
        moris::Cell< moris::mtk::Vertex* >&      aVertexPointers );

    // ----------------------------------------------------------------------------------

    void
    set_integration_cell(
        moris_index                            aCellIndex,
        std::shared_ptr< xtk::Cell_XTK_No_CM > aNewCell );

    // ----------------------------------------------------------------------------------

    void
    add_integration_cell(
        moris_index                            aCellIndex,
        std::shared_ptr< xtk::Cell_XTK_No_CM > aNewCell );

    // ----------------------------------------------------------------------------------

    moris_index
    get_integration_cell_controlled_index(
        moris_index aCellIndex );

    // ----------------------------------------------------------------------------------

    void
    add_cell_to_cell_group(
        moris_index aCellIndex,
        moris_index aCellGroupIndex );

    // ----------------------------------------------------------------------------------

    moris_id
    allocate_entity_ids( moris::size_t aNumIdstoAllocate,
        enum EntityRank                aEntityRank );

    // ----------------------------------------------------------------------------------

    moris_id
    allocate_subphase_ids( moris::size_t aNumIdstoAllocate );

    // ----------------------------------------------------------------------------------

    moris::moris_index
    get_first_available_index( enum EntityRank aEntityRank ) const;

    // ----------------------------------------------------------------------------------

    moris::uint
    get_num_ig_cell_groups();

    // ----------------------------------------------------------------------------------

    std::shared_ptr< IG_Cell_Group >
    get_ig_cell_group( moris_index aGroupIndex );

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const&
    get_ig_cell_group_memberships( moris_index aIgCellIndex );

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell*
    get_ig_cell_group_parent_cell( moris_index aGroupIndex );

    // ----------------------------------------------------------------------------------

    enum CellTopology 
    get_child_element_topology();
    void

    // ----------------------------------------------------------------------------------

    set_child_mesh_subphase(
        moris_index                 aCMIndex,
        moris::Cell< moris_index >& aSubphasesGroups );

    // ----------------------------------------------------------------------------------

    moris::uint
    get_num_subphases();

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< Child_Mesh_Experimental > >&
    get_owned_child_meshes();

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index >&
    get_owned_subphase_indices();

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index >&
    get_not_owned_subphase_indices();

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell*
    get_subphase_parent_cell( moris_index aSubPhaseIndex );

    // ----------------------------------------------------------------------------------

    std::shared_ptr< IG_Cell_Group >
    get_subphase_ig_cells( moris_index aSubPhaseIndex );

    // ----------------------------------------------------------------------------------

    moris_index
    get_subphase_id( moris_index aSubPhaseIndex );

    // ----------------------------------------------------------------------------------

    moris_index
    get_subphase_index( moris_id aSubphaseId );

    // ----------------------------------------------------------------------------------

    moris_index
    get_subphase_bulk_phase( moris_index aSubPhaseIndex );

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const&
    get_parent_cell_subphases( moris_index aParentCellIndex );

    // ----------------------------------------------------------------------------------

    moris_index
    get_ig_cell_subphase_index( moris_index aIgCellIndex );

    // ----------------------------------------------------------------------------------

    bool
    parent_cell_has_children( moris_index aParentCellIndex );

    // ----------------------------------------------------------------------------------

    moris_index
    get_vertex_parent_index( moris_index const& aVertexIndex );

    // ----------------------------------------------------------------------------------

    moris_index
    get_vertex_parent_rank( moris_index const& aVertexIndex );

    // ----------------------------------------------------------------------------------

    void
    finalize_cut_mesh_construction();

    // ----------------------------------------------------------------------------------

    void
    deduce_ig_cell_group_ownership();

    // ----------------------------------------------------------------------------------

    void
    assign_controlled_ig_cell_ids();

    // ----------------------------------------------------------------------------------

    void
    set_face_connectivity( std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity );

    // ----------------------------------------------------------------------------------

    void
    set_face_ancestry( std::shared_ptr< Facet_Based_Ancestry > aFaceAncestry );

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Facet_Based_Connectivity >
    get_face_connectivity();

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Facet_Based_Ancestry >
    get_face_ancestry();

    // ----------------------------------------------------------------------------------

    void
    set_interface_facets( moris::Cell< moris_index >& aInterfaces );

    // ----------------------------------------------------------------------------------

    moris::Cell< moris_index > const&
    get_interface_facets();

    // ----------------------------------------------------------------------------------

    void
    set_bulk_phase_to_bulk_phase_dbl_side_interface( moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > >& aBptoBpDblSideInterfaces );

    // ----------------------------------------------------------------------------------

    moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > const&
    get_bulk_phase_to_bulk_phase_dbl_side_interface();

    // ----------------------------------------------------------------------------------

    void
    set_background_facet_to_child_facet_connectivity( moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > const& aBgtoChildFacet );

    // ----------------------------------------------------------------------------------

    moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > const&
    get_background_facet_to_child_facet_connectivity();

    // ----------------------------------------------------------------------------------

    void
    set_subphase_neighborhood( std::shared_ptr< Subphase_Neighborhood_Connectivity > aSubphaseNeighborhood );

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Subphase_Neighborhood_Connectivity >
    get_subphase_neighborhood();

    // ----------------------------------------------------------------------------------

    void
    setup_glob_to_loc_subphase_map();

    // ----------------------------------------------------------------------------------

    moris_index
    get_cell_bulk_phase( moris_index aCellIndex );

    // ----------------------------------------------------------------------------------

    Cell< moris_index >
    register_side_set_names( moris::Cell< std::string > const& aSideSetNames );

    // ----------------------------------------------------------------------------------

    Cell< moris_index >
    register_block_set_names( moris::Cell< std::string > const& aBlockSetNames,
        enum CellTopology                                       aCellTopo );

    // ----------------------------------------------------------------------------------

    void
    write_mesh( std::string aOutputPath,
        std::string         aOutputFile );

    // ----------------------------------------------------------------------------------

    Cell_Connectivity
    get_background_cell_connectivity( moris_index aBGCellIndex ) const;

    // ----------------------------------------------------------------------------------

    void
    print()
    {
        this->print_cells();
        this->print_vertices();
        this->print_block_sets();
        this->print_groupings();
        this->print_vertex_ancestry();
    }

    // ----------------------------------------------------------------------------------

    void
    print_cells(
        bool        aOmitIndex = false,
        std::string aFile      = "" );

    // ----------------------------------------------------------------------------------

    void
    print_vertices(
        bool        aOmitIndex = false,
        std::string aFile      = "" );

    // ----------------------------------------------------------------------------------

    void
    print_block_sets()
    {
        const char separator = ' ';
        const int  nameWidth = 24;
        const int  numWidth  = 12;

        std::cout << "Num Block Sets: " << this->mBlockSetCellGroup.size() << std::endl;

        std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << "Name";
        std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << "Ordinal";
        std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << "Num Cells";
        std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << "Cells";
        std::cout << std::endl;

        for ( uint iBS = 0; iBS < mBlockSetCellGroup.size(); iBS++ )
        {
            std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << mBlockSetNames( iBS );
            std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << mBlockSetLabelToOrd.find( mBlockSetNames( iBS ) )->second;
            std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << mBlockSetCellGroup( iBS )->mIgCellGroup.size();

            //
            for ( moris::uint iIgCell = 0; iIgCell < mBlockSetCellGroup( iBS )->mIgCellGroup.size(); iIgCell++ )
            {
                std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << mBlockSetCellGroup( iBS )->mIgCellGroup( iIgCell )->get_id();
            }

            std::cout << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    print_groupings( std::string aFile = "" );

    // ----------------------------------------------------------------------------------

    void
    print_vertex_ancestry() const
    {
        moris::print( mIgVertexParentEntityIndex, "Vertex Parent Entity Index" );
        moris::print( mIgVertexParentEntityRank, "Vertex Parent Entity Rank" );

        // moris::Cell< moris::moris_index > mIgVertexParentEntityIndex;
        // moris::Cell< moris::moris_index > mIgVertexParentEntityRank;
    }

    // ----------------------------------------------------------------------------------

    void
    trim_data();

    // ----------------------------------------------------------------------------------

  private:

    // ----------------------------------------------------------------------------------

    void
    create_base_cell_blocks();

    // ----------------------------------------------------------------------------------

    void
    setup_comm_map();

    // ----------------------------------------------------------------------------------

    void
    prepare_child_element_identifier_requests(
        moris::Cell< moris::Cell< moris_id > >&   aNotOwnedChildMeshesToProcs,
        moris::Cell< moris::Matrix< IdMat > >&    aOwnedParentCellId,
        moris::Cell< moris::Matrix< IdMat > >&    aNumOwnedCellIdsOffsets,
        moris::Cell< uint >&                      aProcRanks,
        std::unordered_map< moris_id, moris_id >& aProcRankToDataIndex );

    // ----------------------------------------------------------------------------------

    void
    prepare_child_cell_id_answers(
        Cell< Matrix< IndexMat > >& aReceivedParentCellIds,
        Cell< Matrix< IndexMat > >& aReceivedParentCellNumChildren,
        Cell< Matrix< IndexMat > >& aChildCellIdOffset );

    // ----------------------------------------------------------------------------------

    void
    handle_received_child_cell_id_request_answers(
        Cell< Cell< moris_index > > const& aChildMeshesInInNotOwned,
        Cell< Matrix< IndexMat > > const&  aReceivedChildCellIdOffset );

    // ----------------------------------------------------------------------------------

};

class Child_Mesh_Experimental
{
    friend class Cut_Integration_Mesh;
    friend class Integration_Mesh_Generator;

  public:
    moris::mtk::Cell*                  mParentCell;
    moris::moris_index                 mChildMeshIndex;
    std::shared_ptr< IG_Cell_Group >   mIgCells;
    std::shared_ptr< IG_Vertex_Group > mIgVerts;

    // subphases
    moris::Cell< std::shared_ptr< IG_Cell_Group > > mSubphaseCellGroups;

  public:
    Child_Mesh_Experimental()
    {
    }

    moris_index
    get_parent_element_index()
    {
        return this->get_parent_cell()->get_index();
    }

    moris::mtk::Cell*
    get_parent_cell()
    {
        return mParentCell;
    }

    moris::moris_index
    get_child_mesh_index()
    {
        return mChildMeshIndex;
    }

    void
    set_subphase_groups( moris::Cell< std::shared_ptr< IG_Cell_Group > >& aSubphasesGroups )
    {
        mSubphaseCellGroups = aSubphasesGroups;
    }

    uint
    get_num_subphase_cell_groups() const
    {
        return mSubphaseCellGroups.size();
    }

    std::shared_ptr< IG_Cell_Group >
    get_subphase_cell_group(moris_index aLocalSpIndex)
    {
        return mSubphaseCellGroups(aLocalSpIndex);
    }

};

}// namespace xtk

#endif
