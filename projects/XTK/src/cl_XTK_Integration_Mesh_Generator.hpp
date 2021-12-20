#ifndef MORIS_CL_XTK_Integration_Mesh_Generator_HPP_
#define MORIS_CL_XTK_Integration_Mesh_Generator_HPP_

#include "cl_XTK_Model.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_GEN_Geometric_Query_Interface.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_Tracer.hpp"

namespace xtk
{

struct Integration_Mesh_Generation_Data
{
    // number of child meshes
    moris_index mNumChildMeshes = 0;

    // outer cell geometry index (if the geometry is inactive the cell is empty)
    // inner cell active child mesh index by cut mesh index
    moris::Cell< moris::Cell< moris_index > > mIntersectedBackgroundCellIndex;

    // All intersected backgroun cells (uniques removed from the concatenated version of mIntersectedBackgroundCellIndex)
    moris::Cell< moris_index > mAllIntersectedBgCellInds;


    // // this maps from the background cell index to the child mesh index
    // std::unordered_map< moris_index, moris_index > mIntersectedBackgroundCellIndexToChildMeshIndex;
};

class Geometric_Query_XTK : public moris::ge::Geometric_Query_Interface
{
  private:
    // tell the geometry engine what you are interested
    enum moris::ge::Query_Type mQueryType;

    // associated with a given child mesh
    enum EntityRank    mQueryEntityRank;
    Matrix< IndexMat > mQueryEntityToVertices;
    Matrix< IndexMat > mQueryEntityParamCoords;

    // parent cell is the query object when checking for intersection
    moris::mtk::Cell *mQueryParentCell;

    // for edge based questions
    moris_index                                        mCurrentEdgeIndex;
    std::shared_ptr< Edge_Based_Connectivity >         mEdgeConnectivity;
    moris::Cell< moris::mtk::Cell * > *                mAssociatedBgCellForEdge;
    moris::Cell< std::shared_ptr< IG_Vertex_Group > > *mAssociatedVertexGroup;
    Cut_Integration_Mesh *                             mCutIntegrationMesh;

    // parent info
    moris::mtk::Cell *mParentCell;

    // row-based and indexed based on the indexing in mQueryEntityToVertices.
    moris::Cell< std::shared_ptr< moris::Matrix< moris::DDRMat > > > *mQueryEntityIndexedCoordinates;

    // geometric index
    moris_index mGeometricIndex;

  public:
    Geometric_Query_XTK() :
        mQueryParentCell( nullptr ), mParentCell( nullptr ), mGeometricIndex( MORIS_INDEX_MAX )
    {
    }
    ~Geometric_Query_XTK() {}

    void
    set_query_type( enum moris::ge::Query_Type aQueryType )
    {
        mQueryType = aQueryType;
    }

    void
    set_query_entity_rank( enum moris::EntityRank aEntityRank )
    {
        mQueryEntityRank = aEntityRank;
    }

    void
    set_edge_connectivity( std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity )
    {
        mEdgeConnectivity = aEdgeConnectivity;
    }

    void
    set_edge_associated_background_cell( moris::Cell< moris::mtk::Cell * > *aAssociatedBgCellsForEdge )
    {
        mAssociatedBgCellForEdge = aAssociatedBgCellsForEdge;
    }

    void
    set_associated_vertex_group( moris::Cell< std::shared_ptr< IG_Vertex_Group > > *aAssociatedVertexGroup )
    {
        mAssociatedVertexGroup = aAssociatedVertexGroup;
    }

    void
    set_cut_integration_mesh( Cut_Integration_Mesh *aCutIntegrationMesh )
    {
        mCutIntegrationMesh = aCutIntegrationMesh;
    }

    void
    set_current_edge_index( moris_index aCurrentIndex )
    {
        mCurrentEdgeIndex = aCurrentIndex;
        mQueryEntityToVertices.resize( 1, 2 );
        mQueryEntityToVertices( 0 ) = mEdgeConnectivity->mEdgeVertices( mCurrentEdgeIndex )( 0 )->get_index();
        mQueryEntityToVertices( 1 ) = mEdgeConnectivity->mEdgeVertices( mCurrentEdgeIndex )( 1 )->get_index();
    }

    void
    set_parametric_coordinate_size( moris_index aParametricSize )
    {
        // hard coded for edges right now
        mQueryEntityParamCoords.resize( 2, aParametricSize );
    }

    void
    set_parent_cell( moris::mtk::Cell *aParentCell )
    {
        mParentCell = aParentCell;
    }

    void
    set_query_cell( moris::mtk::Cell *aQueryParentCell )
    {
        mQueryParentCell       = aQueryParentCell;
        mQueryEntityToVertices = aQueryParentCell->get_vertex_inds();
    }

    void
    set_query_edge_connectivity( std::shared_ptr< Edge_Based_Connectivity > aEdgeBasedConnectivity )
    {
    }

    void
    set_coordinates_matrix( moris::Cell< std::shared_ptr< moris::Matrix< moris::DDRMat > > > *aCoordinates )
    {
        mQueryEntityIndexedCoordinates = aCoordinates;
    }

    void
    set_geometric_index( moris_index aGeometricIndex )
    {
        mGeometricIndex = aGeometricIndex;
    }

    enum moris::ge::Query_Type
    get_query_type() const
    {
        return mQueryType;
    }

    moris_index
    get_geomtric_index() const
    {
        return mGeometricIndex;
    }

    enum EntityRank
    get_query_entity_rank() const
    {
        return mQueryEntityRank;
    }

    Matrix< IndexMat > const &
    get_query_entity_to_vertex_connectivity() const
    {
        return mQueryEntityToVertices;
    }

    moris::Cell< std::shared_ptr< moris::Matrix< moris::DDRMat > > > *
    get_query_indexed_coordinates() const
    {
        return mQueryEntityIndexedCoordinates;
    }

    Matrix< DDRMat >
    get_vertex_local_coord_wrt_parent_entity( moris_index aVertexIndex ) const
    {
        std::shared_ptr< IG_Vertex_Group > tCurrentVertexGroup = ( *mAssociatedVertexGroup )( mCurrentEdgeIndex );
        return *tCurrentVertexGroup->get_vertex_local_coords( aVertexIndex );
    }

    enum EntityRank
    get_query_parent_entity_rank() const
    {
        return EntityRank::ELEMENT;
    }

    Matrix< IndexMat >
    get_query_parent_entity_connectivity() const
    {
        return mParentCell->get_vertex_inds();
    }

    Matrix< DDRMat >
    get_query_parent_coordinates() const
    {
        return mParentCell->get_vertex_coords();
    }

    moris_index
    get_query_parent_entity_id() const
    {
        return mParentCell->get_id();
    }

    moris_index
    max_query_entity_intersection() const
    {
        return 1;
    }
};

class Decomposition_Algorithm;
class Integration_Mesh_Generator
{
  private:
    Model *                                mXTKModel;
    moris::ge::Geometry_Engine *           mGeometryEngine;
    moris::Matrix< moris::IndexMat >       mActiveGeometries;
    moris::Cell< enum Subdivision_Method > mSubdivisionMethods;

    bool mOutputCutIgMesh = true;

  public:
    Integration_Mesh_Generator();
    Integration_Mesh_Generator(
        xtk::Model *                     aXTKModelPtr,
        Cell< enum Subdivision_Method >  aMethods,
        moris::Matrix< moris::IndexMat > aActiveGeometries );
    ~Integration_Mesh_Generator();

    std::shared_ptr< Cut_Integration_Mesh >
    perform();

    moris::Matrix< moris::IndexMat > const *
    get_active_geometries();

    moris::ge::Geometry_Engine *
    get_geom_engine();

    uint
    get_spatial_dim();

    uint
    get_ig_mesh_order();

    enum Subdivision_Method
    determine_order_elevation_template();

    void 
    set_cut_IG_mesh_output( bool aOutputCutIgMesh )
    {
        mOutputCutIgMesh = aOutputCutIgMesh;
    }

    bool
    get_cut_IG_mesh_output()
    {
        return mOutputCutIgMesh;
    }

    /**
     * @brief checks whether all intersected background cells are on the same level. The resultant bool is populated
     * in the cut integration mesh mChildMeshSameLevel. ultimately This triggers a remeshing of HMR in the background mesh
     * 
     * @param aMeshGenerationData 
     * @param aBackgroundMesh 
     */

    void
    check_intersected_background_cell_levels(
        Integration_Mesh_Generation_Data &aMeshGenerationData,
        Cut_Integration_Mesh*             aCutIntegrationMesh,
        moris::mtk::Mesh*                 aBackgroundMesh );

    bool
    determine_intersected_background_cells(
        Integration_Mesh_Generation_Data &aMeshGenerationData,
        Cut_Integration_Mesh *            aCutIntegrationMesh,
        moris::mtk::Mesh *                aBackgroundMesh );


    void
    deduce_interfaces(
        Cut_Integration_Mesh *                      aCutIntegrationMesh,
        std::shared_ptr< Facet_Based_Connectivity > aFacetConnectivity,
        moris::Cell< moris_index > &                aInterfaces );

    void
    construct_bulk_phase_to_bulk_phase_interface(
        Cut_Integration_Mesh *                                               aCutIntegrationMesh,
        moris::Cell< moris_index > &                                         aInterfaces,
        moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > > &aInterfaceBulkPhaseToBulk );

    void
    construct_bulk_phase_to_bulk_phase_dbl_side_interface(
        Cut_Integration_Mesh *                                                      aCutIntegrationMesh,
        moris::Cell< moris_index > &                                                aInterfaces,
        moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > &aDblSideInterfaceBulkPhaseToBulk );

    std::string
    get_interface_side_set_name(
        moris_index aGeomIndex,
        moris_index aBulkPhaseIndex0,
        moris_index aBulkPhaseIndex1 );

    void
    construct_interface_sets(
        Cut_Integration_Mesh *                                               aCutIntegrationMesh,
        moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > > &aInterfaceBulkPhaseToBulk );

    void
    construct_bulk_phase_cell_groups(
        Cut_Integration_Mesh *                           aCutIntegrationMesh,
        moris::Cell< std::shared_ptr< IG_Cell_Group > > &aBulkPhaseCellGroups );

    void
    construct_bulk_phase_blocks(
        Cut_Integration_Mesh *                           aCutIntegrationMesh,
        moris::Cell< std::shared_ptr< IG_Cell_Group > > &aBulkPhaseCellGroups );

    bool
    allocate_child_meshes( Integration_Mesh_Generation_Data &aMeshGenerationData,
        Cut_Integration_Mesh *                               aCutIntegrationMesh,
        moris::mtk::Mesh *                                   aBackgroundMesh );

    void
    collect_vertex_groups_for_background_cells(
        Integration_Mesh_Generation_Data *                 aMeshGenerationData,
        Cut_Integration_Mesh *                             aCutIntegrationMesh,
        moris::Cell< moris::mtk::Cell * > *                aBackgroundCells,
        moris::Cell< std::shared_ptr< IG_Vertex_Group > > *aVertexGroups );

    /** 
     * this function checks if first vertex in first map (first input) 
     * exists in last two inputs and returns MORIS_MAX_INDEX if not
     */
    moris_index
    edge_exists(
        moris::Cell< moris::mtk::Vertex * > &               aVerticesOnEdge,
        std::unordered_map< moris_index, moris_index > &    aLocaLVertexMap,
        moris::Cell< moris::Cell< uint > > &                aVertexToEdge,
        moris::Cell< moris::Cell< moris::mtk::Vertex * > > &aFullEdgeVertices );

    moris_index
    facet_exists(
        moris::Cell< moris::mtk::Vertex * > &               aVerticesOnFacet,
        std::unordered_map< moris_index, moris_index > &    aLocaLVertexMap,
        moris::Cell< moris::Cell< uint > > &                aVertexToFacet,
        moris::Cell< moris::Cell< moris::mtk::Vertex * > > &aFullFacetVertices );

    void
    construct_ig_cell_groups_facet_conn(
        moris::Cell< std::shared_ptr< IG_Cell_Group > > &    aActiveIgCellGroups,
        Cell< std::shared_ptr< Facet_Based_Connectivity > > &aGroupFaceConn );

    void
    generate_cell_neighborhood(
        moris::Cell< moris::mtk::Cell * > &               aCells,
        std::shared_ptr< Facet_Based_Connectivity >       aFaceConnectivity,
        std::shared_ptr< Cell_Neighborhood_Connectivity > aNeighborhood );

    void
    create_facet_from_element_to_node(
        moris::Cell< moris::mtk::Cell * > &         aCells,
        std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity );

    void
    create_edges_from_element_to_node(
        moris::Cell< moris::mtk::Cell * >          aCells,
        std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity );

    void
    select_background_cell_for_edge(
        std::shared_ptr< Edge_Based_Connectivity > aEdgeBasedConnectivity,
        Cut_Integration_Mesh *                     aCutIntegrationMesh,
        moris::Cell< moris::mtk::Cell * > &        aBackgroundCellForEdge );

    void
    select_background_cell_for_facet(
        std::shared_ptr< Facet_Based_Connectivity > aEdgeBasedConnectivity,
        Cut_Integration_Mesh *                      aCutIntegrationMesh,
        moris::Cell< moris::mtk::Cell * > &         aBackgroundCellForEdge );

    void
    deduce_facet_ancestry(
        Cut_Integration_Mesh *                      aCutIntegrationMesh,
        moris::mtk::Mesh *                          aBackgroundMesh,
        std::shared_ptr< Facet_Based_Connectivity > aIgCellGroupFacetConnectivity,
        moris::Cell< moris::mtk::Cell * > const &   aParentCellForDeduction,
        std::shared_ptr< Facet_Based_Ancestry >     aIgFacetAncestry );

    void
    compute_bg_facet_to_child_facet_connectivity(
        Cut_Integration_Mesh *                                               aCutIntegrationMesh,
        moris::mtk::Mesh *                                                   aBackgroundMesh,
        std::shared_ptr< Facet_Based_Connectivity >                          aIgCellGroupFacetConnectivity,
        std::shared_ptr< Facet_Based_Ancestry >                              aIgFacetAncestry,
        moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > &aBgFacetToIgFacet );


    void
    deduce_edge_ancestry(
        Cut_Integration_Mesh *                     aCutIntegrationMesh,
        moris::mtk::Mesh *                         aBackgroundMesh,
        std::shared_ptr< Edge_Based_Connectivity > aIgCellGroupEdgeConnectivity,
        moris::Cell< moris::mtk::Cell * > const &  aParentCellForDeduction,
        std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry );

    void
    commit_new_ig_vertices_to_cut_mesh(
        Integration_Mesh_Generation_Data *aMeshGenerationData,
        Decomposition_Data *              aDecompositionData,
        Cut_Integration_Mesh *            aCutIntegrationMesh,
        moris::mtk::Mesh *                aBackgroundMesh,
        Decomposition_Algorithm *         aDecompositionAlgorithm );

    void
    link_new_vertices_to_geometry_engine(
        Decomposition_Data *     aDecompositionData,
        Decomposition_Algorithm *aDecompAlg );

    void
    commit_new_ig_cells_to_cut_mesh(
        Integration_Mesh_Generation_Data *aMeshGenerationData,
        Decomposition_Data *              aDecompositionData,
        Cut_Integration_Mesh *            aCutIntegrationMesh,
        moris::mtk::Mesh *                aBackgroundMesh,
        Decomposition_Algorithm *         aDecompositionAlgorithm );

    void
    identify_and_construct_subphases(
        Integration_Mesh_Generation_Data *                aMeshGenerationData,
        Cut_Integration_Mesh *                            aCutIntegrationMesh,
        moris::mtk::Mesh *                                aBackgroundMesh,
        std::shared_ptr< Cell_Neighborhood_Connectivity > aCutNeighborhood );


    void
    construct_subphase_neighborhood(
        Cut_Integration_Mesh *                                               aCutIntegrationMesh,
        moris::mtk::Mesh *                                                   aBackgroundMesh,
        std::shared_ptr< Facet_Based_Connectivity >                          aFacetConnectivity,
        moris::Cell< std::shared_ptr< moris::Cell< moris::moris_index > > > *aBgFacetToChildFacet,
        std::shared_ptr< Subphase_Neighborhood_Connectivity >                aSubphaseNeighborhood );


    void
    collect_subphases_attached_to_facet_on_cell(
        Cut_Integration_Mesh *                               aCutIntegrationMesh,
        moris::mtk::Cell const *                             aBGCell,
        moris::moris_index                                   aFacetOrdinal,
        moris::moris_index                                   aSharedFacetIndex,
        std::shared_ptr< Facet_Based_Connectivity >          aFacetConnectivity,
        std::shared_ptr< moris::Cell< moris::moris_index > > aBgFacetToChildrenFacets,
        Cell< moris::moris_index > &                         aSubphaseIndices,
        Cell< moris::moris_index > &                         aRepresentativeIgCells,
        Cell< moris::moris_index > &                         aRepresentativeIgCellsOrdinal );

    void
    collect_ig_cells_and_side_ords_on_bg_facet(
        Cut_Integration_Mesh *             aCutIntegrationMesh,
        moris::moris_index                 aBackgroundFacetIndex,
        moris::Cell< moris::mtk::Cell * > &aIgCell,
        moris::Cell< moris_index > &       aIgCellSideOrds );

    void
    assign_subphase_glob_ids(
        Cut_Integration_Mesh *aCutIntegrationMesh,
        moris::mtk::Mesh *    aBackgroundMesh );


    void
    prepare_subphase_identifier_requests(
        Cut_Integration_Mesh *                    aCutIntegrationMesh,
        Cell< Cell< moris_id > > &                aNotOwnedSubphasesToProcs,
        Cell< moris::Matrix< IdMat > > &          aParentCellIds,
        Cell< moris::Matrix< IdMat > > &          aChildCellIds,
        Cell< moris::Matrix< IdMat > > &          aNumChildCellsInSubphase,
        Cell< uint > &                            aProcRanks,
        std::unordered_map< moris_id, moris_id > &aProcRankToDataIndex );

    void
    prepare_subphase_id_answers(
        Cut_Integration_Mesh *      aCutIntegrationMesh,
        moris::mtk::Mesh *          aBackgroundMesh,
        Cell< Matrix< IndexMat > > &aReceivedParentCellIds,
        Cell< Matrix< IndexMat > > &aFirstChildCellIds,
        Cell< Matrix< IndexMat > > &aReceivedNumChildCellsInSubphase,
        Cell< Matrix< IndexMat > > &aSubphaseIds );

    void
    handle_received_subphase_id_request_answers(
        Cut_Integration_Mesh *             aCutIntegrationMesh,
        Cell< Cell< moris_index > > const &aSubphaseIndices,
        Cell< Matrix< IndexMat > > const & aReceivedSubphaseIds );


    moris::Matrix< moris::IndexMat >
    flood_fill_ig_cell_group(
        Cut_Integration_Mesh *                            aCutIntegrationMesh,
        std::shared_ptr< Cell_Neighborhood_Connectivity > aCutNeighborhood,
        std::shared_ptr< IG_Cell_Group >                  aIgCellGroup,
        moris_index &                                     aMaxValueAssigned );


    void
    extract_cells_from_cell_groups(
        moris::Cell< std::shared_ptr< IG_Cell_Group > > const &aCellGroups,
        moris::Cell< moris::mtk::Cell * > &                    aCellsInGroups );

    void
    compute_ig_cell_bulk_phase(
        Cut_Integration_Mesh *aCutIntegrationMesh );

    moris_index
    deduce_ig_cell_bulk_phase_index( moris::mtk::Cell const *aCell );

    moris_index
    get_max_index( moris::Cell< moris::mtk::Cell * > &aCells );


    void
    assign_node_requests_identifiers(
        Decomposition_Data &  aDecompData,
        Cut_Integration_Mesh *aCutIntegrationMesh,
        moris::mtk::Mesh *    aBackgroundMesh,
        moris::moris_index    aMPITag );

    void
    sort_new_node_requests_by_owned_and_not_owned(
        Decomposition_Data &                      tDecompData,
        Cut_Integration_Mesh *                    aCutIntegrationMesh,
        moris::mtk::Mesh *                        aBackgroundMesh,
        Cell< uint > &                            aOwnedRequests,
        Cell< Cell< uint > > &                    aNotOwnedRequests,
        Cell< uint > &                            aProcRanks,
        std::unordered_map< moris_id, moris_id > &aProcRankToIndexInData );

    void
    assign_owned_request_id(
        Decomposition_Data &aDecompData,
        Cell< uint > const &aOwnedRequest,
        moris::moris_id &   aNodeId );

    void
    setup_outward_requests(
        Decomposition_Data const &                aDecompData,
        moris::mtk::Mesh *                        aBackgroundMesh,
        Cell< Cell< uint > > const &              aNotOwnedRequests,
        Cell< uint > const &                      aProcRanks,
        std::unordered_map< moris_id, moris_id > &aProcRankToIndexInData,
        Cell< Matrix< IndexMat > > &              aOutwardRequests );

    void
    prepare_request_answers(
        Decomposition_Data &              aDecompData,
        moris::mtk::Mesh *                aBackgroundMesh,
        Cell< Matrix< IndexMat > > const &aReceiveData,
        Cell< Matrix< IndexMat > > &      aRequestAnswers );
    void
    handle_received_request_answers(
        Decomposition_Data &              aDecompData,
        moris::mtk::Mesh *                aBackgroundMesh,
        Cell< Matrix< IndexMat > > const &aRequests,
        Cell< Matrix< IndexMat > > const &aRequestAnswers,
        moris::moris_id &                 aNodeId );


    moris::uint
    verbosity_level();
    
    /**
     * @brief Removes subphases and deletes contained ig cells
     * 
     * @param aSubphasesToRemove subphase indexes to delete
     */

    void
    remove_subphases_from_cut_mesh(moris::Cell<moris_index> const & aSubphasesToRemove);
};

}// namespace xtk

#endif