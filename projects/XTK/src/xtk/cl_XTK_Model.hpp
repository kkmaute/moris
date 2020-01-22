/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jul 2, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_MODEL_HPP_
#define SRC_XTK_CL_XTK_MODEL_HPP_

// Standard Include
#include <limits>
#include <unordered_map>
#include <mpi.h>
#include <ctime>

// XTKL: Mesh Includes
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Sets_Info.hpp"
#include "cl_MTK_Fields_Info.hpp"
#include "cl_Mesh_Factory.hpp"

// XTKL: Geometry Engine Includes
#include "cl_MGE_Geometry_Engine.hpp"
#include "cl_MGE_Geometry_Object.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Tools includes
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"

#include "cl_XTK_Output_Options.hpp"

// Linalg Includes
#include "cl_Matrix.hpp"
#include "typedefs.hpp"
#include "fn_print.hpp"

#include "cl_Communication_Tools.hpp"

// Topology
//TODO: MOVE THESE WITH CUTTING METHODS SOMEWHERE ELSE
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Quad_4_Topology.hpp"
#include "cl_XTK_Hexahedron_8_Topology.hpp"

#include "fn_tet_volume.hpp"

// general geometry engine class
#include "../projects/GEN/GEN_MAIN/src/geomeng/cl_GEN_Geometry_Engine.hpp"  // FIXME

namespace xtk
{
class Enrichment;
class Enrichment_Parameters;
class Ghost_Stabilization;
class Enriched_Interpolation_Mesh;
class Enriched_Integration_Mesh;
}


using namespace moris;

namespace xtk
{
class Model
{
public:
    // Public member functions/data
    bool mVerbose = false;

    // Forward declare the maximum value of moris::size_t and moris::real
//    moris::real REAL_MAX          = MORIS_REAL_MAX;
//    moris::moris_index INTEGER_MAX = MORIS_INDEX_MAX;

    // friend classes
    friend class Enrichment;
    friend class Enriched_Interpolation_Mesh;
    friend class Enriched_Integration_Mesh;

    //--------------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------------
    Model(){};

    /*
     * using the general geometry engine
     */
    Model(uint aModelDimension,
          moris::mtk::Interpolation_Mesh* aMeshData,
          moris::ge::GEN_Geometry_Engine & aGeometryEngine,
          bool aLinkGeometryOnConstruction = true);

    /**
     * Primary constructor (this constructor is used for all cases except when testing something)
     */
//    Model(uint aModelDimension,
//          moris::mtk::Interpolation_Mesh* aMeshData,
//          Geometry_Engine & aGeometryEngine,
//          bool aLinkGeometryOnConstruction = true);

    // Indicates the background mesh and the geometry are the same thing
    bool mSameMesh;

    ~Model();

    /*!
     * Initializes link between background mesh and the geometry
     */
    void
    link_background_mesh_to_geometry_objects();

    //--------------------------------------------------------------------------------
    // Operations
    //--------------------------------------------------------------------------------
    /*!
     * Decomposes a mesh to conform to a geometry
     * @param aMethods - specify which type of subdivision method to use (this could be changed to command line parsing or XML reading)
     * @param aSetPhase - tell it to set phase information
     */
    void decompose(Cell<enum Subdivision_Method> aMethods,
                   bool                          aSetPhase  = true);

    /*!
    * Uses sub-phase information within a child mesh to construct one interpolation element for each sub-phase cluster
    */
    void
    unzip_child_mesh();

    /*!
     * Unzipps the interface (must be called after decompose but prior
     * to enrichment)
     */
    void
    unzip_interface();


    /*!
     * Compute sensitivity. Must be called after decompose.
     */
    void
    compute_sensitivity();

    /*!
     * Perform the enrichment. The interp ordinal is needed for B-spline interpolation where
     * more than one b-spline of the same order can exist
     */
    void
    perform_basis_enrichment(enum EntityRank  aBasisRank,
                             moris::uint      aInterpIndex = 0);

    /*!
     * returns the basis enrichment class constructed from call to perform basis enrichment
     */
    Enrichment const &
    get_basis_enrichment();

    // ----------------------------------------------------------------------------------

    Enriched_Interpolation_Mesh &
    get_enriched_interp_mesh(moris::moris_index aIndex = 0);

    // ----------------------------------------------------------------------------------

    Enriched_Integration_Mesh &
    get_enriched_integ_mesh(moris::moris_index  aIndex = 0);

    // ----------------------------------------------------------------------------------

    /*!
     * Constructs the face oriented ghost penalization
     */
    void
    construct_face_oriented_ghost_penalization_cells();

    // ----------------------------------------------------------------------------------

    /*!
     * Convert Tet4 elements to Tet10 Elements
     */
    void
    convert_mesh_tet4_to_tet10();

    //--------------------------------------------------------------------------------
    // Member data access functions
    //--------------------------------------------------------------------------------

    Cut_Mesh &              get_cut_mesh()              { return mCutMesh; }
    Cut_Mesh const &        get_cut_mesh() const        { return mCutMesh; }
    Background_Mesh &       get_background_mesh()       { return mBackgroundMesh; }
    Background_Mesh const & get_background_mesh() const { return mBackgroundMesh; }
//    Geometry_Engine &       get_geom_engine()           { return mGeometryEngine; }
    moris::ge::GEN_Geometry_Engine &       get_geom_engine()           { return mGeometryEngine; }

    // ----------------------------------------------------------------------------------
    // Outputting functions
    //-----------------------------------------------------------------------------------
    /*!
     * Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
     */
    moris::mtk::Integration_Mesh* get_output_mesh(Output_Options const & aOutputOptions = Output_Options());


    /*!
     * returns the XTK model as an mtk mesh
     */
    moris::mtk::Mesh* get_xtk_as_mtk();

    /*!
     * Extracts the surface mesh with respect to the provided geometry index.
     * XTK needs to be provided the side sets names
     */
    void
    extract_surface_mesh_to_obj(std::string                      aOutputFile,
                                size_t                           aPhaseIndex,
                                moris::Cell<std::string> const & aBoundingSideSets);

    //--------------------------------------------------------------------------------
    // Cell Neighborhood creation and access
    //--------------------------------------------------------------------------------
    void
    construct_neighborhood();

    void
    construct_subphase_neighborhood();

    void
    construct_cut_mesh_simple_neighborhood();

    void
    construct_cut_mesh_to_uncut_mesh_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace);

    void
    construct_complex_neighborhood();

    void
    construct_uncut_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace);

    void
    delete_neighborhood(){ mElementToElement.resize(0); };

    //--------------------------------------------------------------------------------
    // Data access functions
    //--------------------------------------------------------------------------------

    /*
     * get spatial dimension of model
     */
    moris::uint
    get_spatial_dim() const;

    /*!
     * returns the number of elements in the entire model
     * includes all child elements and all background elements (combination)
     */
    moris::uint
    get_num_elements_total();

    /*!
     * returns the number of elements in unzipped mesh (total - num child meshes)
     */
    moris::uint
    get_num_elements_unzipped();

    /*!
     * Returns the element to element connecivity
     */
    moris::Cell<moris::Cell<moris::mtk::Cell*>> const &
    get_element_to_element(){ return mElementToElement; };

    moris_index
    get_cell_xtk_index(moris_id aCellId);


    moris::Cell<moris::Cell<moris_index>>  const &
    get_subphase_to_subphase(){ return mSubphaseToSubPhase; };

    bool
    subphase_is_in_child_mesh(moris_index aSubphaseIndex);

    /*
     * Get bulk phase of subphase by index
     */
    uint
    get_subphase_bulk_index(moris_index aSubPhaseIndex);


    moris::Matrix<moris::IndexMat>
    get_element_to_subphase();

    moris_id
    get_subphase_id(moris_id aSubphaseIndex);

    moris_index
    get_subphase_index(moris_id aSubphaseId);


    //--------------------------------------------------------------------------------

    // multi grid stuff

    void   perform_multilevel_enrichment_internal();
    //--------------------------------------------------------------------------------
    // FIXME  only temporary
    void set_HMR_mesh_ptr( std::shared_ptr< moris::mtk::Mesh > aMesh )
    {
        mHMRMesh = aMesh;
    };

    std::shared_ptr< moris::mtk::Mesh > mHMRMesh = nullptr;

    //--------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------
    // Printing Functions
    //--------------------------------------------------------------------------------
    void
    print_cells();

    void
    print_vertex_geometry();

    void
    print_neighborhood();

    void
    print_subphase_neighborhood();

    void
    print_interface_vertices();

protected:
    uint                               mModelDimension;
    Background_Mesh                    mBackgroundMesh;
    Cut_Mesh                           mCutMesh;
//    Geometry_Engine                    mGeometryEngine;

    moris::ge::GEN_Geometry_Engine     mGeometryEngine;

    Enrichment*                        mEnrichment;
    Ghost_Stabilization*               mGhostStabilization;
    Cell<Enriched_Interpolation_Mesh*> mEnrichedInterpMesh;
    Cell<Enriched_Integration_Mesh*>   mEnrichedIntegMesh;

private:

    // XTK Model State Flags
    bool mLinkedBackground  = false; // Model background mesh linked to geometry model
    bool mDecomposed        = false; // Model has been decomposed
    bool mSensitivity       = false; // Model has computed sensitivity
    bool mConvertedToTet10s = false; // Model has been converted from tet4's to tet10s
    bool mEnriched          = false; // Model has been enriched
    bool mUnzipped          = false; // Model has been unzipped
    bool mGhost             = false; // Model has setup ghost stabilization

    // cell map
    std::map< moris_id, moris_index> mCellGlbToLocalMap;

    // The midside nodes are stored here currently but this may change
    moris::Matrix< moris::IndexMat > mMidsideElementToNode;

    // element to element neighborhood
    moris::Cell<moris::Cell<moris::mtk::Cell*>> mElementToElement;
    moris::Cell<moris::Cell<moris_index>> mSubphaseToSubPhase;

    // local to global subphase map
    std::unordered_map<moris::moris_id,moris::moris_index> mGlobalToLocalSubphaseMap;


    // Private Functions
private:
    // Internal Decomposition Functions------------------------------------------------------
    /*!
     * formulates node requests in the geometry objects. Dependent on the type of decomposition
     * @param[in] aReqType- specifies which template mesh is going to be used later on
     */
    void decompose_internal(enum Subdivision_Method    const & aSubdivisionMethod,
                            moris::uint                        aGeomIndex,
                            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                            bool const & aFirstSubdivision = true,
                            bool const & aSetIds = false);


    void
    decompose_internal_reg_sub_hex8_make_requests(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                             	  moris::Matrix< moris::IndexMat > & tNewPairBool,
                                             	  Decomposition_Data 			   & tDecompData);

    void
    decompose_internal_reg_sub_quad4_make_requests(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                             	   moris::Matrix< moris::IndexMat > & tNewPairBool,
						   Decomposition_Data 				& tDecompData);

    void
    decompose_internal_set_new_nodes_in_child_mesh_reg_sub(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                           moris::Matrix< moris::IndexMat > & tNewPairBool,
                                                           moris::real                        tNumParamCoords,
                                                           Decomposition_Data &               tDecompData);

    void
    decompose_internal_set_new_nodes_in_child_mesh_nh(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                      Decomposition_Data &               tDecompData);
    void
    create_new_node_association_with_geometry(Decomposition_Data & tDecompData);

    /*
     * Parallel assignment of node request identifiers
     */
    void
    assign_node_requests_identifiers(Decomposition_Data & aDecompData,
                                     moris::moris_index   aMPITag);

    void
    sort_new_node_requests_by_owned_and_not_owned(Decomposition_Data                    & tDecompData,
                                                  Cell<uint>                            & aOwnedRequests,
                                                  Cell<Cell<uint>>                      & aNotOwnedRequests,
                                                  Cell<uint>                            & aProcRanks,
                                                  std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData);

    void
    assign_owned_request_identifiers(Decomposition_Data & aDecompData,
                                     Cell<uint> const &       aOwnedRequest,
                                     moris::moris_id &        aNodeInd,
                                     moris::moris_id &        aNodeId);


    /*!
     * Sets up the decomposition data request for node identifiers from owning proc
     *
     */
    void
    setup_outward_requests(Decomposition_Data              const & aDecompData,
                           Cell<Cell<uint>>                const & aNotOwnedRequests,
                           Cell<uint>                      const & aProcRanks,
                           std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData,
                           Cell<Matrix<IndexMat>>                & aSentRequests);

    void
    send_outward_requests(moris_index            const & aMPITag,
                          Cell<uint>             const & aProcRanks,
                          Cell<Matrix<IndexMat>> & aOutwardRequests);

//    void
//    assign_node_requests_owned_identifiers_and_setup_send(Decomposition_Data & aDecompData,
//                                                          Cell<uint> const &       aOwnedRequest,
//                                                          Cell<Matrix<IndexMat>> & aSendData,
//                                                          moris::moris_id &        aNodeInd,
//                                                          moris::moris_id &        aNodeId);

    void
    inward_receive_requests(moris_index            const & aMPITag,
                            moris::uint                    aNumRows,
                            Cell<Matrix<IndexMat>> &       aReceivedData,
                            Cell<uint>             &       aProcRanksReceivedFrom);

    void
    prepare_request_answers( Decomposition_Data & aDecompData,
                             Cell<Matrix<IndexMat>> const & aReceiveData,
                             Cell<Matrix<IndexMat>>       & aReceivedRequestAnswers);

    void
    return_request_answers(  moris_index const & aMPITag,
                             Cell<Matrix<IndexMat>> const & aRequestAnswers,
                             Cell<uint>              const & aProcRanks);

    void
    inward_receive_request_answers(moris_index            const & aMPITag,
                                   moris::uint            const & aNumRows,
                                   Cell<uint>             const & aProcRanks,
                                   Cell<Matrix<IndexMat>> &       aReceivedData);

    void
    handle_received_request_answers(Decomposition_Data & aDecompData,
                                    Cell<Matrix<IndexMat>> const & aRequests,
                                    Cell<Matrix<IndexMat>> const & aRequestAnswers,
                                    moris::moris_id &    aNodeInd,
                                    moris::moris_id &    aNodeId);

    /*
     * Perform all tasks needed to finalize the decomposition process, such that the model is ready for enrichment, conversion to tet10 etc.
     * Tasks performed here:
     *  - Assign all child elements global ids and processor local indices
     *  - Store phase indices of non-intersected parent elements
     *
     */
    void
    finalize_decomp_in_xtk_mesh(bool aSetPhase);

    /*!
     * assign child element identifiers
     */
    void
    assign_child_element_identifiers();

    void
    prepare_child_element_identifier_requests(Cell<Cell<moris_id>>       & aNotOwnedChildMeshesToProcs,
                                              Cell<moris::Matrix<IdMat>> & aOwnedParentCellId,
                                              Cell<moris::Matrix<IdMat>> & aNumOwnedCellIdsOffsets,
                                              Cell<uint >          & aProcRanks,
                                              std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex);

    void
    prepare_child_cell_id_answers(Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
                                  Cell<Matrix<IndexMat>> & aReceivedParentCellNumChildren,
                                  Cell<Matrix<IndexMat>> & aChildCellIdOffset);

    void
    handle_received_child_cell_id_request_answers( Cell<Cell<moris_index>> const & aChildMeshesInInNotOwned,
                                                   Cell<Matrix<IndexMat>>  const & aReceivedChildCellIdOffset,
                                                   moris::moris_id               & aCellId);

    /*!
     * Add children elements to local to global map
     */
    void
    add_child_elements_to_local_to_global_map();

    /*!
     * Constructs the mtk cell interface for all child elements created during the
     * decomposition process
     */
    void
    create_child_element_mtk_cells();

    /*!
    * Add the vertex pointer to child meshes
    */
    void
    add_vertices_to_child_meshes();

    /*!
     * setup cell id to index map
     */
    void
    setup_cell_glb_to_local_map();

    /*
    *
    * Identifies local sub-phase clusters within a single child mesh.
    * The sub-phase data (result of floodfill) is stored as a member
    * variable in each child mesh as sub-phase bins
    */
    void
    identify_local_subphase_clusters_in_child_meshes();

    void
    prepare_subphase_identifier_requests(Cell<Cell<moris_id>>       & aNotOwnedSubphasesToProcs,
                                         Cell<Cell<moris_id>>       & aSubphaseCMIndices,
                                         Cell<moris::Matrix<IdMat>> & aParentCellIds,
                                         Cell<moris::Matrix<IdMat>> & aChildCellIds,
                                         Cell<uint>                 & aProcRanks,
                                         std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex);

    void
    prepare_subphase_id_answers(Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
                                Cell<Matrix<IndexMat>> & aFirstChildCellIds,
                                Cell<Matrix<IndexMat>> & aSubphaseIds);

    void
    handle_received_subphase_id_request_answers( Cell<Cell<moris_index>>    const & aChildMeshesInNotOwned,
                                                 Cell<Cell<moris_index>>    const & aCMSubphaseIndices,
                                                 Cell<Matrix<IndexMat>>     const & aReceivedSubphaseIds);

    void
    assign_subphase_glob_ids();

    void
    setup_glob_to_loc_subphase_map();


    void
    unzip_child_mesh_internal();

    /**
     * Take the interface faces and create collapsed prisms
     */
    void
    unzip_interface_internal();

    void
    unzip_interface_internal_assign_node_identifiers(moris::uint aNumNodes,
                                                     moris::Matrix<moris::IdMat> & aUnzippedNodeIndices,
                                                     moris::Matrix<moris::IdMat> & aUnzippedNodeIds);

    void
    unzip_interface_internal_modify_child_mesh(moris::uint                         aGeometryIndex,
                                               moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                                               moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                                               moris::Matrix<moris::IdMat> const & aUnzippedNodeIds);

    moris::Matrix< moris::IndexMat >
    unzip_interface_internal_assign_which_element_uses_unzipped_nodes( moris::moris_index aGeometryIndex,
                                                                       moris::Matrix< moris::IndexMat > const & aInterfaceElementPairs );



    moris::Cell<moris::Cell< moris::moris_index >>
    unzip_interface_internal_collect_child_mesh_to_interface_node(moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                                                                  moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                                                                  moris::Matrix<moris::IdMat> const & aUnzippedNodeIds);

    void
    unzip_interface_construct_interface_elements(moris::uint aGeometryIndex,
                                                 moris::Matrix< moris::IndexMat > const & aElementPairs,
                                                 moris::Matrix< moris::IndexMat > const & aSideOrdinalPairs);
    void
    unzip_interface_assign_element_identifiers();




    // Sensitivity computation functions -----------------------------------------------
    void
    compute_interface_sensitivity_internal();


    void
    extract_interface_sensitivity_sparse(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput,
                                         moris::Cell<moris::Matrix<DDRMat>> & adxdpData,
                                         moris::Cell<std::string>           & adxdpNames,
                                         moris::Cell<moris::Matrix<DDRMat>> & aDesVars,
                                         moris::Cell<std::string>           & aDesVarsName,
                                         moris::Matrix<moris::DDRMat>       & aNumDesVars,
                                         std::string                        & aNumDesVarsName) const;

    void
    extract_interface_sensitivity_dense(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput,
                                        moris::Cell<moris::Matrix<DDRMat>>   & adxdpData,
                                        moris::Cell<std::string>             & adxdpNames) const;

    // Enrichment computation functions -----------------------------------------------

    void
    perform_basis_enrichment_internal(enum EntityRank  aBasisRank,
                                      moris::uint      aInterpOrd);

    /*!
     * Links the vertex enrichment to the mtk implementation of the vertex
     */
    void
    link_vertex_enrichment_to_vertex_interpolation();


    // internal ghost functions -------------------------------------------------------


    /*
     * Constructs child mesh groups (owned, owned shared, not owned shared)
     */
    void
    sort_children_meshes_into_groups();

    /*
     * For nodes that are created during the decomposition process, tell
     * the XTK mesh about where they live in child meshes.
     */
    void
    associate_nodes_created_during_decomp_to_child_meshes();

    /*
     * Set element phase index
     */
    void
    set_element_phases();

    /*
     * Tells the XTK mesh about where it's children live in the cut mesh
     */
    void
    set_downward_inheritance();


    /*
     * This algorithm sets up the active child mesh indices and registers new pairs in the downward inheritance
     */

    void
    run_first_cut_routine(enum TemplateType const &          aTemplateType,
                          moris::uint                        aGeomIndex,
                          moris::size_t const &              aNumNodesPerElement,
                          moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                          moris::Matrix< moris::IndexMat > & aNewPairBool);

    /*!
     * Constructs the output mesh using provided Output_Options
     */
    moris::mtk::Integration_Mesh*
    construct_output_mesh( Output_Options const & aOutputOptions );

    /*!
     * Setup interface single sided side sets
     */
    void
    setup_interface_single_side_sets(Output_Options const & aOutputOptions,
                                     Cell<moris::Matrix<moris::IdMat>> & aCellIdsAndSideOrds,
                                     Cell<std::string> &                 aInterfaceSetNames);

    /*!
     * extract surface internal function
     */
    void
    extract_surface_mesh_to_obj_internal(std::string                      aOutputFile,
                                         size_t                           aPhaseIndex,
                                         moris::Cell<std::string> const & aBoundingSideSets);

    moris::Cell< moris::Matrix < moris::DDRMat > >
    assemble_geometry_data_as_mesh_field(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput);

    /*!
     * Takes the whole node local to global map and removes the nodes
     * which are not part of the phase being output
     */
    moris::Matrix<moris::IndexMat>
    get_node_map_restricted_to_output_phases(Output_Options const &           aOutputOptions,
                                             moris::Matrix<moris::IndexMat> & aOutputtedNodeInds);

    moris::Cell<std::string>
    assign_geometry_data_names();


    moris::Cell < enum moris::EntityRank >
    assign_geometry_data_field_ranks();



    /*!
     * Sets up background node sets for mesh output. Propogates the node set from
     * the background mesh to the output mesh
     */
    moris::Cell<moris::mtk::MtkNodeSetInfo>
    propogate_background_node_sets(moris::Cell<moris::Matrix<IndexMat>> & aNodeSetData,
                                   Output_Options const & aOutputOptions);

    /*!
     * Sets up background side sets for mesh output. Propogates the side set from
     * the background mesh to the output mesh
     */
    moris::Cell<moris::mtk::MtkSideSetInfo>
    propogate_background_side_sets(moris::Cell<moris::Matrix<IndexMat>> & aSideSetData,
                                   Output_Options const & aOutputOptions);

    /*!
     * Add a single side set from the background mesh
     */
    void
    propogate_background_side_set( std::string             const &             aSideSetName,
                                   moris::moris_index                          aNoChildIndex,
                                   moris::moris_index                          aChildIndex,
                                   moris::Cell<moris::Matrix<IndexMat>>      & aElementIdsAndSideOrd,
                                   moris::Cell<moris::mtk::MtkSideSetInfo>   & aSideSetData,
                                   Output_Options          const             & aOutputOptions,
                                   bool                                        aOutputIndices);

    /*!
     * This function checks for the side sets which appear in a mesh that comes
     * from a SEACAS generated string. For example:
     * "generated:1x1x1|sideset:xXyYzZ"
     */
    moris::Cell<std::string>
    check_for_and_remove_internal_seacas_side_sets(moris::Cell<std::string> & aSideSetNames);

    /*!
     * Combine interface and non-interface blocks
     */
    Cell<moris::Matrix<moris::IdMat>>
    combine_interface_and_non_interface_blocks(Cell<moris::Matrix<moris::IdMat>> & tChildElementsByPhase,
                                               Cell<moris::Matrix<moris::IdMat>> & tNoChildElementsByPhase);

    /*!
     * Pack the ghost stabilization cells as a side set
     */
    Cell<Matrix<IdMat>>
    pack_ghost_as_side_set();


    // Internal Aura Construction ------------------------------------------------------

    /*
     * Returns the number of phases to output given the output options
     */
    uint
    get_num_phases_to_output(Output_Options const & aOutputOptions);

    /*!
     * Setup clustering data
     */
    void
    setup_cell_clusters_for_output(moris::mtk::Cell_Cluster_Input & aCellClusterInput,
                                   Output_Options const & aOutputOptions,
                                   moris::Cell<Matrix<IdMat>> & tCellIds);

    void
    setup_interface_side_cluster(std::string                      aInterfaceSideLabelBase,
                                 moris::mtk::Side_Cluster_Input & aCellClusterInput,
                                 Output_Options const        & aOutputOptions,
                                 moris::Cell<Matrix<IdMat>>  & tCellIdsandSideOrds,
                                 moris::Cell<Matrix<DDRMat>> & tParametricCoordinates);

    bool
    output_node(moris::moris_index aNodeIndex,
                Output_Options const & aOutputOptions);

    /*
     * Prints the method of decomposition, type of background mesh,
     */
    void
    print_decompsition_preamble(Cell<enum Subdivision_Method> aMethods);


    moris::size_t
    determine_element_phase_index(moris::size_t aRowIndex,
                                  moris::Matrix< moris::IndexMat > const & aElementToNodeIndex);


    void
    collect_subphases_attached_to_facet_on_cell(moris::moris_index aCellIndex,
                                                moris::moris_index aFacetIndex,
                                                Cell<moris::moris_index> & aCellSubphaseIndices,
                                                Cell<moris::moris_index> & aCellSubphaseBulkIndices);


};
}

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
