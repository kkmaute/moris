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
#include <mpi.h>
#include <ctime>

// XTKL: Mesh Includes
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Sets_Info.hpp"
#include "cl_MTK_Fields_Info.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_XTK_Impl.hpp"

// XTKL: Geometry Engine Includes
#include "cl_MGE_Geometry_Engine.hpp"
#include "cl_MGE_Geometry_Object.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Tools includes
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Sensitivity.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"

#include "cl_XTK_Output_Options.hpp"
#include "cl_XTK_Active_Process_Manager.hpp"

// Linalg Includes
#include "cl_Matrix.hpp"
#include "fn_print.hpp"

#include "cl_Communication_Tools.hpp"

// Topology
//TODO: MOVE THESE WITH CUTTING METHODS SOMEWHERE ELSE
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Quad_4_Topology.hpp"
#include "cl_XTK_Hexahedron_8_Topology.hpp"

#include "fn_tet_volume.hpp"

namespace xtk
{
class Enrichment;
class Enrichment_Parameters;
class Ghost_Stabilization;
class Enriched_Interpolation_Mesh;
}




namespace xtk
{
class Model
{
public:
    // Public member functions/data
    bool mVerbose = false;

    // Forward declare the maximum value of moris::size_t and moris::real
    moris::real REAL_MAX          = MORIS_REAL_MAX;
    moris::moris_index INTEGER_MAX = MORIS_INDEX_MAX;

    // friend class
    friend class Enrichment;
    friend class Enriched_Interpolation_Mesh;

    //--------------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------------
    Model(){};

    /**
     * Primary constructor (this constructor is used for all cases except when testing something)
     */
    Model(uint aModelDimension,
          moris::mtk::Interpolation_Mesh* aMeshData,
          Geometry_Engine & aGeometryEngine,
          bool aLinkGeometryOnConstruction = true);

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
     * Perform the generalized heaviside enrichment
     */
    void
    perform_basis_enrichment();

    /*!
     * returns the basis enrichment class constructed from call to perform basis enrichment
     */
    Enrichment const &
    get_basis_enrichment();

    // ----------------------------------------------------------------------------------

    Enriched_Interpolation_Mesh const &
    get_enriched_interp_mesh();


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

    /*!
     * Returns the Cut Mesh
     */
    Cut_Mesh &
    get_cut_mesh()
    {
        return mCutMesh;
    }

    // ----------------------------------------------------------------------------------

    /*!
     * Returns the Cut Mesh
     */
    Cut_Mesh const &
    get_cut_mesh() const
    {
        return mCutMesh;
    }

    // ----------------------------------------------------------------------------------

    /*!
     * Returns the Background Mesh
     */
    Background_Mesh &
    get_background_mesh()
    {
        return mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------------

    /*!
     * Returns the Xtk Mesh
     */
    Background_Mesh const &
    get_background_mesh() const
    {
        return mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------------

    /*!
     * Get geometry engine
     */
    Geometry_Engine &
    get_geom_engine()
    {
        return mGeometryEngine;
    }

    // ----------------------------------------------------------------------------------
    // Outputting functions
    //-----------------------------------------------------------------------------------
    /*!
     * Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
     */
    moris::mtk::Integration_Mesh*
    get_output_mesh(Output_Options const & aOutputOptions = Output_Options());


    /*!
     * returns the XTK model as an mtk mesh
     */
    moris::mtk::Mesh*
    get_xtk_as_mtk();

    /*!
     * Extracts the surface mesh with respect to the provided geometry index.
     * XTK needs to be provided the side sets names
     */
    void
    extract_surface_mesh_to_obj(std::string                      aOutputFile,
                                size_t                           aPhaseIndex,
                                moris::Cell<std::string> const & aBoundingSideSets);


    //--------------------------------------------------------------------------------
    // Data access functions
    //--------------------------------------------------------------------------------

    /*
     * get spatial dimension of model
     */
    moris::uint
    get_spatial_dim();

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


    //--------------------------------------------------------------------------------

    // multi grid stuff

    void
    perform_multilevel_enrichment_internal();
    //--------------------------------------------------------------------------------
    // FIXME  only temporary
    void set_HMR_mesh_ptr( std::shared_ptr< moris::mtk::Mesh > aMesh )
    {
        mHMRMesh = aMesh;
    };

    std::shared_ptr< moris::mtk::Mesh > mHMRMesh = nullptr;

    //--------------------------------------------------------------------------------

protected:
    uint                         mModelDimension;
    Background_Mesh              mBackgroundMesh;
    Cut_Mesh                     mCutMesh;
    Geometry_Engine              mGeometryEngine;
    Enrichment*                  mEnrichment;
    Ghost_Stabilization*         mGhostStabilization;
    Enriched_Interpolation_Mesh* mEnrichedInterpMesh;

private:

    // XTK Model State Flags
    bool mLinkedBackground  = false; // Model background mesh linked to geometry model
    bool mDecomposed        = false; // Model has been decomposed
    bool mSensitivity       = false; // Model has computed sensitivity
    bool mConvertedToTet10s = false; // Model has been converted from tet4's to tet10s
    bool mEnriched          = false; // Model has been enriched
    bool mUnzipped          = false; // Model has been unzipped
    bool mGhost             = false; // Model has setup ghost stabilization

    // The midside nodes are stored here currently but this may change
    moris::Matrix< moris::IndexMat > mMidsideElementToNode;


    // Private Functions
private:
    // Internal Decomposition Functions------------------------------------------------------
    /*!
     * formulates node requests in the geometry objects. Dependent on the type of decomposition
     * @param[in] aReqType- specifies which template mesh is going to be used later on
     */
    void decompose_internal(enum Subdivision_Method    const & aSubdivisionMethod,
                            moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                            bool const & aFirstSubdivision = true,
                            bool const & aSetIds = false);


    void
    decompose_internal_reg_sub_make_requests(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                             moris::Matrix< moris::IndexMat > & tNewPairBool,
                                             Decomposition_Data & tDecompData);

    void
    decompose_internal_set_new_nodes_in_child_mesh_reg_sub(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                           moris::Matrix< moris::IndexMat > & tNewPairBool,
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
    assign_node_requests_identifiers(Decomposition_Data & tDecompData,
                                     moris::moris_index   aMPITag);

    void
    sort_new_node_requests_by_owned_and_not_owned(Decomposition_Data & tDecompData,
                                                  Cell<uint>         & aOwnedRequests,
                                                  Cell<Cell<uint>>   & aSharedRequestFrom);

    void
    assign_node_requests_owned_identifiers_and_setup_send(Decomposition_Data & aDecompData,
                                                          Cell<uint> const &       aOwnedRequest,
                                                          Cell<Matrix<IndexMat>> & aSendData,
                                                          moris::moris_id &        aNodeInd,
                                                          moris::moris_id &        aNodeId);

    void
    outward_communicate_node_requests(Cell<Matrix<IndexMat>> & aSendData,
                                      moris_index              aMPITag);

    void
    inward_receive_node_requests(Cell<Matrix<IndexMat>> & aReceiveData,
                                 moris_index              aMPITag);

    void
    set_received_node_ids(Cell<Matrix<IndexMat>> & aReceiveData,
                          Decomposition_Data &     aDecompData,
                          moris::moris_id &        aNodeInd);

    /*!
     * During the regular subdivision at the boundary of a processor mesh a hanging node can appear on a face. This happens if an element connected
     * to a face on processor boundary is intersected but the element across the processor boundary is not. This function handles this case.
     */
    void
    handle_hanging_nodes_in_reg_sub( Decomposition_Data & aDecompData,
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

    /*
    *
    * Identifies local sub-phase clusters within a single child mesh.
    * The sub-phase data (result of floodfill) is stored as a member
    * variable in each child mesh as sub-phase bins
    */
    void
    identify_local_subphase_clusters_in_child_meshes();


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

    // Enrichment computation functions -----------------------------------------------

    void
    perform_basis_enrichment_internal();

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
    set_element_phases()
    {
        // Set element phase indices
         mBackgroundMesh.initialize_element_phase_indices(this->get_num_elements_total());

        moris::size_t tNumElem = mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);

         for(moris::size_t i = 0; i<tNumElem; i++)
         {
             if(mBackgroundMesh.entity_has_children(i,EntityRank::ELEMENT))
             {
                 moris::size_t tChildMeshIndex = mBackgroundMesh.child_mesh_index(i,EntityRank::ELEMENT);

                 Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                 moris::Matrix< moris::IndexMat > tElemToNode = tChildMesh.get_element_to_node();

                 moris::Matrix< moris::IndexMat > const & tElemInds  = tChildMesh.get_element_inds();

                 tChildMesh.initialize_element_phase_mat();

                 moris::size_t tNumElem = tChildMesh.get_num_entities(EntityRank::ELEMENT);

                 for( moris::size_t j = 0; j<tNumElem; j++)
                 {
                     moris::size_t tElemPhaseIndex = determine_element_phase_index(j,tElemToNode);
                     mBackgroundMesh.set_element_phase_index(tElemInds(0,j),tElemPhaseIndex);
                     tChildMesh.set_element_phase_index(j,tElemPhaseIndex);
                 }
             }

             else
             {
                 moris::Matrix< moris::IndexMat > tElementNodes = mBackgroundMesh.get_mesh_data().get_entity_connected_to_entity_loc_inds(i,moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

                 moris::size_t tElemPhaseIndex = determine_element_phase_index(0,tElementNodes);

                 mBackgroundMesh.set_element_phase_index(i,tElemPhaseIndex);
             }


         }
    }

    /*
     * Tells the XTK mesh about where it's children live in the cut mesh
     */
    void set_downward_inheritance()
    {
        moris::size_t tNumChildMesh = mCutMesh.get_num_child_meshes();
        Cell<std::pair<moris::moris_index,moris::moris_index>> tXTKElementToCutMeshPairs(tNumChildMesh);

        for(moris::size_t iMesh = 0; iMesh<tNumChildMesh; iMesh++)
        {
            tXTKElementToCutMeshPairs(iMesh) = std::pair<moris::moris_index,moris::moris_index> (mCutMesh.get_parent_element_index(iMesh),iMesh);
        }

        mBackgroundMesh.register_new_downward_inheritance(tXTKElementToCutMeshPairs);
    }


    /*
     * This algorithm sets up the active child mesh indices and registers new pairs in the downward inheritance
     */

    void  run_first_cut_routine(enum TemplateType const &          aTemplateType,
                                moris::size_t const &              aNumNodesPerElement,
                                moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                moris::Matrix< moris::IndexMat > & aNewPairBool)
    {
        // Note this method is independent of node ids for this reason Background_Mesh is not given the node Ids during this subdivision
        moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

        // Package up node to element connectivity
        moris::moris_index tParentElementIndex = INTEGER_MAX;
        moris::size_t tNumElements = mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);

        moris::Matrix< moris::IndexMat > tElementToNodeConnInd (tNumElements, aNumNodesPerElement);
        moris::Matrix< moris::IndexMat > tElementToNodeVector (1, aNumNodesPerElement);

        for (moris::size_t i = 0; i < tNumElements; i++)
        {
            tElementToNodeVector = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

            tElementToNodeConnInd.set_row(i, tElementToNodeVector);
        }

        // Get the Node Coordinates
        moris::Matrix< moris::DDRMat > tAllNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

        // Intersected elements are flagged via the Geometry_Engine
        Cell<Geometry_Object> tGeoObjects;
        mGeometryEngine.is_intersected(tAllNodeCoords, tElementToNodeConnInd, 0,tGeoObjects);

        // Count number intersected
        moris::size_t tIntersectedCount = tGeoObjects.size();

        // Loop over and determine how many new meshes that need to be registered (Avoids dynamic allocation in the child mesh)
        // Also register active mesh pairs
        Cell<std::pair<moris::moris_index,moris::moris_index>> tNewChildElementPair;
        aNewPairBool = moris::Matrix< moris::IndexMat >(1,tIntersectedCount,0);
        tNewChildElementPair.reserve(tIntersectedCount);

        moris::size_t tNumNewChildMeshes = 0;
        moris::moris_index tNewIndex = 0;
        aActiveChildMeshIndices.resize(1,tIntersectedCount);
        for (moris::size_t j = 0; j < tIntersectedCount; j++)
        {
            tParentElementIndex = tGeoObjects(j).get_parent_entity_index();
            if(!mBackgroundMesh.entity_has_children(tParentElementIndex,EntityRank::ELEMENT))
            {
                tNewIndex = tNumNewChildMeshes+mCutMesh.get_num_child_meshes();
                tNewChildElementPair.push_back( std::pair<moris::moris_index,moris::moris_index>(tParentElementIndex, tNewIndex));
                aActiveChildMeshIndices(0,j) = tNewIndex;
                tNumNewChildMeshes++;
            }

            else
            {
                aActiveChildMeshIndices(0,j) = mBackgroundMesh.child_mesh_index(tParentElementIndex,EntityRank::ELEMENT);
                aNewPairBool(0,j) = 1;
            }
        }


        // Add the downward pair to the mesh for all the newly created element pairs
        mBackgroundMesh.register_new_downward_inheritance(tNewChildElementPair);


        // Allocate space for more simple meshes in XTK mesh
        mCutMesh.inititalize_new_child_meshes(tNumNewChildMeshes, mModelDimension);


        moris::Matrix< moris::IndexMat > tPlaceHolder(1, 1);
        for (moris::size_t j = 0; j < tIntersectedCount; j++)
        {
            if(aNewPairBool(0,j) == 0)
            {
                tParentElementIndex = tGeoObjects(j).get_parent_entity_index();

                // Get information to provide ancestry
                // This could be replaced with a proper topology implementation that knows faces, edges based on parent element nodes
                Matrix< IndexMat > tNodetoElemConnVec = tElementToNodeConnInd.get_row(tParentElementIndex);
                Matrix< IndexMat > tEdgetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, moris::EntityRank::ELEMENT, moris::EntityRank::EDGE);
                Matrix< IndexMat > tFacetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, moris::EntityRank::ELEMENT, moris::EntityRank::FACE);
                Matrix< IndexMat > tElementMat        = {{tParentElementIndex}};

                for(moris::uint i = 0; i < tEdgetoElemConnInd.numel(); i++)
                {
                    moris_index tEdgeIndex = tEdgetoElemConnInd(i);
                    Matrix<IndexMat> tEdgeNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tEdgeIndex,EntityRank::EDGE,EntityRank::NODE);
                }

                // Set parent element, nodes, and entity ancestry
                moris::Matrix< moris::IndexMat > tElemToNodeIndices(tNodetoElemConnVec);

                Cell<moris::Matrix< moris::IndexMat >> tAncestorInformation = {tPlaceHolder, tEdgetoElemConnInd, tFacetoElemConnInd, tElementMat};
                mCutMesh.initialize_new_mesh_from_parent_element(aActiveChildMeshIndices(0,j), aTemplateType, tNodetoElemConnVec, tAncestorInformation);

                // add node ids
                mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,tNodetoElemConnVec);

                // get child mesh
                moris::moris_index tCMIndex = mBackgroundMesh.child_mesh_index(tParentElementIndex,EntityRank::ELEMENT);

                // get child mesh
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);

                // add node ids
                tChildMesh.add_node_ids(tNodetoElemConnVec);

            }
        }
    }

    /*!
     * Constructs the output mesh using provided Output_Options
     */
    moris::mtk::Integration_Mesh*
    construct_output_mesh( Output_Options const & aOutputOptions );

    /*!
     * extract surface internal function
     */
    void
    extract_surface_mesh_to_obj_internal(std::string                      aOutputFile,
                                         size_t                           aPhaseIndex,
                                         moris::Cell<std::string> const & aBoundingSideSets);

public:
    moris::Cell< moris::Matrix < moris::DDRMat > >
    assemble_geometry_data_as_mesh_field(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput);

    /*!
     * Takes the whole node local to global map and removes the nodes
     * which are not part of the phase being output
     */
    moris::Matrix<moris::IndexMat>
    get_node_map_restricted_to_output_phases(Output_Options const &           aOutputOptions,
                                             moris::Matrix<moris::IndexMat> & aOutputtedNodeInds);
private:

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
    /*!
     * Package up the elements in the aura
     */
    void
    package_aura_block(Output_Options              const & aOutputOptions,
                       Cell<std::string>                 & aAuraChildrenBlockNames,
                       Cell<moris::Matrix<moris::IdMat>> & aAuraChildrenCellIdsByPhase,
                       Cell<std::string>                 & aAuraNoChildrenBlockNames,
                       Cell<moris::Matrix<moris::IdMat>> & aAuraNoChildrenCellIdsByPhase);

    /*
     * setup the block names
     */
    void
    setup_aura_block_names(Output_Options              const & aOutputOptions,
                           Cell<std::string>                 & aAuraChildrenBlockNames,
                           Cell<std::string>                 & aAuraNoChildrenBlockNames);

    /*
     * place cells into aura blocks
     */
    void
    setup_aura_cells_into_blocks(Output_Options              const & aOutputOptions,
                                 Cell<std::string>                 & aAuraChildrenBlockNames,
                                 Cell<moris::Matrix<moris::IdMat>> & aAuraChildrenCellIdsByPhase,
                                 Cell<std::string>                 & aAuraNoChildrenBlockNames,
                                 Cell<moris::Matrix<moris::IdMat>> & aAuraNoChildrenCellIdsByPhase);

    moris_index
    get_aura_block_index(Output_Options const & aOutputOptions,
                         moris_index            aBulkPhase,
                         moris_index            aProcessorRank,
                         moris_index            aNumBulkPhases);




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
    print_decompsition_preamble(Cell<enum Subdivision_Method> aMethods)
    {
        // Only process with rank 0 prints the preamble


        if(moris::par_rank() == 0 && mVerbose)
        {
            std::cout<<"XTK: Specified Decomposition Routines: ";

            for(moris::size_t i = 0 ; i<aMethods.size(); i++)
            {
                std::cout<<"["<<get_enum_str(aMethods(i))<<  "] ";
            }

            std::cout<<std::endl;
        }
    }


    moris::size_t
    determine_element_phase_index(moris::size_t aRowIndex,
                                  moris::Matrix< moris::IndexMat > const & aElementToNodeIndex)
    {
        moris::size_t tNumGeom = mGeometryEngine.get_num_geometries();
        moris::size_t tNumNodesPerElem = aElementToNodeIndex.n_cols();
        moris::Matrix< moris::IndexMat > tNodalPhaseVals(1,tNumGeom,INTEGER_MAX);

        for(moris::size_t i = 0; i<tNumGeom; i++)
        {
            bool tFoundNonInterfaceNode = false;
            for( moris::size_t j = 0; j<tNumNodesPerElem; j++)
            {
                if(!mBackgroundMesh.is_interface_node(aElementToNodeIndex(aRowIndex,j),i))
                {
                    tNodalPhaseVals(0,i) = mGeometryEngine.get_node_phase_index_wrt_a_geometry(aElementToNodeIndex(aRowIndex,j),i);
                    tFoundNonInterfaceNode = true;

                }
            }

            if(!tFoundNonInterfaceNode)
            {
                std::cout<<"Did not find a non-interface node for this element"<<std::endl;
                tNodalPhaseVals(0,i) = 1;
            }
        }


        moris::moris_index tElemPhaseVal = mGeometryEngine.get_elem_phase_index(tNodalPhaseVals);

        return tElemPhaseVal;
    }


};
}

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
