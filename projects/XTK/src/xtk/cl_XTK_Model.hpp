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
#include "cl_MTK_Mesh.hpp"
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


#include "cl_XTK_Request_Handler.hpp"
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
}




namespace xtk
{
class Model
{
public:
    // Public member functions
    bool mVerbose = false;


    // Forward declare the maximum value of moris::size_t and moris::real
    moris::real REAL_MAX          = MORIS_REAL_MAX;
    moris::moris_index INTEGER_MAX = MORIS_INDEX_MAX;

    Model(){};

    /**
     * Primary constructor (this constructor is used for all cases except when testing something)
     */
    Model(uint aModelDimension,
          moris::mtk::Mesh* aMeshData,
          Geometry_Engine & aGeometryEngine,
          bool aLinkGeometryOnConstruction = true);

    // Indicates the background mesh and the geometry are the same thing
    bool mSameMesh;

    ~Model();

    // ----------------------------------------------------------------------------------

    /*!
     * Initializes link between background mesh and the geometry
     */
    void
    link_background_mesh_to_geometry_objects();

    // ----------------------------------------------------------------------------------

    /*!
     * Decomposes a mesh to conform to a geometry
     * @param aMethods - specify which type of subdivision method to use (this could be changed to command line parsing or XML reading)
     * @param aSetPhase - tell it to set phase information
     */
    void decompose(Cell<enum Subdivision_Method> aMethods,
                   bool                          aSetPhase  = true);

    // ----------------------------------------------------------------------------------

    /*!
     * Unzipps the interface (must be called after decompose but prior
     * to enrichment)
     */
    void
    unzip_interface();

    // ----------------------------------------------------------------------------------

    /*!
     * Compute sensitivity. Must be called after decompose.
     */
    void
    compute_sensitivity();

    // ----------------------------------------------------------------------------------

    /*!
     * Perform the generalized heaviside enrichment
     */
    void
    perform_basis_enrichment();

    // ----------------------------------------------------------------------------------

    /*!
     * returns the basis enrichment class constructed from call to perform basis enrichment
     */
    Enrichment const &
    get_basis_enrichment();

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

    // ----------------------------------------------------------------------------------


    // Accessor functions
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
     * Get geomtry engine
     */
    Geometry_Engine &
    get_geom_engine()
    {
        return mGeometryEngine;
    }

    // ----------------------------------------------------------------------------------

    /*!
     * Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
     */
    moris::mtk::Mesh*
    get_output_mesh(Output_Options const & aOutputOptions = Output_Options());


    /*!
     * retuns the XTK model as an mtk mesh
     */
    moris::mtk::Mesh*
    get_xtk_as_mtk();

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

    void
    perform_multilevel_enrichment_internal();
    //--------------------------------------------------------------------------------
    // FIXME  only temporary
    void set_HMR_mesh_ptr( std::shared_ptr< moris::mtk::Mesh > aMesh )
    {
        mHMRMesh = aMesh;
    };

    std::shared_ptr< moris::mtk::Mesh > mHMRMesh = nullptr;


private:
    uint                 mModelDimension;
    Background_Mesh      mBackgroundMesh;
    Cut_Mesh             mCutMesh;
    Geometry_Engine      mGeometryEngine;
    Enrichment*          mEnrichment;
    Ghost_Stabilization* mGhostStabilization;

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
                                     moris::moris_index   aMPITag)
    {
        // asserts
        MORIS_ASSERT(tDecompData.tNewNodeId.size() == tDecompData.tNewNodeIndex.size(),      "Dimension mismatch in assign_node_requests_identifiers");
        MORIS_ASSERT(tDecompData.tNewNodeId.size() == tDecompData.tNewNodeParentRank.size(), "Dimension mismatch in assign_node_requests_identifiers");
        MORIS_ASSERT(tDecompData.tNewNodeId.size() == tDecompData.tNewNodeParentIndex.size(),"Dimension mismatch in assign_node_requests_identifiers");

        // owned requests and shared requests sorted by owning proc
        Cell<uint> tOwnedRequest;
        Cell<Cell<uint>> tSharedRequestFrom;
        sort_new_node_requests_by_owned_and_not_owned(tDecompData,tOwnedRequest,tSharedRequestFrom);

        // allocate ids for nodes I own
        moris::moris_id tNodeId  = mBackgroundMesh.allocate_entity_ids(tDecompData.tNewNodeId.size(), EntityRank::NODE);

        // get first available index
        moris::moris_id tNodeInd = mBackgroundMesh.get_first_available_index(EntityRank::NODE);

        // setup the sending data
        Cell<Matrix<IndexMat>> tSendData;
        assign_node_requests_owned_identifiers_and_setup_send(tDecompData, tOwnedRequest,tSendData,tNodeInd,tNodeId);

        // send data to other processors
        outward_communicate_node_requests(tSendData,aMPITag);
        barrier();

        //recieve
        Cell<Matrix<IndexMat>> tReceivedData(par_size());
        inward_receive_node_requests(tReceivedData,aMPITag);

        // set received node ids
        set_received_node_ids(tReceivedData,tDecompData,tNodeInd);

        barrier();

        // it is possible during the regular subdivision method that a hanging node appears at the processor boundary.
        // this handles that case
        if(tDecompData.mSubdivisionMethod ==Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 && tDecompData.mNumNewNodesWithIds != tDecompData.tNewNodeId.size())
        {
            handle_hanging_nodes_in_reg_sub(tDecompData, tNodeInd, tNodeId);
        }

        // Return the node indexes to the background mesh
        mBackgroundMesh.update_first_available_index(tNodeInd,EntityRank::NODE);

        MORIS_ASSERT(tDecompData.mNumNewNodesWithIds == tDecompData.tNewNodeId.size(),"Not all ids have been set at the end of the assign_node_requests_identifiers function");

    }

    void
    sort_new_node_requests_by_owned_and_not_owned(Decomposition_Data & tDecompData,
                                                  Cell<uint>         & aOwnedRequests,
                                                  Cell<Cell<uint>>   & aSharedRequestFrom)
    {
        aSharedRequestFrom = Cell<Cell<uint>>(par_size());

        // access the mesh data behind the background mesh
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        // number of new nodes
        moris::uint tNumNewNodes = tDecompData.tNewNodeParentIndex.size();

        // Par rank
        moris::moris_index tParRank = par_rank();

        // iterate through each node request and figure out the owner
        for(moris::uint i = 0; i <tNumNewNodes; i++)
        {
            // Parent Rank
            enum EntityRank    tParentRank  = tDecompData.tNewNodeParentRank(i);
            moris::moris_index tParentIndex = tDecompData.tNewNodeParentIndex(i);

            // get the owner processor
            moris::moris_index tOwnerProc = tMeshData.get_entity_owner(tParentIndex,tParentRank);

            // If i own the request keep track of the index
            if(tOwnerProc == tParRank)
            {
                aOwnedRequests.push_back(i);
            }
            else
            {
                aSharedRequestFrom(tOwnerProc).push_back(i);
            }
        }
    }

    void
    assign_node_requests_owned_identifiers_and_setup_send(Decomposition_Data & aDecompData,
                                                          Cell<uint> const &       aOwnedRequest,
                                                          Cell<Matrix<IndexMat>> & aSendData,
                                                          moris::moris_id &        aNodeInd,
                                                          moris::moris_id &        aNodeId)
    {
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();
        moris::moris_index tParSize = par_size();
        moris::moris_index tParRank = par_rank();

        // Allocate data which needs to be sent
        // outer cell proc index
        aSendData = Cell<Matrix<IndexMat>>(tParSize);
        Cell<uint> tSendDataSizes(tParSize,0);

        // allocate send data
        for(moris::uint i =0; i <(uint)tParSize; i++)
        {
            // in a given columns
            //   r0 - parent entity index
            //   r1 - parent entity rank
            //   r2 - Secondary id
            //   r3 - node id
            aSendData(i) = moris::Matrix<IndexMat>(4,aOwnedRequest.size());
        }

        // iterate through node request owned by this proc and assign an Id
        for(moris::uint iOwn = 0; iOwn<aOwnedRequest.size(); iOwn++)
        {
            // location in the requests
            moris::moris_index tRequestIndex = aOwnedRequest(iOwn);

            // set the new node index
            aDecompData.tNewNodeIndex(tRequestIndex) = aNodeInd;
            aNodeInd++;

            // set the new node id
            aDecompData.tNewNodeId(tRequestIndex) = aNodeId;
            aNodeId++;

            // increment number of new nodes with set ids (for assertion purposes)
            aDecompData.mNumNewNodesWithIds++;

            // determine who to tell about this node
            enum EntityRank tParentRank  = aDecompData.tNewNodeParentRank(tRequestIndex);
            moris_index tParentIndex     = aDecompData.tNewNodeParentIndex(tRequestIndex);
            Matrix< IdMat > tSharedProcs;
            tMeshData.get_processors_whom_share_entity(tParentIndex,tParentRank,tSharedProcs);


            // iterate through the shared processors
            for(moris::uint iShare = 0; iShare<tSharedProcs.numel(); iShare++)
            {
                // shared processor id
                moris::moris_id tProcId = tSharedProcs(iShare);

                if(tProcId != tParRank)
                {
                    // Count
                    moris::uint tCount = tSendDataSizes(tProcId);

                    //   r0 - parent entity index
                    aSendData(tProcId)(0,tCount) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tParentIndex,tParentRank);
                    //   r1 - parent entity rank (cast to an index)
                    aSendData(tProcId)(1,tCount) = (moris_index) tParentRank;
                    //   r2 - Secondary id
                    aSendData(tProcId)(2,tCount) = aDecompData.tSecondaryIdentifiers(tRequestIndex);
                    //   r3 - node id
                    aSendData(tProcId)(3,tCount) = aDecompData.tNewNodeId(tRequestIndex);

                    // increment the size
                    tSendDataSizes(tProcId) ++;
                }
            }

        }

        // resize matrices in send data cell
        for(moris::uint  i =0; i<aSendData.size(); i++)
        {
            aSendData(i).resize(4,tSendDataSizes(i));
        }

    }

    void
    outward_communicate_node_requests(Cell<Matrix<IndexMat>> & aSendData,
                                      moris_index              aMPITag)
    {
        for(moris::uint i = 0; i <aSendData.size(); i++)
        {
            if(aSendData(i).numel()!=0)
            {
                nonblocking_send(aSendData(i),aSendData(i).n_rows(),aSendData(i).n_cols(),i,aMPITag);
            }
        }
    }

    void
    inward_receive_node_requests(Cell<Matrix<IndexMat>> & aReceiveData,
                                 moris_index              aMPITag)
    {
        moris::moris_index tNumRow = 4;
        moris::moris_index tParRank = par_rank();

        for(moris::uint i = 0; i<(moris::uint)par_size(); i++)
        {
            if((moris_index)i != tParRank)
            {
                // if there is a sent message from a processor go receive it
                if(sent_message_exists(i,aMPITag))
                {
                    aReceiveData(i).resize(1,1);
                    receive(aReceiveData(i),tNumRow, i,aMPITag);
                }
            }
        }
    }

    void
    set_received_node_ids(Cell<Matrix<IndexMat>> & aReceiveData,
                          Decomposition_Data &     aDecompData,
                          moris::moris_id &        aNodeInd)
    {
        moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

        moris::uint tNumCells = aReceiveData.size();

        for(moris::uint i = 0; i <tNumCells; i++)
        {
            moris::uint tNumNodeIdsReceived = aReceiveData(i).n_cols();

            for(moris::uint iN =0; iN<tNumNodeIdsReceived; iN++)
            {
                moris::moris_index tParentRank      = aReceiveData(i)(1,iN);
                moris::moris_index tParentEntityInd = tMeshData.get_loc_entity_ind_from_entity_glb_id(aReceiveData(i)(0,iN),(enum EntityRank)tParentRank);
                moris::moris_index tSecondaryId     = aReceiveData(i)(2,iN);
                moris::moris_index tNodeId          = aReceiveData(i)(3,iN);
                moris::moris_index tRequestIndex    = MORIS_INDEX_MAX;
                bool               tRequestExists   = false;

                if(aDecompData.mHasSecondaryIdentifier)
                {
                    tRequestExists = aDecompData.request_exists(tParentEntityInd,tSecondaryId,(EntityRank)tParentRank,tRequestIndex);
                }

                else
                {
                    tRequestExists = aDecompData.request_exists(tParentEntityInd,(EntityRank)tParentRank,tRequestIndex);
                }

                if(tRequestExists)
                {
                    // set the new node index
                    aDecompData.tNewNodeIndex(tRequestIndex) = aNodeInd;
                    aNodeInd++;

                    // set the new node id
                    aDecompData.tNewNodeId(tRequestIndex) = tNodeId;

                    aDecompData.mNumNewNodesWithIds++;
                }
#ifdef DEBUG
                else
                {
                    MORIS_ASSERT(aDecompData.mSubdivisionMethod == Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8," The only known case where hanging nodes are allowed to show up is on a face during NC_REGULAR_SUBDIVISION_HEX8");
                }
#endif
            }


        }
    }

    /*!
     * During the regular subdivision at the boundary of a processor mesh a hanging node can appear on a face. This happens if an element connected
     * to a face on processor boundary is intersected but the element across the processor boundary is not. This function handles this case.
     */
    void
    handle_hanging_nodes_in_reg_sub( Decomposition_Data & aDecompData,
                                     moris::moris_id &    aNodeInd,
                                     moris::moris_id &    aNodeId)
    {
        // iterate through node requests and assign node ids to nodes without one
        for(moris::uint i = 0; i < aDecompData.tNewNodeId.size(); i++)
        {
            if(aDecompData.tNewNodeId(i) == MORIS_ID_MAX)
            {
                // TODO: Check that this face is on the boundary of the aura
                MORIS_ASSERT(aDecompData.tNewNodeIndex(i)      == MORIS_INDEX_MAX,  "New Node index assigned when the node id has not been set");
                MORIS_ASSERT(aDecompData.tNewNodeParentRank(i) == EntityRank::FACE, "Hanging node in regular subdivision can only be on a face");

                // assign node index
                aDecompData.tNewNodeIndex(i) = aNodeInd;
                aNodeInd++;

                //assign node id
                // set the new node id
                aDecompData.tNewNodeId(i) = aNodeId;
                aNodeId++;

                // increment tally
                aDecompData.mNumNewNodesWithIds++;
            }
        }
    }

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
     * Constructs the mtk cell interface for all child elements created during the
     * decomposition process
     */
    void
    create_child_element_mtk_cells();

    // Internal interface functions

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
     * Convert the child meshes into tet4's
     */

//    void convert_mesh_tet4_to_tet10_internal()
//    {
//
//        mesh::Mesh_Data<moris::real, moris::size_t, moris::DDRMat, moris::DDSTMat> & tXTKMeshData = mXTKMesh.get_mesh_data();
//
//        // Make sure there is at lease one element in the mesh
//        moris::size_t tNumElems = mXTKMesh.get_num_entities(EntityRank::ELEMENT);
//
//        MORIS_ASSERT(tNumElems>0," There needs to be at lease one element in the background mesh to convert to tet 10s");
//
//        // Check background element type
//        // If I start with a tet mesh, need to convert all children and background elements to tet10s
//        // If I start with a hex mesh only convert children elements to tet10
//        bool tTetStart = true;
//        moris::Matrix< moris::DDSTMat > tElem1Nodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(0,EntityRank::ELEMENT,EntityRank::NODE);
//        moris::size_t tNumNodesPerElem = tElem1Nodes.n_cols();
//        MORIS_ASSERT((tNumNodesPerElem == 4 || tNumNodesPerElem == 8), "Background needs to be tet4s or hex8s");
//
//        if(tNumNodesPerElem == 8)
//        {
//            tTetStart = false;
//        }
//
//
//        // Number of Children Meshes
//        moris::size_t tNumChildrenMesh = mCutMesh.get_num_simple_meshes();
//        moris::size_t tTotalNumChildrenElem = mCutMesh.get_num_entities(EntityRank::ELEMENT);
//        moris::Matrix< moris::DDRMat >tNonInterfaceDxDp(3,1,0.0);
//
//        // Initialize request list for faces and elements
//        // Number of children allowed on a parent mesh entity
//        moris::size_t tEdgeChildren = 160;
//        moris::size_t tFaceChildren = 4*4*10;
//        moris::size_t tElemChildren = 4*24*10;
//
//        // Of all the edges in a regular subdivision template 12 live on parent mesh edges, 24 live on parent mesh faces and 14 live in the parent mesh element
//        // the number of parent edges and faces were done by hand and are hardcoded here.
//        // TODO: Ask XTKMesh how many parents of each rank it has (This would eliminate the dependency on the regular subdivision coming first)
//        moris::size_t tNumParentEdge = 4*12;
//        moris::size_t tNumParentFace = 4*24;
//        moris::size_t tNumParentElem = 4*14;
//
//        moris::size_t tNumEdges = mXTKMesh.get_num_entities(EntityRank::EDGE);
//
//        Request_Handler tEdgeRequests(tNumChildrenMesh * tNumParentEdge + tNumEdges, tEdgeChildren, EntityRank::EDGE,    EntityRank::NODE, tXTKMeshData, mCutMesh);
//        Request_Handler tFaceRequests(tNumChildrenMesh * tNumParentFace,             tEdgeChildren, EntityRank::FACE,    EntityRank::NODE, tXTKMeshData, mCutMesh);
//        Request_Handler tElemRequests(tNumChildrenMesh * tNumParentElem,             tEdgeChildren, EntityRank::ELEMENT, EntityRank::NODE, tXTKMeshData, mCutMesh);
//
//        // Tell cut mesh to initialize auxiliary connectivity
//        mCutMesh.init_intersect_connectivity();
//
//        // Make all the children elements request nodes to convert them to
//
//        //TODO: RESTRUCTURING MOVED ABILITY TO ASK CHILD MESH IF A NODE IS ON THE INTERFACE NEED THIS FUNCTION REENABLED
////        make_tet4_to_tet10_child_mesh_node_requests(tEdgeRequests,tFaceRequests,tElemRequests);
//
//
//
//        // Make requests for all non-intersected elements
//        Cell<Cell<moris::size_t*>> tPendingNodes;
//        if(tTetStart)
//        {
//            // Initialize midside node member variables
//            mMidsideElementToNode(tNumElems,6);
//            tPendingNodes =  make_tet4_to_tet10_unintersected_parent_mesh_node_requests(tEdgeRequests);
//        }
//
//        // Assign unique node ids to all requests
//        bool tCoordinateFlag = false;
//
//        tEdgeRequests.handle_requests(tCoordinateFlag, mSameMesh, false, mGeometryEngines);
//
//        tFaceRequests.handle_requests(tCoordinateFlag, mSameMesh, false, mGeometryEngines);
//
//        tElemRequests.handle_requests(tCoordinateFlag, mSameMesh, false, mGeometryEngines);
//
//        // Tell XTK mesh to grab the pending node indices and that these nodes are on the interface and have sensitivity
//        mCutMesh.retrieve_pending_node_inds();
//
//        // Tell the cut mesh to have all of its children meshes generated tet10s from their tet4s
//        mCutMesh.convert_cut_mesh_to_tet10s();
//
//        // Commit midside node connectivity of unintersected parent elements
//        if(tTetStart)
//        {
//            retrieve_midside_node_inds_unintersected(tPendingNodes);
//        }
//
//        // Set new node ids
//        moris::Matrix< moris::DDSTMat > tNodeIds;
//        for (moris::size_t j = 0; j < mCutMesh.get_num_simple_meshes(); j++)
//        {
//            moris::Matrix< moris::DDSTMat > const & tNodeIndices = mCutMesh.get_node_indices(j);
//            tNodeIds = mesh::Mesh_Helper::get_glb_entity_id_from_entity_loc_index_range(tXTKMeshData, tNodeIndices, EntityRank::NODE);
//
//            mCutMesh.set_node_ids(j, tNodeIds);
//        }
//
//    }


    /*
     * Makes node requests at the midside for non-intersected elements, Note only edge requests are needed here because all edges
     * in the parent mesh cannot be from faces or elements
//     */
//    Cell<Cell<moris::size_t*>>
//    make_tet4_to_tet10_unintersected_parent_mesh_node_requests(Request_Handler & aEdgeRequests)
//    {
//        mesh::Mesh_Data<moris::real, moris::size_t, moris::DDRMat, moris::DDSTMat> & tXTKMeshData = mXTKMesh.get_mesh_data();
//        moris::size_t tNumElem = mXTKMesh.get_num_entities(EntityRank::ELEMENT);
//        moris::size_t tNumEdgePerElem  = 6;
//        moris::size_t tNumNodesPerElem = 10;
//        moris::size_t tNumChildMeshs = mCutMesh.get_num_simple_meshes();
//
//        //
//        moris::Matrix< moris::DDRMat >       tEdgeCoords(2, 3, REAL_MAX);
//        moris::Matrix< moris::DDRMat >       tLocalCoord(1, 1, 0.0);
//        moris::Matrix< moris::DDRMat >       tNodeGlobCoord(1, 3, REAL_MAX);
//        moris::Matrix< moris::DDRMat >       tCoordSwapper(1, 3, REAL_MAX);
//        moris::Matrix< moris::DDRMat >       tNodeCoord(1,3,REAL_MAX);
//        moris::Matrix< moris::DDSTMat >  tEdgeNodes(1,2,INTEGER_MAX);
//        moris::Matrix< moris::DDSTMat >  tElementEdges(1,2,INTEGER_MAX);
//        Edge_Topology tEdgeTopology;
//
//        // cell of pending indices
//        Cell<Cell<moris::size_t*>> tNodeInds(tNumElem,6);
//
//        for(moris::size_t i = 0; i<tNumElem; i++)
//        {
//            // If this entity does not have children
//            if(!mXTKMesh.entity_has_children(i,EntityRank::ELEMENT))
//            {
//
//                // Get edges attached to element
//                tElementEdges = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::EDGE);
//
//                // Iterate through element edges and make a request
//                for(moris::size_t j = 0; j <tNumEdgePerElem; j++)
//                {
//                    // Get nodes attached to edge
//                    tEdgeNodes =  tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElementEdges(0,j),EntityRank::EDGE,EntityRank::NODE);
//
//                    // Compute Midisde Node Coordinate
//                    tEdgeCoords = mXTKMesh.get_selected_node_coordinates_loc_inds(tEdgeNodes);
//
//                    // Interpolate coordinate on edge using node coordinates
//                    // Must happen using local XTK index
//                    //TODO: MOVE ALL INTERPOLATION TO GEOMETRY ENGINE
//                    Interpolation::linear_interpolation_location(tEdgeCoords, tLocalCoord,tNodeGlobCoord);
//
//                    // tEdge Nodes is now the processor local index (notice the 1 versus above)
//                    tEdgeTopology.set_node_indices(tEdgeNodes);
//
//                    // Convert to global id using mesh
//                    (*tEdgeNodes)(0, 0) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
//                    (*tEdgeNodes)(0, 1) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);
//
//                    // Create a cantor pairing (since this is in conjunction with the child mesh a secondary identifier is needed)
//                    // Although it does not do anything in this case
//                    moris::size_t tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
//
//                    tNodeInds(i)(j) = aEdgeRequests.set_request_info(tElementEdges(0,j),tSecondaryId, tEdgeTopology, tNodeGlobCoord, tLocalCoord);
//
//                }
//            }
//        }
//
//        return tNodeInds;

//    }


    /*
     * After all edges have been assigned a midside node id, this function stores those values
     */
    void
    retrieve_midside_node_inds_unintersected(Cell<Cell<moris::size_t*>> aNodePointers)
    {
//        mesh::Mesh_Data<moris::real, moris::size_t, moris::DDRMat, moris::DDSTMat> & tXTKMeshData = mXTKMesh.get_mesh_data();
//        moris::size_t tNumElems = mMidsideElementToNode.n_rows();
//        moris::size_t tNumMidSideNodes = mMidsideElementToNode.n_cols();
//
//        for(moris::size_t i = 0; i<tNumElems; i++)
//        {
//            if(!mXTKMesh.entity_has_children(i,EntityRank::ELEMENT))
//            {
//                for(moris::size_t j = 0; j<tNumMidSideNodes; j++)
//                {
//                    mMidsideElementToNode(i,j) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(*(aNodePointers(i)(j)),EntityRank::NODE);
//                }
//            }
//        }
    }


    /*
     */
//    void
//    make_tet4_to_tet10_child_mesh_node_requests(Request_Handler & aEdgeRequests,
//                                                Request_Handler & aFaceRequests,
//                                                Request_Handler & aElementRequests)
//    {
//
//        mesh::Mesh_Data<moris::real, moris::size_t, moris::DDRMat, moris::DDSTMat> & tXTKMeshData = mXTKMesh.get_mesh_data();
//
//        // Number of Children Meshes
//        moris::size_t tNumChildrenMesh = mCutMesh.get_num_simple_meshes();
//        moris::size_t tNumChildrenElem = 0;
//        moris::size_t tNumFacePerElem  = 4;
//        moris::size_t tTotalNumChildrenElem = mCutMesh.get_num_entities(EntityRank::ELEMENT);
//        moris::Matrix< moris::DDSTMat >     tElemToNodeConn     = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        moris::Matrix< moris::DDSTMat >     tElemToEdgeConn     = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        moris::Matrix< moris::DDSTMat >     tEdgeToNodeIndices  = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        moris::Matrix< moris::DDSTMat >     tEdgeToNodeChildLoc = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        std::shared_ptr<Matrix_Base<moris::real, moris::DDRMat>>           tCoordSwapper       = mMatrixFactory.create_real_type_matrix_base(1, 3, REAL_MAX);
//        std::shared_ptr<Matrix_Base<moris::real, moris::DDRMat>>           tEdgeCoords         = mMatrixFactory.create_real_type_matrix_base(2, 3, REAL_MAX);
//        std::shared_ptr<Matrix_Base<moris::real, moris::DDRMat>>           tLocalCoord         = mMatrixFactory.create_real_type_matrix_base(1, 1, 0.0);
//        std::shared_ptr<Matrix_Base<moris::real, moris::DDRMat>>           tNodeGlobCoord      = mMatrixFactory.create_real_type_matrix_base(1, 3, REAL_MAX);
//        Edge_Topology tEdgeTopology;
//
//
//
//        //Number of elements in a children mesh
//        for(moris::size_t i = 0; i<tNumChildrenMesh; i++)
//        {
//
//            // If it the first subdivision we need to find the intersected before placing the conformal nodes
//            // Intersected elements are flagged via the Geometry_Engine
//
//            // Initialize
//            moris::size_t tEdgeInd = INTEGER_MAX;
//            moris::size_t tCount = 0;
//
//            // Parent Information
//            moris::Matrix< moris::DDSTMat > tParentInfo = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//
//            // Get full edge to element connectivity from XTK Mesh (probably slow)
//            // 0 specifies XTK local indices if an analytic geometry
//            // Otherwise this needs to be the processor local index
//            // or even the ID
//            moris::size_t tNumEdges = mCutMesh.get_num_entities(i,EntityRank::EDGE);
//
//            moris::Matrix< moris::DDSTMat > tNodeIndices = mCutMesh.get_all_node_inds(i);
//            Matrix<moris::real, moris::DDRMat>tNodeCoords        = tXTKMeshData.get_selected_node_coordinates_loc_inds(*tNodeIndices);
//
//            // Initialize node index pointers based on number of elements and a tet4 having 6 edges
//            tNumChildrenElem = tElemToNodeConn->n_rows();
//
//            moris::size_t tTotal = tNumChildrenElem*tNumFacePerElem;
//
//            std::shared_ptr<Matrix_Base<moris::real, moris::DDRMat>> tHalf = mMatrixFactory.create_real_type_matrix_base(1,1,0.5);
//
//
//
//
//            Cell<moris::size_t*> tNodeInds(tNumEdges);
//
//            // Loop over all edges
//            for (moris::size_t ed = 0; ed < tNumEdges; ed++)
//            {
//
//                bool tInterfaceFlag = false;
//                // Edge to Node
//                // Child Local
//                tEdgeToNodeChildLoc = mCutMesh.get_entities_connected_to_entity(i,EntityRank::EDGE,EntityRank::NODE,ed,0);
//
//                // Processor Indices
//                tEdgeToNodeIndices = mCutMesh.get_entities_connected_to_entity(i,EntityRank::EDGE,EntityRank::NODE,ed,1);
//
//                // Column 1 parent rank, column 2 parent index
//                tParentInfo = mCutMesh.get_parent_entity(i,EntityRank::EDGE,ed);
//
//                // Get Node 1s Sensitivity (node 1 relative to edge)
//                // If it is an interface node I can expect there to be sensitivity information, otherwise I need to assume it is 0
//                bool tNode1Interface = tNodes((*tEdgeToNodeChildLoc)(0,0)).is_interface_node();
//
//                // Get Node 2s Sensitivity
//                // If it is an interface node I can expect there to be sensitivity information, otherwise I need to assume it is 0
//                bool tNode2Interface = tNodes((*tEdgeToNodeChildLoc)(0,0)).is_interface_node();
//
//                // If both are an interface node, they have sensitivity and we need to average them
//                std::shared_ptr<Matrix_Base<moris::real, moris::DDRMat>>      tDxDp = mMatrixFactory.create_real_type_matrix_base(1,1);
//                moris::Matrix< moris::DDSTMat > tNewNodeADVS;
//
//                if( tNode1Interface && tNode2Interface)
//                {
//                    if(mGeometryEngines.mComputeDxDp)
//                    {
//                        Matrix_Base<moris::real,moris::DDRMat> * tNode1DxDp       = tNodes((*tEdgeToNodeChildLoc)(0,0)).get_dx_dp();
//                        Matrix_Base<moris::real,moris::DDRMat> * tNode2DxDp       = tNodes((*tEdgeToNodeChildLoc)(0,1)).get_dx_dp();
//                        moris::Matrix< moris::DDSTMat > * tNode1ADVS = tNodes((*tEdgeToNodeChildLoc)(0,0)).get_adv_indices();
//                        tNewNodeADVS = tNode1ADVS->copy();
//                        tDxDp->matrix_data() = 0.5*(tNode1DxDp->matrix_data() + tNode2DxDp->matrix_data());
//
//                    }
//                    tInterfaceFlag = true;
//                }
//
//                // If node 1 is an interface node and node 2 is not we take half of node 1's sensitivity
//                else if( tNode1Interface && !tNode2Interface )
//                {
//                    if(mGeometryEngines.mComputeDxDp)
//                    {
//                    Matrix_Base<moris::real,moris::DDRMat> * tNode1DxDp =  tNodes((*tEdgeToNodeChildLoc)(0,0)).get_dx_dp();
//
//                    tDxDp->matrix_data() = 0.5*(tNode1DxDp->matrix_data());
//                    moris::Matrix< moris::DDSTMat > * tNode1ADVS = tNodes((*tEdgeToNodeChildLoc)(0,0)).get_adv_indices();
//
//                    tNewNodeADVS = tNode1ADVS->copy();
//                    }
//
//                    tInterfaceFlag = true;
//                }
//
//                // If node 2 is an interface node and node 1 is not we take half of node 1 sensitivity
//                else if( !tNode1Interface && tNode2Interface )
//                {
//                    if(mGeometryEngines.mComputeDxDp)
//                    {
//                        Matrix_Base<moris::real,moris::DDRMat> * tNode2DxDp =  tNodes((*tEdgeToNodeChildLoc)(0,1)).get_dx_dp();
//                        tDxDp->matrix_data() = 0.5*(tNode2DxDp->matrix_data());
//                        moris::Matrix< moris::DDSTMat > * tNode2ADVS = tNodes((*tEdgeToNodeChildLoc)(0,1)).get_adv_indices();
//                        tNewNodeADVS = tNode2ADVS->copy();
//                    }
//                    tInterfaceFlag = true;
//                }
//
//
//                mCutMesh.add_entity_to_aux_connectivity(i, tCount, ed, 0);
//                tCount++;
//                // Add edge to the entity auxiliary connectivity
//
//                tNodeCoords->get_row((*tEdgeToNodeChildLoc)(0, 0), *tCoordSwapper);
//                tEdgeCoords->set_row(0, *tCoordSwapper);
//                tNodeCoords->get_row((*tEdgeToNodeChildLoc)(0, 1), *tCoordSwapper);
//                tEdgeCoords->set_row(1, *tCoordSwapper);
//
//                // Interpolate coordinate on edge using node coordinates
//                // Must happen using local XTK index
//                //TODO: MOVE ALL INTERPOLATION TO GEOMETRY ENGINE
//                Interpolation::linear_interpolation_location(*tEdgeCoords, *tLocalCoord,*tNodeGlobCoord);
//
//                // Get information about parent in stk mesh from XTK mesh
//                moris::size_t tParentRank = (*tParentInfo)(0, 0);
//
//                // tEdge Nodes is now the processor local index (notice the 1 versus above)
//                tEdgeTopology.set_node_indices(*tEdgeToNodeIndices);
//                // Convert to global id using mesh
//                (*tEdgeToNodeIndices)(0, 0) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index((*tEdgeToNodeIndices)(0, 0), EntityRank::NODE);
//                (*tEdgeToNodeIndices)(0, 1) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index((*tEdgeToNodeIndices)(0, 1), EntityRank::NODE);
//
//                // Order the nodes in ascending order
//                if((*tEdgeToNodeIndices)(0, 1) < (*tEdgeToNodeIndices)(0, 0))
//                {
//                    moris::size_t tSwap = (*tEdgeToNodeIndices)(0, 0);
//                    (*tEdgeToNodeIndices)(0, 0) = (*tEdgeToNodeIndices)(0, 1);
//                    (*tEdgeToNodeIndices)(0, 1) = tSwap;
//                }
//
//                if (tParentRank == 1)
//                {
//
//                    // Intersected edge is an existing stk edge
//                    // Make request in edge requests
//                    // This does not require a supplemental identifier
//
//                    moris::size_t tSecondaryId = xtk::cantor_pairing((*tEdgeToNodeIndices)(0, 0),(*tEdgeToNodeIndices)(0, 1));
//
//                    tNodeInds(ed) = aEdgeRequests.set_request_info((*tParentInfo)(0, 1),tSecondaryId, tEdgeTopology, *tNodeGlobCoord, *tLocalCoord);
//                }
//
//
//                else if (tParentRank == 2)
//                {
//                    // Intersected edge was built on an stk face
//                    // Make request in face requests
//                    // This requires a supplemental identifier
//                    moris::size_t tSecondaryId = xtk::cantor_pairing((*tEdgeToNodeIndices)(0, 0),(*tEdgeToNodeIndices)(0, 1));
//
//                    tNodeInds(ed) = aFaceRequests.set_request_info((*tParentInfo)(0, 1), tSecondaryId, tEdgeTopology, *tNodeGlobCoord, *tLocalCoord);
//                }
//                //
//                else if (tParentRank == 3)
//                {
//                    // Intersected edge was built in stk element
//                    // Make request in element requests
//                    // This requires a supplemental identifier
//                    moris::size_t tSecondaryId = xtk::cantor_pairing((*tEdgeToNodeIndices)(0, 0),(*tEdgeToNodeIndices)(0, 1));
//                    tNodeInds(ed) = aElementRequests.set_request_info((*tParentInfo)(0, 1), tSecondaryId, tEdgeTopology, *tNodeGlobCoord, *tLocalCoord);
//                }
//
//                else
//                {
//                    std::cout << "Invalid ancestry returned from XTK Mesh";
//                }
//
//
//                if(!tInterfaceFlag)
//                {
//                    mCutMesh.set_pending_node_index_pointers(i, tNodeInds(ed));
//                }
//                else
//                {
//
//                    Matrix_Base<moris::real,moris::DDRMat> *       tDxDpLocation;
//                    moris::Matrix< moris::DDSTMat > * tNodeADVLocation;
//                    mGeometryEngines.store_dx_dp(tDxDp,tNewNodeADVS,tDxDpLocation,tNodeADVLocation);
//
//
//
//                    mCutMesh.set_pending_node_index_pointers_with_dx_dp(i, tNodeInds(ed),tDxDpLocation,tNodeADVLocation);
//                }
//
//                tInterfaceFlag = false;
//            }
//
//
//            // Add node types to cell
//
//        }
//        bool tCoordinateFlag = false;
//    }

    /*
     * For nodes that are created during the decomposition process, tell
     * the XTK mesh about where they live in child meshes.
     */
    void
    associate_nodes_created_during_decomp_to_child_meshes()
    {
        // Initialize the data in the XTK mesh
        mBackgroundMesh.allocate_external_node_to_child_mesh_associations();

        // Number of children meshes
        size_t tNumCM = mCutMesh.get_num_child_meshes();
        for(size_t i = 0 ; i < tNumCM; i++)
        {
            // Get reference to the child mesh
            Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

            // Get reference to the nods in the child mesh node indices
            moris::Matrix<moris::IndexMat> const & tNodeIndices = tChildMesh.get_node_indices();

            // Associate these node indices with their child mesh index
            mBackgroundMesh.associate_external_nodes_to_child_mesh(i,tNodeIndices);
        }
    }

    /*
     * Set element phase index
     */
    void
    set_element_phases(moris::size_t aElementIndexOffset)
    {
        // Set element phase indices
         mBackgroundMesh.initialize_element_phase_indices(aElementIndexOffset);

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
    moris::mtk::Mesh*
    construct_output_mesh( Output_Options const & aOutputOptions );

public:
    moris::Cell< moris::Matrix < moris::DDRMat > >
    assemble_geometry_data_as_mesh_field(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput);
private:

    moris::Cell<std::string>
    assign_geometry_data_names();


    moris::Cell < enum moris::EntityRank >
    assign_geometry_data_field_ranks();


    /*!
     * Takes the whole node local to global map and removes the nodes
     * which are not part of the phase being output
     */
    moris::Matrix<moris::IndexMat>
    get_node_map_restricted_to_output_phases(Output_Options const &           aOutputOptions,
                                             moris::Matrix<moris::IndexMat> & aOutputtedNodeInds);


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
