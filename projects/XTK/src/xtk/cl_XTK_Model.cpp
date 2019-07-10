/*
 * cl_XTK_Model.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "fn_all_true.hpp"
#include "fn_unique.hpp"
#include "op_equal_equal.hpp"
#include "fn_sort.hpp"
#include "HDF5_Tools.hpp"
#include "cl_MTK_Visualization_STK.hpp"
//#include "cl_XTK_Enrichment.hpp"
namespace xtk
{
// ----------------------------------------------------------------------------------
// Constructor/Deconstructor Source code
// ----------------------------------------------------------------------------------
Model::~Model()
{
    if(mEnriched)
    {
        delete mEnrichment;
    }

    if(mGhost)
    {
        delete mGhostStabilization;
    }
}

Model::Model(uint aModelDimension,
             moris::mtk::Interpolation_Mesh* aMeshData,
             Geometry_Engine & aGeometryEngine,
             bool aLinkGeometryOnConstruction) :
                                  mSameMesh(false),
                                  mModelDimension(aModelDimension),
                                  mBackgroundMesh(aMeshData,aGeometryEngine),
                                  mCutMesh(mModelDimension),
                                  mGeometryEngine(aGeometryEngine),
                                  mConvertedToTet10s(false)
{
    MORIS_ASSERT(mModelDimension == 3,"Currently, XTK only supports 3D model dimensions");

    if(aLinkGeometryOnConstruction == true)
    {
        link_background_mesh_to_geometry_objects();
    }

    mBackgroundMesh.initialize_interface_node_flags(mBackgroundMesh.get_num_entities(EntityRank::NODE),mGeometryEngine.get_num_geometries());
}


void
Model::link_background_mesh_to_geometry_objects()
{
    // initialize geometry objects
    moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

    mGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes(tNodeCoords.n_rows());
    mGeometryEngine.initialize_geometry_object_phase_values(tNodeCoords);
    mLinkedBackground = true;
}

// ----------------------------------------------------------------------------------
// Decomposition Source code
// ----------------------------------------------------------------------------------
void
Model::decompose(Cell<enum Subdivision_Method> aMethods,
                 bool                          aSetPhase)
{
    // Start clock
    std::clock_t tTotalTime = std::clock();

    // Assert that there has been a link between geometry model and background mesh
    MORIS_ERROR(mLinkedBackground, "Geometry model and background mesh have not been linked via call to link_background_mesh_to_geometry_objects");

    // Process for a decomposition
    uint tNumDecompositions = aMethods.size();
    uint tNumGeometries = mGeometryEngine.get_num_geometries();

    print_decompsition_preamble(aMethods);

    // Tell the subdivision to assign node Ids if it is the only subdivision method (critical for outputting)
    // This is usually only going to happen in test cases
    // Note: the Conformal subdivision methods dependent on node ids for subdivision routine, the node Ids are set regardless of the below boolean

    bool tNonConformingMeshFlag = false;
    if(aMethods.size() == 1)
    {
        tNonConformingMeshFlag = true;
    }

    // Loop over each geometry and have an active child mesh indices list for each
    for(moris::size_t iGeom = 0; iGeom<tNumGeometries; iGeom++)
    {
        bool tFirstSubdivisionFlag = true;
        moris::Matrix< moris::IndexMat > tActiveChildMeshIndices(1,1,0);

        for (moris::size_t iDecomp = 0; iDecomp < tNumDecompositions; iDecomp++)
        {
            // start timing on this decomposition
            std::clock_t start = std::clock();

            // Perform subdivision
            this->decompose_internal(aMethods(iDecomp), tActiveChildMeshIndices, tFirstSubdivisionFlag, tNonConformingMeshFlag);

            // Change the first subdivision flag as false
            tFirstSubdivisionFlag = false;

            // print timing
            if(moris::par_rank() == 0 && mVerbose)
            {
                std::cout<<"XTK: Decomposition "<<get_enum_str(aMethods(iDecomp))<<" for geometry "<<iGeom<< " completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
                std::cout<<"XTK: Decomposition "<<get_enum_str(aMethods(iDecomp))<<" for geometry "<<iGeom<< " had "<<  tActiveChildMeshIndices.numel()<<" intersected background elements."<<std::endl;
            }
        }
        // If it's not the last geometry tell the geometry engine we're moving on
        if(iGeom!= tNumGeometries-1)
        {
            mGeometryEngine.advance_geometry_index();
        }
    }

    // Tell the xtk mesh to set all necessary information to finalize decomposition allowing
    // i.e set element ids, indices for children elements
    finalize_decomp_in_xtk_mesh(aSetPhase);

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Decomposition completed in " <<(std::clock() - tTotalTime) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

void
Model::decompose_internal(enum Subdivision_Method    const & aSubdivisionMethod,
                          moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                          bool const &                       aFirstSubdivision,
                          bool const &                       aSetIds)
{
    switch (aSubdivisionMethod)
    {
        case (Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8):
            {
            //            MORIS_ASSERT(tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, moris::EntityRank::ELEMENT, moris::EntityRank::NODE).numel() == 8, "NC_REGULAR_SUBDIVISION_HEX8 is for HEX8 meshes only.");
            MORIS_ASSERT(aFirstSubdivision,"NC_REGULAR_SUBDIVISION_HEX8 needs to be the first subdivision routine for each geometry");
            MORIS_ASSERT(mModelDimension == 3,"NC_REGULAR_SUBDIVISION_HEX8 needs to be done on a 3D mesh");

            // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones dont
            moris::Matrix< moris::IndexMat > tNewPairBool;
            run_first_cut_routine(TemplateType::HEX_8, 8,  aActiveChildMeshIndices,tNewPairBool);


            // initialize a struct of all the data we are keeping track of in this decomposition
            // intended to reduce the clutter of function inputs etc
            Decomposition_Data tDecompData;
            tDecompData.mSubdivisionMethod = Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;

            // number of intersected elements
            moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

            // make node requests for each intersected element
            this->decompose_internal_reg_sub_make_requests(aActiveChildMeshIndices,tNewPairBool,tDecompData);

            // specify a dummy secondary id (not really needed for this type of decomposition)
            tDecompData.tSecondaryIdentifiers = Cell<moris_index>(tDecompData.tNewNodeParentIndex.size(), MORIS_INDEX_MAX);

            moris_index tMessageTag = 60000; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine.get_num_geometries());

            // add nodes to child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_reg_sub(aActiveChildMeshIndices,tNewPairBool,tDecompData);

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex, tDecompData.tNewNodeOwner,tDecompData.tNewNodeSharing,tDecompData.tNewNodeCoordinate);

            // associate new nodes with geometry objects
            create_new_node_association_with_geometry(tDecompData);

            for(moris::size_t i = 0; i< tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {
                    mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(i),TemplateType::REGULAR_SUBDIVISION_HEX8);
                }
            }

            break;
            }
        case (Subdivision_Method::C_HIERARCHY_TET4):
            {

            // If it the first subdivision we need to find the intersected before placing the conformal nodes
            // Intersected elements are flagged via the Geometry_Engine
            if(aFirstSubdivision)
            {
                moris::Matrix< moris::IndexMat > tNewPairBool;
                run_first_cut_routine(TemplateType::TET_4, 4, aActiveChildMeshIndices,tNewPairBool);

                for(moris::size_t i = 0; i<aActiveChildMeshIndices.n_cols(); i++)
                {
                    Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,i));
                    tChildMesh.generate_connectivities(true,true,true);
                }

            }

            // For hex background meshes we have a three dimension parametric coordinate
            moris::size_t tDimParamCoord = 3;

            // For tet background meshes we have a 4-d parametric coordinate
            if(aFirstSubdivision)
            {
                tDimParamCoord = 4;
            }

            // initialize a struct of all the data we are keeping track of in this decomposition
            // intended to reduce the clutter of function inputs etc
            Decomposition_Data tDecompData;
            tDecompData.mSubdivisionMethod = Subdivision_Method::C_HIERARCHY_TET4;
            tDecompData.mConformalDecomp = true;
            tDecompData.mHasSecondaryIdentifier = true;

            // Initialize
            moris::size_t tEdgeInd = INTEGER_MAX;

            // Initialize topologies used in this method (all local coordinates are with respect to an edge)
            Edge_Topology tEdgeTopology;

            // initialize a couple of commonly used matrices in this method
            moris::Matrix< moris::DDRMat > tLocalCoordRelativeToEdge(1,1, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tGlobalCoord(1,3, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates(2,tDimParamCoord); // parametric coordinate of end nodes wrt parent element

            // Check type specified as conformal (could change this to enum)
            moris::size_t tCheckType = 1;
            moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

            // get the underlying background mesh data
            moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

            // resize child mesh to new node information
            tDecompData.tCMNewNodeLoc.resize(aActiveChildMeshIndices.n_cols());
            tDecompData.tCMNewNodeParamCoord.resize(aActiveChildMeshIndices.n_cols());

            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {

                // Initialize geometry objects
                Cell<Geometry_Object> tGeoObjects;

                // Get the child mesh that is active
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));

                // edge to node connectivity from child mesh
                moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                // Ask geometry engine which edges are intersected (Simple mesh local indexed edges)
                mGeometryEngine.is_intersected(tNodeCoords, tEdgeToNode, tCheckType, tGeoObjects);

                // Initialize node index pointers based on number of intersected edges and parametric coordinates
                uint tNumNewNodes = 0;
                Cell<moris::moris_index*>      tNodeInds(tGeoObjects.size());
                moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem(tGeoObjects.size(),tDimParamCoord);

                // get reference to child mesh edge parent information
                moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

                for (moris::size_t k = 0; k < tGeoObjects.size(); k++)
                {
                    if(!tGeoObjects(k).has_parent_nodes_on_interface())
                    {
                        // Local index to XTK Mesh
                        tEdgeInd = tGeoObjects(k).get_parent_entity_index();

                        // get a local coordinate along the intersected edge [-1,1]
                        tLocalCoordRelativeToEdge(0,0) = tGeoObjects(k).get_interface_lcl_coord();

                        // get the interpolated global coordinate
                        tGlobalCoord = tGeoObjects(k).get_interface_glb_coord();

                        // Add edge to the entity intersection connectivity
                        mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), tNumNewNodes, tEdgeInd, 0);

                        // Edge nodes
                        moris::Matrix<moris::IndexMat> tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);

                        // Compute new node parametric coordinate with respect to the current parent element
                        tEdgeNodeParamCoordinates.set_row(0, tChildMesh.get_parametric_coordinates(tEdgeNodes(0)));
                        tEdgeNodeParamCoordinates.set_row(1, tChildMesh.get_parametric_coordinates(tEdgeNodes(1)));
                        moris::Matrix< moris::DDRMat > tParametricCoordsRelativeToParentElem = Interpolation::linear_interpolation_location(tEdgeNodeParamCoordinates,tLocalCoordRelativeToEdge);

                        // Parent edge information
                        moris::size_t      tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                        moris::moris_index tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                        // get the owning processor for an entity
                        moris::moris_index tOwningProc = tMeshData.get_entity_owner(tParentIndex, (enum EntityRank)tParentRank);

                        // parent entity sharing information
                        moris::Matrix<moris::IdMat> tSharedProcs;
                          tMeshData.get_processors_whom_share_entity(tParentIndex, (enum EntityRank)tParentRank,tSharedProcs);

                        // Convert to global id using mesh
                        tEdgeNodes(0, 0) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                        tEdgeNodes(0, 1) = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                        // Order the nodes in ascending order
                        if(tEdgeNodes(0, 1) < tEdgeNodes(0, 0))
                        {
                            moris::size_t tSwap = tEdgeNodes(0, 0);
                            tEdgeNodes(0, 0) = tEdgeNodes(0, 1);
                            tEdgeNodes(0, 1) = tSwap;
                        }

                        // Intersected edge is an existing  edge
                        // Make request in edge requests
                        // This does not require a supplemental identifier
                        // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING!!!!!!
                        moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                        moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;
                        bool tRequestExist = tDecompData.request_exists(tParentIndex,tSecondaryId,(enum EntityRank)tParentRank,tNewNodeIndexInSubdivision);

                        // location for this face in the map
                        if(!tRequestExist)
                        {
                            tNewNodeIndexInSubdivision = tDecompData.register_new_request(tParentIndex,
                                                                                          tSecondaryId,
                                                                                          tOwningProc,
                                                                                          tSharedProcs,
                                                                                          (enum EntityRank)tParentRank,
                                                                                          tGlobalCoord,
                                                                                          new Edge_Topology(tEdgeToNode.get_row(tEdgeInd)), /*Note: this is deleted in the decomp data deconstructor*/
                                                                                          tLocalCoordRelativeToEdge.get_row(0));
                        }

                        // add to pending node pointers for child mesh
                        tDecompData.tCMNewNodeLoc(j).push_back(tNewNodeIndexInSubdivision);

                        // add parametric coordinate to decomp data
                        tDecompData.tCMNewNodeParamCoord(j).push_back(tParametricCoordsRelativeToParentElem);

                        // Creating a new node add 1 to count
                        tNumNewNodes++;
                    }

                    else if(tGeoObjects(k).all_parent_nodes_on_interface())
                    {
                        moris::moris_index tParentIndex = tGeoObjects(k).get_parent_entity_index();

                        // Tell the child mesh this edge is actually on the interface already
                        tChildMesh.mark_edge_as_on_interface(tParentIndex);

                        // Tell the xtk mesh that these edge nodes are interface nodes
                        mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tParentIndex,0),mGeometryEngine.get_active_geometry_index());
                        mBackgroundMesh.mark_node_as_interface_node(tEdgeToNode(tParentIndex,1),mGeometryEngine.get_active_geometry_index());
                    }
                } // geometry object

                tChildMesh.mark_interface_faces_from_interface_coincident_faces();
            } // XTK Mesh loop

            moris_index tMessageTag = 60001; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine.get_num_geometries());

            // add nodes to child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_nh(aActiveChildMeshIndices,tDecompData);

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeOwner, tDecompData.tNewNodeSharing,tDecompData.tNewNodeCoordinate);

            // associate new nodes with geometry objects
            create_new_node_association_with_geometry(tDecompData);

            // mark nodes as interface nodes
            moris_index tGeomIndex = mGeometryEngine.get_active_geometry_index();
            for(moris::uint i = 0; i <tDecompData.tNewNodeId.size(); i++)
            {
                mBackgroundMesh.mark_node_as_interface_node(tDecompData.tNewNodeIndex(i),tGeomIndex);
            }


            // Set Node Ids and tell the child mesh to update
            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {
                moris::Matrix< moris::IndexMat > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));
                moris::Matrix< moris::IdMat > tNodeIds = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index_range(tNodeIndices, EntityRank::NODE);

                mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j), tNodeIds);
                mCutMesh.modify_templated_mesh(aActiveChildMeshIndices(0,j), TemplateType::HIERARCHY_TET4);
            }

            break;
            }

        default:
        {
            moris::size_t breaker = 0;
            MORIS_ERROR(breaker != 0, "formulate_node_request should not enter the default case, check to see if your aCheckType is undefined.");
        }
    }
}


void
Model::decompose_internal_reg_sub_make_requests(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                moris::Matrix< moris::IndexMat > & tNewPairBool,
                                                Decomposition_Data & tDecompData)
{
    // mesh data accessor
    moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

    // number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // allocate child mesh to new node location
    tDecompData.tCMNewNodeLoc.resize(tIntersectedCount,7);
    tDecompData.tCMNewNodeParamCoord.resize(tIntersectedCount);

    // parametric coordinates relative to hex where we put the nodes
    // Parametric coordinates for this subdivision routine
    const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToElem(
            {{ 0.0, -1.0,  0.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        {-1.0,  0.0,  0.0},
        { 0.0,  0.0, -1.0},
        { 0.0,  0.0,  1.0},
        { 0.0,  0.0,  0.0}});

    const moris::Matrix< moris::DDRMat > tParamCoordsRelativeToFace(
            {{ 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0},
        { 0.0,  0.0}});


    // get the underlying background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // setup Child mesh to new node location
    for (moris::size_t i = 0; i < tIntersectedCount; i++)
    {
        if(tNewPairBool(0,i) == 0)
        {

            // Get element index
            moris::moris_index tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

            // Get local index of faces connected to element using local element index
            moris::Matrix<moris::IndexMat> tFaceIndices = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::FACE);

            // Loop over faces (6 in a hex 8) and set a node request.
            // Request will return a pointer to where the created node index will be placed
            for (moris::size_t fi = 0; fi < 6; fi++)
            {

                moris_index tRequestLoc = MORIS_INDEX_MAX;
                bool tRequestExists = tDecompData.request_exists(tFaceIndices(fi),EntityRank::FACE,tRequestLoc);

                // if we haven't created a node on this face then create one
                if(!tRequestExists)
                {
                    // node indices attached to face fi
                    moris::Matrix<moris::IndexMat> tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndices(fi), moris::EntityRank::FACE, moris::EntityRank::NODE);

                    // face owner
                    moris::moris_index tOwningProc = tMeshData.get_entity_owner(tFaceIndices(fi), EntityRank::FACE);

                    // procs whom share this face
                    moris::Matrix<moris::IdMat> tSharedProcs;
                    tMeshData.get_processors_whom_share_entity(tFaceIndices(fi), EntityRank::FACE,tSharedProcs);

                    // coordinates of nodes attached to the nodes of this face
                    moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                    // bilinearly interpolate to the center of this face fi
                    moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                    xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToFace.get_row(fi),tNewNodeCoordinates);

                    // location for this face in the map
                    moris_index tNewNodeIndexInSubdivision = tDecompData.register_new_request(tFaceIndices(fi),
                                                                                              tOwningProc,
                                                                                              tSharedProcs,
                                                                                              EntityRank::FACE,
                                                                                              tNewNodeCoordinates,
                                                                                              new Quad_4_Topology(tFaceNodes), /*Note: this is deleted in the decomp data deconstructor*/
                                                                                              tParamCoordsRelativeToFace.get_row(fi));
                    // add to pending node pointers for child mesh
                    tDecompData.tCMNewNodeLoc(i)(fi) = tNewNodeIndexInSubdivision;

                    // add parametric coordinate to decomp data
                    tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(fi));
                }

                // if debug check the coordinate will be the same
                else
                {
                    tDecompData.tCMNewNodeLoc(i)(fi) = tRequestLoc;
                    tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(fi));
#ifdef DEBUG
                    moris::uint tNewNodeIndexInSubdivision = tRequestLoc;

                    // node indices attached to face fi
                    moris::Matrix<moris::IndexMat> tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndices(fi), moris::EntityRank::FACE, moris::EntityRank::NODE);

                    // coordinates of nodes attached to the nodes of this face
                    moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                    // bilinearly interpolate to the center of this face fi
                    moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                    xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToFace.get_row(fi),tNewNodeCoordinates);

                    // other coordinate
                    moris::Matrix<moris::DDRMat> tExistingNodeCoordinate = tDecompData.tNewNodeCoordinate(tNewNodeIndexInSubdivision);

                    MORIS_ASSERT(all_true(tNewNodeCoordinates == tExistingNodeCoordinate) ,"Node coordinates created on same face do not match");
#endif

                }
            }

            // Place node at center of element
            // get the nodes attached to the element
            moris::Matrix<moris::IndexMat>tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);

            // coordinates of nodes attached to element
            moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodes);

            // trilinearly interpolate to the center of the element
            moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
            xtk::Interpolation::trilinear_interpolation(tCoordinates, tParamCoordsRelativeToElem.get_row(6), tNewNodeCoordinates);

            // add the new node at center of element to the map
            // location for this face in the map
            moris_index tNewNodeIndexInSubdivision = MORIS_INDEX_MAX;

            // owner of element
            moris::moris_index tOwningProc = tMeshData.get_entity_owner(tElemInd, EntityRank::ELEMENT);

            // procs sharing element
            moris::Matrix<moris::IdMat> tSharedProcs;
            tMeshData.get_processors_whom_share_entity(tElemInd, EntityRank::ELEMENT,tSharedProcs);


            MORIS_ASSERT(!tDecompData.request_exists(tElemInd,EntityRank::ELEMENT,tNewNodeIndexInSubdivision),"All element requests should be unique, therefore tNewRequest is expected to be true here");


            tNewNodeIndexInSubdivision = tDecompData.register_new_request(tElemInd,
                                                                          tOwningProc,
                                                                          tSharedProcs,
                                                                          EntityRank::ELEMENT,
                                                                          tNewNodeCoordinates,
                                                                          new Hexahedron_8_Topology(tElementNodes), /*Note: this is deleted in the decomp data deconstructor*/
                                                                          tParamCoordsRelativeToElem.get_row(6));

            // add child mesh new node location and parametric coordinate relative to element
            tDecompData.tCMNewNodeLoc(i)(6) = tNewNodeIndexInSubdivision;
            tDecompData.tCMNewNodeParamCoord(i).push_back(tParamCoordsRelativeToElem.get_row(6));

        }

    }

}


void
Model::decompose_internal_set_new_nodes_in_child_mesh_reg_sub(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                              moris::Matrix< moris::IndexMat > & tNewPairBool,
                                                              Decomposition_Data &               tDecompData)
{
    // number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // iterate through active child mesh indices
    for(moris::uint i = 0 ; i <tIntersectedCount; i++)
    {
        // only regularly subdivide if it hasnt already been regularly subdivided
        if(tNewPairBool(0,i) == 0)
        {

            // number of new nodes for child mesh i
            moris::uint tNumNewNodesForCM = tDecompData.tCMNewNodeLoc(i).size();

            // matrix of new node indices
            moris::Matrix<IndexMat> tCMNewNodeInds(1,tNumNewNodesForCM);

            // matrix of new node ids
            moris::Matrix<IdMat> tCMNewNodeIds(1,tNumNewNodesForCM);

            // iterate through new nodes for child mesh i to collect index and id
            for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
            {
                // location relative to the decomposition data
                moris::moris_index tNodeIndexInRequestVect = tDecompData.tCMNewNodeLoc(i)(iN);

                // retreive node index and id
                tCMNewNodeInds(iN) = tDecompData.tNewNodeIndex(tNodeIndexInRequestVect);
                tCMNewNodeIds(iN)  = tDecompData.tNewNodeId(tNodeIndexInRequestVect);
            }

            // retrieve child mesh
            Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(i));


            // add node indices and ids to child mesh
            tChildMesh.add_node_indices(tCMNewNodeInds);
            tChildMesh.add_node_ids(tCMNewNodeIds);

            // allocate space for parametric coordinates
            tChildMesh.allocate_parametric_coordinates(tNumNewNodesForCM,3);

            // iterate through nods and add parametric coordinate
            for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
            {
                tChildMesh.add_node_parametric_coordinate( tCMNewNodeInds(iN),tDecompData.tCMNewNodeParamCoord(i)(iN));
            }
        }
    }

}



void
Model::decompose_internal_set_new_nodes_in_child_mesh_nh(moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                                         Decomposition_Data &               tDecompData)
{
    // number of intersected elements
    moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();

    // iterate through active child mesh indices
    for(moris::uint i = 0 ; i <tIntersectedCount; i++)
    {
        // number of new nodes for child mesh i
        moris::uint tNumNewNodesForCM = tDecompData.tCMNewNodeLoc(i).size();

        // matrix of new node indices
        moris::Matrix<IndexMat> tCMNewNodeInds(1,tNumNewNodesForCM);

        // matrix of new node ids
        moris::Matrix<IdMat> tCMNewNodeIds(1,tNumNewNodesForCM);

        // iterate through new nodes for child mesh i to collect index and id
        for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
        {
            // location relative to the decomposition data
            moris::moris_index tNodeIndexInRequestVect = tDecompData.tCMNewNodeLoc(i)(iN);

            // retreive node index and id
            tCMNewNodeInds(iN) = tDecompData.tNewNodeIndex(tNodeIndexInRequestVect);
            tCMNewNodeIds(iN)  = tDecompData.tNewNodeId(tNodeIndexInRequestVect);
        }

        // retrieve child mesh
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(i));

        // add node indices and ids to child mesh
        tChildMesh.add_node_indices(tCMNewNodeInds);
        tChildMesh.add_node_ids(tCMNewNodeIds);

        // allocate space for parametric coordinates
        tChildMesh.allocate_parametric_coordinates(tNumNewNodesForCM,3);

        // iterate through nods and add parametric coordinate
        for(moris::uint iN =0; iN< tNumNewNodesForCM; iN++)
        {
            tChildMesh.add_node_parametric_coordinate( tCMNewNodeInds(iN),tDecompData.tCMNewNodeParamCoord(i)(iN));
        }
    }
}

void
Model::create_new_node_association_with_geometry(Decomposition_Data & tDecompData)
{
    // create geometry objects for each node
    mGeometryEngine.create_new_node_geometry_objects(tDecompData.tNewNodeIndex,
                                                     tDecompData.mConformalDecomp,
                                                     tDecompData.tNewNodeParentTopology,
                                                     tDecompData.tParamCoordRelativeToParent,
                                                     tDecompData.tNewNodeCoordinate);
}


void
Model::assign_node_requests_identifiers(Decomposition_Data & tDecompData,
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
Model::sort_new_node_requests_by_owned_and_not_owned(Decomposition_Data & tDecompData,
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
Model::assign_node_requests_owned_identifiers_and_setup_send(Decomposition_Data & aDecompData,
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
Model::outward_communicate_node_requests(Cell<Matrix<IndexMat>> & aSendData,
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
Model::inward_receive_node_requests(Cell<Matrix<IndexMat>> & aReceiveData,
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
Model::set_received_node_ids(Cell<Matrix<IndexMat>> & aReceiveData,
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
            else
            {

                std::cout<<"Parent Id   = "<<aReceiveData(i)(0,iN);
                std::cout<<" | Parent Rank  = " <<tParentRank;
                std::cout<<" | Secondary Id = " <<tSecondaryId;
                std::cout<<" | Proc Rank    = " <<par_rank()<<std::endl;

                MORIS_ASSERT(aDecompData.mSubdivisionMethod == Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 && tParentRank == 2," The only known case where hanging nodes are allowed to show up is on a face during NC_REGULAR_SUBDIVISION_HEX8");

                // nodes attached to face
                moris::Matrix<moris::IndexMat> tNodesOnFace = tMeshData.get_entity_connected_to_entity_loc_inds(tParentEntityInd,EntityRank::FACE,EntityRank::NODE);

                // compute vertex coordinate
                moris::Matrix<moris::DDRMat> tNodeCoordsOnFace(4,3);
                for(moris::uint j = 0; j < 4; j++)
                {
                    tNodeCoordsOnFace.get_row(j) = tMeshData.get_node_coordinate(tNodesOnFace(j)).get_row(0);
                }

                moris::Matrix<moris::DDRMat> tParamCoordCenterQuad( {{ 0.0,  0.0}} );
                moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                xtk::Interpolation::bilinear_interpolation(tNodeCoordsOnFace, tParamCoordCenterQuad, tNewNodeCoordinates);

                // specify a dummy parametric coordinate, this node should not show up in any place where it uses parametric coords on this processor
                moris::Matrix<moris::DDRMat> tDummyParamCoord = {{std::numeric_limits<real>::max(),std::numeric_limits<real>::max(),std::numeric_limits<real>::max()}};

                // process sharing is the proc which we have received the unexpected node from
                moris::Matrix<moris::IdMat> tSharingProc = {{(moris_index)i}};

                // make new request to add to decomp data
                moris_index tHangingNodeRequIndex = aDecompData.register_new_request(tParentEntityInd,
                                                                                     (moris_index)i,
                                                                                     tSharingProc,
                                                                                     EntityRank::FACE,
                                                                                     tNewNodeCoordinates,
                                                                                     new Quad_4_Topology(tNodesOnFace), /*Note: this is deleted in the decomp data deconstructor*/
                                                                                     tDummyParamCoord.get_row(0));

                // give node id and index
                aDecompData.tNewNodeIndex(tHangingNodeRequIndex) = aNodeInd;
                aNodeInd++;

                // set the new node id
                aDecompData.tNewNodeId(tHangingNodeRequIndex) = tNodeId;
                aDecompData.mNumNewNodesWithIds++;

                // add this node to be created
                aDecompData.tNewNodeHangingFlag.push_back(tRequestIndex);
                aDecompData.tNewNodeHangingWRTProcRank.push_back((moris_index)i);

            }
        }


    }
}

void
Model::handle_hanging_nodes_in_reg_sub( Decomposition_Data & aDecompData,
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

void
Model::finalize_decomp_in_xtk_mesh(bool aSetPhase)
{
    // Change XTK model decomposition state flag
    mDecomposed = true;

    // Sort the children meshes into groups
    this->sort_children_meshes_into_groups();

    // give each child cell its id (parallel consistent) and index (not parallel consistent)
    this->assign_child_element_identifiers();

    // add child element to local to global map
    this->add_child_elements_to_local_to_global_map();

    // Associate nodes created during decomposition to their child meshes
    this->associate_nodes_created_during_decomp_to_child_meshes();

    // Compute the child element phase using the geometry engine
    // a case where the phase may not be set is when we only do a
    // non-conformal decomposition
    if(aSetPhase)
    {
        this->set_element_phases();
    }

    // creates mtk cells for all child elements (parent elements are assumed to have mtk cells in the mtk mesh)
    this->create_child_element_mtk_cells();

    // identify local subphases in child mesh
    this->identify_local_subphase_clusters_in_child_meshes();

}

void
Model::assign_child_element_identifiers()
{
    // Set child element ids and indices
    moris::size_t tNumElementsInCutMesh = mCutMesh.get_num_entities(EntityRank::ELEMENT);

    // Allocate global element ids (these need to be give to the children meshes)
    moris_id    tElementIdOffset = mBackgroundMesh.allocate_entity_ids(tNumElementsInCutMesh, moris::EntityRank::ELEMENT);
    moris_index tElementIndOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);

    // set child elements ids in the children meshes which I own and dont share
    Cell<Child_Mesh*> const & tOwnedNotSharedChildMeshes = mCutMesh.get_owned_not_shared_child_meshes();
    for(moris::size_t i = 0; i<tOwnedNotSharedChildMeshes.size(); i++)
    {
        tOwnedNotSharedChildMeshes(i)->set_child_element_ids(tElementIdOffset);
        tOwnedNotSharedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
    }

    // set the child element ids in the children meshe which I own and share
    Cell<Child_Mesh*> const & tOwnedSharedChildMeshes = mCutMesh.get_owned_shared_child_meshes();
    for(moris::size_t i = 0; i<tOwnedSharedChildMeshes.size(); i++)
    {
        tOwnedSharedChildMeshes(i)->set_child_element_ids(tElementIdOffset);
        tOwnedSharedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
    }

    // prepare outward communication of owned not shared child meshes


    // prepare inward receive of not owned shared child meshes

    // set the child element ids in the children meshes which I do not own and share
    Cell<Child_Mesh*> const & tNotOwnedSharedChildMeshes = mCutMesh.get_not_owned_shared_child_meshes();
    for(moris::size_t i = 0; i<tNotOwnedSharedChildMeshes.size(); i++)
    {
        tNotOwnedSharedChildMeshes(i)->set_child_element_ids(tElementIdOffset);
        tNotOwnedSharedChildMeshes(i)->set_child_element_inds(tElementIndOffset);
    }

    // tell the background mesh about the new first available index
    mBackgroundMesh.update_first_available_index(tElementIndOffset,EntityRank::ELEMENT);
}


void
Model::add_child_elements_to_local_to_global_map()
{
    // get all child element ids and indexes
    Matrix<IndexMat> tChildElementInds = mCutMesh.get_all_element_inds();
    Matrix<IndexMat> tChildElementIds  = mCutMesh.get_all_element_ids();

    mBackgroundMesh.add_cells_to_global_to_local_map(tChildElementInds,tChildElementIds);


}

void
Model::sort_children_meshes_into_groups()
{
    // background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // my proc rank
    moris_index tProcRank = par_rank();

    // number of children meshes
    uint tNumChildrenMeshes = mCutMesh.get_num_child_meshes();

    // allocate data
    Cell<Child_Mesh*>   tOwnedNotSharedChildrenMeshes(tNumChildrenMeshes);
    Cell<Child_Mesh*>   tOwnedSharedChildrenMeshes(tNumChildrenMeshes);
    Cell<Child_Mesh*>   tNotOwnedSharedChildrenMeshes(tNumChildrenMeshes);
    Cell<Matrix<IdMat>> tOwnedSharedOtherProcs(tNumChildrenMeshes);
    Cell<moris_id>      tNotOwnedSharedOwningProc(tNumChildrenMeshes);

    // keep track of the number in each group
    uint tOwnedNotSharedCount = 0;
    uint tOwnedSharedCount    = 0;
    uint tNotOwnedSharedCount = 0;

    for(moris::size_t i = 0; i<mCutMesh.get_num_child_meshes(); i++)
    {
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        moris_index tParentCellInd = tChildMesh.get_parent_element_index();

        // get owner of parent cell
        moris_index tOwnerProc = tMeshData.get_entity_owner(tParentCellInd,EntityRank::ELEMENT);

        // if this processor does not own the element add it to the not owned shared list
        if(tOwnerProc != tProcRank)
        {
            tNotOwnedSharedChildrenMeshes(tNotOwnedSharedCount) = &tChildMesh;
            tNotOwnedSharedOwningProc(tNotOwnedSharedCount) = tOwnerProc;
            tNotOwnedSharedCount++;
        }

        else
        {
            Matrix<IdMat> tSharingProcs;
            tMeshData.get_processors_whom_share_entity(tParentCellInd,EntityRank::ELEMENT,tSharingProcs);

            // if there are no sharing procssors add it to owned not shared
            if(tSharingProcs.numel() == 0)
            {
                tOwnedNotSharedChildrenMeshes(tOwnedNotSharedCount) = &tChildMesh;
                tOwnedNotSharedCount++;
            }

            // if there are sharing processors add to owned shared
            else
            {
                tOwnedSharedChildrenMeshes(tOwnedSharedCount) = &tChildMesh;
                tOwnedSharedOtherProcs(tOwnedSharedCount) = tSharingProcs.copy();
                tOwnedSharedCount++;
            }
        }
    }

    // size out extra space
    tNotOwnedSharedChildrenMeshes.resize(tNotOwnedSharedCount);
    tNotOwnedSharedOwningProc.resize(tNotOwnedSharedCount);
    tOwnedNotSharedChildrenMeshes.resize(tOwnedNotSharedCount);
    tOwnedSharedChildrenMeshes.resize(tOwnedSharedCount);
    tOwnedSharedOtherProcs.resize(tOwnedSharedCount);

    // add to cut mesh
    mCutMesh.add_child_mesh_groups( tOwnedNotSharedChildrenMeshes,
                                    tOwnedSharedChildrenMeshes,
                                    tNotOwnedSharedChildrenMeshes,
                                    tOwnedSharedOtherProcs,
                                    tNotOwnedSharedOwningProc);
}

void
Model::associate_nodes_created_during_decomp_to_child_meshes()
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

void
Model::create_child_element_mtk_cells()
{

    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    moris::mtk::Mesh const & tMeshData = mBackgroundMesh.get_mesh_data();

    for(moris::uint i=0; i<tNumChildMeshes; i++)
    {
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        moris::moris_index tOwnerProc = tMeshData.get_entity_owner(tChildMesh.get_parent_element_index(),EntityRank::ELEMENT);

        moris::Matrix< moris::IdMat >    const & tElementIds  = tChildMesh.get_element_ids();
        moris::Matrix< moris::IndexMat > const & tElementInds = tChildMesh.get_element_inds();

        // Iterate over elements
        for(moris::uint j = 0; j<tElementIds.numel(); j++)
        {
            mBackgroundMesh.add_child_element_to_mtk_cells(tElementInds(j),tElementIds(j),tOwnerProc,j, &tChildMesh);
        }

    }
}


void
Model::identify_local_subphase_clusters_in_child_meshes()
{

    // get the number of children meshes
    moris::size_t tNumChildMeshes =  mCutMesh.get_num_child_meshes();

    // iterate over children meshes and perform local flood-fill
    for(moris::size_t i = 0; i<tNumChildMeshes; i++)
    {
        // Get child mesh index
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        // Perform local flood-fill on child mesh to identify subphase
        moris::Matrix< moris::IndexMat > tLocalFloodFill = local_child_mesh_flood_fill(tChildMesh);

        // Set the local floodfill data as the elemental subphase values in the child mesh
        // The child mesh then sorts the elements into bins
        tChildMesh.set_elemental_subphase(tLocalFloodFill);
    }
}

// ----------------------------------------------------------------------------------
// Sensitivity Source code
// ----------------------------------------------------------------------------------
void
Model::compute_sensitivity()
{
    // Start the clock
    std::clock_t start = std::clock();

    // verify the state of the xtk model
    MORIS_ERROR(mDecomposed,"Prior to computing sensitivity, the decomposition process must be called");
    MORIS_ERROR(!mConvertedToTet10s,"Prior to computing sensitivity, the convert tet4 to tet10 process was called");
    MORIS_ERROR(!mSensitivity,"Calling compute interface sensitivity twice is not supported");

    // Compute interface sensitivity
    compute_interface_sensitivity_internal();

    // Change the sensitivity computation state flag
    mSensitivity = true;

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Sensitivity computation completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

// ----------------------------------------------------------------------------------
// Unzipping Child Mesh Source code
// ----------------------------------------------------------------------------------
void
Model::unzip_child_mesh()
{
    // start the clock
    std::clock_t start = std::clock();

    MORIS_ERROR(mDecomposed,"Prior to unzip_child_mesh, the decomposition process must be called");

    // unzip the interface
    unzip_child_mesh_internal();

    mUnzipped = true;
    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Child mesh unzipping completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

void
Model::unzip_child_mesh_internal()
{
    // get the number of children meshes
    moris::size_t tNumChildMeshes =  mCutMesh.get_num_child_meshes();

    for(moris::size_t i = 0; i<tNumChildMeshes; i++)
    {
        // Get child mesh index
        //        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);
    }
}

// ----------------------------------------------------------------------------------
// Unzipping Interface Source code
// ----------------------------------------------------------------------------------
void
Model::unzip_interface()
{
    // start the clock
    std::clock_t start = std::clock();

    // unzip the interface
    unzip_interface_internal();

    mUnzipped = true;
    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Interface unzipping completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

void
Model::unzip_interface_internal()
{
    // Get the number of geometries (we need to unzip each interface)
    uint tNumGeoms = mGeometryEngine.get_num_geometries();

    // Get the interface nodes (wrt all geometries)
    Cell<moris::Matrix<moris::IndexMat>> tAllInterfaceNodeIds = mBackgroundMesh.get_interface_nodes_glob_ids();
    Cell<moris::Matrix<moris::IndexMat>> tAllInterfaceNodeInds = mBackgroundMesh.get_interface_nodes_loc_inds();

    // Keep count of which interface node index in matrices were in
    for(uint iG = 0; iG<tNumGeoms; iG++)
    {
        // Interface node wrt geometry iG
        moris::Matrix<moris::IndexMat> const & tInterfaceNodeIds = tAllInterfaceNodeIds(iG);
        moris::Matrix<moris::IndexMat> const & tInterfaceNodeInds = tAllInterfaceNodeInds(iG);

        // Number of interface nodes wrt geometry iG (and check the sizes of the interface information)
        uint tNumInterfaceNodes = tInterfaceNodeIds.numel();
        MORIS_ASSERT(tInterfaceNodeIds.numel() == tInterfaceNodeInds.numel(), "Interface Ids and Indices dimension mismatch");

        // Assign node ids and indices ( row - 0 Node ids, row 1 - New node ids)
        moris::Matrix<moris::IdMat> tNewUnzippedNodeIds((size_t)tNumInterfaceNodes);
        moris::Matrix<moris::IndexMat> tNewUnzippedNodeInds((size_t)tNumInterfaceNodes);
        this->unzip_interface_internal_assign_node_identifiers(tNumInterfaceNodes,tNewUnzippedNodeInds,tNewUnzippedNodeIds);

        // Add new nodes to the mesh (as a copy of the existing node)
        mBackgroundMesh.batch_create_new_nodes_as_copy_of_other_nodes(tInterfaceNodeInds,tNewUnzippedNodeIds,tNewUnzippedNodeInds);

        // Allocate space in background mesh interface node flags
        mBackgroundMesh.allocate_space_in_interface_node_flags(tNumInterfaceNodes, mGeometryEngine.get_num_geometries());

        // Mark the newly created nodes as interface nodes
        mBackgroundMesh.mark_nodes_as_interface_node_loc_inds(tNewUnzippedNodeInds,iG);

        // Link the new nodes to the geometry object of the node they were copied to
        mGeometryEngine.link_new_nodes_to_existing_geometry_objects(tInterfaceNodeInds,tNewUnzippedNodeInds);

        // unzip_child_mesh_index
        this->unzip_interface_internal_modify_child_mesh(iG,tInterfaceNodeInds,tNewUnzippedNodeInds,tNewUnzippedNodeIds);

    }

}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_internal_assign_node_identifiers(moris::uint aNumNodes,
                                                        moris::Matrix<moris::IdMat> & aUnzippedNodeIndices,
                                                        moris::Matrix<moris::IdMat> & aUnzippedNodeIds)
{
    // Verify sizes
    MORIS_ASSERT(aUnzippedNodeIndices.numel() == aNumNodes, "Size mismatch between aNumNodes and aUnzippedNodeIndices. Please pre-allocate these matrices ");
    MORIS_ASSERT(aUnzippedNodeIds.numel() == aNumNodes, "Size mismatch between aNumNodes and aUnzippedNodeIds.  Please pre-allocate these matrices ");

    // Ask the mesh for new node ids
    moris::moris_index tNodeIndexOffset = mBackgroundMesh.get_first_available_index(EntityRank::NODE);
    moris::moris_id    tNodeIdOffset    = mBackgroundMesh.allocate_entity_ids(aNumNodes,EntityRank::NODE);

    // Iterate new nodes and assign new node ids
    for( uint iN = 0; iN<aNumNodes; iN++)
    {

        // TODO: ADD PARALLEL OWNERSHIP STUFF HERE TO ASSIGN CONSISTENT NODE IDS
        // Give node global ids
        aUnzippedNodeIds(iN) = tNodeIdOffset;

        // increase the id offset
        tNodeIdOffset++;

        // Give nodes processor indices
        aUnzippedNodeIndices(iN) = tNodeIndexOffset;

        // increase the node index offset
        tNodeIndexOffset++;
    }

    // update the first available node index in our background mesh
    mBackgroundMesh.update_first_available_index(tNodeIndexOffset,EntityRank::NODE);

}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_internal_modify_child_mesh(moris::uint                         aGeometryIndex,
                                                  moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                                                  moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                                                  moris::Matrix<moris::IdMat> const & aUnzippedNodeIds)
{

    // from interface node indices, figure out which interface nodes live in which interface
    moris::Cell<moris::Cell< moris::moris_index >> tChildMeshInterfaceNodes = unzip_interface_internal_collect_child_mesh_to_interface_node(aInterfaceNodeIndices,aUnzippedNodeIndices,aUnzippedNodeIds);

    // Flag indicating there is an interface without an element pair
    bool tNoPairFlag = false;

    // Iterate through child meshes and add new unzipped nodes
    uint tCMIndex = 0;
    for(auto iCM = tChildMeshInterfaceNodes.begin(); iCM != tChildMeshInterfaceNodes.end(); ++iCM)
    {
        // number of interface nodes in this child mesh
        uint tNumCMInterfaceNodes = iCM->size();

        // Allocate matrices of interface node indices
        moris::Matrix< moris::IndexMat > tCMInterfaceNodeIndices(1,tNumCMInterfaceNodes);
        moris::Matrix< moris::IndexMat > tCMUnzippedInterfaceNodeIndices(1,tNumCMInterfaceNodes);
        moris::Matrix< moris::IdMat >    tCMUnzippedInterfaceNodeIds(1,tNumCMInterfaceNodes);

        // Collect information on the interface nodes on this child mesh
        for(moris::uint iN = 0; iN<tNumCMInterfaceNodes; iN++)
        {
            // node index local to the numbering scheme in interface nodes
            moris::moris_index tInterfaceLocInd = (*iCM)(iN);

            tCMInterfaceNodeIndices(iN)         = aInterfaceNodeIndices(tInterfaceLocInd);
            tCMUnzippedInterfaceNodeIndices(iN) = aUnzippedNodeIndices(tInterfaceLocInd);
            tCMUnzippedInterfaceNodeIds(iN)     = aUnzippedNodeIds(tInterfaceLocInd);
        }

        // Tell the child mesh to unzip it's interface
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(tCMIndex);

        // initialize the unzipping (which basically allocated node to element connectivity
        // in the child mesh because it is not typically needed)
        tChildMesh.initialize_unzipping();

        // Ask the child mesh to construct interface element pairs
        moris::Matrix<moris::IndexMat> tInterfaceElementPairsCMIndex;
        moris::Matrix<moris::IndexMat> tInterfaceSideOrdinals;

        tChildMesh.unzip_child_mesh_interface_get_interface_element_pairs(aGeometryIndex,tNoPairFlag,tInterfaceElementPairsCMIndex,tInterfaceSideOrdinals);

        // Convert the pairs to processor local indices because we need to be able to access the element phase index
        moris::Matrix<moris::IndexMat> tInterfaceElementPairs = tChildMesh.convert_to_proc_local_elem_inds(tInterfaceElementPairsCMIndex);

        // TODO: Add method to resolve cross child mesh element pairs for when the interface coincides with a parent face
        // NOTE: By using the sign of the geometry value, it really shouldnt take a whole lot of work to accomodated
        MORIS_ERROR(!tNoPairFlag," in unzip_interface_internal_modify_child_mesh, interface detected on a child mesh boundary. Currently, no method is implemented to resolve this");

        // Take the child mesh pairs and determine who gets which id
        // This output is either a 0 or 1, meaning the first or second element of the pair gets the unzipped nodes
        moris::Matrix< moris::IndexMat > tElementWhichKeepsOriginalNodes =
                this->unzip_interface_internal_assign_which_element_uses_unzipped_nodes(aGeometryIndex,tInterfaceElementPairs);

        // Get the elements on the boundary
        tChildMesh.unzip_child_mesh_interface(aGeometryIndex,
                                              tInterfaceElementPairsCMIndex,
                                              tElementWhichKeepsOriginalNodes,
                                              tCMInterfaceNodeIndices,
                                              tCMUnzippedInterfaceNodeIndices,
                                              tCMUnzippedInterfaceNodeIds);

        // Construct interface elements
        unzip_interface_construct_interface_elements(aGeometryIndex,
                                                     tInterfaceElementPairs,
                                                     tInterfaceSideOrdinals);


        tChildMesh.finalize_unzipping();

        tCMIndex++;
    }

    unzip_interface_assign_element_identifiers();


}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat >
Model::unzip_interface_internal_assign_which_element_uses_unzipped_nodes( moris::moris_index aGeometryIndex,
                                                                          moris::Matrix< moris::IndexMat > const & aInterfaceElementPairs )
{

    // specify which geometry sign gets to keep the nodes
    moris::moris_index tValWhichUsesUnzipped = 0;

    // The rule used here is whichever phase gets a 1 from the phase table with respect to the current geometry index
    // gets to keep the original node. The other element changes it's nodes to the unzipped indices.
    moris::Matrix< moris::IndexMat > tElementWhichKeepsUsesUnzippedNodes(aInterfaceElementPairs.n_cols());

    // number of pairs
    moris::uint tNumPairs = aInterfaceElementPairs.n_cols();

    // allocate
    moris::moris_index tElement0           = MORIS_INDEX_MAX;
    moris::moris_index tElement1           = MORIS_INDEX_MAX;
    moris::moris_index tElement0PhaseIndex = MORIS_INDEX_MAX;
    moris::moris_index tElement1PhaseIndex = MORIS_INDEX_MAX;

    // iterate through pairs
    for(moris::uint iP = 0; iP<tNumPairs; iP++)
    {
        // set this back to moris index max so we dont run into issues if there is actually only one element in the pair
        moris::moris_index tElement0GeomSign   = MORIS_INDEX_MAX; // sign of phase of element wrt to aGeometryIndex (0 for negative, 1 for positive)
        moris::moris_index tElement1GeomSign   = MORIS_INDEX_MAX;// sign of phase of element wrt to aGeometryIndex (0 for negative, 1 for positive)

        // Element indices
        tElement0 = aInterfaceElementPairs(0,iP);
        tElement1 = aInterfaceElementPairs(1,iP);

        // Figure out which element in the pair gets to keep the original
        if(tElement0 != MORIS_INDEX_MAX)
        {
            tElement0PhaseIndex = mBackgroundMesh.get_element_phase_index(tElement0);
            tElement0GeomSign = mGeometryEngine.get_phase_sign_of_given_phase_and_geometry(tElement0PhaseIndex,aGeometryIndex);

            if(tElement0GeomSign == tValWhichUsesUnzipped)
            {
                tElementWhichKeepsUsesUnzippedNodes(iP) = 0;
            }
        }

        if(tElement1 != MORIS_INDEX_MAX)
        {
            tElement1PhaseIndex = mBackgroundMesh.get_element_phase_index(tElement1);
            tElement1GeomSign = mGeometryEngine.get_phase_sign_of_given_phase_and_geometry(tElement1PhaseIndex,aGeometryIndex);
            if(tElement1GeomSign == tValWhichUsesUnzipped)
            {
                tElementWhichKeepsUsesUnzippedNodes(iP) = 1;
            }
        }

        // Make sure both don't end up with the same sign
        MORIS_ASSERT(tElement1GeomSign != tElement0GeomSign,"Both elements in an interface pair returned the same phase sign");

        MORIS_ASSERT(tElement0GeomSign != MORIS_INDEX_MAX,"tElement0GeomSign no pair");
        MORIS_ASSERT(tElement1GeomSign != MORIS_INDEX_MAX,"tElement1GeomSign no pair");

    }

    return tElementWhichKeepsUsesUnzippedNodes;

}


// ----------------------------------------------------------------------------------
moris::Cell<moris::Cell< moris::moris_index >>
Model::unzip_interface_internal_collect_child_mesh_to_interface_node(moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                                                                     moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                                                                     moris::Matrix<moris::IdMat> const & aUnzippedNodeIds)
{
    // Allocate cell to keep track of the node indices of each child mesh
    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    moris::Cell<moris::Cell< moris::moris_index >> tChildMeshInterfaceNodes(tNumChildMeshes);

    // Iterate through interface node indices
    for(moris::uint iN = 0; iN < aInterfaceNodeIndices.numel(); iN++)
    {
        // Get the child meshes which have the interface node
        moris::Matrix< moris::IndexMat > tNodeChildMeshIndices = mBackgroundMesh.get_node_child_mesh_assocation(aInterfaceNodeIndices(iN));

        // Iterate through child meshes and mark interface node indices in this child mesh
        for(moris::uint iCM = 0; iCM<tNodeChildMeshIndices.numel(); iCM++)
        {
            moris::moris_index tCMIndex = tNodeChildMeshIndices(iCM);
            tChildMeshInterfaceNodes(tCMIndex).push_back(iN);
        }
    }

    return tChildMeshInterfaceNodes;
}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_construct_interface_elements(moris::uint aGeometryIndex,
                                                    moris::Matrix< moris::IndexMat > const & aElementPairs,
                                                    moris::Matrix< moris::IndexMat > const & aSideOrdinalPairs)
{
    moris::uint tNumPairs = aElementPairs.n_cols();

    moris::Cell<const moris::mtk::Cell*> tPairCells(2);

    for(moris::uint i = 0; i <tNumPairs; i++)
    {
        tPairCells(0) = mBackgroundMesh.get_child_element_mtk_cell_ptr(aElementPairs(0,i));
        tPairCells(1) = mBackgroundMesh.get_child_element_mtk_cell_ptr(aElementPairs(1,i));

        // give ownership information to left element
        moris_index tOwnerOfElem0 = tPairCells(0)->get_owner();

        // construct an interface element
        Interface_Element tInterfaceElement;
        tInterfaceElement.set_element_pair_and_side_ordinal(tPairCells,aSideOrdinalPairs.get_column(i));
        tInterfaceElement.set_element_owner(tOwnerOfElem0);

        // Add interface element to cut mesh
        mCutMesh.add_interface_element(tInterfaceElement);
    }
}
// ----------------------------------------------------------------------------------
void
Model::unzip_interface_assign_element_identifiers()
{
    moris::Cell<Interface_Element> & tInterfaceElements = mCutMesh.get_interface_elements();
    moris::uint tNumInterfaceElements = tInterfaceElements.size();

    // Allocate ids
    moris::moris_index tIdOffset    = mBackgroundMesh.allocate_entity_ids(tNumInterfaceElements,EntityRank::ELEMENT);
    moris::moris_index tIndexOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);

    for(moris::uint i = 0; i <tNumInterfaceElements; i++)
    {
        tInterfaceElements(i).set_element_identifiers(tIndexOffset,tIdOffset);
        tIndexOffset++;
        tIdOffset++;
    }

    mBackgroundMesh.update_first_available_index(tIndexOffset,EntityRank::ELEMENT);

}



// ----------------------------------------------------------------------------------
// Enrichment Source code
// ----------------------------------------------------------------------------------
void
Model::perform_basis_enrichment()
{
    // Start the clock
    std::clock_t start = std::clock();


    MORIS_ERROR(mDecomposed,"Prior to computing basis enrichment, the decomposition process must be called");
    MORIS_ERROR(mUnzipped,"Prior to computing basis enrichment, the interface unzipping process must be called");
    MORIS_ERROR(!mEnriched,"Calling perform_basis_enrichment twice is not supported");
    perform_basis_enrichment_internal();

    // Change the enrichment flag
    mEnriched = true;

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Basis enrichment computation completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}


Enrichment const &
Model::get_basis_enrichment()
{
    MORIS_ASSERT(mEnriched,"Cannot get basis enrichment from an XTK model which has not called perform_basis_enrichment ");
    return *mEnrichment;
}

void
Model::perform_basis_enrichment_internal()
{
    // initialize enrichment (ptr because of circular dependency)
    mEnrichment = new Enrichment(Enrichment_Method::USE_INTERPOLATION_CELL_BASIS,
                                 mGeometryEngine.get_num_phases(),
                                 this,
                                 &mCutMesh,
                                 &mBackgroundMesh);

    // Set verbose flag to match XTK.
    mEnrichment->mVerbose = mVerbose;

    // perform the enrichment
    mEnrichment->perform_enrichment();

    // Link Enrichment Results to Vertex interpolation
    this->link_vertex_enrichment_to_vertex_interpolation();
}

void
Model::link_vertex_enrichment_to_vertex_interpolation()
{
    // Iterate through nodes and add vertex interpolation for each
    moris::uint tNumNodes = mBackgroundMesh.get_num_entities(EntityRank::NODE);
    for(moris::moris_index  iN = 0; iN<(moris::moris_index)tNumNodes; iN++)
    {
        // Vertex Interpolation for node iN
        moris::mtk::Vertex_Interpolation_XTK & tVertexInterp = mBackgroundMesh.get_mtk_vertex_interpolation(iN);

        // Vertex enrichment for node iN
        Vertex_Enrichment & tVertexEnrichment = mEnrichment->get_vertex_enrichment(iN);

        // Vertex enrichment for node iN
        moris::mtk::Vertex_XTK & tVertex = mBackgroundMesh.get_mtk_vertex_xtk(iN);

        // set the vertex enrichme in the vertex interpolation
        tVertexInterp.set_vertex_enrichment(&tVertexEnrichment);

        // Set the vertex interpolation in the vertex
        tVertex.set_vertex_interpolation(&tVertexInterp);

    }

}

void
Model::construct_face_oriented_ghost_penalization_cells()
{
    MORIS_ERROR(mDecomposed,"Mesh needs to be decomposed prior to calling ghost penalization");
    MORIS_ERROR(!mGhost,"Ghost penalization has already been called");

    std::clock_t start = std::clock();

    mGhostStabilization = new Ghost_Stabilization(this);

    mGhostStabilization->setup_ghost_stabilization_facets();

    mGhost = true;

    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Ghost stabilization setup completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }


}


void
Model::perform_multilevel_enrichment_internal()
{
    //        mEnrichment->create_multilevel_enrichments();
}

// ----------------------------------------------------------------------------------
// Tet 10 conversion Source code
// ----------------------------------------------------------------------------------
void
Model::convert_mesh_tet4_to_tet10()
{
    MORIS_ASSERT(0,"not currently working");
    mConvertedToTet10s = true;

    // start timing on this decomposition
    std::clock_t start = std::clock();
    //        convert_mesh_tet4_to_tet10_internal();

    if(moris::par_rank() == 0)
    {
        std::cout<<"Tet4 to Tet10 conversion completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
}

// ----------------------------------------------------------------------------------
// Export mesh Source code
// ----------------------------------------------------------------------------------
moris::mtk::Mesh*
Model::get_xtk_as_mtk()
{
    return new moris::mtk::XTK_Impl(this);
}

void
Model::extract_surface_mesh_to_obj(std::string                      aOutputFile,
                                   size_t                           aPhaseIndex,
                                   moris::Cell<std::string> const & aBoundingSideSets)
{
    // start timing on this decomposition
    std::clock_t start = std::clock();

    // create the output mtk mesh
    extract_surface_mesh_to_obj_internal(aOutputFile,aPhaseIndex,aBoundingSideSets);

    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Extract surface to obj completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        std::cout<<"XTK: OBJ File: "<<aOutputFile<<std::endl;
    }
}

void
Model::extract_surface_mesh_to_obj_internal(std::string                      aOutputFile,
                                            size_t                           aPhaseIndex,
                                            moris::Cell<std::string> const & aBoundingSideSets)
{
    MORIS_ERROR(aBoundingSideSets.size() == 6," There needs to be 6 side sets which skin the mesh to extract the surface");

    // allocate side set data
    moris::Cell<moris::Matrix<IndexMat>> tElementIndAndSideOrds(2*6);
    moris::Cell<moris::mtk::MtkSideSetInfo> tBackgroundSideSets(2*6);

    // background mesh data
    moris::mtk::Mesh const & tMeshData = mBackgroundMesh.get_mesh_data();

    // setup outputting options
    Output_Options tOutputOptions;

    // Specify there are 2 possible phases
    size_t tNumPhases = mGeometryEngine.get_num_bulk_phase();

    // Say I only want to output phase 1
    Cell<size_t> tPhasesToOutput = {aPhaseIndex};
    tOutputOptions.change_phases_to_output(tNumPhases,tPhasesToOutput);

    // Initialize Sets information structure
    moris::mtk::MtkSetsInfo tMtkMeshSets;

    moris::uint tNumChildCells   = 0;
    moris::uint tNumNoChildCells = 0;
    for(moris::uint i = 0; i <aBoundingSideSets.size(); i++)
    {
        moris::uint tNoChildInd = 2*i;
        moris::uint tChildInd = 2*i+1;

        this->propogate_background_side_set(aBoundingSideSets(i),tNoChildInd,tChildInd,tElementIndAndSideOrds,tBackgroundSideSets,tOutputOptions,true);

        tNumNoChildCells = tNumNoChildCells + tElementIndAndSideOrds(tNoChildInd).n_rows();
        tNumChildCells = tNumChildCells + tElementIndAndSideOrds(tChildInd).n_rows();
    }

    // get the interface information
    Cell<moris::Matrix<moris::IndexMat>> tInterfaceNodes = mBackgroundMesh.get_interface_nodes_glob_ids();

    // interface sides
    moris::Matrix<moris::IdMat> tInterfaceElemIndandSideOrd = mCutMesh.pack_interface_sides(true,aPhaseIndex);


    // tri 3s
    moris::Matrix<moris::IdMat> tTri3ElemToNode(tNumChildCells + tInterfaceElemIndandSideOrd.n_rows() + tNumNoChildCells*4,3);

    // keep track of nodes that are in skinned mesh
    moris::uint tNodeCount = 0;
    moris::Matrix<moris::IdMat> tNodesOnBoundary((tNumNoChildCells*5+tNumChildCells*3 + tInterfaceNodes(aPhaseIndex).numel()) ,1);
    tNodesOnBoundary.fill(MORIS_ID_MAX);

    // collect tri nodes
    moris_index tTriCount = 0;
    for(moris::uint  i = 0; i<6; i++)
    {
        moris::uint tChildInd = 2*i+1;


        // iterate through cells in this side set
        for(moris::uint  j = 0; j <tElementIndAndSideOrds(tChildInd).n_rows(); j++)
        {
            // get the mtk cell
            moris::mtk::Cell const & tCell = mBackgroundMesh.get_mtk_cell(tElementIndAndSideOrds(tChildInd)(j,0));

            moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tCell.get_vertices_on_side_ordinal(tElementIndAndSideOrds(tChildInd)(j,1));

            Matrix<IdMat> tVertexIds( 1, tVerticesOnSide.size());
            for(moris::uint k =0; k < 3; k++)
            {
                tVertexIds(k) = tVerticesOnSide(k)->get_id();
                tNodesOnBoundary(tNodeCount) = tVertexIds(k);
                tNodeCount++;
            }

            tTri3ElemToNode.get_row(tTriCount) = tVertexIds.get_row(0);
            tTriCount++;
        }
    }

    // Collect the nodes on the surface and construct a quad 4 / tri 3 mesh
    // quad 4s
    uint tQuadCount = 0;
    moris::Matrix<moris::IdMat> tQuad4ElemToNode(tNumNoChildCells,4);

    for(moris::uint i = 0; i < 6; i++)
    {
        moris::uint tNoChildInd = 2*i;

        for(moris::uint  j = 0; j <tElementIndAndSideOrds(tNoChildInd).n_rows(); j++)
        {
            // get the mtk cell
            moris::mtk::Cell const & tCell = mBackgroundMesh.get_mtk_cell(tElementIndAndSideOrds(tNoChildInd)(j,0));

            moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tCell.get_vertices_on_side_ordinal(tElementIndAndSideOrds(tNoChildInd)(j,1));

            Matrix<IdMat> tVertexIds( 1, tVerticesOnSide.size());

            for(moris::uint k =0; k < 4; k++)
            {
                tVertexIds(k) = tVerticesOnSide(k)->get_id();
                tNodesOnBoundary(tNodeCount) = tVertexIds(k);
                tNodeCount++;
            }
            tQuad4ElemToNode.get_row(tQuadCount) = tVertexIds.get_row(0);
            tQuadCount++;
        }
    }

    // triangulate quad 4s
    //
    // allocate ids for triangulation
    moris_id tNodeId = mBackgroundMesh.allocate_entity_ids(tQuad4ElemToNode.n_rows(),EntityRank::NODE);

    // template for splitting the quad 4 into tris
    moris::Matrix<moris::IdMat> tTriangulatedQuadMap = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};

    moris::Matrix<moris::IdMat> tNodeIdsForTemplate(1,5);

    moris::Matrix<moris::IdMat>  tNewVertexId(tQuad4ElemToNode.n_rows(),1);
    moris::Matrix<moris::DDRMat> tNewVertexCoords(tQuad4ElemToNode.n_rows(),3);
    moris::Matrix<moris::DDRMat> tParamCoordCenterQuad( {{ 0.0,  0.0}} );

    for(moris::uint i = 0; i<tQuad4ElemToNode.n_rows(); i++)
    {
        // assign ids needed for template
        tNodeIdsForTemplate({0,0},{0,3})              = tQuad4ElemToNode.get_row(i);
        tNodeIdsForTemplate(4)                        = tNodeId;
        tNewVertexId(i)                               = tNodeId;
        tNodeId++;

        // turn quad into triangles
        moris::Matrix<moris::IdMat> tTriangulatedQuad = reindex_matrix(tTriangulatedQuadMap,0,tNodeIdsForTemplate);

        // compute vertex coordinate
        moris::Matrix<moris::DDRMat> tNodeCoordsOnFace(4,3);
        for(moris::uint j = 0; j < 4; j++)
        {
            moris_index tNodeIndex = tMeshData.get_loc_entity_ind_from_entity_glb_id(tNodeIdsForTemplate(j),EntityRank::NODE);
            tNodeCoordsOnFace.get_row(j) = tMeshData.get_node_coordinate(tNodeIndex).get_row(0);
        }

        // bilinearly interpolate to the center of this face fi
        moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
        xtk::Interpolation::bilinear_interpolation(tNodeCoordsOnFace, tParamCoordCenterQuad, tNewNodeCoordinates);
        tNewVertexCoords.get_row(i) = tNewNodeCoordinates.get_row(0);

        // add quad to tri 3s
        for(moris::uint  j = 0; j < 4; j++)
        {
            tTri3ElemToNode.get_row(tTriCount) = tTriangulatedQuad.get_row(j);
            tTriCount++;
        }
    }

    // add interface nodes to nodes on boundary
    tNodesOnBoundary({tNodeCount, tNodeCount + tInterfaceNodes(aPhaseIndex).numel()-1},{0,0}) = moris::trans(tInterfaceNodes(aPhaseIndex));
    tNodeCount = tNodeCount + tInterfaceNodes(aPhaseIndex).numel();

    // add interface facets

    // iterate through cells in interface
    for(moris::uint  i = 0; i <tInterfaceElemIndandSideOrd.n_rows(); i++)
    {
        // get the mtk cell

        moris::mtk::Cell const & tCell = mBackgroundMesh.get_mtk_cell(tInterfaceElemIndandSideOrd(i,0));
        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tCell.get_vertices_on_side_ordinal(tInterfaceElemIndandSideOrd(i,1));

        Matrix<IdMat> tVertexIds( 1, tVerticesOnSide.size());
        for(moris::uint k =0; k < 3; k++)
        {
            tVertexIds(k) = tVerticesOnSide(k)->get_id();
        }

        tTri3ElemToNode.get_row(tTriCount) = tVertexIds.get_row(0);
        tTriCount++;
    }

    // assign element ids
    //fixme: Make these ids unique across processors
    moris::Matrix<moris::IdMat> tChildFaceElementIds(tTriCount,1);
    for(moris::uint  i = 0; i <tChildFaceElementIds.numel(); i++)
    {
        tChildFaceElementIds(i) = (moris_id)i+1;
    }


    // Interface elements
    moris::mtk::MtkBlockSetInfo tTri3Block;
    tTri3Block.mCellIdsInSet = &tChildFaceElementIds;
    tTri3Block.mBlockSetName = "tri_surf";
    tTri3Block.mBlockSetTopo = CellTopology::TRI3;
    tMtkMeshSets.add_block_set(&tTri3Block);

    // Convert to vertex indices
    tNodesOnBoundary.resize(tNodeCount,1);
    moris::Matrix<moris::IdMat> tNodeMap;

    moris::unique( tNodesOnBoundary,tNodeMap);
    moris::Matrix<moris::IndexMat> tNodeIndices(tNodeMap.numel(),1);

    for(moris::uint i = 0; i <tNodeMap.numel(); i++)
    {
        tNodeIndices(i) = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
    }

    // Get the node coordinates
    moris::Matrix<moris::DDRMat> tNodeCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tNodeIndices);

    // add nodes created during quad 4
    uint tNumNodeNoMidside = tNodeCoordinates.n_rows();
    uint tNumMidsideNodes = tNewVertexId.numel();
    uint tTotalNumNodes  =  tNumNodeNoMidside + tNumMidsideNodes;
    tNodeCoordinates.resize(tNumNodeNoMidside + tNumMidsideNodes,3);
    tNodeCoordinates({tNumNodeNoMidside,tTotalNumNodes-1},{0,2}) = tNewVertexCoords.matrix_data();

    tNodeMap.resize(tTotalNumNodes, 1);
    tNodeMap({tNumNodeNoMidside, tTotalNumNodes-1},{0,0}) = tNewVertexId.matrix_data();

    // make vertices consecutive for obj output
    std::unordered_map<moris_index, moris_index> tObjIndex;
    for(moris::uint i = 0; i < tNodeMap.numel(); i++)
    {
        tObjIndex[tNodeMap(i)] = i+1;
    }

    // renumber vertex
    for(moris::uint i  = 0; i <tTri3ElemToNode.n_rows(); i++)
    {
        for(moris::uint  j = 0; j<tTri3ElemToNode.n_cols(); j++)
        {
            auto tIter = tObjIndex.find(tTri3ElemToNode(i,j));
            tTri3ElemToNode(i,j) = tIter->second;
        }
    }

     // write to an obj file
    std::string tParObjPath =  moris::make_path_parallel( aOutputFile );

    std::ofstream tOFS (tParObjPath, std::ofstream::out);

     tOFS << std::setprecision (15);

     tOFS << "# XTK OBJ EXTRACTION"<<std::endl;

     // add vertices
     for(moris::uint i = 0 ; i < tNodeCoordinates.n_rows(); i++)
     {
         tOFS << std::scientific << "v "<< tNodeCoordinates(i,0) <<" "<< tNodeCoordinates(i,1) << " "<< tNodeCoordinates(i,2) <<std::endl;
     }

     // add facets
     for(moris::uint i = 0; i < tTri3ElemToNode.n_rows(); i++)
     {
         tOFS << "f "<< tTri3ElemToNode(i,0) <<" "<< tTri3ElemToNode(i,1)<< " "<< tTri3ElemToNode(i,2) <<std::endl;
     }

     tOFS.close();
}

moris::mtk::Integration_Mesh*
Model::get_output_mesh(Output_Options const & aOutputOptions)

{
    // start timing on this decomposition
    std::clock_t start = std::clock();

    // create the output mtk mesh
    moris::mtk::Integration_Mesh* tOutputMesh = construct_output_mesh(aOutputOptions);

    if(moris::par_rank() == 0  && mVerbose)
    {
        std::cout<<"XTK: Mesh output completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }
    return tOutputMesh;
}

//------------------------------------------------------------------------------

moris::uint
Model::get_num_elements_total()
{
    MORIS_ASSERT(mDecomposed,"Prior to using get_num_elements, the decomposition process must be finished");

    return mBackgroundMesh.get_num_entities(EntityRank::ELEMENT) + mCutMesh.get_num_entities(EntityRank::ELEMENT);
}

moris::uint
Model::get_num_elements_unzipped()
{
    return this->get_num_elements_total()- mCutMesh.get_num_child_meshes();
}

//------------------------------------------------------------------------------

moris::mtk::Integration_Mesh*
Model::construct_output_mesh( Output_Options const & aOutputOptions )
{

    // start timing on this decomposition
    std::clock_t start = std::clock();

    // Get mesh information ready for outputting
    // Package element to Node Connectivity
    moris::uint tSpatialDim = mBackgroundMesh.get_mesh_data().get_spatial_dim();

    // Children element nodes connected to elements
    moris::Cell<moris::Matrix<moris::IdMat>>  tElementToNodeChildrenByPhase = mCutMesh.get_full_element_to_node_by_phase_glob_ids(mGeometryEngine.get_num_bulk_phase(),mBackgroundMesh.get_mesh_data());

    // Child element ids
    moris::Cell<moris::Matrix<moris::IdMat>>  tChildElementsByPhase = mCutMesh.get_child_elements_by_phase(mGeometryEngine.get_num_bulk_phase(),mBackgroundMesh.get_mesh_data());

    // Parent elements without children
    Cell<moris::Matrix<moris::IdMat>>  tElementNoChildrenIdsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(mGeometryEngine.get_num_bulk_phase());

    // Connectivity of parent elements without children
    Cell<moris::Matrix<moris::IdMat>>  tElementToNodeNoChildrenByPhase = mBackgroundMesh.get_non_intersected_element_to_node_by_phase(mGeometryEngine.get_num_bulk_phase());

    // Node map  of nodes in the phase we are interested in
    moris::Matrix<moris::IndexMat> tOutputtedNodeInds;
    moris::Matrix<moris::IdMat>  tLocalToGlobalNodeMap = this->get_node_map_restricted_to_output_phases(aOutputOptions,tOutputtedNodeInds);

    // All node coordinates
    moris::Matrix<moris::DDRMat> tNodeCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tOutputtedNodeInds);


    // Number of bulk phases
    uint tNumBulkPhases = mGeometryEngine.get_num_phases();

    // Get non-interescted parent elements by phase
    Cell<moris::Matrix<moris::IdMat>> tNoChildElementsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(tNumBulkPhases);

    // combination of the elements by phase (if specified)
    Cell<moris::Matrix<moris::IdMat>> tCombinedElementsByPhase(tNoChildElementsByPhase.size());
    if(!aOutputOptions.mSeparateInterfaceBlock)
    {
        MORIS_ASSERT(mBackgroundMesh.get_XTK_mesh_element_topology() == CellTopology::TET4," Combining the interface block w/ non-interface block only valid on tet background mesh");

        tCombinedElementsByPhase = combine_interface_and_non_interface_blocks(tChildElementsByPhase,tNoChildElementsByPhase);

    }

    // Interface nodes
    Cell<moris::Matrix<moris::IndexMat>> tInterfaceNodes = mBackgroundMesh.get_interface_nodes_glob_ids();

    // Assemble geometry data as field for mesh output
    moris::Cell< moris::Matrix < moris::DDRMat > > tGeometryFieldData = assemble_geometry_data_as_mesh_field(tOutputtedNodeInds);

    // Give the geometry data a name
    moris::Cell<std::string> tGeometryFieldNames = assign_geometry_data_names();

    // Get rank of the geometry data field
    moris::Cell < enum moris::EntityRank > tFieldRanks =  assign_geometry_data_field_ranks();

    // Get the packaged interface side sets from the cut mesh
    moris::Matrix<moris::IdMat> tInterfaceElemIdandSideOrd = mCutMesh.pack_interface_sides();

    // number of phases being output
    moris::uint tNumPhasesOutput = get_num_phases_to_output(aOutputOptions);

    // Set up field data structure for MTK input
    moris::mtk::MtkFieldsInfo tFieldsInfo;
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tGeometryFields(tGeometryFieldData.size());

    for(uint i = 0; i <tGeometryFieldData.size(); i++)
    {
        tGeometryFields(i).set_field_name(tGeometryFieldNames(i));
        tGeometryFields(i).set_field_entity_rank(moris::EntityRank::NODE);
        tGeometryFields(i).add_field_data(&tLocalToGlobalNodeMap, &tGeometryFieldData(i));
        tFieldsInfo.mRealScalarFields.push_back(&tGeometryFields(i));
    }

    // External Fields - real cell fields
    uint tNumExtRealCellScalarFields = aOutputOptions.mRealElementExternalFieldNames.size();
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealCellScalarFields(tNumExtRealCellScalarFields);
    for(uint i = 0; i<tNumExtRealCellScalarFields; i++)
    {
        tExternalRealCellScalarFields(i).set_field_name(aOutputOptions.mRealElementExternalFieldNames(i));
        tExternalRealCellScalarFields(i).set_field_entity_rank(moris::EntityRank::ELEMENT);
        add_field_for_mesh_input(&tExternalRealCellScalarFields(i),tFieldsInfo);
    }

    // External Fields - real vertex fields
    uint tNumExtRealVertexScalarFields = aOutputOptions.mRealNodeExternalFieldNames.size();
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealVertexScalarFields(tNumExtRealVertexScalarFields);
    for(uint i = 0; i<tNumExtRealVertexScalarFields; i++)
    {
        tExternalRealVertexScalarFields(i).set_field_name(aOutputOptions.mRealNodeExternalFieldNames(i));
        tExternalRealVertexScalarFields(i).set_field_entity_rank(moris::EntityRank::NODE);
        add_field_for_mesh_input(&tExternalRealVertexScalarFields(i),tFieldsInfo);
    }


    // sensitivity fields
    moris::Cell<moris::Matrix<DDRMat>> adxdpData;
    moris::Cell<std::string>           adxdpNames;
    moris::Cell<moris::Matrix<DDRMat>> aDesVars;
    moris::Cell<std::string>           aDesVarsName;
    moris::Matrix<moris::DDRMat>       aNumDesVars;
    std::string                        aNumDesVarsName;

    if(aOutputOptions.mPackageDxDpSparsely && mSensitivity)
    {
        this->extract_interface_sensitivity_sparse(tOutputtedNodeInds,adxdpData,adxdpNames,aDesVars,aDesVarsName,aNumDesVars,aNumDesVarsName);
    }

    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tdxdpDataFields(adxdpData.size());
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tDesVarFields(aDesVars.size());
    moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tNumDesVarsField(1);

    if(aOutputOptions.mPackageDxDpSparsely && mSensitivity)
    {
        // place into a field
        for(moris::uint  i = 0; i <tdxdpDataFields.size(); i++)
        {
            tdxdpDataFields(i).set_field_name(adxdpNames(i));
            tdxdpDataFields(i).set_field_entity_rank(moris::EntityRank::NODE);
            tdxdpDataFields(i).add_field_data( &tLocalToGlobalNodeMap, &adxdpData(i));
            add_field_for_mesh_input(&tdxdpDataFields(i),tFieldsInfo);
        }

        for(moris::uint  i = 0; i <tDesVarFields.size(); i++)
        {
            tDesVarFields(i).set_field_name(aDesVarsName(i));
            tDesVarFields(i).set_field_entity_rank(moris::EntityRank::NODE);
            tDesVarFields(i).add_field_data( &tLocalToGlobalNodeMap, &aDesVars(i));
            add_field_for_mesh_input(&tDesVarFields(i),tFieldsInfo);
        }

        tNumDesVarsField(0).set_field_name(aNumDesVarsName);
        tNumDesVarsField(0).set_field_entity_rank(moris::EntityRank::NODE);
        tNumDesVarsField(0).add_field_data( &tLocalToGlobalNodeMap, &aNumDesVars);
        add_field_for_mesh_input(&tNumDesVarsField(0),tFieldsInfo);
    }

    //TODO: implement node owner (currently set to owned by this proc)
//    moris::Matrix<moris::IdMat> tNodeOwner(1,tOutputtedNodeInds.numel(),moris::par_rank());

    moris::Matrix<moris::IdMat> tNodeOwner = mBackgroundMesh.get_vertices_owner(tOutputtedNodeInds);

    moris::Matrix<moris::IdMat> tNodeSharing = mBackgroundMesh.get_vertices_sharing(tOutputtedNodeInds);

    // Set up mesh sets
    // Initialize Sets information structure
    moris::mtk::MtkSetsInfo tMtkMeshSets;

    //
    moris::uint tNumBlocksPerPhase = 2;
    if(!aOutputOptions.mSeparateInterfaceBlock)
    {
        MORIS_ASSERT(mBackgroundMesh.get_XTK_mesh_element_topology() == CellTopology::TET4," Combining the interface block w/ non-interface block only valid on tet background mesh");
        tNumBlocksPerPhase = 1;
    }

    // Setup block sets
    Cell<moris::mtk::MtkBlockSetInfo> tBlockSets(tNumPhasesOutput*tNumBlocksPerPhase);
    uint tCount= 0;

    for(uint i = 0; i <tNumBulkPhases; i++)
    {
        if(aOutputOptions.output_phase(i) && aOutputOptions.mSeparateInterfaceBlock)
        {
            // Children of material phase i
            tBlockSets(tCount).mCellIdsInSet = &tChildElementsByPhase(i);
            tBlockSets(tCount).mBlockSetName = "child_"+std::to_string(i);
            tBlockSets(tCount).mBlockSetTopo = CellTopology::TET4;

            tMtkMeshSets.add_block_set(&tBlockSets(tCount));
            tCount++;

            // Children of material phase i
            tBlockSets(tCount).mCellIdsInSet = &tNoChildElementsByPhase(i);
            tBlockSets(tCount).mBlockSetName = "parent_"+std::to_string(i);
            tBlockSets(tCount).mBlockSetTopo = mBackgroundMesh.get_XTK_mesh_element_topology();

            tMtkMeshSets.add_block_set(&tBlockSets(tCount));
            tCount++;
        }

        else if(aOutputOptions.output_phase(i) && !aOutputOptions.mSeparateInterfaceBlock)
        {
            // Children of material phase i
            tBlockSets(tCount).mCellIdsInSet = &tCombinedElementsByPhase(i);
            tBlockSets(tCount).mBlockSetName = "phase_"+std::to_string(i);
            tBlockSets(tCount).mBlockSetTopo = CellTopology::TET4;

            tMtkMeshSets.add_block_set(&tBlockSets(tCount));
            tCount++;
        }
    }

    // Interface elements
    moris::Matrix<moris::IndexMat> tInterfaceElements(0,0);
    moris::Matrix<moris::IndexMat> tInterfaceElementIds(0,0);
    moris::mtk::MtkBlockSetInfo tUnzippedInterfaceBlockSet;
    if(mUnzipped && aOutputOptions.mHaveInterface)
    {
        // get the interface elements local node index element connectivity
        tInterfaceElements = mCutMesh.get_extracted_interface_elements_loc_inds();

        mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,tInterfaceElements);

        tInterfaceElementIds = mCutMesh.get_interface_element_ids();

        tUnzippedInterfaceBlockSet.mCellIdsInSet = &tInterfaceElementIds;
        tUnzippedInterfaceBlockSet.mBlockSetName = "interface";
        tUnzippedInterfaceBlockSet.mBlockSetTopo = CellTopology::PRISM6;
        tMtkMeshSets.add_block_set(&tUnzippedInterfaceBlockSet);
    }

    // propogate background mesh node sets
    moris::Cell<moris::Matrix<IndexMat>> tBackgroundNodeSetData;
    moris::Cell<moris::mtk::MtkNodeSetInfo> tBackgroundNodeSets;
    if(aOutputOptions.mAddNodeSets)
    {
        tBackgroundNodeSets = propogate_background_node_sets(tBackgroundNodeSetData,aOutputOptions);

        for(moris::uint i = 0; i<tBackgroundNodeSets.size(); i++)
        {
            tMtkMeshSets.add_node_set(&tBackgroundNodeSets(i));
        }
    }


    moris::Cell<moris::mtk::MtkNodeSetInfo> tInterfaceNodeSets(tInterfaceNodes.size());
    if(aOutputOptions.mHaveInterface)
    {


        for(uint i = 0; i<tInterfaceNodes.size(); i++)
        {
            tInterfaceNodeSets(i).mNodeIds     = &tInterfaceNodes(i);
            tInterfaceNodeSets(i).mNodeSetName = "inodes_" +std::to_string(i) ;
            tMtkMeshSets.add_node_set(&tInterfaceNodeSets(i));
        }
    }

    moris::mtk::MtkSideSetInfo tInterfaceSideSet;
    if(aOutputOptions.mHaveInterface)
    {
        tInterfaceSideSet.mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd;
        tInterfaceSideSet.mSideSetName        = "iside_0" ;

        // Add side side set to mesh sets
        tMtkMeshSets.add_side_set(&tInterfaceSideSet);
    }


    Cell<moris::mtk::MtkSideSetInfo> tGhostSideSets;
    Cell<Matrix<IdMat>> tGhostElementsAndSideOrds;
    if(mGhost)
    {
        tGhostSideSets.resize(mGeometryEngine.get_num_bulk_phase());

        tGhostElementsAndSideOrds = this->pack_ghost_as_side_set();

        // set up ghost names
        for(moris::uint i = 0; i <tGhostElementsAndSideOrds.size(); i++)
        {
            tGhostSideSets(i).mElemIdsAndSideOrds = &tGhostElementsAndSideOrds(i);
            tGhostSideSets(i).mSideSetName = "ghost_" + std::to_string(i);
            tMtkMeshSets.add_side_set(&tGhostSideSets(i));
        }

    }

    // propogate side sets from background mesh
    moris::Cell<moris::Matrix<IndexMat>> tSideSetData;
    moris::Cell<moris::mtk::MtkSideSetInfo> tBackgroundSideSets;
    if(aOutputOptions.mAddSideSets)
    {
        // collect information about background side set
        tBackgroundSideSets = this->propogate_background_side_sets(tSideSetData,aOutputOptions);

        // add to mesh input structure
        for(moris::uint i = 0; i<tBackgroundSideSets.size(); i++)
        {
            tMtkMeshSets.add_side_set(&tBackgroundSideSets(i));
        }
    }

    // Mesh data input structure (with multi element type mesh)
    moris::mtk::MtkMeshData tMeshDataInput;

    moris::uint tNumElemTypes = tNumPhasesOutput*2;
    if(mUnzipped)
    {
        tNumElemTypes = tNumElemTypes + 1;
    }

    tMeshDataInput.ElemConn             = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
    tMeshDataInput.LocaltoGlobalElemMap = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
    tMeshDataInput.CellTopology         = moris::Cell<enum CellTopology>(tNumElemTypes,CellTopology::INVALID);

    tMeshDataInput.Verbose                 = mVerbose;
    tMeshDataInput.SpatialDim              = &tSpatialDim;

    tCount = 0;
    for(moris::uint  i = 0 ; i <mGeometryEngine.get_num_bulk_phase(); i++)
    {
        if(aOutputOptions.output_phase(i))
        {
            tMeshDataInput.ElemConn(tCount)             = &tElementToNodeChildrenByPhase(i);
            tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tChildElementsByPhase(i);
            tMeshDataInput.CellTopology(tCount)         = CellTopology::TET4;
            tCount++;
            tMeshDataInput.ElemConn(tCount)             = &tElementToNodeNoChildrenByPhase(i);
            tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tElementNoChildrenIdsByPhase(i);
            tMeshDataInput.CellTopology(tCount)         = mBackgroundMesh.get_XTK_mesh_element_topology();
            tCount++;
        }
    }

    tMeshDataInput.NodeCoords              = &tNodeCoordinates;


    Cell<std::string>                 aAuraChildrenBlockNames;
    Cell<moris::Matrix<moris::IdMat>> aAuraChildrenCellIdsByPhase;
    Cell<std::string>                 aAuraNoChildrenBlockNames;
    Cell<moris::Matrix<moris::IdMat>> aAuraNoChildrenCellIdsByPhase;
    Cell<moris::mtk::MtkBlockSetInfo> tAuraNoChildrenBlocks;
    Cell<moris::mtk::MtkBlockSetInfo> tAuraChildrenBlocks;

    if(par_size()>1)
    {
        tMeshDataInput.NodeProcOwner           = &tNodeOwner;

        //FIXME: THIS CAUSES SOME REAL ISSUES IN STK
//        tMeshDataInput.NodeProcsShared         = &tNodeSharing;

        // construct aura part
        this->package_aura_block( aOutputOptions, aAuraChildrenBlockNames,aAuraChildrenCellIdsByPhase,aAuraNoChildrenBlockNames,aAuraNoChildrenCellIdsByPhase);

        tAuraNoChildrenBlocks.resize(aAuraNoChildrenBlockNames.size());
        tAuraChildrenBlocks.resize(aAuraChildrenBlockNames.size());

        MORIS_ASSERT(aAuraNoChildrenBlockNames.size() == aAuraChildrenBlockNames.size()," There should be an aura block for each material phase on each proc. And there should be one for children and no children cells");

        for(moris::uint  iA = 0; iA<aAuraNoChildrenBlockNames.size(); iA++)
        {

            // no children
            tAuraNoChildrenBlocks(iA).mBlockSetName           =   aAuraNoChildrenBlockNames(iA);
            tAuraNoChildrenBlocks(iA).mCellIdsInSet           = & aAuraNoChildrenCellIdsByPhase(iA);
            tAuraNoChildrenBlocks(iA).mBlockSetTopo           =   mBackgroundMesh.get_XTK_mesh_element_topology();
            tAuraNoChildrenBlocks(iA).mParallelConsistencyReq =   true;
            tMtkMeshSets.add_block_set(&tAuraNoChildrenBlocks(iA));

            // children
            tAuraChildrenBlocks(iA).mBlockSetName           =   aAuraChildrenBlockNames(iA);
            tAuraChildrenBlocks(iA).mCellIdsInSet           = & aAuraChildrenCellIdsByPhase(iA);
            tAuraChildrenBlocks(iA).mBlockSetTopo           =   CellTopology::TET4;
            tAuraChildrenBlocks(iA).mParallelConsistencyReq =   true;
            tMtkMeshSets.add_block_set(&tAuraChildrenBlocks(iA));
        }
    }


    tMeshDataInput.FieldsInfo              = &tFieldsInfo;
    tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
    tMeshDataInput.SetsInfo                = &tMtkMeshSets;
    tMeshDataInput.MarkNoBlockForIO        = false;
    tMeshDataInput.CreateAllEdgesAndFaces  = true;
    tMeshDataInput.AutoAuraOptionInSTK     = false;

    //Add clustering information
    moris::mtk::Cell_Cluster_Input tCellClusterInput;                // Cell clusters
    moris::Cell<Matrix<IdMat>>     tClusterCellIds;                  // Cell cluster Ids
    moris::mtk::Side_Cluster_Input tInterfaceSideClusterInput;       // Side clusters
    moris::Cell<Matrix<IdMat>>     tInterfaceCellIdsandSideOrds;     // side cluster ids and side ordinals
    moris::Cell<Matrix<DDRMat>>    tInterfaceSideClusterParamCoords; // side cluster vertex parametric coordinates

    if(aOutputOptions.mAddClusters)
    {
        // cell clustering
        this->setup_cell_clusters_for_output(tCellClusterInput,aOutputOptions,tClusterCellIds);

        tMeshDataInput.CellClusterInput = &tCellClusterInput;


        setup_interface_side_cluster(tInterfaceSideSet.mSideSetName, tInterfaceSideClusterInput, aOutputOptions, tInterfaceCellIdsandSideOrds, tInterfaceSideClusterParamCoords);

        tMeshDataInput.SideClusterInput = &tInterfaceSideClusterInput;
    }

    // Interface elements
    if(mUnzipped)
    {
        tMeshDataInput.ElemConn(tNumElemTypes-1)             = &tInterfaceElements;
        tMeshDataInput.LocaltoGlobalElemMap(tNumElemTypes-1) = &tInterfaceElementIds;
        tMeshDataInput.CellTopology(tNumElemTypes-1)         = CellTopology::PRISM6;
    }

    // add parallel information
    moris::mtk::Visualization_STK tVizTool;
    if(aOutputOptions.mAddParallelFields && par_size()>1)
    {
        moris::mtk::MtkFieldsInfo* tParFields = tVizTool.setup_parallel_cell_fields_for_declaration();
        tFieldsInfo.combine_fields_info(*tParFields);
    }


    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Mesh data setup completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }


    start = std::clock();

//    tMeshDataInput.print_details();

    // cast background mesh to an interpolation mesh and pass in
    moris::mtk::Integration_Mesh* tMeshData = nullptr;
    if(aOutputOptions.mAddClusters)
    {
        moris::mtk::Interpolation_Mesh* tInterpMesh = dynamic_cast<moris::mtk::Interpolation_Mesh*>(&mBackgroundMesh.get_mesh_data());

        tMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput, tInterpMesh );
    }
    else
    {
        tMeshData = moris::mtk::create_integration_mesh( MeshType::STK, tMeshDataInput);
    }

    if(aOutputOptions.mAddParallelFields && par_size()>1)
    {
        tVizTool.populate_parallel_cell_fields_on_mesh(tMeshData);
    }

    if(moris::par_rank() == 0 && mVerbose)
    {
        std::cout<<"XTK: Write to MTK completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
    }

    return tMeshData;

}

//------------------------------------------------------------------------------

moris::Matrix<moris::IndexMat>
Model::get_node_map_restricted_to_output_phases(Output_Options const & aOutputOptions,
                                                moris::Matrix<moris::IndexMat> & aOutputtedNodeInds)
{
    moris::Matrix<moris::IndexMat> tNodeMap = mBackgroundMesh.get_local_to_global_map(EntityRank::NODE);

    moris_index tMyProcRank = par_rank();

    aOutputtedNodeInds.resize(tNodeMap.n_rows(),tNodeMap.n_cols());
    moris::uint tNumNodes = tNodeMap.numel();
    // if we are returning all phases there is no reason to restrict the map
    if(aOutputOptions.output_all_phases())
    {

        for(moris::uint i = 0; i <tNumNodes; i++)
        {
            aOutputtedNodeInds(i) = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
        }

        return tNodeMap;
    }

    else
    {
        moris::Matrix<moris::IndexMat> tRestrictedNodeMap(tNodeMap.n_rows(),tNodeMap.n_cols());

        moris::uint tCount = 0;
        for(moris::uint i = 0; i <tNumNodes; i++)
        {
            if(output_node(i,aOutputOptions) && mBackgroundMesh.get_vertex_owner(i) == tMyProcRank)
            {

                moris::size_t tPhaseIndex = 0;
                moris_index   tVertexIndex = mBackgroundMesh.get_loc_entity_ind_from_entity_glb_id(tNodeMap(i),EntityRank::NODE);
                mGeometryEngine.get_phase_index(tVertexIndex,tPhaseIndex);
                aOutputtedNodeInds(tCount) = tVertexIndex;
                tRestrictedNodeMap(tCount) = tNodeMap(i);

                tCount++;
            }
        }


        tRestrictedNodeMap.resize(1,tCount);
        aOutputtedNodeInds.resize(1,tCount);

        return tRestrictedNodeMap;
    }

}

//------------------------------------------------------------------------------

moris::Cell<moris::mtk::MtkNodeSetInfo>
Model::propogate_background_node_sets(moris::Cell<moris::Matrix<IndexMat>> & aNodeSetData,
                                      Output_Options const & aOutputOptions)
{
    // access background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // get all node set names in background mesh
    moris::Cell<std::string> tSetNames = tMeshData.get_set_names(EntityRank::NODE);

    // allocate output
    aNodeSetData = moris::Cell<moris::Matrix<IndexMat>>(tSetNames.size());
    moris::Cell<moris::mtk::MtkNodeSetInfo> tNodeSetInfo(tSetNames.size());
    for(moris::uint i = 0; i <tSetNames.size(); i++)
    {
        moris::uint tCount = 0;
        moris::Cell<moris::mtk::Vertex const *> tNodesInSetInds = tMeshData.get_vertices_in_vertex_set_no_aura(tSetNames(i));

        aNodeSetData(i) = moris::Matrix<moris::IndexMat>(tNodesInSetInds.size(),1);

        for(moris::uint iNode =0; iNode<tNodesInSetInds.size(); iNode++)
        {
            moris_index tVertexInd = tNodesInSetInds(iNode)->get_index();

            if(this->output_node(tVertexInd,aOutputOptions))
            {
                aNodeSetData(i)(tCount) = tNodesInSetInds(iNode)->get_index();
                tCount++;
            }
        }


        aNodeSetData(i).resize(tCount,1);

        mBackgroundMesh.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,aNodeSetData(i));

        tNodeSetInfo(i).mNodeIds = &aNodeSetData(i);
        tNodeSetInfo(i).mNodeSetName = tSetNames(i);

    }
    return tNodeSetInfo;

}

//------------------------------------------------------------------------------

moris::Cell<moris::mtk::MtkSideSetInfo>
Model::propogate_background_side_sets(moris::Cell<moris::Matrix<IndexMat>> & aSideSetElemIdsAndOrds,
                                      Output_Options const & aOutputOptions)
{
    // access background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // get all side set names in background mesh
    moris::Cell<std::string> tSetNames = tMeshData.get_set_names(EntityRank::FACE);

    // remove internal side sets which show up with a generated string
    tSetNames = check_for_and_remove_internal_seacas_side_sets(tSetNames);

    // allocate output side sets
    moris::Cell<moris::mtk::MtkSideSetInfo> tSideSets(2*tSetNames.size());
    aSideSetElemIdsAndOrds = moris::Cell<moris::Matrix<IndexMat>>(2*tSetNames.size());


    for(moris::uint i = 0; i <tSetNames.size(); i++)
    {
        moris::uint tNoChildInd = 2*i;
        moris::uint tChildInd = 2*i+1;

        this->propogate_background_side_set(tSetNames(i),tNoChildInd,tChildInd,aSideSetElemIdsAndOrds,tSideSets,aOutputOptions,false);
    }


    return tSideSets;
}

void
Model::propogate_background_side_set( std::string             const &             aSideSetName,
                                      moris::moris_index                          aNoChildIndex,
                                      moris::moris_index                          aChildIndex,
                                      moris::Cell<moris::Matrix<IndexMat>>      & aElementIdsAndSideOrd,
                                      moris::Cell<moris::mtk::MtkSideSetInfo>   & aSideSetData,
                                      Output_Options          const             & aOutputOptions,
                                      bool                                        aOutputIndices)
{

    // access background mesh data
    moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

    // declare matrices used through
    moris::Matrix< moris::IndexMat > tElementsAttachedToFace(1,1);
    moris::Matrix< moris::IdMat >    tChildElemsIdsOnFace;
    moris::Matrix< moris::IndexMat > tChildElemsCMIndOnFace;
    moris::Matrix< moris::IndexMat > tChildElemOnFaceOrdinal;

    moris::uint tElementIndex          = 0;
    moris::uint tPhaseIndex            = 0;
    moris::uint tFaceOrdinal           = 0;
    moris::moris_index tChildMeshIndex = 0;
    moris::moris_id    tElementId      = 0;
    moris::moris_index tMyProcRank     = par_rank();
    bool        tHasChildren = false;


    // get cells and sides in side set
    moris::Cell< mtk::Cell const * > tCellsInSideSet;
    moris::Matrix< moris::IndexMat > tSideSetOrdinals;
    tMeshData.get_sideset_cells_and_ords(aSideSetName, tCellsInSideSet, tSideSetOrdinals );

    // side set data non-intersected
    aElementIdsAndSideOrd(aNoChildIndex)   = Matrix<IndexMat>(tCellsInSideSet.size()*2,2);

    // intersected data
    //TODO: FIGURE OUT MAXIMUM VALUE
    aElementIdsAndSideOrd(aChildIndex) = Matrix<IndexMat>(tCellsInSideSet.size()*2*10,2);

    // keep count
    moris::Cell<moris::uint> tCount(2,0);

    // iterate through sides in set i
    for(moris::uint iSide= 0; iSide<tCellsInSideSet.size(); iSide++)
    {
        tElementIndex = tCellsInSideSet(iSide)->get_index();

        // sides attached to cell
        moris::Matrix<moris::IdMat> tElementFaces = tMeshData.get_entity_connected_to_entity_loc_inds(tElementIndex,EntityRank::ELEMENT,EntityRank::FACE);

        moris_index tSideIndex = tElementFaces(tSideSetOrdinals(iSide));

        if(tMeshData.get_entity_owner(tElementIndex, EntityRank::ELEMENT) == tMyProcRank)
        {
            tHasChildren = mBackgroundMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT);
            // get the faces from the child mesh
            if(tHasChildren)
            {
                tChildMeshIndex = mBackgroundMesh.child_mesh_index(tElementIndex,EntityRank::ELEMENT);

                Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                tChildMesh.get_child_elements_connected_to_parent_face(tSideIndex,
                                                                       tChildElemsIdsOnFace,
                                                                       tChildElemsCMIndOnFace,
                                                                       tChildElemOnFaceOrdinal);

                moris::Matrix< moris::IndexMat > const & tChildElementPhaseIndices = tChildMesh.get_element_phase_indices();
                moris::Matrix< moris::IndexMat > const & tChildElementIndices = tChildMesh.get_element_inds();
                moris::Matrix< moris::IndexMat > const & tElementIds = tChildMesh.get_element_ids();

                for(moris::moris_index iCElem  = 0; iCElem < (moris::moris_index)tChildElemsCMIndOnFace.numel(); iCElem++)
                {
                    tPhaseIndex = tChildElementPhaseIndices(0,tChildElemsCMIndOnFace(0,iCElem));
                    if(aOutputOptions.output_phase(tPhaseIndex))
                    {
                        // Child Element Id
                        if(!aOutputIndices)
                        {
                            tElementId = tElementIds(tChildElemsCMIndOnFace(iCElem));
                            tFaceOrdinal   = tChildElemOnFaceOrdinal(iCElem);
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),0) = tElementId;
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),1) = tFaceOrdinal;
                            tCount(1)++;
                        }
                        else
                        {
                            tFaceOrdinal   = tChildElemOnFaceOrdinal(iCElem);
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),0) = tChildElementIndices(tChildElemsCMIndOnFace(iCElem));
                            aElementIdsAndSideOrd(aChildIndex)(tCount(1),1) = tFaceOrdinal;
                            tCount(1)++;
                        }
                    }
                }
            }

            else
            {

                tFaceOrdinal = tSideSetOrdinals(iSide);
                tElementId   = tCellsInSideSet(iSide)->get_id();

                if(aOutputOptions.output_phase(mBackgroundMesh.get_element_phase_index(tElementIndex)))
                {
                    if(!aOutputIndices)
                    {
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),0) = tElementId;
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),1) = tFaceOrdinal;
                        tCount(0)++;
                    }
                    else
                    {
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),0) = tElementIndex;
                        aElementIdsAndSideOrd(aNoChildIndex)(tCount(0),1) = tFaceOrdinal;
                        tCount(0)++;
                    }
                }

            }
        }

    }

    // resize data
    aElementIdsAndSideOrd(aChildIndex).resize(tCount(1),2);
    aElementIdsAndSideOrd(aNoChildIndex).resize(tCount(0),2);

    // Add data to side set info
    // no child
    aSideSetData(aNoChildIndex).mElemIdsAndSideOrds = &aElementIdsAndSideOrd(aNoChildIndex);
    aSideSetData(aNoChildIndex).mSideSetName        = aSideSetName;
    aSideSetData(aNoChildIndex).mSideTopology       = CellTopology::QUAD4;
    aSideSetData(aChildIndex).mElemIdsAndSideOrds   = &aElementIdsAndSideOrd(aChildIndex);
    aSideSetData(aChildIndex).mSideSetName          = aSideSetName + "_i";
    aSideSetData(aChildIndex).mSideTopology         = CellTopology::TRI3;


}

moris::Cell<std::string>
Model::check_for_and_remove_internal_seacas_side_sets(moris::Cell<std::string> & aSideSetNames)
{
    for(std::vector<std::string>::iterator iSet = aSideSetNames.begin(); iSet != aSideSetNames.end(); ++iSet)
    {
        if(iSet->compare("surface_1_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }

        else if(iSet->compare("surface_2_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_3_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_4_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_5_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_6_quad4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad_1") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad_2") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_1") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_2") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_3") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
        else if(iSet->compare("surface_hex8_quad4_4") == 0)
        {
            aSideSetNames.data().erase(iSet--);
        }
    }
    return aSideSetNames;

}

//------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IdMat>>
Model::combine_interface_and_non_interface_blocks(Cell<moris::Matrix<moris::IdMat>> & aChildElementsByPhase,
                                                  Cell<moris::Matrix<moris::IdMat>> & aNoChildElementsByPhase)
{
    moris::uint tNumPhase = aChildElementsByPhase.size();

    Cell<moris::Matrix<moris::IdMat>> tCombinedElementsByPhase(tNumPhase);

    for(moris::uint i =0; i<tNumPhase; i++)
    {
        moris::uint tNumChildElems   = aChildElementsByPhase(i).numel();
        moris::uint tNumNoChildElems = aNoChildElementsByPhase(i).numel();

        tCombinedElementsByPhase(i) = moris::Matrix<moris::IdMat>(1,tNumChildElems + tNumNoChildElems);

        tCombinedElementsByPhase(i)({0,0},{0,tNumChildElems-1}) = aChildElementsByPhase(i).get_row(0);
        tCombinedElementsByPhase(i)({0,0},{tNumChildElems,tNumChildElems + tNumNoChildElems -1}) = aNoChildElementsByPhase(i).get_row(0);
    }

    return tCombinedElementsByPhase;

}

//------------------------------------------------------------------------------

Cell<Matrix<IdMat>>
Model::pack_ghost_as_side_set()
{
    MORIS_ASSERT(mGhost,"Trying to pack ghost side set on a mesh without ghost setup");

    // number of bulk phases
    moris::uint tNumBulkPhase = mGeometryEngine.get_num_bulk_phase();

    Cell<Matrix<IdMat>> tGhostCellIdsAndSideOrdByBulkPhase(tNumBulkPhase);

    // iterate through bulk phases
    for(moris::uint i = 0 ; i<tNumBulkPhase; i++)
    {
        // get the bulk phase ghost elements
        Cell<Ghost_Cell> const & tBulkPhaseGhostCells = mGhostStabilization->get_ghost_cells_by_bulk_phase(i);

        tGhostCellIdsAndSideOrdByBulkPhase(i).resize(tBulkPhaseGhostCells.size()*2,2);

        // iterate through cells and extract the side set information
        uint tCount = 0;
        for(moris::uint j = 0; j <tBulkPhaseGhostCells.size(); j++)
        {
            // Ghost cell reference
            Ghost_Cell const & tGhostCell = tBulkPhaseGhostCells(j);

            // figure out left cell topology

            // left cell element id and side ordinal

            tGhostCellIdsAndSideOrdByBulkPhase(i)(tCount,0) = tGhostCell.get_left_cell().get_id();
            tGhostCellIdsAndSideOrdByBulkPhase(i)(tCount,1) = tGhostCell.get_left_cell_side_ordinal();
            tCount++;

            // right cell element id and side ordinal
            tGhostCellIdsAndSideOrdByBulkPhase(i)(tCount,0) = tGhostCell.get_right_cell().get_id();
            tGhostCellIdsAndSideOrdByBulkPhase(i)(tCount,1) = tGhostCell.get_right_cell_side_ordinal();
            tCount++;
        }

    }

    return tGhostCellIdsAndSideOrdByBulkPhase;
}

void
Model::package_aura_block(Output_Options              const & aOutputOptions,
                          Cell<std::string>                 & aAuraChildrenBlockNames,
                          Cell<moris::Matrix<moris::IdMat>> & aAuraChildrenCellIdsByPhase,
                          Cell<std::string>                 & aAuraNoChildrenBlockNames,
                          Cell<moris::Matrix<moris::IdMat>> & aAuraNoChildrenCellIdsByPhase)
{
    // setup the aura names
    this->setup_aura_block_names(aOutputOptions, aAuraChildrenBlockNames, aAuraNoChildrenBlockNames);

    // setup the blocks
    this->setup_aura_cells_into_blocks(aOutputOptions, aAuraChildrenBlockNames, aAuraChildrenCellIdsByPhase, aAuraNoChildrenBlockNames, aAuraNoChildrenCellIdsByPhase);

}

//------------------------------------------------------------------------------

void
Model::setup_aura_block_names(Output_Options              const & aOutputOptions,
                              Cell<std::string>                 & aAuraChildrenBlockNames,
                              Cell<std::string>                 & aAuraNoChildrenBlockNames)
{
    // my process rank
    moris::moris_index tParSize    = par_size();

    // number of bulk phases
    moris::uint tNumBulkPhase = mGeometryEngine.get_num_bulk_phase();

    // matrix which keeps track of the processor combo I need to declare a block set for
    aAuraNoChildrenBlockNames = Cell<std::string>(tParSize*tNumBulkPhase);
    aAuraChildrenBlockNames   = Cell<std::string>(tParSize*tNumBulkPhase);

    // Setup aura labels
    uint tCount = 0;
    for(moris::moris_index i = 0; i <tParSize; i++)
    {
        std::string tBaseStringNoChild = "aura_nc_pr_"+std::to_string(i);
        std::string tBaseStringChild   = "aura_c_pr_"+std::to_string(i);

        for(moris::uint iBulk = 0; iBulk < tNumBulkPhase; iBulk++)
        {
            aAuraNoChildrenBlockNames(tCount) = tBaseStringNoChild + "_p_" + std::to_string(iBulk);
            aAuraChildrenBlockNames(tCount)   = tBaseStringChild + "_p_" + std::to_string(iBulk);
            tCount++;
        }
    }
}

//------------------------------------------------------------------------------


void
Model::setup_aura_cells_into_blocks(Output_Options              const & aOutputOptions,
                                    Cell<std::string>                 & aAuraChildrenBlockNames,
                                    Cell<moris::Matrix<moris::IdMat>> & aAuraChildrenCellIdsByPhase,
                                    Cell<std::string>                 & aAuraNoChildrenBlockNames,
                                    Cell<moris::Matrix<moris::IdMat>> & aAuraNoChildrenCellIdsByPhase)
{
    // access background mesh data
     moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

     // my process rank
     moris::moris_index tParSize    = par_size();
     moris::moris_index tParRank    = par_rank();

     // number of bulk phases
     moris::uint tNumBulkPhase = mGeometryEngine.get_num_bulk_phase();

     // size block cells
     aAuraNoChildrenCellIdsByPhase =  Cell<moris::Matrix<moris::IdMat>>(tParSize*tNumBulkPhase);
     aAuraChildrenCellIdsByPhase   =  Cell<moris::Matrix<moris::IdMat>>(tParSize*tNumBulkPhase);

     Cell<Cell<Matrix< IdMat >>> tChildElementIdsByPhase(tParSize*tNumBulkPhase);
     Cell<Cell<moris_id>> tNoChildElementIdsByPhase(tParSize*tNumBulkPhase);

     // place entities into the cell of cells above, this will then be concatenated into a single matrix
     for(moris::uint iEl = 0; iEl < tMeshData.get_num_entities(EntityRank::ELEMENT); iEl++)
     {
         moris_index tOwningProc = tMeshData.get_entity_owner((moris_index)iEl,EntityRank::ELEMENT);

         Matrix<moris::IndexMat> tSharingProcessors;
         tMeshData.get_processors_whom_share_entity((moris_index)iEl,EntityRank::ELEMENT,tSharingProcessors);


         // if we have children, add to count of children block
         if(mBackgroundMesh.entity_has_children((moris_index)iEl,EntityRank::ELEMENT))
         {
             // child mesh index
             moris::moris_index tChildMeshIndex = mBackgroundMesh.child_mesh_index((moris_index)iEl, EntityRank::ELEMENT);

             // Get child mesh
             Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

             // pack the element ids into phase grouping
             Cell<Matrix< IdMat >> tSingleCMElementIds;
             Cell<Matrix< IdMat >> tCMElementInds;
             tChildMesh.pack_child_mesh_by_phase(mGeometryEngine.get_num_bulk_phase(), tSingleCMElementIds, tCMElementInds);

             bool tAddToMyAura = false;

             // iterate through phases
             for(moris::uint iP = 0; iP<tSingleCMElementIds.size(); iP++)
             {
                 // iterate through shared processors
                 for(moris::uint iShare = 0; iShare<tSharingProcessors.numel(); iShare++)
                 {
                     // if we share this entity with another process, mark it to also be added to my aura
                     if(tSharingProcessors(iShare) != tParRank)
                     {
                         tAddToMyAura = true;

                         // if the other proc owns the parent element then it does not belong in their aura
                         if(tSharingProcessors(iShare) != tOwningProc)
                         {
                         // get the index in the cells with the sharing processor
                         moris_index tIndex = get_aura_block_index(aOutputOptions,iP,tSharingProcessors(iShare),tNumBulkPhase);

                         // add ids to the correct process aura block
                         tChildElementIdsByPhase(tIndex).push_back(tSingleCMElementIds(iP));
                         }

                     }
                 }
             }

             // add also to my aura
             if(tAddToMyAura && tOwningProc != tParRank)
             {
                 for(moris::uint iP = 0; iP<tSingleCMElementIds.size(); iP++)
                 {

                     // get the index in the cells with the sharing processor
                     moris_index tIndex = get_aura_block_index(aOutputOptions,iP,tParRank,tNumBulkPhase);

                     // add ids to the correct process aura block
                     tChildElementIdsByPhase(tIndex).push_back(tSingleCMElementIds(iP));

                 }
             }
         }

         // entity does not have children
         else
         {
             // get phase
             moris_index tPhaseIndex = mBackgroundMesh.get_element_phase_index((moris_index)iEl);

             bool tAddToMyAura = false;

             // iterate through shared processors
             for(moris::uint iShare = 0; iShare<tSharingProcessors.numel(); iShare++)
             {
                 // if we share this entity with another process, mark it to also be added to my aura
                 if(tSharingProcessors(iShare) != tParRank)
                 {
                     tAddToMyAura = true;

                     // if the other proc owns the parent element then it does not belong in their aura
                     if(tSharingProcessors(iShare) != tOwningProc)
                     {
                         // get the index in the cells with the sharing processor
                         moris_index tIndex = this->get_aura_block_index(aOutputOptions,tPhaseIndex,tSharingProcessors(iShare),tNumBulkPhase);

                         // add ids to the correct process aura block
                         tNoChildElementIdsByPhase(tIndex).push_back(tMeshData.get_glb_entity_id_from_entity_loc_index((moris_index)iEl,EntityRank::ELEMENT));
                     }
                 }
             }

             if(tAddToMyAura && tOwningProc != tParRank)
             {
                 // get the index in the cells with the sharing processor
                 moris_index tIndex = this->get_aura_block_index(aOutputOptions,tPhaseIndex,tParRank,tNumBulkPhase);

                 // add ids to the correct process aura block
                 tNoChildElementIdsByPhase(tIndex).push_back(tMeshData.get_glb_entity_id_from_entity_loc_index((moris_index)iEl,EntityRank::ELEMENT));
             }

         }
     }

     // take the cell of cell of matrices and turn the interior cell into a single matrix
     // for children cells
     for(moris::uint i = 0; i <tChildElementIdsByPhase.size(); i++)
     {
         // count total size
         moris::uint tSize = 0;

         // iterate through matrices in interior cell and tally the size
         for(moris::uint j = 0; j <tChildElementIdsByPhase(i).size(); j++)
         {
             tSize = tSize + tChildElementIdsByPhase(i)(j).numel();
         }

         // allocate a matrix
         aAuraChildrenCellIdsByPhase(i) = Matrix<IdMat>(1,tSize);

         // place data into matrix
         uint tStart = 0;
         for(moris::uint j = 0; j <tChildElementIdsByPhase(i).size(); j++)
         {

             if(tChildElementIdsByPhase(i)(j).numel() != 0)
             {
                 uint tEnd = tStart + tChildElementIdsByPhase(i)(j).numel() - 1;

                 aAuraChildrenCellIdsByPhase(i)({0,0},{tStart,tEnd}) = tChildElementIdsByPhase(i)(j).get_row(0);

                 tStart = tEnd+1;
             }
         }
     }

     // take the cells of cell of ids and turn the interior cell into a single matrix
     for(moris::uint i = 0; i <tNoChildElementIdsByPhase.size(); i++)
     {
         // count total size
         moris::uint tSize = tNoChildElementIdsByPhase(i).size();

         // allocate a matrix
         aAuraNoChildrenCellIdsByPhase(i) = Matrix<IdMat>(1,tSize);

         // place data into matrix
         uint tLoc = 0;
         for(moris::uint j = 0; j <tSize; j++)
         {

             aAuraNoChildrenCellIdsByPhase(i)(tLoc) = tNoChildElementIdsByPhase(i)(j);

             tLoc++;
         }
     }

}

moris_index
Model::get_aura_block_index(Output_Options const & aOutputOptions,
                            moris_index            aBulkPhase,
                            moris_index            aProcessorRank,
                            moris_index            aNumBulkPhases)
{
    MORIS_ASSERT(aOutputOptions.output_all_phases() == true, "This function has not been setup to figure out index of output which is not outputting all material phases");

    moris_index tIndex = aProcessorRank*aNumBulkPhases+aBulkPhase;

    return tIndex;
}

//------------------------------------------------------------------------------

uint
Model::get_num_phases_to_output(Output_Options const & aOutputOptions)
{
    uint tNumPhasesOutput = 0;
    if(aOutputOptions.output_all_phases())
    {
        tNumPhasesOutput = mGeometryEngine.get_num_bulk_phase();
    }
    else
    {
        tNumPhasesOutput = aOutputOptions.num_phases_to_output();
    }

    return tNumPhasesOutput;
}

//------------------------------------------------------------------------------

void
Model::setup_cell_clusters_for_output(moris::mtk::Cell_Cluster_Input & aCellClusterInput,
                                      Output_Options const & aOutputOptions,
                                      moris::Cell<Matrix<IdMat>> & aCellIds)
{
    // iterate through child meshes and construct cells
    uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    for(moris::uint i = 0; i < tNumChildMeshes; i ++)
    {
        // Get child mesh
        Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

        // pack the element ids into phase grouping
        Cell<moris::Matrix< moris::IdMat >> tElementIds;
        Cell<moris::Matrix< moris::IdMat >> tCMElementInds;
        tChildMesh.pack_child_mesh_by_phase(mGeometryEngine.get_num_bulk_phase(), tElementIds, tCMElementInds);

        // add them to cell to keep in scope
        aCellIds.push_back(tElementIds(0));
        aCellIds.push_back(tElementIds(1));
    }

    for(moris::uint i = 0; i < tNumChildMeshes; i ++)
    {
        // Get child mesh
        Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(i);

        // primary index
        moris_index tPrimaryCellIndex = 2*i;
        moris_index tVoidCellIndex    = 2*i+1;

        // parent index
        moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

        // access the parent element from the background mesh
        moris::mtk::Cell* tInterpCell = &mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex);

        // add to cluster
        aCellClusterInput.add_cluster_data(tInterpCell,&aCellIds(tPrimaryCellIndex),&aCellIds(tVoidCellIndex),&tChildMesh.get_node_ids(),&tChildMesh.get_parametric_coordinates());

    }
}

void
Model::setup_interface_side_cluster(std::string                      aInterfaceSideLabelBase,
                                    moris::mtk::Side_Cluster_Input & aSideClusterInput,
                                    Output_Options           const & aOutputOptions,
                                    moris::Cell<Matrix<IdMat>>     & aCellIdsandSideOrds,
                                    moris::Cell<Matrix<DDRMat>>    & aParametricCoordinates)
{
    moris::uint tNumPhases = mGeometryEngine.get_num_bulk_phase();

    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    for(moris::uint  iP = 0; iP<tNumPhases; iP++)
    {
        // if we are outputting this phase
//        if(aOutputOptions.output_phase((size_t)iP))
        if(iP == 0)
        {
            // add side set to output
            std::string tSetName = aInterfaceSideLabelBase;

            //iterate through children meshes
            for(moris::uint iC = 0; iC < tNumChildMeshes; iC ++)
            {
                // Get child mesh
                Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(iC);

                // package this child element by bulk phase
                moris::Matrix< moris::IdMat > tInterfaceElementIdsAndSideOrd = tChildMesh.pack_interface_sides( false, iP );

                // add to data which will stay in scope
                aCellIdsandSideOrds.push_back(tInterfaceElementIdsAndSideOrd);

            }
        }
    }

    uint tCount = 0;
    for(moris::uint  iP = 0; iP<tNumPhases; iP++)
    {
        // if we are outputting this phase
//        if(aOutputOptions.output_phase((size_t)iP))
        if(iP == 0)
        {
            // add side set to output
            std::string tSetName = aInterfaceSideLabelBase;
           moris_index tSideSetOrd = aSideClusterInput.add_side_set_label(tSetName);

            //iterate through children meshes
            for(moris::uint iC = 0; iC < tNumChildMeshes; iC ++)
            {
                // Get child mesh
                Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(iC);

                // parent cell index
                moris_index tParentCellIndex = tChildMesh.get_parent_element_index();

                // access the parent element from the background mesh
                moris::mtk::Cell* tInterpCell = &mBackgroundMesh.get_mesh_data().get_mtk_cell(tParentCellIndex);

                // add to cluster input data
                //fixme: Add only vertex indices on the interface to cluster. Adding all.
                aSideClusterInput.add_cluster_data(false,tSideSetOrd,tInterpCell,&aCellIdsandSideOrds(tCount),&tChildMesh.get_node_ids(),&tChildMesh.get_parametric_coordinates());

                tCount++;
            }
        }
    }
}
//------------------------------------------------------------------------------


bool
Model::output_node(moris::moris_index aNodeIndex,
                   Output_Options const & aOutputOptions)
{
    bool tIsInterface = mBackgroundMesh.is_interface_node(aNodeIndex,0);
    moris::size_t tPhaseIndex = 0;
    mGeometryEngine.get_phase_index(aNodeIndex,tPhaseIndex);


    if(aOutputOptions.output_phase(tPhaseIndex) && !tIsInterface)
    {
        return true;
    }
    else if(tIsInterface)
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
Model::compute_interface_sensitivity_internal()
{
    // Number of geometries in the geometry engine (we need to compute sensitivity wrt each)
    uint tNumGeoms = mGeometryEngine.get_num_geometries();

    // Node coordinates
    moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

    for(uint iGeo = 0; iGeo <tNumGeoms; iGeo++)
    {
        // Get interface nodes
        moris::Matrix< moris::IndexMat > tInterfaceNodes = mBackgroundMesh.get_interface_nodes_loc_inds(iGeo);

        // Compute interface sensitivity
        mGeometryEngine.compute_interface_sensitivity(tInterfaceNodes,tNodeCoords,iGeo);
    }
}


void
Model::extract_interface_sensitivity_sparse(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput,
                                            moris::Cell<moris::Matrix<DDRMat>>   & adxdpData,
                                            moris::Cell<std::string>             & adxdpNames,
                                            moris::Cell<moris::Matrix<DDRMat>>   & aDesVars,
                                            moris::Cell<std::string>             & aDesVarsName,
                                            moris::Matrix<moris::DDRMat>         & aNumDesVars,
                                            std::string                          & aNumDesVarsName) const
{
    // names of sparsely packaged fields
    moris::uint tNumFields = 6;
    adxdpNames = moris::Cell<std::string>({{"dx0dp0"},
        {"dx1dp0"},
        {"dx2dp0"},
        {"dx0dp1"},
        {"dx1dp1"},
        {"dx2dp1"}});



    moris::uint tNumNodes = aNodeIndsToOutput.numel();
    adxdpData = moris::Cell<moris::Matrix<moris::DDRMat>>(tNumFields,moris::Matrix<moris::DDRMat>(tNumNodes,1,0.0));

    //TODO: hardcoded to 2
    tNumFields = 2;
    aDesVarsName = moris::Cell<std::string>({{"DesVar0"},{"DesVar1"}});
    aDesVars = moris::Cell<moris::Matrix<moris::DDRMat>>(tNumFields,moris::Matrix<moris::DDRMat>(tNumNodes,1));

    aNumDesVarsName = "NumDesVar";
    aNumDesVars = moris::Matrix<moris::DDRMat>(1,tNumNodes,0);

    for(moris::uint iNode = 0; iNode<tNumNodes; iNode++)
    {
        moris::moris_index tNodeIndex = aNodeIndsToOutput(iNode);

        if(mBackgroundMesh.is_interface_node(tNodeIndex,0))
        {
            Geometry_Object const & tNodeGeoObj = mGeometryEngine.get_geometry_object(tNodeIndex);

            moris::Matrix< moris::DDRMat > const & tdxdp = tNodeGeoObj.get_sensitivity_dx_dp();

            MORIS_ASSERT(tdxdp.n_rows() == 2,"Invalid dxdp size for sparse packing, This function only works on tet meshes with discrete fields at the moment");
            MORIS_ASSERT(tdxdp.n_cols() == 3,"Invalid dxdp size for sparse packing, This function only works on tet meshes with discrete fields at the moment");

            adxdpData(0)(iNode) = tdxdp(0,0);
            adxdpData(1)(iNode) = tdxdp(0,1);
            adxdpData(2)(iNode) = tdxdp(0,2);
            adxdpData(3)(iNode) = tdxdp(1,0);
            adxdpData(4)(iNode) = tdxdp(1,1);
            adxdpData(5)(iNode) = tdxdp(1,2);

            moris::Matrix< moris::IndexMat > const & tDesVarInds = tNodeGeoObj.get_node_adv_indices();

            aDesVars(0)(iNode) = (moris::real)mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tDesVarInds(0),EntityRank::NODE);
            aDesVars(1)(iNode) = (moris::real)mBackgroundMesh.get_mesh_data().get_glb_entity_id_from_entity_loc_index(tDesVarInds(1),EntityRank::NODE);

            aNumDesVars(iNode) = tDesVarInds.numel();
        }

        else
        {

        }
    }
}

//------------------------------------------------------------------------------

moris::Cell< moris::Matrix < moris::DDRMat > >
Model::assemble_geometry_data_as_mesh_field(moris::Matrix<moris::IndexMat> const & aNodeIndsToOutput)
{
    uint tNumGeometries = mGeometryEngine.get_num_geometries();
    uint tNumNodes      = aNodeIndsToOutput.numel();


    // Allocate output data
    moris::Cell< moris::Matrix < moris::DDRMat > > tGeometryData(tNumGeometries, moris::Matrix<moris::DDRMat>(tNumNodes,1));

    //Iterate through geometries
    for(uint iG = 0; iG <tNumGeometries; iG++)
    {
        // Iterate through nodes
        for(uint iN = 0; iN<tNumNodes; iN++)
        {
            tGeometryData(iG)(iN) = mGeometryEngine.get_entity_phase_val(aNodeIndsToOutput(iN),iG);
        }
    }

    return tGeometryData;
}

//------------------------------------------------------------------------------

moris::Cell<std::string>
Model::assign_geometry_data_names()
{
    uint tNumGeometries = mGeometryEngine.get_num_geometries();

    // base string of geometry data
    std::string tBaseName = "gd_";

    // Allocate output
    moris::Cell<std::string> tGeometryFieldName(tNumGeometries);

    //Iterate through geometries
    for(uint iG = 0; iG <tNumGeometries; iG++)
    {
        tGeometryFieldName(iG) = tBaseName+std::to_string(iG);
    }

    return tGeometryFieldName;
}

//------------------------------------------------------------------------------

moris::Cell < enum moris::EntityRank >
Model::assign_geometry_data_field_ranks()
{
    uint tNumGeometries = mGeometryEngine.get_num_geometries();

    // base string of geometry data
    std::string tBaseName = "gd_";

    // Allocate output
    // Note: for now this is a nodal field always
    moris::Cell<enum moris::EntityRank> tGeometryFieldRank(tNumGeometries,moris::EntityRank::NODE);

    return tGeometryFieldRank;
}


}
