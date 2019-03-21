/*
 * cl_XTK_Model.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
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
}

Model::Model(uint aModelDimension,
             moris::mtk::Mesh* aMeshData,
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
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeCoordinate);

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

                        // Intersected edge is an existing stk edge
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


            } // XTK Mesh loop

            moris_index tMessageTag = 60001; /*arbitrary tag for regular subdivision*/
            assign_node_requests_identifiers(tDecompData,tMessageTag);

            // Allocate interface flag space in XTK mesh even though these are not interface nodes
            mBackgroundMesh.allocate_space_in_interface_node_flags(tDecompData.tNewNodeIndex.size(),mGeometryEngine.get_num_geometries());

            // add nodes to child mesh
            this->decompose_internal_set_new_nodes_in_child_mesh_nh(aActiveChildMeshIndices,tDecompData);

            // add nodes to the background mesh
            mBackgroundMesh.batch_create_new_nodes(tDecompData.tNewNodeId,tDecompData.tNewNodeIndex,tDecompData.tNewNodeCoordinate);

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

                // coordinates of nodes attached to the nodes of this face
                moris::Matrix<moris::DDRMat> tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);

                // bilinearly interpolate to the center of this face fi
                moris::Matrix<moris::DDRMat> tNewNodeCoordinates;
                xtk::Interpolation::bilinear_interpolation(tCoordinates, tParamCoordsRelativeToFace.get_row(fi),tNewNodeCoordinates);

                // location for this face in the map
                moris_index tNewNodeIndexInSubdivision = tDecompData.register_new_request(tFaceIndices(fi),
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

        MORIS_ASSERT(!tDecompData.request_exists(tElemInd,EntityRank::ELEMENT,tNewNodeIndexInSubdivision),"All element requests should be unique, therefore tNewRequest is expected to be true here");


        tNewNodeIndexInSubdivision = tDecompData.register_new_request(tElemInd,
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
Model::finalize_decomp_in_xtk_mesh(bool aSetPhase)
{

    // Set child element ids and indices
    moris::size_t tNumElementsInCutMesh = mCutMesh.get_num_entities(EntityRank::ELEMENT);

    // Allocate global element ids (these need to be give to the children meshes)
    moris_id    tElementIdOffset = mBackgroundMesh.allocate_entity_ids(tNumElementsInCutMesh, moris::EntityRank::ELEMENT);
    moris_index tElementIndOffset = mBackgroundMesh.get_first_available_index(EntityRank::ELEMENT);
    for(moris::size_t i = 0; i<mCutMesh.get_num_child_meshes(); i++)
    {
        mCutMesh.set_child_element_ids(i,tElementIdOffset);
        mCutMesh.set_child_element_inds(i,tElementIndOffset);
    }

    mBackgroundMesh.update_first_available_index(tElementIndOffset,EntityRank::ELEMENT);

    // Associate nodes created during decomposition to their child meshes
    associate_nodes_created_during_decomp_to_child_meshes();

    // Compute the child element phase using the geometry engine
    if(aSetPhase)
    {
        this->set_element_phases(tElementIndOffset);
    }

    // creates mtk cells for all child elements (parent elements are assumed to have mtk cells in the mtk mesh)
    create_child_element_mtk_cells();


    // Change XTK model decomposition state flag
    mDecomposed = true;
}


void
Model::create_child_element_mtk_cells()
{

    moris::uint tNumChildMeshes = mCutMesh.get_num_child_meshes();

    for(moris::uint i=0; i<tNumChildMeshes; i++)
    {
        Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(i);

        moris::Matrix< moris::IdMat >    const & tElementIds  = tChildMesh.get_element_ids();
        moris::Matrix< moris::IndexMat > const & tElementInds = tChildMesh.get_element_inds();

        // Iterate over elements
        for(moris::uint j = 0; j<tElementIds.numel(); j++)
        {
            mBackgroundMesh.add_child_element_to_mtk_cells(tElementInds(j),tElementIds(j),j, &tChildMesh);
        }

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
// Unzipping Source code
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

        // construct an interface element
        Interface_Element tInterfaceElement;
        tInterfaceElement.set_element_pair_and_side_ordinal(tPairCells,aSideOrdinalPairs.get_column(i));

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
        MORIS_ERROR(mDecomposed,"Prior to computing basis enrichment, the decomposition process must be called");
        MORIS_ERROR(mUnzipped,"Prior to computing basis enrichment, the interface unzipping process must be called");
        MORIS_ERROR(!mEnriched,"Calling perform_basis_enrichment twice is not supported");
        perform_basis_enrichment_internal();

        // Change the enrichment flag
        mEnriched = true;
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
        mEnrichment = new Enrichment(mGeometryEngine.get_num_phases(),
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
    Model::perform_multilevel_enrichment_internal()
    {
        mEnrichment->create_multilevel_enrichments();
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

    moris::mtk::Mesh*
    Model::get_output_mesh(Output_Options const & aOutputOptions)

    {
        // start timing on this decomposition
        std::clock_t start = std::clock();

        // create the output mtk mesh
        moris::mtk::Mesh* tOutputMesh = construct_output_mesh(aOutputOptions);

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

    moris::mtk::Mesh*
    Model::construct_output_mesh( Output_Options const & aOutputOptions )
        {

            // start timing on this decomposition
            std::clock_t start = std::clock();

            // Get mesh information ready for outputting
            // Package element to Node Connectivity
            moris::uint tSpatialDim = mBackgroundMesh.get_mesh_data().get_spatial_dim();

            // Children element nodes connected to elements
            moris::Cell<moris::Matrix<moris::IdMat>>  tElementToNodeChildrenByPhase = mCutMesh.get_full_element_to_node_by_phase_glob_ids(mGeometryEngine.get_num_bulk_phase());

            // Child element ids
            moris::Cell<moris::Matrix<moris::IdMat>>  tElementChildrenIdsByPhase = mCutMesh.get_child_elements_by_phase(mGeometryEngine.get_num_bulk_phase());

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

            // Get child elements sorted by phase
            Cell<moris::Matrix<moris::IdMat>> tChildElementsByPhase = mCutMesh.get_child_elements_by_phase(tNumBulkPhases);

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
            moris::uint tNumPhasesOutput = 0;
            if(aOutputOptions.output_all_phases())
            {
                tNumPhasesOutput = mGeometryEngine.get_num_bulk_phase();
            }
            else
            {
                tNumPhasesOutput = aOutputOptions.num_phases_to_output();
            }

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

            // External Fields
            uint tNumExtRealScalarFields = aOutputOptions.mRealElementExternalFieldNames.size();
            moris::Cell<moris::mtk::Scalar_Field_Info<DDRMat>> tExternalRealScalarFields(tNumExtRealScalarFields);
            for(uint i = 0; i<tNumExtRealScalarFields; i++)
            {
                tExternalRealScalarFields(i).set_field_name(aOutputOptions.mRealElementExternalFieldNames(i));
                tExternalRealScalarFields(i).set_field_entity_rank(moris::EntityRank::ELEMENT);
                add_field_for_mesh_input(&tExternalRealScalarFields(i),tFieldsInfo);
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
            moris::Matrix<moris::IdMat> tNodeOwner(1,tOutputtedNodeInds.numel(),moris::par_rank());

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
                 moris::mtk::MtkSideSetInfo tInterfaceSideSet;
                 tInterfaceSideSet.mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd;
                 tInterfaceSideSet.mSideSetName        = "iside" ;

                 // Add side side set to mesh sets
                 tMtkMeshSets.add_side_set(&tInterfaceSideSet);
             }

             // propogate side sets from background mesh
             moris::Cell<moris::Matrix<IndexMat>> tSideSetData;
             moris::Cell<moris::mtk::MtkSideSetInfo> tBackgroundSideSets;
             if(aOutputOptions.mAddSideSets)
             {
                 tBackgroundSideSets = this->propogate_background_side_sets(tSideSetData,aOutputOptions);

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

            tMeshDataInput.CreateAllEdgesAndFaces  = false;
            tMeshDataInput.Verbose                 = true;
            tMeshDataInput.SpatialDim              = &tSpatialDim;

            tCount = 0;
            for(moris::uint  i = 0 ; i <mGeometryEngine.get_num_bulk_phase(); i++)
            {
                if(aOutputOptions.output_phase(i))
                {
                tMeshDataInput.ElemConn(tCount)             = &tElementToNodeChildrenByPhase(i);
                tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tElementChildrenIdsByPhase(i);
                tCount++;
                tMeshDataInput.ElemConn(tCount)             = &tElementToNodeNoChildrenByPhase(i);
                tMeshDataInput.LocaltoGlobalElemMap(tCount) = &tElementNoChildrenIdsByPhase(i);
                tCount++;
                }
            }

            tMeshDataInput.NodeCoords              = &tNodeCoordinates;
            tMeshDataInput.NodeProcOwner           = &tNodeOwner;
            tMeshDataInput.FieldsInfo              = &tFieldsInfo;
            tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
            tMeshDataInput.SetsInfo                = &tMtkMeshSets;
            tMeshDataInput.MarkNoBlockForIO        = false;

            // Interface elements
            if(mUnzipped)
            {
                tMeshDataInput.ElemConn(tNumElemTypes-1)             = &tInterfaceElements;
                tMeshDataInput.LocaltoGlobalElemMap(tNumElemTypes-1) = &tInterfaceElementIds;
            }

            if(moris::par_rank() == 0 && mVerbose)
            {
                std::cout<<"XTK: Mesh data setup completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
            }


            start = std::clock();

            moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshDataInput );

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

        aOutputtedNodeInds.resize(tNodeMap.n_rows(),tNodeMap.n_cols());
        moris::uint tNumNodes = tNodeMap.numel();
        // if we are returning all phases there is no reason to restrict the map
        if(aOutputOptions.output_all_phases())
        {

            for(moris::uint i = 0; i <tNumNodes; i++)
            {
                aOutputtedNodeInds(i) = i;
            }

            return tNodeMap;
        }

        else
        {
            moris::Matrix<moris::IndexMat> tRestrictedNodeMap(tNodeMap.n_rows(),tNodeMap.n_cols());

            moris::uint tCount = 0;
            for(moris::uint i = 0; i <tNumNodes; i++)
            {
                if(output_node(i,aOutputOptions))
                {

                    moris::size_t tPhaseIndex = 0;
                    mGeometryEngine.get_phase_index(i,tPhaseIndex);
                    aOutputtedNodeInds(tCount) = i;
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
             moris::Matrix<moris::IndexMat> tNodesInSetInds = tMeshData.get_set_entity_loc_inds(EntityRank::NODE, tSetNames(i));

             aNodeSetData(i) = moris::Matrix<moris::IndexMat>(tNodesInSetInds.n_rows(),tNodesInSetInds.n_cols());

             for(moris::uint iNode =0; iNode<tNodesInSetInds.numel(); iNode++)
             {
                 if(this->output_node(tNodesInSetInds(iNode),aOutputOptions))
                 {
                     aNodeSetData(i)(tCount) = tNodesInSetInds(iNode);
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
    Model::propogate_background_side_sets(moris::Cell<moris::Matrix<IndexMat>> & aSideSetData,
                                       Output_Options const & aOutputOptions)
        {
            // access background mesh data
            moris::mtk::Mesh & tMeshData = mBackgroundMesh.get_mesh_data();

            // get all side set names in background mesh
            moris::Cell<std::string> tSetNames = tMeshData.get_set_names(EntityRank::FACE);

            // allocate output side sets
            moris::Cell<moris::mtk::MtkSideSetInfo> tSideSets(2*tSetNames.size());
            aSideSetData = moris::Cell<moris::Matrix<IndexMat>>(2*tSetNames.size());

            // declare matrices used through
            moris::Matrix< IndexMat > tElementsAttachedToFace(1,1);
            moris::Matrix< moris::IdMat >    tChildElemsIdsOnFace;
            moris::Matrix< moris::IndexMat > tChildElemsCMIndOnFace;
            moris::Matrix< moris::IndexMat > tChildElemOnFaceOrdinal;

            moris::uint tNumElementsAttachedToFace = 0;
            moris::uint tElementIndex = 0;
            moris::uint tPhaseIndex = 0;
            moris::uint tFaceOrdinal =0;
            moris::moris_index tChildMeshIndex = 0;
            moris::moris_id    tElementId =0;
            bool        tHasChildren = false;

            for(moris::uint i = 0; i <tSetNames.size(); i++)
            {
                // get sides in set i
                Matrix< IndexMat > tSidesInSet = tMeshData.get_set_entity_loc_inds(EntityRank::FACE, tSetNames(i));

                moris::uint tNoChildInd = 2*i;
                moris::uint tChildInd = 2*i+1;

                // side set data non-intersected
                aSideSetData(tNoChildInd)   = Matrix<IndexMat>(tSidesInSet.numel()*2,2);

                // intersected data
                //TODO: FIGURE OUT MAXIMUM VALUE
                aSideSetData(tChildInd) = Matrix<IndexMat>(tSidesInSet.numel()*2*10,2);

                // keep count
                moris::Cell<moris::uint> tCount(2,0);

                // iterate through sides in set i
                for(moris::uint iSide= 0; iSide<tSidesInSet.numel(); iSide++)
                {
                    moris::uint tSideIndex = tSidesInSet(iSide);

                    tElementsAttachedToFace = tMeshData.get_entity_connected_to_entity_loc_inds(tSideIndex,EntityRank::FACE,EntityRank::ELEMENT);

                    tNumElementsAttachedToFace = tElementsAttachedToFace.numel();
                    for( moris::uint  iElem = 0; iElem<tNumElementsAttachedToFace; iElem++)
                    {
                        tElementIndex = tElementsAttachedToFace(iElem);
                        tHasChildren = mBackgroundMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT);
                        // get the faces from the child mesh
                        if(tHasChildren)
                        {
                            tChildMeshIndex = mBackgroundMesh.child_mesh_index(tElementsAttachedToFace(iElem),EntityRank::ELEMENT);

                            Child_Mesh const & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                            tChildMesh.get_child_elements_connected_to_parent_face(tSideIndex,
                                                                                   tChildElemsIdsOnFace,
                                                                                   tChildElemsCMIndOnFace,
                                                                                   tChildElemOnFaceOrdinal);

                            moris::Matrix< moris::IndexMat > const & tChildElementPhaseIndices = tChildMesh.get_element_phase_indices();
                            moris::Matrix< moris::IndexMat > const & tElementIds = tChildMesh.get_element_ids();

                            for(moris::moris_index iCElem  = 0; iCElem < (moris::moris_index)tChildElemsCMIndOnFace.numel(); iCElem++)
                            {
                                tPhaseIndex = tChildElementPhaseIndices(0,tChildElemsCMIndOnFace(0,iCElem));
                                if(aOutputOptions.output_phase(tPhaseIndex))
                                {
                                    // Child Element Id
                                    tElementId = tElementIds(tChildElemsCMIndOnFace(iCElem));
                                    tFaceOrdinal   = tChildElemOnFaceOrdinal(iCElem);
                                    aSideSetData(tChildInd)(tCount(1),0) = tElementId;
                                    aSideSetData(tChildInd)(tCount(1),1) = tFaceOrdinal;
                                    tCount(1)++;
                                }
                            }
                        }

                        else
                        {

                            tFaceOrdinal = tMeshData.get_facet_ordinal_from_cell_and_facet_loc_inds(tSideIndex,tElementsAttachedToFace(0,iElem));
                            tElementId = tMeshData.get_glb_entity_id_from_entity_loc_index(tElementsAttachedToFace(0,iElem),EntityRank::ELEMENT);

                            if(aOutputOptions.output_phase(mBackgroundMesh.get_element_phase_index(tElementsAttachedToFace(0,iElem))))
                            {
                                aSideSetData(tNoChildInd)(tCount(0),0) = tElementId;
                                aSideSetData(tNoChildInd)(tCount(0),1) = tFaceOrdinal;
                                tCount(0)++;
                            }

                        }
                    }
                }

                // resize data
                aSideSetData(tChildInd).resize(tCount(1),2);
                aSideSetData(tNoChildInd).resize(tCount(0),2);

                // Add data to side set info
                // no child
                tSideSets(tNoChildInd).mElemIdsAndSideOrds = &aSideSetData(tNoChildInd);
                tSideSets(tNoChildInd).mSideSetName = tSetNames(i);
                tSideSets(tChildInd).mElemIdsAndSideOrds = &aSideSetData(tChildInd);
                tSideSets(tChildInd).mSideSetName = tSetNames(i) + "_i";
            }


            return tSideSets;
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
