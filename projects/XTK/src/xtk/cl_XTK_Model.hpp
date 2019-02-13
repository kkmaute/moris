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
#include "cl_XTK_Background_Mesh.hpp"


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
class Model
{
public:
    // Public member functions
    bool mVerbose = false;


    // Forward declare the maximum value of moris::size_t and moris::real
    moris::real REAL_MAX    = std::numeric_limits<moris::real>::max();
    moris::real INTEGER_MAX = MORIS_INDEX_MAX;

    Model(){};

    /**
     * Primary constructor (this constructor is used for all cases except when testing something)
     */
    Model(uint aModelDimension,
          moris::mtk::Mesh* aMeshData,
          Geometry_Engine & aGeometryEngine,
          bool aLinkGeometryOnConstruction = true) :
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

    bool mSameMesh;

    ~Model(){}


    void
    link_background_mesh_to_geometry_objects()
    {
        // initialize geometry objects
        moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

        mGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes(tNodeCoords.n_rows());
        mGeometryEngine.initialize_geometry_object_phase_values(tNodeCoords);
        mLinkedBackground = true;
    }

    /*!
     * Decomposes a mesh using a geometry engine (split up across processors)
     * @param aMethods - specify which type of subdivision method to use (this could be changed to command line parsing or XML reading)
     * @param aSetPhase - tell it to set phase information
     */
    void decompose(Cell<enum Subdivision_Method> aMethods,
                   bool                          aSetPhase  = true)
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
                this->subdivide(aMethods(iDecomp), tActiveChildMeshIndices, tFirstSubdivisionFlag, tNonConformingMeshFlag);

                // Change the first subdivision flag as false
                tFirstSubdivisionFlag = false;

                // print timing
                if(moris::par_rank() == 0 && mVerbose)
                {
                    std::cout<<"XTK: Decomposition "<<get_enum_str(aMethods(iDecomp))<<" for geometry "<<iGeom<< " completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
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
    unzip_interface()
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

    /*!
     * Compute sensitivity. Must be called after decompose.
     */
    void
    compute_sensitivity()
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

    /*!
     * Perform the generalized heaviside enrichment
     */
    void
    perform_basis_enrichment()
    {
        MORIS_ERROR(mDecomposed,"Prior to computing sensitivity, the decomposition process must be called");
        MORIS_ERROR(!mEnriched,"Calling perform_basis_enrichment twice is not supported");
        perform_basis_enrichment_internal();

        mEnriched = true;
    }

    /*!
     * Convert Tet4 elements to Tet10 Elements
     */
    void
    convert_mesh_tet4_to_tet10()
    {

        mConvertedToTet10s = true;

        // start timing on this decomposition
        std::clock_t start = std::clock();
//        convert_mesh_tet4_to_tet10_internal();

        if(moris::par_rank() == 0)
        {
            std::cout<<"Tet4 to Tet10 conversion completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }

    /*!
     * Returns the Cut Mesh
     */
    Cut_Mesh &
    get_cut_mesh()
    {
        return mCutMesh;
    }

    /*!
     * Returns the Cut Mesh
     */
    Cut_Mesh const &
    get_cut_mesh() const
    {
        return mCutMesh;
    }

    /*!
     * Returns the Background Mesh
     */
    Background_Mesh &
    get_background_mesh()
    {
        return mBackgroundMesh;
    }

    /*!
     * Returns the Xtk Mesh
     */
    Background_Mesh const &
    get_background_mesh() const
    {
        return mBackgroundMesh;
    }


    /*!
     * Get geomtry engine
     */

    Geometry_Engine &
    get_geom_engine()
    {
        return mGeometryEngine;
    }

    /*!
     * Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
     */
    moris::mtk::Mesh*
    get_output_mesh(Output_Options const & aOutputOptions = Output_Options())

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


private:
    uint            mModelDimension;
    Background_Mesh mBackgroundMesh;
    Cut_Mesh        mCutMesh;
    Geometry_Engine mGeometryEngine;
    Enrichment      mEnrichment;

    // XTK Model State Flags
    bool mLinkedBackground  = false; // Model background mesh linked to geometry model
    bool mDecomposed        = false; // Model has been decomposed
    bool mSensitivity       = false; // Model has computed sensitivity
    bool mConvertedToTet10s = false; // Model has been converted from tet4's to tet10s
    bool mEnriched          = false; // Model has been enriched
    bool mUnzipped          = false; // Model has been enriched

    // The midside nodes are stored here currently but this may change
    moris::Matrix< moris::IndexMat > mMidsideElementToNode;



    // Private Functions
private:
    // Decomposition Functions------------------------------------------------------
    /*!
     * formulates node requests in the geometry objects. Dependent on the type of decomposition
     * @param[in] aReqType- specifies which template mesh is going to be used later on
     */
    void subdivide(enum Subdivision_Method    const & aSubdivisionMethod,
                   moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                   bool const & aFirstSubdivision = true,
                   bool const & aSetIds = false)
    {
        moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();
        switch (aSubdivisionMethod)
        {
        case (Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8):
        {
            bool tCoordFlag = false;
            MORIS_ASSERT(tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, moris::EntityRank::ELEMENT, moris::EntityRank::NODE).numel() == 8, "NC_REGULAR_SUBDIVISION_HEX8 is for HEX8 meshes only.");
            MORIS_ASSERT(aFirstSubdivision,"NC_REGULAR_SUBDIVISION_HEX8 needs to be the first subdivision routine for each geometry");


            // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones dont
            moris::Matrix< moris::IndexMat > tNewPairBool;
            run_first_cut_routine(TemplateType::HEX_8, 8,  aActiveChildMeshIndices,tNewPairBool);

            // Initialize request list for faces and elements
            moris::uint tIntersectedCount = aActiveChildMeshIndices.n_cols();
            moris::uint tNumFaceRequests = 6; // (one per face on each element)
            moris::uint tNumElemRequests = 1; // (one internal element request)
            moris::uint tNumChildAllow = 1; // (one child allowed per parent entity for this subdivision)
            moris::uint tNumNewNodes = tNumFaceRequests + tNumElemRequests; //(7)
            Request_Handler tFaceRequests(tIntersectedCount * tNumFaceRequests, tNumChildAllow, EntityRank::FACE, EntityRank::NODE, mBackgroundMesh, mCutMesh); // 6 face requests per element
            Request_Handler tElemRequests(tIntersectedCount * tNumElemRequests, tNumChildAllow, EntityRank::ELEMENT, EntityRank::NODE, mBackgroundMesh, mCutMesh); // 1 face request per element

            // Initialize topologies used in this method
            Edge_Topology         tEdgeTopology;
            Quad_4_Topology       tFaceTopology;
            Hexahedron_8_Topology tElementTopology;

            // Initialize a cell of pointers to future node index
            Cell<moris::moris_index*> tNodeInds(tNumNewNodes);
            moris::Matrix< moris::IndexMat > tFaceNodes(1,4);
            moris::Matrix< moris::IndexMat > tElementNodes(1,8);
            moris::Matrix< moris::IndexMat > tFaceIndices(1,6);
            moris::Matrix< moris::DDRMat > tCoordinates(0,0);
            moris::Matrix< moris::DDRMat > tNewNodeCoordinates(1,mModelDimension,0);
            moris::Matrix< moris::DDRMat > tCenterFaceLocCoordinate(1,2,0.0);
            moris::Matrix< moris::DDRMat > tCenterElementLocCoordinate(1,3,0.0);

            // Parametric coordinates for this subdivision routine
            const moris::Matrix< moris::DDRMat > tParamCoords(
                    {{ 0.0, -1.0,  0.0},
                     { 1.0,  0.0,  0.0},
                     { 0.0,  1.0,  0.0},
                     {-1.0,  0.0,  0.0},
                     { 0.0,  0.0, -1.0},
                     { 0.0,  0.0,  1.0},
                     { 0.0,  0.0,  0.0}});

            // Loop over xtk meshes and place a node request on each face and at center of element volume
            for (moris::size_t i = 0; i < tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {

                // Get element index
                moris::moris_index tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

                // Get local index of faces connected to element using local element index
                tFaceIndices = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::FACE);

                // Loop over faces (6 in a hex 8) and set a node request.
                // Request will return a pointer to where the created node index will be placed
                for (moris::size_t fi = 0; fi < 6; fi++)
                {
                    tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds((moris_index)tFaceIndices(fi), moris::EntityRank::FACE, moris::EntityRank::NODE);
                    tFaceTopology.set_node_indices(tFaceNodes);

                    tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tFaceNodes);
                    xtk::Interpolation::bilinear_interpolation(tCoordinates, tCenterFaceLocCoordinate,tNewNodeCoordinates);
                    tNodeInds(fi) = tFaceRequests.set_request_info(tFaceIndices(fi), tFaceTopology, tNewNodeCoordinates, tCenterFaceLocCoordinate);
                }

                // Place node at center of element
                tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);
                tElementTopology.set_node_indices(tElementNodes);
                tCoordinates = mBackgroundMesh.get_selected_node_coordinates_loc_inds(tElementNodes);
                xtk::Interpolation::trilinear_interpolation(tCoordinates, tCenterElementLocCoordinate, tNewNodeCoordinates);
                tNodeInds(6) = tElemRequests.set_request_info(tElemInd, tElementTopology, tNewNodeCoordinates, tCenterElementLocCoordinate);



                // Give XTK Mesh pointers to where its node indices will be located
                mCutMesh.set_pending_node_index_pointers(aActiveChildMeshIndices(i), tNodeInds, tParamCoords);
                }
            }
            tFaceRequests.handle_requests(tCoordFlag, mSameMesh, false, mGeometryEngine);
            tElemRequests.handle_requests(tCoordFlag, mSameMesh, false, mGeometryEngine);

            // Allocate interface flag space in XTK mesh
            mBackgroundMesh.allocate_space_in_interface_node_flags(tFaceRequests.get_num_requests(),mGeometryEngine.get_num_geometries());
            mBackgroundMesh.allocate_space_in_interface_node_flags(tElemRequests.get_num_requests(),mGeometryEngine.get_num_geometries());
            // Tell XTK mesh to retrieve pending node indices and generate tabulated meshes
            mCutMesh.retrieve_pending_node_inds();
            for(moris::size_t i = 0; i< tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {
                mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(i),TemplateType::REGULAR_SUBDIVISION_HEX8);
                }
            }
            if(aSetIds)
            {
                moris::Matrix< moris::IdMat > tNodeIds;
                for (moris::size_t j = 0; j < tIntersectedCount; j++)
                    if(tNewPairBool(0,j) == 0)
                {
                    {
                        moris::Matrix< moris::IndexMat > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(j));
                        tNodeIds = mBackgroundMesh.get_glb_entity_id_from_entity_loc_index_range(tNodeIndices, EntityRank::NODE);
                        mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j),tNodeIds);
                    }
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

            // Initialize
            moris::size_t tEdgeInd = INTEGER_MAX;
            moris::size_t tNumMesh = aActiveChildMeshIndices.n_cols();

            // Initialize topologies used in this method (all local coordinates are with respect to an edge)
            Edge_Topology tEdgeTopology;

            moris::Matrix< moris::DDRMat > tLocalCoord(1,1, 0); // ALong an edge
            moris::Matrix< moris::DDRMat > tEdgeNodeParamCoordinates(2,tDimParamCoord); // parametric coordinate of end nodes wrt parent element
            moris::Matrix< moris::DDRMat > tNewNodeParamCoord(1,tDimParamCoord); // new node parametric coordinate wrt parent element
            moris::Matrix< moris::DDRMat > tEdgeCoords(2, 3, REAL_MAX);
            moris::Matrix< moris::DDRMat > tGlobalCoord(1, 3, REAL_MAX);

            moris::Matrix< moris::IndexMat > tEdgeNodes(1, 2, INTEGER_MAX);
            moris::Matrix< moris::IndexMat >  tParentInfo(1, 2, INTEGER_MAX);

            // Initialize request list for faces and elements
            // Number of children allowed on a parent mesh entity
            moris::size_t tEdgeChildren = 4*2*10;
            moris::size_t tFaceChildren = 4*4*10;
            moris::size_t tElemChildren = 4*24*10;

            // Of all the edges in a regular subdivision template 12 live on parent mesh edges, 24 live on parent mesh faces and 14 live in the parent mesh element
            // the number of parent edges and faces were done by hand and are hardcoded here.
            // TODO: Ask XTKMesh how many parents of each rank it has (This would eliminate the dependency on the regular subdivision coming first)
            moris::size_t tNumParentEdge = 2*12;
            moris::size_t tNumParentFace = 2*24;
            moris::size_t tNumParentElem = 2*14;

            //TODO: MAKE Request_Handler MORE MEMORY EFFICIENT!!!!!
            Request_Handler tEdgeRequests(tNumMesh * tNumParentEdge, tEdgeChildren, EntityRank::EDGE, EntityRank::NODE, mBackgroundMesh, mCutMesh);
            Request_Handler tFaceRequests(tNumMesh * tNumParentFace, tFaceChildren, EntityRank::FACE, EntityRank::NODE, mBackgroundMesh, mCutMesh);
            Request_Handler tElemRequests(tNumMesh * tNumParentElem, tElemChildren, EntityRank::ELEMENT, EntityRank::NODE, mBackgroundMesh, mCutMesh);
            // Tell XTKMesh to initialize intersection connectivity
            mCutMesh.init_intersect_connectivity(aActiveChildMeshIndices);

            // Check type specified as conformal (could change this to enum)
            moris::size_t tCheckType = 1;
            moris::Matrix< moris::DDRMat > tNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();
            // Ask the geometry engine whether it has sensitivity information
            bool tHasDxDp = mGeometryEngine.mComputeDxDp;

            for (moris::size_t j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {
                // Initialize geometry objects
                Cell<Geometry_Object> tGeoObjects;

                // Get the child mesh that is active
                Child_Mesh & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));


                // Get full edge to element connectivity from XTK Mesh (probably slow)
                // 0 specifies XTK local indices if an analytic geometry
                // Otherwise this needs to be the processor local index
                // or even the ID
                moris::Matrix< moris::IndexMat > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                // Ask geometry engine which edges are intersected (Simple mesh local indexed edges)
                mGeometryEngine.is_intersected(tNodeCoords, tEdgeToNode, tCheckType, tGeoObjects);

                // Initialize node index pointers based on number of intersected edges and parametric coordinates
                uint tNumNewNodes = 0;
                Cell<moris::moris_index*> tNodeInds(tGeoObjects.size());
                moris::Matrix< moris::DDRMat > tParametricCoords(tGeoObjects.size(),tDimParamCoord);

                // get reference to child mesh edge parent information
                moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                moris::Matrix< moris::DDSTMat > const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();
                for (moris::size_t k = 0; k < tGeoObjects.size(); k++)
                {

                    if(!tGeoObjects(k).has_parent_nodes_on_interface())
                    {
                        // Local index to XTK Mesh
                        tEdgeInd = tGeoObjects(k).get_parent_entity_index();

                        // get a local coordinate along the intersected edge [-1,1]
                        tLocalCoord(0,0) = tGeoObjects(k).get_interface_lcl_coord();

                        // get the interpolated global coordinate
                        tGlobalCoord = tGeoObjects(k).get_interface_glb_coord();

                        // Add edge to the entity intersection connectivity
                        mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), tNumNewNodes, tEdgeInd, 0);

                        // Edge nodes
                        tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);

                        // Compute new node parametric coordinate with respect to the current parent element
                        tEdgeNodeParamCoordinates.set_row(0, tChildMesh.get_parametric_coordinates(tEdgeNodes(0)));
                        tEdgeNodeParamCoordinates.set_row(1, tChildMesh.get_parametric_coordinates(tEdgeNodes(1)));
                        tParametricCoords.set_row(tNumNewNodes,Interpolation::linear_interpolation_location(tEdgeNodeParamCoordinates,tLocalCoord));

                        // Parent edge information
                        moris::size_t tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                        moris::moris_index tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                        // Set parent edge topology information with process local indices
                        // prior to converting to global ids for secondary id calculation.
                        tEdgeTopology.set_node_indices(tEdgeNodes);

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

                        if (tParentRank == 1)
                        {

                            // Intersected edge is an existing stk edge
                            // Make request in edge requests
                            // This does not require a supplemental identifier
                            // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING!!!!!!
                            moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));

                            tNodeInds(tNumNewNodes) = tEdgeRequests.set_request_info(tParentIndex,
                                                                          tSecondaryId,
                                                                          tEdgeTopology,
                                                                          tGlobalCoord,
                                                                          tLocalCoord,
                                                                          tGeoObjects(k).get_sensitivity_dx_dp(),
                                                                          tGeoObjects(k).get_node_adv_indices(),
                                                                          tHasDxDp,
                                                                          tHasDxDp);
                        }


                        else if (tParentRank == 2)
                        {
                            // Intersected edge was built on an stk face
                            // Make request in face requests
                            // This requires a supplemental identifier
                            moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));

                            tNodeInds(tNumNewNodes) = tFaceRequests.set_request_info(tParentIndex,
                                                                          tSecondaryId,
                                                                          tEdgeTopology,
                                                                          tGlobalCoord,
                                                                          tLocalCoord,
                                                                          tGeoObjects(k).get_sensitivity_dx_dp(),
                                                                          tGeoObjects(k).get_node_adv_indices(),
                                                                          tHasDxDp,
                                                                          tHasDxDp);
                        }
                        //
                        else if (tParentRank == 3)
                        {
                            // Intersected edge was built in stk element
                            // Make request in element requests
                            // This requires a supplemental identifier
                            moris::moris_index tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                            tNodeInds(tNumNewNodes) = tElemRequests.set_request_info(tParentIndex,
                                                                          tSecondaryId,
                                                                          tEdgeTopology,
                                                                          tGlobalCoord,
                                                                          tLocalCoord,
                                                                          tGeoObjects(k).get_sensitivity_dx_dp(),
                                                                          tGeoObjects(k).get_node_adv_indices(),
                                                                          tHasDxDp,
                                                                          tHasDxDp);
                        }

                        else
                        {
                            std::cout << "Invalid ancestry returned from XTK Mesh";
                        }

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

                tNodeInds.resize(tNumNewNodes,NULL);
                tParametricCoords.resize(tNumNewNodes,3);

                mCutMesh.set_pending_node_index_pointers(aActiveChildMeshIndices(0,j), tNodeInds,tParametricCoords);

                // Coincident edges have been marked, use this to create interface elements with side ordinal
                tChildMesh.mark_interface_faces_from_interface_coincident_faces();

            } // XTK Mesh loop

            bool tCoordinateFlag = false;

            // handle requests
            tEdgeRequests.handle_requests(tCoordinateFlag, mSameMesh, true, mGeometryEngine);
            tFaceRequests.handle_requests(tCoordinateFlag, mSameMesh, true, mGeometryEngine);
            tElemRequests.handle_requests(tCoordinateFlag, mSameMesh, true, mGeometryEngine);

            // Allocate interface flag space in XTK mesh
            mBackgroundMesh.allocate_space_in_interface_node_flags(tEdgeRequests.get_num_requests(),mGeometryEngine.get_num_geometries());
            mBackgroundMesh.allocate_space_in_interface_node_flags(tFaceRequests.get_num_requests(),mGeometryEngine.get_num_geometries());
            mBackgroundMesh.allocate_space_in_interface_node_flags(tElemRequests.get_num_requests(),mGeometryEngine.get_num_geometries());

            // Mark these nodes as interface nodes with respect to the active geometry
            tEdgeRequests.mark_pending_nodes_as_interface_nodes(mBackgroundMesh,mGeometryEngine.get_active_geometry_index());
            tFaceRequests.mark_pending_nodes_as_interface_nodes(mBackgroundMesh,mGeometryEngine.get_active_geometry_index());
            tElemRequests.mark_pending_nodes_as_interface_nodes(mBackgroundMesh,mGeometryEngine.get_active_geometry_index());

            // Tell XTK mesh to grab the pending node indices and that these nodes are on the interface and have sensitivity
            mCutMesh.retrieve_pending_node_inds();

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

    /*
     * Perform all tasks needed to finalize the decomposition process, such that the model is ready for enrichment, conversion to tet10 etc.
     * Tasks performed here:
     *  - Assign all child elements global ids and processor local indices
     *  - Store phase indices of non-intersected parent elements
     *
     */
    void
    finalize_decomp_in_xtk_mesh(bool aSetPhase)
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

        create_child_element_mtk_cells();


        // Change XTK model decomposition state flag
        mDecomposed = true;
    }

    /*!
     * Constructs the mtk cell interface for all child elements created during the
     * decomposition process
     */
    void
    create_child_element_mtk_cells()
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
    /**
     * Take the interface faces and create collapsed prisms
     */
    void
    unzip_interface_internal()
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

    void
    unzip_interface_internal_assign_node_identifiers(moris::uint aNumNodes,
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

    void
    unzip_interface_internal_modify_child_mesh(moris::uint                         aGeometryIndex,
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

    moris::Matrix< moris::IndexMat >
    unzip_interface_internal_assign_which_element_uses_unzipped_nodes( moris::moris_index aGeometryIndex,
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



    moris::Cell<moris::Cell< moris::moris_index >>
    unzip_interface_internal_collect_child_mesh_to_interface_node(moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
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

    void
    unzip_interface_construct_interface_elements(moris::uint aGeometryIndex,
                                                 moris::Matrix< moris::IndexMat > const & aElementPairs,
                                                 moris::Matrix< moris::IndexMat > const & aSideOrdinalPairs)
    {
        moris::uint tNumPairs = aElementPairs.n_cols();

        moris::Cell<const moris::mtk::Cell*> tPairCells(2);

        for(moris::uint i = 0; i <tNumPairs; i++)
        {
            tPairCells(0) = mBackgroundMesh.get_child_element_mtk_cell(aElementPairs(0,i));
            tPairCells(1) = mBackgroundMesh.get_child_element_mtk_cell(aElementPairs(1,i));

            // construct an interface element
            Interface_Element tInterfaceElement;
            tInterfaceElement.set_element_pair_and_side_ordinal(tPairCells,aSideOrdinalPairs.get_column(i));

            // Add interface element to cut mesh
            mCutMesh.add_interface_element(tInterfaceElement);
        }
    }

    void
    unzip_interface_assign_element_identifiers()
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


    // Sensitivity computation functions -----------------------------------------------
    void
    compute_interface_sensitivity_internal()
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

    // Enrichment computation functions -----------------------------------------------

    void
    perform_basis_enrichment_internal()
    {
        // initialize enrichment
        mEnrichment = Enrichment(mGeometryEngine.get_num_phases(),
                                 &this->get_cut_mesh(),
                                 &this->get_background_mesh());

        // Set verbose flag to match XTK.
        mEnrichment.mVerbose = mVerbose;

        // perform the enrichment
        mEnrichment.perform_enrichment();
    }

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

    void  run_first_cut_routine(enum TemplateType const & aTemplateType,
                                moris::size_t const & tNumNodesPerElement,
                                moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                                moris::Matrix< moris::IndexMat > & aNewPairBool)
    {
        // Note this method is independent of node ids for this reason Background_Mesh is not given the node Ids during this subdivision
        moris::mtk::Mesh & tXTKMeshData = mBackgroundMesh.get_mesh_data();

        // Package up node to element connectivity
        moris::moris_index tParentElementIndex = INTEGER_MAX;
        moris::size_t tNumElements = mBackgroundMesh.get_num_entities(EntityRank::ELEMENT);
        moris::Matrix< moris::IndexMat > tNodetoElemConnInd (tNumElements, tNumNodesPerElement);
        moris::Matrix< moris::IndexMat > tNodetoElemConnVec (1, tNumNodesPerElement);
        moris::Matrix< moris::IndexMat > tEdgetoElemConnInd (1, 1);
        moris::Matrix< moris::IndexMat > tFacetoElemConnInd (1, 1);
        moris::Matrix< moris::IndexMat > tElementMat(1, 1);
        moris::Matrix< moris::IndexMat > tPlaceHolder(1, 1);

        for (moris::size_t i = 0; i < tNumElements; i++)
        {
            tNodetoElemConnVec = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i, moris::EntityRank::ELEMENT, moris::EntityRank::NODE);
            tNodetoElemConnInd.set_row(i, tNodetoElemConnVec);
        }

        // Get the Node Coordinates
        moris::Matrix< moris::DDRMat > tAllNodeCoords = mBackgroundMesh.get_all_node_coordinates_loc_inds();

        // Intersected elements are flagged via the Geometry_Engine
        Cell<Geometry_Object> tGeoObjects;
        mGeometryEngine.is_intersected(tAllNodeCoords, tNodetoElemConnInd, 0,tGeoObjects);

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

        for (moris::size_t j = 0; j < tIntersectedCount; j++)
        {
            if(aNewPairBool(0,j) == 0)
            {
                tParentElementIndex = tGeoObjects(j).get_parent_entity_index();

                // Get information to provide ancestry
                // This could be replaced with a proper topology implementation that knows faces, edges based on parent element nodes
                tNodetoElemConnVec = tNodetoElemConnInd.get_row(tParentElementIndex);
                tEdgetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, moris::EntityRank::ELEMENT, moris::EntityRank::EDGE);
                tFacetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, moris::EntityRank::ELEMENT, moris::EntityRank::FACE);
                tElementMat(0,0) = tParentElementIndex;

                // Set parent element, nodes, and entity ancestry
                moris::Matrix< moris::IndexMat > tElemToNodeMat(tNodetoElemConnVec);

                Cell<moris::Matrix< moris::IndexMat >> tAncestorInformation = {tPlaceHolder, tEdgetoElemConnInd, tFacetoElemConnInd, tElementMat};
                mCutMesh.initialize_new_mesh_from_parent_element(aActiveChildMeshIndices(0,j), aTemplateType, tNodetoElemConnVec, tAncestorInformation);
            }
        }
    }


    moris::mtk::Mesh*
    construct_output_mesh( Output_Options const & aOutputOptions )
    {

        // start timing on this decomposition
        std::clock_t start = std::clock();

        // Get mesh information ready for outputting
        // Package element to Node Connectivity
        moris::uint                         tSpatialDim         = mBackgroundMesh.get_mesh_data().get_spatial_dim();

        // Children element nodes connected to elements
        moris::Matrix<moris::IdMat>  tElementToNodeChildren     = mCutMesh.get_full_element_to_node_glob_ids();

        // Child element ids
        moris::Matrix<moris::IdMat>  tElementChildrenIds        = mCutMesh.get_all_element_ids();

        // Parent elements without children
        moris::Matrix<moris::IdMat>  tElementNoChildrenIds      = mBackgroundMesh.get_all_non_intersected_elements();

        // Connectivity of parent elements without children
        moris::Matrix<moris::IdMat>  tElementToNodeNoChildren   = mBackgroundMesh.get_full_non_intersected_node_to_element_glob_ids();

        // Node map
        moris::Matrix<moris::IdMat>  tLocalToGlobalNodeMap      = mBackgroundMesh.get_local_to_global_map(EntityRank::NODE);

        // All node coordinates
        moris::Matrix<moris::DDRMat> tNodeCoordinates           = mBackgroundMesh.get_all_node_coordinates_loc_inds();

        // Number of bulk phases
        uint tNumBulkPhases = mGeometryEngine.get_num_phases();

        // Get child elements sorted by phase
        Cell<moris::Matrix<moris::IdMat>> tChildElementsByPhase   = mCutMesh.get_child_elements_by_phase(tNumBulkPhases);

        // Get non-interescted parent elements by phase
        Cell<moris::Matrix<moris::IdMat>> tNoChildElementsByPhase = mBackgroundMesh.get_all_non_intersected_elements_by_phase(tNumBulkPhases);

        // Interface nodes
        Cell<moris::Matrix<moris::IndexMat>> tInterfaceNodes = mBackgroundMesh.get_interface_nodes_glob_ids();

        // Assemble geometry data as field for mesh output
        moris::Cell< moris::Matrix < moris::DDRMat > > tGeometryFieldData = assemble_geometry_data_as_mesh_field();

        // Give the geometry data a name
        moris::Cell<std::string> tGeometryFieldNames = assign_geometry_data_names();

        // Get rank of the geometry data field
        moris::Cell < enum moris::EntityRank > tFieldRanks =  assign_geometry_data_field_ranks();

        // Get the packaged interface side sets from the cut mesh
        moris::Matrix<moris::IdMat> tInterfaceElemIdandSideOrd = mCutMesh.pack_interface_sides();


        // Set up field data structure for MTK input.
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


        //TODO: implement node owner (currently set to owned by this proc)
        moris::Matrix<moris::IdMat> tNodeOwner(1,mBackgroundMesh.get_num_entities(EntityRank::NODE),moris::par_rank());

        // Set up mesh sets
        // Initialize Sets information structure
         moris::mtk::MtkSetsInfo tMtkMeshSets;

         // Setup block sets
         Cell<moris::mtk::MtkBlockSetInfo> tBlockSets(tNumBulkPhases*2);
         uint tCount= 0;

         for(uint i = 0; i <tNumBulkPhases; i++)
         {
             // Children of material phase i
             tBlockSets(tCount).mCellIdsInSet = &tChildElementsByPhase(i);
             tBlockSets(tCount).mBlockSetName = "child_"+std::to_string(i);
             tBlockSets(tCount).mBlockSetTopo = CellTopology::TET4;

             tMtkMeshSets.add_block_set(&tBlockSets(tCount));
             tCount++;

             // Children of material phase i
             tBlockSets(tCount).mCellIdsInSet = &tNoChildElementsByPhase(i);
             tBlockSets(tCount).mBlockSetName = "no_child_"+std::to_string(i);
             tBlockSets(tCount).mBlockSetTopo = CellTopology::HEX8;

             tMtkMeshSets.add_block_set(&tBlockSets(tCount));
             tCount++;
         }

         // Interface elements
         moris::Matrix<moris::IndexMat> tInterfaceElements(0,0);
         moris::Matrix<moris::IndexMat> tInterfaceElementIds(0,0);
         moris::mtk::MtkBlockSetInfo tUnzippedInterfaceBlockSet;
         if(mUnzipped)
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



         moris::Cell<moris::mtk::MtkNodeSetInfo> tInterfaceNodeSets(tInterfaceNodes.size());
         for(uint i = 0; i<tInterfaceNodes.size(); i++)
         {
             tInterfaceNodeSets(i).mNodeIds     = &tInterfaceNodes(i);
             tInterfaceNodeSets(i).mNodeSetName = "inodes_" +std::to_string(i) ;
             tMtkMeshSets.add_node_set(&tInterfaceNodeSets(i));
         }

         moris::mtk::MtkSideSetInfo tInterfaceSideSet;
         tInterfaceSideSet.mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd;
         tInterfaceSideSet.mSideSetName        = "iside" ;

         // Add side side set to mesh sets
         tMtkMeshSets.add_side_set(&tInterfaceSideSet);


        // Mesh data input structure (with multi element type mesh)
        moris::mtk::MtkMeshData tMeshDataInput;

        moris::uint tNumElemTypes = 2;
        if(mUnzipped)
        {
            tNumElemTypes = 3;
        }


        tMeshDataInput.ElemConn             = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);
        tMeshDataInput.LocaltoGlobalElemMap = moris::Cell<moris::Matrix<IdMat>*>(tNumElemTypes);

        tMeshDataInput.CreateAllEdgesAndFaces  = false;
        tMeshDataInput.Verbose                 = true;
        tMeshDataInput.SpatialDim              = &tSpatialDim;
        tMeshDataInput.ElemConn(0)             = &tElementToNodeChildren;
        tMeshDataInput.ElemConn(1)             = &tElementToNodeNoChildren;
        tMeshDataInput.LocaltoGlobalElemMap(0) = &tElementChildrenIds;
        tMeshDataInput.LocaltoGlobalElemMap(1) = &tElementNoChildrenIds;
        tMeshDataInput.NodeCoords              = &tNodeCoordinates;
        tMeshDataInput.NodeProcOwner           = &tNodeOwner;
        tMeshDataInput.FieldsInfo              = &tFieldsInfo;
//        tMeshDataInput.LocaltoGlobalElemMap    = &aElemLocaltoGlobal;
        tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
        tMeshDataInput.SetsInfo                = &tMtkMeshSets;

        // Interface elements
        if(mUnzipped)
        {
            tMeshDataInput.ElemConn(2)             = &tInterfaceElements;
            tMeshDataInput.LocaltoGlobalElemMap(2) = &tInterfaceElementIds;
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

public:
    moris::Cell< moris::Matrix < moris::DDRMat > >
    assemble_geometry_data_as_mesh_field()
    {
        uint tNumGeometries = mGeometryEngine.get_num_geometries();
        uint tNumNodes      = mBackgroundMesh.get_num_entities(EntityRank::NODE);

        // Allocate output data
        moris::Cell< moris::Matrix < moris::DDRMat > > tGeometryData(tNumGeometries, moris::Matrix<moris::DDRMat>(tNumNodes,1));

        //Iterate through geometries
        for(uint iG = 0; iG <tNumGeometries; iG++)
        {
            // Iterate through nodes
            for(uint iN = 0; iN<tNumNodes; iN++)
            {
                tGeometryData(iG)(iN) = mGeometryEngine.get_entity_phase_val(iN,iG);
            }
        }

        return tGeometryData;
    }
private:

    moris::Cell<std::string>
    assign_geometry_data_names()
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


    moris::Cell < enum moris::EntityRank >
    assign_geometry_data_field_ranks()
    {
        uint tNumGeometries = mGeometryEngine.get_num_geometries();

        // base string of geometry data
        std::string tBaseName = "gd_";

        // Allocate output
        // Note: for now this is a nodal field always
        moris::Cell<enum moris::EntityRank> tGeometryFieldRank(tNumGeometries,moris::EntityRank::NODE);

        return tGeometryFieldRank;
    }


//        /*
//         * Initialize and Allocate
//         */
//        moris::mtk::Mesh const & tXTKMeshData = mXTKMesh.get_mesh_data();
//        enum EntityRank tElementRank          = EntityRank::ELEMENT;
//        moris::size_t tNumElementsInXTKMesh         = mXTKMesh.get_num_entities(tElementRank);
//        moris::size_t tNumChildrenMeshes            = mCutMesh.get_num_simple_meshes();
//        moris::size_t tNumElementsInCutMesh         = mCutMesh.get_num_entities(tElementRank);
//
//        moris::size_t tNumPhases = mGeometryEngines.get_num_phases();
//
//
//        // Package Elements into buckets
//        Cell<mesh::Bucket< moris::size_t,moris::DDSTMat >> tXTKElementBuckets;
//        // Not reserved assuming all children mesh is tetrahedrons and parent mesh has hex 8
//        // Also worst case would be there is a different field name for every element
//        tXTKElementBuckets.reserve((tNumElementsInXTKMesh-tNumChildrenMeshes)*8+tNumElementsInCutMesh*4);
//
//        package_elements_into_buckets( aOutputOptions, tXTKElementBuckets);
//
//
//        Cell<mesh::Side_Set_Input<moris::size_t,moris::DDSTMat>> tXTKSideSets;
//
//        // -----------------------------------------------------------------------------------
//        // Package Side Sets from Parent Mesh
//
//        Cell<std::string> tFacePartNames;
//        Cell<std::string> tAppendedFacePartNames;
//        if(aOutputOptions.mAddSideSets)
//        {
//            package_parent_side_sets_for_mesh_input(tXTKSideSets,aOutputOptions);
//
//            // Append parent face partspart names
//            tXTKMeshData.get_all_part_names(EntityRank::FACE, tFacePartNames);
//            tFacePartNames.resize(tFacePartNames.size()/2,"");
//            this->append_face_parts(aOutputOptions,tNumPhases, tFacePartNames,tAppendedFacePartNames);
//        }
//
//        // -----------------------------------------------------------------------------------
//        // Package Interface Side Sets
//
//        // Append interface name to all side sets
//        std::string       tIntefaceSideName;
//        Cell<std::string> tInterfaceSideNames;
//        for(moris::size_t iPhase = 0; iPhase<tNumPhases; iPhase++)
//        {
//            if(aOutputOptions.output_phase(iPhase))
//            {
//                tIntefaceSideName = aOutputOptions.mMaterialAppendix + aOutputOptions.mInterfaceAppendix + std::to_string(iPhase);
//                tInterfaceSideNames.push_back(tIntefaceSideName);
//            }
//        }
//
//        // Add to face part names
//        tAppendedFacePartNames.append(tInterfaceSideNames);
//
//        package_interface_side_sets_for_mesh_input(aOutputOptions, tXTKSideSets, tInterfaceSideNames);
//
//        Cell<mesh::Node_Set_Input<moris::real, moris::size_t, moris::DDRMat, moris::DDSTMat>> tXTKNodeSets;
//
//        if(aOutputOptions.mAddNodeSets)
//        {
//            package_node_sets_for_mesh_input(tXTKNodeSets);
//        }
//
//
//        /*
//         * Package Interface Node Set
//         */
//        std::string tNodeBaseStr = "nodes";
//        xtk::Sensitivity<moris::real, moris::size_t, moris::DDRMat, moris::DDSTMat> tSensitivity;
//        package_interface_node_sets_for_mesh_input(tXTKNodeSets,aOutputOptions.mInterfaceAppendix,tNodeBaseStr,tSensitivity,aOutputOptions);
//        tSensitivity.commit_sensitivities();
//
//        /*
//         * Initialize and get all part names which have a primary rank of element
//         */
//        Cell<std::string> tElementPartNames;
//        Cell<std::string> tAppendedElementPartNames;
//        Cell<enum CellTopology> tPartTopologies;
//        tXTKMeshData.get_all_part_names(tElementRank, tElementPartNames);
//
//        for(moris::size_t i = 0; i<tElementPartNames.size(); i++)
//        {
//            if(tElementPartNames(i).compare("unnamed_block_id:_1_type:_hex8") == 0)
//            {
//
//                std::cout<<"Match check 2"<<std::endl;
//                tElementPartNames(i) = "hex_bl_1";
//            }
//        }
//
//        enum CellTopology tParentTopology = mXTKMesh.get_XTK_mesh_element_topology();
//        if(tParentTopology == CellTopology::TET4 && mConvertedToTet10s)
//        {
//            tParentTopology = EntityTopology::TET_10;
//        }
//
//        append_element_part_names_and_assign_topology(tParentTopology,aOutputOptions,tNumPhases,tElementPartNames,tAppendedElementPartNames,tPartTopologies);
//
//
//        /*
//         * Append and store all face rank mesh parts
//         */
//
//
//        /*
//         * Get node parts and add interface part
//         */
//        Cell<std::string> tNodePartNames;
//        tXTKMeshData.get_all_part_names(EntityRank::NODE, tNodePartNames);
//
//
//        Cell<std::string> tInterfaceNodeSetNames(mGeometryEngines.get_num_geometries());
//        for(moris::size_t i = 0; i<mGeometryEngines.get_num_geometries(); i++)
//        {
//            tInterfaceNodeSetNames(i) = tNodeBaseStr+aOutputOptions.mInterfaceAppendix+"_" + std::to_string(i);
//        }
//
//
//        Cell<std::string> tInterfaceSideSetNames;
//
//
//        // Scalar field names
//        Cell<std::string> tScalarNodeFieldNames = aOutputOptions.mRealNodeExternalFieldNames;
//        Cell<std::string> tIntElementScalarFieldNames = aOutputOptions.mIntElementExternalFieldNames;
//        Cell<std::string> tVectorFieldNames = {};
//
//        // Get nodal coordinates
//        moris::Matrix< moris::DDRMat >  tNodeCoords = mXTKMesh.get_all_node_coordinates_loc_inds();
//
//        // Get node to global map
//        moris::Matrix< moris::DDSTMat > tNodeLocaltoGlobal = mXTKMesh.get_local_to_global_map(EntityRank::NODE);
////
////        moris::mtk::Mesh*
////        tMeshData = aMeshBuilder.build_mesh_from_data(mModelDimension,
////                                                      tXTKElementBuckets,
////                                                      tXTKSideSets,
////                                                      tXTKNodeSets,
////                                                      tNodeCoords,
////                                                      tNodeLocaltoGlobal,
////                                                      tAppendedElementPartNames,
////                                                      tPartTopologies,
////                                                      tAppendedFacePartNames,
////                                                      tNodePartNames,
////                                                      tInterfaceNodeSetNames,
////                                                      tInterfaceSideSetNames,
////                                                      tScalarNodeFieldNames,
////                                                      tVectorFieldNames,
////                                                      tIntElementScalarFieldNames,
////                                                      tSensitivity,
////                                                      aOutputOptions.mInternalUseFlag);

//    }


//    void package_elements_into_buckets(Output_Options<moris::size_t> const & aOutputOptions,
//                                       Cell<mesh::Bucket< moris::size_t,moris::DDSTMat >> & aXTKBuckets)
//    {
//
//        moris::mtk::Mesh const & tXTKMeshData = mXTKMesh.get_mesh_data();
//        enum EntityRank tElementRank  = EntityRank::ELEMENT;
//        moris::size_t tElementId;
//        moris::size_t tElementIndex;
//        moris::size_t tChildMeshIndex;
//        Cell<moris::size_t>   tPartOrdinals;
//        Cell<std::string> tPartNames;
//        Cell<std::string> tAppendedPartNames;
//        Cell<Cell<std::string>> tAppendedNonInterfacePartNames;
//        Cell<Cell<std::string>> tAppendedInterfacePartNames;
//
//        moris::size_t tNumPhases  = mGeometryEngines.get_num_phases();
//        moris::size_t tNumBuckets = tXTKMeshData.get_num_buckets(tElementRank);
//
//       moris::Matrix< moris::DDSTMat > tPhaseIndex(1,1);
//       moris::Matrix< moris::DDSTMat > tNodeIndiceMat(1,1);
//       moris::Matrix< moris::DDSTMat > tElementsInBucket(1,1);
//
//           /*
//           * Iterate over element buckets and sort into XTK Buckets
//           */
//          for(moris::size_t iBuck = 0; iBuck<tNumBuckets; iBuck++)
//          {
//
//              tElementsInBucket = tXTKMeshData.get_entities_in_bucket_loc_index(iBuck,tElementRank);
//              moris::size_t tNumElementsInBuckets = tElementsInBucket.n_cols();
//              tXTKMeshData.get_entity_part_membership_ordinals(tElementsInBucket(0,0),tElementRank,tPartOrdinals);
//              tXTKMeshData.get_part_name_from_part_ordinals(tPartOrdinals,tPartNames);
//
//              // TODO: REMOVE THIS HACK (NOTE the unnamed_block... causes a truncation issue in STK)
//              if(tPartNames(0).compare("unnamed_block_id:_1_type:_hex8") == 0)
//              {
//
//                  std::cout<<"Match"<<std::endl;
//                  tPartNames(0) = "hex_bl_1";
//              }
//
//              /*
//               * Add appendixes to the part names
//               */
//              append_interface_element_parts(tNumPhases,aOutputOptions,tPartNames,tAppendedInterfacePartNames);
//              append_non_interface_parts(tNumPhases,aOutputOptions,tPartNames,tAppendedNonInterfacePartNames);
//
//
//
//              /*
//               * Setup a temporary bucket which takes the elements that do not have chidlren and sorts them into different phases
//               */
//              Cell<mesh::Bucket< moris::size_t,moris::DDSTMat >> tTemporaryBucketForNoChildrenElements(tNumPhases,mesh::Bucket< moris::size_t,moris::DDSTMat >(tNumElementsInBuckets,8,tAppendedNonInterfacePartNames.size()*tNumPhases,20));
//              for(moris::size_t iPhase = 0; iPhase<tNumPhases; iPhase++)
//              {
//                  if(aOutputOptions.output_phase(iPhase))
//                  {
//                      tTemporaryBucketForNoChildrenElements(iPhase).add_part_names(tAppendedNonInterfacePartNames(iPhase));
//                  }
//              }
//              /*
//               * Setup a temporary bucket for elements which do have children
//               * Assuming all are tet 4s
//               */
//
//              Cell<mesh::Bucket< moris::size_t,moris::DDSTMat >> tTemporaryBucketForChildrenElements(tNumPhases,mesh::Bucket< moris::size_t,moris::DDSTMat >(tNumElementsInBuckets*24*6,8,tAppendedInterfacePartNames.size()*tNumPhases,20));
//              /**
//               * Setup a temporary bucket which takes the elements that do not have children and sorts them into different phases
//               */
//              for(moris::size_t iPhase = 0; iPhase<tNumPhases; iPhase++)
//              {
//                  if(aOutputOptions.output_phase(iPhase))
//                  {
//                      tTemporaryBucketForChildrenElements(iPhase).add_part_names(tAppendedInterfacePartNames(iPhase));
//                  }
//              }
//
//
//
//              for(moris::size_t iElem = 0; iElem < tElementsInBucket.n_cols(); iElem++)
//              {
//                  tElementIndex = tElementsInBucket(0,iElem);
//
//                  if(mXTKMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT))
//                  {
//
//                      tChildMeshIndex = mXTKMesh.child_mesh_index(tElementIndex,EntityRank::ELEMENT);
//
//                      Cell<moris::Matrix< moris::DDSTMat >> tElementCMInds;
//                      Cell<moris::Matrix< moris::DDSTMat >> tElementIds;
//                      Cell<std::string> tPhaseNames;
//                      mCutMesh.pack_cut_mesh_by_phase(tChildMeshIndex, tNumPhases, tElementCMInds,tElementIds);
//
//                      moris::Matrix< moris::DDSTMat > tElementToNode = mCutMesh.get_child_mesh(tChildMeshIndex).get_element_to_node_global();
//                      for(moris::size_t iPhase = 0; iPhase<tElementCMInds.size(); iPhase++ )
//                      {
//                          if(aOutputOptions.output_phase(iPhase))
//                          {
//                              tTemporaryBucketForChildrenElements(iPhase).add_entity_ids(tElementIds(iPhase));
//                              for(moris::size_t iE = 0; iE<tElementCMInds(iPhase).n_cols(); iE++)
//                              {
//                                  moris::Matrix< moris::DDSTMat > tSingleElemToNode = tElementToNode.get_row(tElementCMInds(iPhase)(0,iE));
//                                  tTemporaryBucketForChildrenElements(iPhase).add_entity(tSingleElemToNode);
//                              }
//                          }
//                      }
//
//                  }
//
//                  else
//                  {
//                      moris::Matrix< moris::DDSTMat > tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElementIndex,EntityRank::ELEMENT,EntityRank::NODE);
//
//                      tNodeIndiceMat(0,0) = tElementNodes(0,0);
//
//                      tElementId = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tElementIndex,EntityRank::ELEMENT);
//
//                      mGeometryEngines.get_phase_index(tNodeIndiceMat,tPhaseIndex);
//
//                      tElementNodes = tXTKMeshData.get_entity_connected_to_entity_glb_ids(tElementIndex,EntityRank::ELEMENT,EntityRank::NODE);
//
//                      if(tElementNodes.n_cols() == 4 && mConvertedToTet10s)
//                      {
//                          tElementNodes.resize(1,10);
//                          for(moris::size_t iMN = 4; iMN<10; iMN++)
//                          {
//                              tElementNodes(0,iMN) = mMidsideElementToNode(tElementIndex,iMN-4);
//                          }
//
//                      }
//
//
//                      if(aOutputOptions.output_phase(tPhaseIndex(0,0)))
//                      {
//                          tTemporaryBucketForNoChildrenElements(tPhaseIndex(0,0)).add_entity(tElementNodes);
//                          tTemporaryBucketForNoChildrenElements(tPhaseIndex(0,0)).add_entity_ids(tElementId);
//                      }
//                  }
//
//              }
//              /*
//               * Add temporary buckets to full bucket
//               */
//
//              aXTKBuckets.append(tTemporaryBucketForNoChildrenElements);
//              aXTKBuckets.append(tTemporaryBucketForChildrenElements);
//          }
//    }

    void append_face_parts(Output_Options    const & aOutputOptions,
                           moris::size_t     const & aNumPhases,
                           Cell<std::string> const & aPartNames,
                           Cell<std::string>       & aAppendedPartNames)
    {
        for(moris::size_t iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(moris::size_t iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedPartNames.push_back(aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase));
                    }
                }
        }
    }

    void append_element_part_names_and_assign_topology(enum CellTopology const & aParentElementTopology,
                                                       Output_Options const & aOutputOptions,
                                                       moris::size_t const & aNumPhases,
                                                       Cell<std::string> const & aPartNames,
                                                       Cell<std::string> & aAppendedPartNames,
                                                       Cell<enum CellTopology> & aPartTopologies)
    {
        aAppendedPartNames.reserve(aPartNames.size()*aNumPhases*2*20);
        aPartTopologies.reserve(aPartNames.size()*aNumPhases*2);

        /*
         * Append the non interface part names
         */
        for(moris::size_t iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(moris::size_t iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedPartNames.push_back(aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase));
                        aPartTopologies.push_back(aParentElementTopology);
                    }
                }
        }

        enum CellTopology tChildTopo = mCutMesh.get_child_element_topology();
        /*
         * Append the interface part names
         */
        for(moris::size_t iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(moris::size_t iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedPartNames.push_back(aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase)+aOutputOptions.mInterfaceAppendix);
                        aPartTopologies.push_back(tChildTopo);
                    }
                }
        }


        /*
         * Append the interface names
         */

        aAppendedPartNames.shrink_to_fit();
        aPartTopologies.shrink_to_fit();
    }

    void append_non_interface_parts(moris::size_t     const & aNumPhases,
                                    Output_Options    const & aOutputOptions,
                                    Cell<std::string> const & aPartNames,
                                    Cell<Cell<std::string>> & aAppendedNonInterfacePartNames)
    {
        aAppendedNonInterfacePartNames = Cell<Cell<std::string>>(aNumPhases, moris::Cell<std::string>(aPartNames.size()));
        for(moris::size_t iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(moris::size_t iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedNonInterfacePartNames(iPhase)(iNames) = aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase);
                    }
                }
        }
    }

    void append_interface_element_parts(moris::size_t     const & aNumPhases,
                                        Output_Options    const & aOutputOptions,
                                        Cell<std::string> const & aPartNames,
                                        Cell<Cell<std::string>> & aAppendedInterfacePartNames)
    {

        aAppendedInterfacePartNames.resize(aNumPhases, moris::Cell<std::string>(aPartNames.size()));
        for(moris::size_t iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(moris::size_t iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    aAppendedInterfacePartNames(iPhase)(iNames) = aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase)+aOutputOptions.mInterfaceAppendix;
                }
        }
    }

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
