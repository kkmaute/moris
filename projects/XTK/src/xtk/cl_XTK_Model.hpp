/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jul 2, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_MODEL_HPP_
#define SRC_XTK_CL_XTK_MODEL_HPP_

// Standard Include
#include <memory>
#include <limits>
#include <mpi.h>
#include <ctime>

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Tools.hpp"
#include "mesh/cl_Mesh_Bucket.hpp"
#include "mesh/cl_Mesh_Builder.hpp"
#include "mesh/cl_Mesh_Side_Set_Input.hpp"
#include "mesh/cl_Mesh_Node_Set_Input.hpp"
#include "mesh/cl_Mesh_Tools.hpp"

// XTKL: Geometry Engine Includes
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "geomeng/cl_MGE_Geometry_Object.hpp"

#include "assert/fn_xtk_assert.hpp"
#include "ios/cl_Logger.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Tools includes
#include "xtk/cl_XTK_Mesh.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Request_Handler.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "xtk/cl_XTK_Active_Process_Manager.hpp"
#include "xtk/cl_XTK_Sensitivity.hpp"

// XTKL: Linalg Includes
#include "linalg/cl_XTK_Matrix.hpp"

// Topology
//TODO: MOVE THESE WITH CUTTING METHODS SOMEWHERE ELSE
#include "topology/cl_XTK_Topology.hpp"
#include "topology/cl_XTK_Edge_Topology.hpp"
#include "topology/cl_XTK_Quad_4_Topology.hpp"
#include "topology/cl_XTK_Hexahedron_8_Topology.hpp"

#include "tools/fn_tet_volume.hpp"
#include "tools/cl_MPI_Tools.hpp"


namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Model
{
public:
    // Forward declare the maximum value of Integer and Real
    Real REAL_MAX = std::numeric_limits<Real>::max();
    Integer INTEGER_MAX = std::numeric_limits<Integer>::max();

    Model(){};

    /**
     * Primary constructor (this constructor is used for all cases except when testing something)
     */
    Model(Integer aModelDimension,
          std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> aMeshData,
          Geometry_Engine<Real, Integer, Real_Matrix, Integer_Matrix> & aGeometryEngine) :
           mSameMesh(false),
           mModelDimension(aModelDimension),
           mXTKMesh(aMeshData,aGeometryEngine),
           mCutMesh(mModelDimension),
           mGeometryEngines(aGeometryEngine),
           mConvertedToTet10s(false)
    {
        moris::Matrix< Real_Matrix > tNodeCoords = aMeshData->get_all_node_coordinates_loc_inds();
        mGeometryEngines.create_geometry_objects_for_background_mesh_nodes(tNodeCoords);
        mXTKMesh.initialize_interface_node_flags(aMeshData->get_num_entities(EntityRank::NODE),mGeometryEngines.get_num_geometries());
    }

    bool mSameMesh;

    ~Model(){}

    /**
     * @brief Decomposes a mesh using a geometry engine (split up across processors)
     * @param[in] aMethods - specify which type of subdivision method to use (this could be changed to command line parsing or XML reading)
     */

    void decompose(Cell<enum Subdivision_Method> aMethods,
                   bool aSetPhase  = true)
    {

        std::clock_t tTotalTime = std::clock();
        // Process for a decomposition
        Integer tNumDecompositions = aMethods.size();
        Integer tNumGeometries = mGeometryEngines.get_num_geometries();

        print_decompsition_preamble(aMethods);

        // Tell the subdivision to assign node Ids if it is the only subdivision method (critical for outputting)
        // This is usually only going to happen in test cases
        // Note: the Conformal subdivision methods dependent on node ids for subdivision routine, the node Ids are set regardless of the below boolean

        bool tNonconformingSubdivisionAssignNodeIds = false;
        if(aMethods.size() == 1)
        {
            tNonconformingSubdivisionAssignNodeIds = true;
        }

        // Loop over each geometry and have an active child mesh indices list for each
        for(Integer iGeom = 0; iGeom<tNumGeometries; iGeom++)
        {
            bool tFirstSubdivisionFlag = true;
            moris::Matrix< Integer_Matrix > tActiveChildMeshIndices(1,1,0);

            for (Integer iDecomp = 0; iDecomp < tNumDecompositions; iDecomp++)
            {
                // start timing on this decomposition
                std::clock_t start = std::clock();

                // Perform subdivision
                this->subdivide(aMethods(iDecomp),tActiveChildMeshIndices,tFirstSubdivisionFlag,tNonconformingSubdivisionAssignNodeIds);
                tFirstSubdivisionFlag = false;

                // print timing
                if(get_rank(get_comm()) == 0)
                {
                    std::cout<<"    "<<get_enum_str(aMethods(iDecomp))<<" completed in " <<(std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
                }
            }
            // If it's not the last geometry tell the geometry engine we're moving on
            if(iGeom!= tNumGeometries-1)
            {
                mGeometryEngines.advance_geometry_index();
            }
        }

        // Tell the xtk mesh to set all necessary information to finalize decomposition allowing
        // i.e set element ids, indices for children elements
        finalize_decomp_in_xtk_mesh(aSetPhase);

        if(get_rank(get_comm()) == 0)
        {
            std::cout<<"Decomposition completed in " <<(std::clock() - tTotalTime) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }

    /*
     * Convert Tet4 elements to Tet10 Elements
     */
    void
    convert_mesh_tet4_to_tet10()
    {
        mConvertedToTet10s = true;

        // start timing on this decomposition
        std::clock_t start = std::clock();
        convert_mesh_tet4_to_tet10_internal();

        if(get_rank(get_comm()) == 0)
        {
            std::cout<<"Tet4 to Tet10 conversion completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
    }

    /*
     * Returns the Cut Mesh
     */
    Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &
    get_cut_mesh()
    {
        return mCutMesh;
    }

    /*
     * Returns the Xtk Mesh
     */
    XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> &
    get_xtk_mesh()
    {
        return mXTKMesh;
    }


    /*
     * Get geomtry engine
     */

    Geometry_Engine<Real, Integer, Real_Matrix, Integer_Matrix> &
    get_geom_engine()
    {
        return mGeometryEngines;
    }

    /*
     * Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
     */
    std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>>
    get_output_mesh(mesh::Mesh_Builder<Real,Integer,Real_Matrix,Integer_Matrix> const & aMeshBuilder,
                    Output_Options<Integer> const &                                     aOutputOptions = Output_Options<Integer>())

    {
        // start timing on this decomposition
        std::clock_t start = std::clock();
        std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> tOutputMesh = construct_output_mesh(aMeshBuilder,aOutputOptions);

        if(get_rank(get_comm()) == 0)
        {
        std::cout<<"Mesh output completed in "<< (std::clock() - start) / (double)(CLOCKS_PER_SEC)<<" s."<<std::endl;
        }
        return tOutputMesh;
    }


private:
    Integer mModelDimension;
    XTK_Mesh<Real,Integer,Real_Matrix,Integer_Matrix>           mXTKMesh;
    Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix>        mCutMesh;
    Geometry_Engine<Real, Integer, Real_Matrix, Integer_Matrix> mGeometryEngines;



    // Member variables indicating the mesh has been converted to tet10s. The midside nodes are stored here currently but this may change
    bool mConvertedToTet10s;
    moris::Matrix< Integer_Matrix > mMidsideElementToNode;


    // Private Functions
private:
    // Decomposition Functions------------------------------------------------------
    /**
     * formulates node requests in the geometry objects. Dependent on the type of decomposition
     * @param[in] aReqType- specifies which template mesh is going to be used later on
     *
     *see cl_model.cpp for template on how to code a new request type
     */
    void subdivide(enum Subdivision_Method const & aRequestType,
                   moris::Matrix< Integer_Matrix > & aActiveChildMeshIndices,
                   bool const & aFirstSubdivision = true,
                   bool const & aSetIds = false)
    {
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = mXTKMesh.get_mesh_data();
        switch (aRequestType)
        {
        case (Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8):
        {
            bool tCoordFlag = false;
            XTK_ASSERT(tXTKMeshData.get_entity_connected_to_entity_loc_inds(0, EntityRank::ELEMENT, EntityRank::NODE).n_cols() == 8, "NC_REGULAR_SUBDIVISION_HEX8 is for HEX8 meshes only.");
            XTK_ASSERT(aFirstSubdivision,"NC_REGULAR_SUBDIVISION_HEX8 needs to be the first subdivision routine for each geometry");


            // Runs the first cut routine to get the new active child mesh indices and indicate which are new and need to be regularly subdivided and which ones dont
            moris::Matrix< Integer_Matrix > tNewPairBool;
            run_first_cut_routine(TemplateType::REGULAR_SUBDIVISION_HEX8, 8,  aActiveChildMeshIndices,tNewPairBool);


            // Initialize request list for faces and elements
            Integer tIntersectedCount = aActiveChildMeshIndices.n_cols();
            Integer tNumFaceRequests = 6; // (one per face on each element)
            Integer tNumElemRequests = 1; // (one internal element request)
            Integer tNumChildAllow = 1; // (one child allowed per parent entity for this subdivision)
            Integer tNumNewNodes = tNumFaceRequests + tNumElemRequests; //(7)
            Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tFaceRequests(tIntersectedCount * tNumFaceRequests, tNumChildAllow, EntityRank::FACE, EntityRank::NODE, tXTKMeshData, mCutMesh); // 6 face requests per element
            Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tElemRequests(tIntersectedCount * tNumElemRequests, tNumChildAllow, EntityRank::ELEMENT, EntityRank::NODE, tXTKMeshData, mCutMesh); // 1 face request per element

            // Initialize topologies used in this method
            Edge_Topology<Real, Integer, Real_Matrix, Integer_Matrix> tEdgeTopology;
            Quad_4_Topology<Real, Integer, Real_Matrix, Integer_Matrix> tFaceTopology;
            Hexahedron_8_Topology<Real, Integer, Real_Matrix, Integer_Matrix> tElementTopology;

            // Initialize a cell of pointers to future node index
            Cell<Integer*> tNodeInds(tNumNewNodes);
            moris::Matrix< Integer_Matrix > tFaceNodes(1,4);
            moris::Matrix< Integer_Matrix > tElementNodes(1,8);
            moris::Matrix< Integer_Matrix > tFaceIndices(1,6);
            moris::Matrix< Real_Matrix > tCoordinates(0,0);
            moris::Matrix< Real_Matrix > tNewNodeCoordinates(1,mModelDimension,0);
            moris::Matrix< Real_Matrix > tCenterFaceLocCoordinate(1,2,0.0);
            moris::Matrix< Real_Matrix > tCenterElementLocCoordinate(1,3,0.0);


            // Loop over xtk meshes and place a node request on each face and at center of element volume
            for (Integer i = 0; i < tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {

                // Get element index
                Integer tElemInd = mCutMesh.get_parent_element_index(aActiveChildMeshIndices(0,i));

                // Get local index of faces connected to element using local element index
                tFaceIndices = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, EntityRank::ELEMENT, EntityRank::FACE);

                // Loop over faces (6 in a hex 8) and set a node request.
                // Request will return a pointer to where the created node index will be placed
                for (Integer fi = 0; fi < 6; fi++)
                {
                    tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndices(0, fi), EntityRank::FACE, EntityRank::NODE);
                    tFaceTopology.set_node_indices(tFaceNodes);

                    tCoordinates = tXTKMeshData.get_selected_node_coordinates_loc_inds(tFaceNodes);
                    xtk::Interpolation::bilinear_interpolation(tCoordinates, tCenterFaceLocCoordinate,tNewNodeCoordinates);
                    tNodeInds(fi) = tFaceRequests.set_request_info(tFaceIndices(0, fi), tFaceTopology, tNewNodeCoordinates, tCenterFaceLocCoordinate);
                }

                // Place node at center of element
                tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElemInd, EntityRank::ELEMENT, EntityRank::NODE);
                tElementTopology.set_node_indices(tElementNodes);
                tCoordinates = tXTKMeshData.get_selected_node_coordinates_loc_inds(tElementNodes);
                xtk::Interpolation::trilinear_interpolation(tCoordinates, tCenterElementLocCoordinate, tNewNodeCoordinates);
                tNodeInds(6) = tElemRequests.set_request_info(tElemInd, tElementTopology, tNewNodeCoordinates, tCenterElementLocCoordinate);


                // Give XTK Mesh pointers to where its node indices will be located
                mCutMesh.set_pending_node_index_pointers(aActiveChildMeshIndices(0,i), tNodeInds);
                }
            }

            // Handle the requests (MPI communication here, happens in mesh communicate_mesh_info)
            tFaceRequests.handle_requests(tCoordFlag, mSameMesh, false, mGeometryEngines);
            tElemRequests.handle_requests(tCoordFlag, mSameMesh, false, mGeometryEngines);

            // Allocate interface flag space in XTK mesh
            mXTKMesh.allocate_space_in_interface_node_flags(tFaceRequests.get_num_requests(),mGeometryEngines.get_num_geometries());
            mXTKMesh.allocate_space_in_interface_node_flags(tElemRequests.get_num_requests(),mGeometryEngines.get_num_geometries());

            // Tell XTK mesh to retrieve pending node indices and generate tabulated meshes
            mCutMesh.retrieve_pending_node_inds();

            for(Integer i = 0; i< tIntersectedCount; i++)
            {
                if(tNewPairBool(0,i) == 0)
                {

                mCutMesh.generate_templated_mesh(aActiveChildMeshIndices(0,i),TemplateType::REGULAR_SUBDIVISION_HEX8);
                }
            }


            if(aSetIds)
            {
                moris::Matrix< Integer_Matrix > tNodeIds;
                for (Integer j = 0; j < tIntersectedCount; j++)
                    if(tNewPairBool(0,j) == 0)
                {
                    {
                        moris::Matrix< Integer_Matrix > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));
                        tNodeIds = mesh::Mesh_Helper::get_glb_entity_id_from_entity_loc_index_range(tXTKMeshData, tNodeIndices, EntityRank::NODE);
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
                moris::Matrix< Integer_Matrix > tNewPairBool;
                run_first_cut_routine(TemplateType::TET_4, 4, aActiveChildMeshIndices,tNewPairBool);

                for(Integer i = 0; i<aActiveChildMeshIndices.n_cols(); i++)
                {
                    Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,i));
                    tChildMesh.generate_connectivities(true,true,true);
                }

            }


            // Initialize
            Integer tEdgeInd = INTEGER_MAX;
            Integer tNumMesh = aActiveChildMeshIndices.n_cols();


            // Initialize topologies used in this method (all local coordinates are with respect to an edge)
            Edge_Topology<Real, Integer, Real_Matrix, Integer_Matrix> tEdgeTopology;

            moris::Matrix< Real_Matrix > tLocalCoord(1,1, 0);
            moris::Matrix< Real_Matrix > tEdgeCoords(2, 3, REAL_MAX);
            moris::Matrix< Real_Matrix > tGlobalCoord(1, 3, REAL_MAX);
            moris::Matrix< Real_Matrix > tCoordSwapper(1, 3, REAL_MAX);

            moris::Matrix< Integer_Matrix > tEdgeNodes(1, 2, INTEGER_MAX);
            moris::Matrix< Integer_Matrix > tParentInfo(1, 2, INTEGER_MAX);

            Cell<Geometry_Object<Real, Integer, Real_Matrix,Integer_Matrix>> tGeoObjects;

            // Initialize request list for faces and elements
            // Number of children allowed on a parent mesh entity
            Integer tEdgeChildren = 4*2*10;
            Integer tFaceChildren = 4*4*10;
            Integer tElemChildren = 4*24*10;

            // Of all the edges in a regular subdivision template 12 live on parent mesh edges, 24 live on parent mesh faces and 14 live in the parent mesh element
            // the number of parent edges and faces were done by hand and are hardcoded here.
            // TODO: Ask XTKMesh how many parents of each rank it has (This would eliminate the dependency on the regular subdivision coming first)
            Integer tNumParentEdge = 2*12;
            Integer tNumParentFace = 2*24;
            Integer tNumParentElem = 2*14;


            Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tEdgeRequests(tNumMesh * tNumParentEdge, tEdgeChildren, EntityRank::EDGE, EntityRank::NODE, tXTKMeshData, mCutMesh);
            Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tFaceRequests(tNumMesh * tNumParentFace, tFaceChildren, EntityRank::FACE, EntityRank::NODE, tXTKMeshData, mCutMesh);
            Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tElemRequests(tNumMesh * tNumParentElem, tElemChildren, EntityRank::ELEMENT, EntityRank::NODE, tXTKMeshData, mCutMesh);

            // Tell XTKMesh to initialize intersection connectivity
            mCutMesh.init_intersect_connectivity(aActiveChildMeshIndices);

            // Check type specified as conformal (could change this to enum)
            Integer tCheckType = 1;

            moris::Matrix< Real_Matrix > tNodeCoords = tXTKMeshData.get_all_node_coordinates_loc_inds();


            // Ask the geometry engine whether it has sensitivity information
            bool tHasDxDp = mGeometryEngines.mComputeDxDp;

            for (Integer j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {

                // Get the child mesh that is active
                Child_Mesh_Test<Real, Integer, Real_Matrix,Integer_Matrix> & tChildMesh = mCutMesh.get_child_mesh(aActiveChildMeshIndices(0,j));


                // Get full edge to element connectivity from XTK Mesh (probably slow)
                // 0 specifies XTK local indices if an analytic geometry
                // Otherwise this needs to be the processor local index
                // or even the ID
                moris::Matrix< Integer_Matrix > const & tEdgeToNode = tChildMesh.get_edge_to_node();

                // Ask geometry engine which edges are intersected (Simple mesh local indexed edges)
                mGeometryEngines.is_intersected(tNodeCoords, tEdgeToNode, tCheckType, tGeoObjects);

                // Initialize node index pointers based on number of intersected edges
                Cell<Integer*> tNodeInds(tGeoObjects.size());

                // get reference to child mesh edge parent information
                moris::Matrix< Integer_Matrix > const & tEdgeParentIndices = tChildMesh.get_edge_parent_inds();
                moris::Matrix< Integer_Matrix > const & tEdgeParentRanks   = tChildMesh.get_edge_parent_ranks();

                for (Integer k = 0; k < tGeoObjects.size(); k++)
                {
                    // Local index to XTK Mesh
                    tEdgeInd = tGeoObjects(k).get_parent_entity_index();

                    // get a local coordinate [-1,1]
                    tLocalCoord(0,0) = tGeoObjects(k).get_interface_lcl_coord();

                    tGlobalCoord = tGeoObjects(k).get_interface_glb_coord();

                    // Add edge to the entity auxiliary connectivity
                    mCutMesh.add_entity_to_intersect_connectivity(aActiveChildMeshIndices(0,j), k, tEdgeInd, 0);

                    tEdgeNodes = tEdgeToNode.get_row(tEdgeInd);
                    tCoordSwapper = tNodeCoords.get_row(tEdgeNodes(0, 0));
                    tEdgeCoords.set_row(0, tCoordSwapper);
                    tCoordSwapper = tNodeCoords.get_row(tEdgeNodes(0, 1));
                    tEdgeCoords.set_row(1, tCoordSwapper);


                    Integer tParentRank  = tEdgeParentRanks(0, tEdgeInd);
                    Integer tParentIndex = tEdgeParentIndices(0, tEdgeInd);

                    tEdgeTopology.set_node_indices(tEdgeNodes);

                    // Convert to global id using mesh
                    tEdgeNodes(0, 0) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                    tEdgeNodes(0, 1) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                    // Order the nodes in ascending order
                    if(tEdgeNodes(0, 1) < tEdgeNodes(0, 0))
                    {
                        Integer tSwap = tEdgeNodes(0, 0);
                        tEdgeNodes(0, 0) = tEdgeNodes(0, 1);
                        tEdgeNodes(0, 1) = tSwap;
                    }

                    if (tParentRank == 1)
                    {

                        // Intersected edge is an existing stk edge
                        // Make request in edge requests
                        // This does not require a supplemental identifier
                        // TODO: ADD OVERFLOW CHECK IN CANTOR PAIRING!!!!!!
                        Integer tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));

                        tNodeInds(k) = tEdgeRequests.set_request_info(tParentIndex,
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
                        Integer tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));

                        tNodeInds(k) = tFaceRequests.set_request_info(tParentIndex,
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
                        Integer tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));
                        tNodeInds(k) = tElemRequests.set_request_info(tParentIndex,
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
                        XTK_ERROR << "Invalid ancestry returned from XTK Mesh";
                    }
                } // geometry object

                mCutMesh.set_pending_node_index_pointers(aActiveChildMeshIndices(0,j), tNodeInds);

            } // XTK Mesh loop

            bool tCoordinateFlag = false;


            // handle requests
            tEdgeRequests.handle_requests(tCoordinateFlag, mSameMesh, true, mGeometryEngines);
            tFaceRequests.handle_requests(tCoordinateFlag, mSameMesh, true, mGeometryEngines);
            tElemRequests.handle_requests(tCoordinateFlag, mSameMesh, true, mGeometryEngines);

            // Allocate interface flag space in XTK mesh
            mXTKMesh.allocate_space_in_interface_node_flags(tEdgeRequests.get_num_requests(),mGeometryEngines.get_num_geometries());
            mXTKMesh.allocate_space_in_interface_node_flags(tFaceRequests.get_num_requests(),mGeometryEngines.get_num_geometries());
            mXTKMesh.allocate_space_in_interface_node_flags(tElemRequests.get_num_requests(),mGeometryEngines.get_num_geometries());

            // Mark these nodes as interface nodes with respect to the active geometry
            tEdgeRequests.mark_pending_nodes_as_interface_nodes(mXTKMesh,mGeometryEngines.get_active_geometry_index());
            tFaceRequests.mark_pending_nodes_as_interface_nodes(mXTKMesh,mGeometryEngines.get_active_geometry_index());
            tElemRequests.mark_pending_nodes_as_interface_nodes(mXTKMesh,mGeometryEngines.get_active_geometry_index());

            // Tell XTK mesh to grab the pending node indices and that these nodes are on the interface and have sensitivity
            mCutMesh.retrieve_pending_node_inds();

            // Set Node Ids and tell the child mesh to update
            for (Integer j = 0; j < aActiveChildMeshIndices.n_cols(); j++)
            {
                moris::Matrix< Integer_Matrix > const & tNodeIndices = mCutMesh.get_node_indices(aActiveChildMeshIndices(0,j));
                moris::Matrix< Integer_Matrix > tNodeIds = mesh::Mesh_Helper::get_glb_entity_id_from_entity_loc_index_range(tXTKMeshData, tNodeIndices, EntityRank::NODE);

                mCutMesh.set_node_ids(aActiveChildMeshIndices(0,j), tNodeIds);
                mCutMesh.modify_templated_mesh(aActiveChildMeshIndices(0,j), TemplateType::HIERARCHY_TET4);
            }



            break;
        }

        default:
        {
            Integer breaker = 0;
            XTK_ASSERT(breaker != 0, "formulate_node_request should not enter the default case, check to see if your aCheckType is undefined.");
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

        mesh::Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > const & tXTKMeshData = mXTKMesh.get_mesh_data();

        // Set child element ids and indices
        Integer tNumElementsInCutMesh = mCutMesh.get_num_entities(EntityRank::ELEMENT);

        // Allocate global element ids (these need to be give to the children meshes)
        Integer tElementIdOffset = tXTKMeshData.allocate_entity_ids(tNumElementsInCutMesh, EntityRank::ELEMENT);
        Integer tElementIndOffset = tXTKMeshData.get_num_entities(EntityRank::ELEMENT);
        for(Integer i = 0; i<mCutMesh.get_num_simple_meshes(); i++)
        {
            mCutMesh.set_child_element_ids(i,tElementIdOffset);
            mCutMesh.set_child_element_inds(i,tElementIndOffset);
        }


        // Compute the child element phase using the geometry engine
        if(aSetPhase)
        {
            // Set element phase indices
            mXTKMesh.initialize_element_phase_indices(tElementIndOffset);



            Integer tNumElem = tXTKMeshData.get_num_entities(EntityRank::ELEMENT);

            for(Integer i = 0; i<tNumElem; i++)
            {
                if(mXTKMesh.entity_has_children(i,EntityRank::ELEMENT))
                {
                    Integer tChildMeshIndex = mXTKMesh.child_mesh_index(i,EntityRank::ELEMENT);
                    Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                    moris::Matrix< Integer_Matrix > tElemToNode = tChildMesh.get_element_to_node();

                    moris::Matrix< Integer_Matrix > const & tElemInds  = tChildMesh.get_element_inds();


                    tChildMesh.initialize_element_phase_mat();

                    Integer tNumElem = tChildMesh.get_num_entities(EntityRank::ELEMENT);

                    for( Integer j = 0; j<tNumElem; j++)
                    {
                        Integer tElemPhaseIndex = determine_element_phase_index(j,tElemToNode);
                        mXTKMesh.set_element_phase_index(tElemInds(0,j),tElemPhaseIndex);
                        tChildMesh.set_element_phase_index(j,tElemPhaseIndex);
                    }

                }

                else
                {
                    moris::Matrix< Integer_Matrix > tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::NODE);

                    Integer tElemPhaseIndex = determine_element_phase_index(0,tElementNodes);

                    mXTKMesh.set_element_phase_index(i,tElemPhaseIndex);
                }


            }
        }

    }

    /*
     * Convert the child meshes into tet4's
     */

    void convert_mesh_tet4_to_tet10_internal()
    {

        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = mXTKMesh.get_mesh_data();

        // Make sure there is at lease one element in the mesh
        Integer tNumElems = tXTKMeshData.get_num_entities(EntityRank::ELEMENT);

        XTK_ASSERT(tNumElems>0," There needs to be at lease one element in the background mesh to convert to tet 10s");

        // Check background element type
        // If I start with a tet mesh, need to convert all children and background elements to tet10s
        // If I start with a hex mesh only convert children elements to tet10
        bool tTetStart = true;
        moris::Matrix< Integer_Matrix > tElem1Nodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(0,EntityRank::ELEMENT,EntityRank::NODE);
        Integer tNumNodesPerElem = tElem1Nodes.n_cols();
        XTK_ASSERT((tNumNodesPerElem == 4 || tNumNodesPerElem == 8), "Background needs to be tet4s or hex8s");

        if(tNumNodesPerElem == 8)
        {
            tTetStart = false;
        }


        // Number of Children Meshes
        Integer tNumChildrenMesh = mCutMesh.get_num_simple_meshes();
        Integer tTotalNumChildrenElem = mCutMesh.get_num_entities(EntityRank::ELEMENT);
        moris::Matrix< Real_Matrix >tNonInterfaceDxDp(3,1,0.0);

        // Initialize request list for faces and elements
        // Number of children allowed on a parent mesh entity
        Integer tEdgeChildren = 160;
        Integer tFaceChildren = 4*4*10;
        Integer tElemChildren = 4*24*10;

        // Of all the edges in a regular subdivision template 12 live on parent mesh edges, 24 live on parent mesh faces and 14 live in the parent mesh element
        // the number of parent edges and faces were done by hand and are hardcoded here.
        // TODO: Ask XTKMesh how many parents of each rank it has (This would eliminate the dependency on the regular subdivision coming first)
        Integer tNumParentEdge = 4*12;
        Integer tNumParentFace = 4*24;
        Integer tNumParentElem = 4*14;

        Integer tNumEdges = tXTKMeshData.get_num_entities(EntityRank::EDGE);

        Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tEdgeRequests(tNumChildrenMesh * tNumParentEdge + tNumEdges, tEdgeChildren, EntityRank::EDGE,    EntityRank::NODE, tXTKMeshData, mCutMesh);
        Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tFaceRequests(tNumChildrenMesh * tNumParentFace,             tEdgeChildren, EntityRank::FACE,    EntityRank::NODE, tXTKMeshData, mCutMesh);
        Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> tElemRequests(tNumChildrenMesh * tNumParentElem,             tEdgeChildren, EntityRank::ELEMENT, EntityRank::NODE, tXTKMeshData, mCutMesh);

        // Tell cut mesh to initialize auxiliary connectivity
        mCutMesh.init_intersect_connectivity();

        // Make all the children elements request nodes to convert them to

        //TODO: RESTRUCTURING MOVED ABILITY TO ASK CHILD MESH IF A NODE IS ON THE INTERFACE NEED THIS FUNCTION REENABLED
//        make_tet4_to_tet10_child_mesh_node_requests(tEdgeRequests,tFaceRequests,tElemRequests);



        // Make requests for all non-intersected elements
        Cell<Cell<Integer*>> tPendingNodes;
        if(tTetStart)
        {
            // Initialize midside node member variables
            mMidsideElementToNode(tNumElems,6);
            tPendingNodes =  make_tet4_to_tet10_unintersected_parent_mesh_node_requests(tEdgeRequests);
        }

        // Assign unique node ids to all requests
        bool tCoordinateFlag = false;

        tEdgeRequests.handle_requests(tCoordinateFlag, mSameMesh, false, mGeometryEngines);

        tFaceRequests.handle_requests(tCoordinateFlag, mSameMesh, false, mGeometryEngines);

        tElemRequests.handle_requests(tCoordinateFlag, mSameMesh, false, mGeometryEngines);

        // Tell XTK mesh to grab the pending node indices and that these nodes are on the interface and have sensitivity
        mCutMesh.retrieve_pending_node_inds();

        // Tell the cut mesh to have all of its children meshes generated tet10s from their tet4s
        mCutMesh.convert_cut_mesh_to_tet10s();

        // Commit midside node connectivity of unintersected parent elements
        if(tTetStart)
        {
            retrieve_midside_node_inds_unintersected(tPendingNodes);
        }

        // Set new node ids
        moris::Matrix< Integer_Matrix > tNodeIds;
        for (Integer j = 0; j < mCutMesh.get_num_simple_meshes(); j++)
        {
            moris::Matrix< Integer_Matrix > const & tNodeIndices = mCutMesh.get_node_indices(j);
            tNodeIds = mesh::Mesh_Helper::get_glb_entity_id_from_entity_loc_index_range(tXTKMeshData, tNodeIndices, EntityRank::NODE);

            mCutMesh.set_node_ids(j, tNodeIds);
        }

    }


    /*
     * Makes node requests at the midside for non-intersected elements, Note only edge requests are needed here because all edges
     * in the parent mesh cannot be from faces or elements
     */
    Cell<Cell<Integer*>>
    make_tet4_to_tet10_unintersected_parent_mesh_node_requests(Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> & aEdgeRequests)
    {
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = mXTKMesh.get_mesh_data();
        Integer tNumElem = tXTKMeshData.get_num_entities(EntityRank::ELEMENT);
        Integer tNumEdgePerElem  = 6;
        Integer tNumNodesPerElem = 10;
        Integer tNumChildMeshs = mCutMesh.get_num_simple_meshes();

        //
        moris::Matrix< Real_Matrix >       tEdgeCoords(2, 3, REAL_MAX);
        moris::Matrix< Real_Matrix >       tLocalCoord(1, 1, 0.0);
        moris::Matrix< Real_Matrix >       tNodeGlobCoord(1, 3, REAL_MAX);
        moris::Matrix< Real_Matrix >       tCoordSwapper(1, 3, REAL_MAX);
        moris::Matrix< Real_Matrix >       tNodeCoord(1,3,REAL_MAX);
        moris::Matrix< Integer_Matrix >  tEdgeNodes(1,2,INTEGER_MAX);
        moris::Matrix< Integer_Matrix >  tElementEdges(1,2,INTEGER_MAX);
        Edge_Topology<Real, Integer, Real_Matrix, Integer_Matrix> tEdgeTopology;

        // cell of pending indices
        Cell<Cell<Integer*>> tNodeInds(tNumElem,6);

        for(Integer i = 0; i<tNumElem; i++)
        {
            // If this entity does not have children
            if(!mXTKMesh.entity_has_children(i,EntityRank::ELEMENT))
            {

                // Get edges attached to element
                tElementEdges = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::EDGE);

                // Iterate through element edges and make a request
                for(Integer j = 0; j <tNumEdgePerElem; j++)
                {
                    // Get nodes attached to edge
                    tEdgeNodes =  tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElementEdges(0,j),EntityRank::EDGE,EntityRank::NODE);

                    // Compute Midisde Node Coordinate
                    tEdgeCoords = tXTKMeshData.get_selected_node_coordinates_loc_inds(tEdgeNodes);

                    // Interpolate coordinate on edge using node coordinates
                    // Must happen using local XTK index
                    //TODO: MOVE ALL INTERPOLATION TO GEOMETRY ENGINE
                    Interpolation::linear_interpolation_location(tEdgeCoords, tLocalCoord,tNodeGlobCoord);

                    // tEdge Nodes is now the processor local index (notice the 1 versus above)
                    tEdgeTopology.set_node_indices(tEdgeNodes);

                    // Convert to global id using mesh
                    (*tEdgeNodes)(0, 0) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 0), EntityRank::NODE);
                    (*tEdgeNodes)(0, 1) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tEdgeNodes(0, 1), EntityRank::NODE);

                    // Create a cantor pairing (since this is in conjunction with the child mesh a secondary identifier is needed)
                    // Although it does not do anything in this case
                    Integer tSecondaryId = xtk::cantor_pairing(tEdgeNodes(0, 0),tEdgeNodes(0, 1));

                    tNodeInds(i)(j) = aEdgeRequests.set_request_info(tElementEdges(0,j),tSecondaryId, tEdgeTopology, tNodeGlobCoord, tLocalCoord);

                }
            }
        }

        return tNodeInds;

    }

    /*
     * After all edges have been assigned a midside node id, this function stores those values
     */
    void
    retrieve_midside_node_inds_unintersected(Cell<Cell<Integer*>> aNodePointers)
    {
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = mXTKMesh.get_mesh_data();
        Integer tNumElems = mMidsideElementToNode.n_rows();
        Integer tNumMidSideNodes = mMidsideElementToNode.n_cols();

        for(Integer i = 0; i<tNumElems; i++)
        {
            if(!mXTKMesh.entity_has_children(i,EntityRank::ELEMENT))
            {
                for(Integer j = 0; j<tNumMidSideNodes; j++)
                {
                    mMidsideElementToNode(i,j) = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(*(aNodePointers(i)(j)),EntityRank::NODE);
                }
            }
        }
    }


    /*
     */
//    void
//    make_tet4_to_tet10_child_mesh_node_requests(Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> & aEdgeRequests,
//                                                Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> & aFaceRequests,
//                                                Request_Handler<Real, Integer, Real_Matrix, Integer_Matrix> & aElementRequests)
//    {
//
//        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = mXTKMesh.get_mesh_data();
//
//        // Number of Children Meshes
//        Integer tNumChildrenMesh = mCutMesh.get_num_simple_meshes();
//        Integer tNumChildrenElem = 0;
//        Integer tNumFacePerElem  = 4;
//        Integer tTotalNumChildrenElem = mCutMesh.get_num_entities(EntityRank::ELEMENT);
//        moris::Matrix< Integer_Matrix >     tElemToNodeConn     = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        moris::Matrix< Integer_Matrix >     tElemToEdgeConn     = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        moris::Matrix< Integer_Matrix >     tEdgeToNodeIndices  = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        moris::Matrix< Integer_Matrix >     tEdgeToNodeChildLoc = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//        std::shared_ptr<Matrix_Base<Real, Real_Matrix>>           tCoordSwapper       = mMatrixFactory.create_real_type_matrix_base(1, 3, REAL_MAX);
//        std::shared_ptr<Matrix_Base<Real, Real_Matrix>>           tEdgeCoords         = mMatrixFactory.create_real_type_matrix_base(2, 3, REAL_MAX);
//        std::shared_ptr<Matrix_Base<Real, Real_Matrix>>           tLocalCoord         = mMatrixFactory.create_real_type_matrix_base(1, 1, 0.0);
//        std::shared_ptr<Matrix_Base<Real, Real_Matrix>>           tNodeGlobCoord      = mMatrixFactory.create_real_type_matrix_base(1, 3, REAL_MAX);
//        Edge_Topology<Real, Integer, Real_Matrix, Integer_Matrix> tEdgeTopology;
//
//
//
//        //Number of elements in a children mesh
//        for(Integer i = 0; i<tNumChildrenMesh; i++)
//        {
//
//            // If it the first subdivision we need to find the intersected before placing the conformal nodes
//            // Intersected elements are flagged via the Geometry_Engine
//
//            // Initialize
//            Integer tEdgeInd = INTEGER_MAX;
//            Integer tCount = 0;
//
//            // Parent Information
//            moris::Matrix< Integer_Matrix > tParentInfo = mMatrixFactory.create_integer_type_matrix_base(1, 2, INTEGER_MAX);
//
//            // Get full edge to element connectivity from XTK Mesh (probably slow)
//            // 0 specifies XTK local indices if an analytic geometry
//            // Otherwise this needs to be the processor local index
//            // or even the ID
//            Integer tNumEdges = mCutMesh.get_num_entities(i,EntityRank::EDGE);
//
//            moris::Matrix< Integer_Matrix > tNodeIndices = mCutMesh.get_all_node_inds(i);
//            Matrix<Real, Real_Matrix>tNodeCoords        = tXTKMeshData.get_selected_node_coordinates_loc_inds(*tNodeIndices);
//
//            // Initialize node index pointers based on number of elements and a tet4 having 6 edges
//            tNumChildrenElem = tElemToNodeConn->n_rows();
//
//            Integer tTotal = tNumChildrenElem*tNumFacePerElem;
//
//            std::shared_ptr<Matrix_Base<Real, Real_Matrix>> tHalf = mMatrixFactory.create_real_type_matrix_base(1,1,0.5);
//
//
//
//
//            Cell<Integer*> tNodeInds(tNumEdges);
//
//            // Loop over all edges
//            for (Integer ed = 0; ed < tNumEdges; ed++)
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
//                std::shared_ptr<Matrix_Base<Real, Real_Matrix>>      tDxDp = mMatrixFactory.create_real_type_matrix_base(1,1);
//                moris::Matrix< Integer_Matrix > tNewNodeADVS;
//
//                if( tNode1Interface && tNode2Interface)
//                {
//                    if(mGeometryEngines.mComputeDxDp)
//                    {
//                        Matrix_Base<Real,Real_Matrix> * tNode1DxDp       = tNodes((*tEdgeToNodeChildLoc)(0,0)).get_dx_dp();
//                        Matrix_Base<Real,Real_Matrix> * tNode2DxDp       = tNodes((*tEdgeToNodeChildLoc)(0,1)).get_dx_dp();
//                        moris::Matrix< Integer_Matrix > * tNode1ADVS = tNodes((*tEdgeToNodeChildLoc)(0,0)).get_adv_indices();
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
//                    Matrix_Base<Real,Real_Matrix> * tNode1DxDp =  tNodes((*tEdgeToNodeChildLoc)(0,0)).get_dx_dp();
//
//                    tDxDp->matrix_data() = 0.5*(tNode1DxDp->matrix_data());
//                    moris::Matrix< Integer_Matrix > * tNode1ADVS = tNodes((*tEdgeToNodeChildLoc)(0,0)).get_adv_indices();
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
//                        Matrix_Base<Real,Real_Matrix> * tNode2DxDp =  tNodes((*tEdgeToNodeChildLoc)(0,1)).get_dx_dp();
//                        tDxDp->matrix_data() = 0.5*(tNode2DxDp->matrix_data());
//                        moris::Matrix< Integer_Matrix > * tNode2ADVS = tNodes((*tEdgeToNodeChildLoc)(0,1)).get_adv_indices();
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
//                Integer tParentRank = (*tParentInfo)(0, 0);
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
//                    Integer tSwap = (*tEdgeToNodeIndices)(0, 0);
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
//                    Integer tSecondaryId = xtk::cantor_pairing((*tEdgeToNodeIndices)(0, 0),(*tEdgeToNodeIndices)(0, 1));
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
//                    Integer tSecondaryId = xtk::cantor_pairing((*tEdgeToNodeIndices)(0, 0),(*tEdgeToNodeIndices)(0, 1));
//
//                    tNodeInds(ed) = aFaceRequests.set_request_info((*tParentInfo)(0, 1), tSecondaryId, tEdgeTopology, *tNodeGlobCoord, *tLocalCoord);
//                }
//                //
//                else if (tParentRank == 3)
//                {
//                    // Intersected edge was built in stk element
//                    // Make request in element requests
//                    // This requires a supplemental identifier
//                    Integer tSecondaryId = xtk::cantor_pairing((*tEdgeToNodeIndices)(0, 0),(*tEdgeToNodeIndices)(0, 1));
//                    tNodeInds(ed) = aElementRequests.set_request_info((*tParentInfo)(0, 1), tSecondaryId, tEdgeTopology, *tNodeGlobCoord, *tLocalCoord);
//                }
//
//                else
//                {
//                    XTK_ERROR << "Invalid ancestry returned from XTK Mesh";
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
//                    Matrix_Base<Real,Real_Matrix> *       tDxDpLocation;
//                    moris::Matrix< Integer_Matrix > * tNodeADVLocation;
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
     * Tells the XTK mesh about where it's children live in the cut mesh
     */
    void set_downward_inheritance()
    {
        Integer tNumChildMesh = mCutMesh.get_num_simple_meshes();
        Cell<std::pair<Integer,Integer>> tXTKElementToCutMeshPairs(tNumChildMesh);

        for(Integer iMesh = 0; iMesh<tNumChildMesh; iMesh++)
        {
            tXTKElementToCutMeshPairs(iMesh) = std::pair<Integer,Integer> (mCutMesh.get_parent_element_index(iMesh),iMesh);
        }

        mXTKMesh.register_new_downward_inheritance(tXTKElementToCutMeshPairs);
    }


    /*
     * This algorithm sets up the active child mesh indices and registers new pairs in the downward inheritance
     */

    void  run_first_cut_routine(enum TemplateType const & aTemplateType,
                                Integer const & tNumNodesPerElement,
                                moris::Matrix< Integer_Matrix > & aActiveChildMeshIndices,
                                moris::Matrix< Integer_Matrix > & aNewPairBool)
    {
        // Note this method is independent of node ids for this reason XTK_Mesh is not given the node Ids during this subdivision
        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & tXTKMeshData = mXTKMesh.get_mesh_data();

        // Package up node to element connectivity
        Integer tParentElementIndex = INTEGER_MAX;
        Integer tNumElements = tXTKMeshData.get_num_entities(EntityRank::ELEMENT);
        moris::Matrix< Integer_Matrix > tNodetoElemConnInd (tNumElements, tNumNodesPerElement);
        moris::Matrix< Integer_Matrix > tNodetoElemConnRow (1, tNumNodesPerElement, INTEGER_MAX);
        moris::Matrix< Integer_Matrix > tEdgetoElemConnInd (1, 1, INTEGER_MAX);
        moris::Matrix< Integer_Matrix > tFacetoElemConnInd (1, 1, INTEGER_MAX);
        moris::Matrix< Integer_Matrix > tElementMat(1, 1, INTEGER_MAX);
        moris::Matrix< Integer_Matrix > tPlaceHolder(1, 1, INTEGER_MAX);

        // TODO: Nest this in a mesh function
        for (Integer i = 0; i < tNumElements; i++)
        {
            tNodetoElemConnRow = tXTKMeshData.get_entity_connected_to_entity_loc_inds(i, EntityRank::ELEMENT, EntityRank::NODE);
            tNodetoElemConnInd.set_row(i, tNodetoElemConnRow);
        }

        // Get the Node Coordinates
        moris::Matrix< Real_Matrix > tAllNodeCoords = tXTKMeshData.get_all_node_coordinates_loc_inds();

        // Intersected elements are flagged via the Geometry_Engine
        Cell<Geometry_Object<Real, Integer, Real_Matrix,Integer_Matrix>> tGeoObjects;
        mGeometryEngines.is_intersected(tAllNodeCoords, tNodetoElemConnInd, 0,tGeoObjects);

        // Count number intersected
        Integer tIntersectedCount = tGeoObjects.size();

        // Loop over and determine how many new meshes that need to be registered (Avoids dynamic allocation in the child mesh)
        // Also register active mesh pairs
        Cell<std::pair<Integer,Integer>> tNewChildElementPair;
        aNewPairBool = moris::Matrix< Integer_Matrix >(1,tIntersectedCount,0);
        tNewChildElementPair.reserve(tIntersectedCount);

        Integer tNumNewChildMeshes = 0;
        Integer tNewIndex = 0;
        aActiveChildMeshIndices.resize(1,tIntersectedCount);
        for (Integer j = 0; j < tIntersectedCount; j++)
        {
            tParentElementIndex = tGeoObjects(j).get_parent_entity_index();

            if(!mXTKMesh.entity_has_children(tParentElementIndex,EntityRank::ELEMENT))
            {
                tNewIndex = tNumNewChildMeshes+mCutMesh.get_num_simple_meshes();
                tNewChildElementPair.push_back( std::pair<Integer,Integer>(tParentElementIndex, tNewIndex));
                aActiveChildMeshIndices(0,j) = tNewIndex;
                tNumNewChildMeshes++;
            }

            else
            {
                aActiveChildMeshIndices(0,j) = mXTKMesh.child_mesh_index(tParentElementIndex,EntityRank::ELEMENT);
                aNewPairBool(0,j) = 1;
            }
        }


        // Add the downward pair to the mesh for all the newly created element pairs
        mXTKMesh.register_new_downward_inheritance(tNewChildElementPair);


        // Allocate space for more simple meshes in XTK mesh
        mCutMesh.inititalize_new_child_meshes(tNumNewChildMeshes, mModelDimension);

        for (Integer j = 0; j < tIntersectedCount; j++)
        {
            if(aNewPairBool(0,j) == 0)
            {
                tParentElementIndex = tGeoObjects(j).get_parent_entity_index();

                // Get information to provide ancestry
                // This could be replaced with a proper topology implementation that knows faces, edges based on parent element nodes
                tNodetoElemConnRow = tNodetoElemConnInd.get_row(tParentElementIndex);
                tEdgetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, EntityRank::ELEMENT, EntityRank::EDGE);
                tFacetoElemConnInd = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tParentElementIndex, EntityRank::ELEMENT, EntityRank::FACE);
                tElementMat(0,0) = tParentElementIndex;

                // Set parent element, nodes, and entity ancestry
                moris::Matrix< Integer_Matrix > tElemToNodeMat(tNodetoElemConnRow);
//                mCutMesh.set_node_index(aActiveChildMeshIndices(0,j), tElemToNodeMat);          // Set the node indexes (only the node indexes from parent element)

                Cell<moris::Matrix< Integer_Matrix >> tAncestorInformation = {tPlaceHolder, tEdgetoElemConnInd, tFacetoElemConnInd, tElementMat};
                mCutMesh.initialize_new_mesh_from_parent_element(aActiveChildMeshIndices(0,j), aTemplateType, tElemToNodeMat, tAncestorInformation);
            }
        }
    }


    std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>>
    construct_output_mesh( mesh::Mesh_Builder<Real,Integer,Real_Matrix,Integer_Matrix> const & aMeshBuilder,
                           Output_Options<Integer> const & aOutputOptions)
    {


        /*
         * Initialize and Allocate
         */
        mesh::Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > const & tXTKMeshData = mXTKMesh.get_mesh_data();
        enum EntityRank tElementRank  = EntityRank::ELEMENT;
        Integer tNumElementsInXTKMesh = tXTKMeshData.get_num_entities(tElementRank);
        Integer tNumChildrenMeshes    = mCutMesh.get_num_simple_meshes();
        Integer tNumElementsInCutMesh = mCutMesh.get_num_entities(tElementRank);

        Integer tNumPhases = mGeometryEngines.get_num_phases();


        // Package Elements into buckets
        Cell<mesh::Bucket< Integer,Integer_Matrix >> tXTKElementBuckets;
        // Not reserved assuming all children mesh is tetrahedrons and parent mesh has hex 8
        // Also worst case would be there is a different field name for every element
        tXTKElementBuckets.reserve((tNumElementsInXTKMesh-tNumChildrenMeshes)*8+tNumElementsInCutMesh*4);

        package_elements_into_buckets( aOutputOptions, tXTKElementBuckets);


        Cell<mesh::Side_Set_Input<Integer,Integer_Matrix>> tXTKSideSets;

        // -----------------------------------------------------------------------------------
        // Package Side Sets from Parent Mesh

        Cell<std::string> tFacePartNames;
        Cell<std::string> tAppendedFacePartNames;
        if(aOutputOptions.mAddSideSets)
        {
            package_parent_side_sets_for_mesh_input(tXTKSideSets,aOutputOptions);

            // Append parent face partspart names
            tXTKMeshData.get_all_part_names(EntityRank::FACE, tFacePartNames);
            tFacePartNames.resize(tFacePartNames.size()/2,"");
            this->append_face_parts(aOutputOptions,tNumPhases, tFacePartNames,tAppendedFacePartNames);
        }

        // -----------------------------------------------------------------------------------
        // Package Interface Side Sets

        // Append interface name to all side sets
        std::string       tIntefaceSideName;
        Cell<std::string> tInterfaceSideNames;
        for(Integer iPhase = 0; iPhase<tNumPhases; iPhase++)
        {
            if(aOutputOptions.output_phase(iPhase))
            {
                tIntefaceSideName = aOutputOptions.mMaterialAppendix + aOutputOptions.mInterfaceAppendix + std::to_string(iPhase);
                tInterfaceSideNames.push_back(tIntefaceSideName);
            }
        }

        // Add to face part names
        tAppendedFacePartNames.append(tInterfaceSideNames);


        package_interface_side_sets_for_mesh_input(aOutputOptions, tXTKSideSets, tInterfaceSideNames);

        Cell<mesh::Node_Set_Input<Real, Integer, Real_Matrix, Integer_Matrix>> tXTKNodeSets;

        if(aOutputOptions.mAddNodeSets)
        {
            package_node_sets_for_mesh_input(tXTKNodeSets);
        }


        /*
         * Package Interface Node Set
         */
        std::string tNodeBaseStr = "nodes";
        xtk::Sensitivity<Real, Integer, Real_Matrix, Integer_Matrix> tSensitivity;
        package_interface_node_sets_for_mesh_input(tXTKNodeSets,aOutputOptions.mInterfaceAppendix,tNodeBaseStr,tSensitivity,aOutputOptions);
        tSensitivity.commit_sensitivities();


        /*
         * Initialize and get all part names which have a primary rank of element
         */
        Cell<std::string> tElementPartNames;
        Cell<std::string> tAppendedElementPartNames;
        Cell<enum EntityTopology> tPartTopologies;
        tXTKMeshData.get_all_part_names(tElementRank, tElementPartNames);

        for(Integer i = 0; i<tElementPartNames.size(); i++)
        {
            if(tElementPartNames(i).compare("unnamed_block_id:_1_type:_hex8") == 0)
            {

                std::cout<<"Match check 2"<<std::endl;
                tElementPartNames(i) = "hex_bl_1";
            }
        }

        enum EntityTopology tParentTopology = mXTKMesh.get_XTK_mesh_element_topology();
        if(tParentTopology == EntityTopology::TET_4 && mConvertedToTet10s)
        {
            tParentTopology = EntityTopology::TET_10;
        }

        append_element_part_names_and_assign_topology(tParentTopology,aOutputOptions,tNumPhases,tElementPartNames,tAppendedElementPartNames,tPartTopologies);


        /*
         * Append and store all face rank mesh parts
         */


        /*
         * Get node parts and add interface part
         */
        Cell<std::string> tNodePartNames;
        tXTKMeshData.get_all_part_names(EntityRank::NODE, tNodePartNames);


        Cell<std::string> tInterfaceNodeSetNames(mGeometryEngines.get_num_geometries());
        for(Integer i = 0; i<mGeometryEngines.get_num_geometries(); i++)
        {
            tInterfaceNodeSetNames(i) = tNodeBaseStr+aOutputOptions.mInterfaceAppendix+"_" + std::to_string(i);
        }


        Cell<std::string> tInterfaceSideSetNames;


        // Scalar field names
        Cell<std::string> tScalarNodeFieldNames = aOutputOptions.mRealNodeExternalFieldNames;
        Cell<std::string> tIntElementScalarFieldNames = aOutputOptions.mIntElementExternalFieldNames;
        Cell<std::string> tVectorFieldNames = {};

        // Get nodal coordinates
        moris::Matrix< Real_Matrix >  tNodeCoords = tXTKMeshData.get_all_node_coordinates_loc_inds();

        // Get node to global map
        moris::Matrix< Integer_Matrix > const & tNodeLocaltoGlobal = tXTKMeshData.get_local_to_global_map(EntityRank::NODE);


        std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>>
        tMeshData = aMeshBuilder.build_mesh_from_data(mModelDimension,
                                                      tXTKElementBuckets,
                                                      tXTKSideSets,
                                                      tXTKNodeSets,
                                                      tNodeCoords,
                                                      tNodeLocaltoGlobal,
                                                      tAppendedElementPartNames,
                                                      tPartTopologies,
                                                      tAppendedFacePartNames,
                                                      tNodePartNames,
                                                      tInterfaceNodeSetNames,
                                                      tInterfaceSideSetNames,
                                                      tScalarNodeFieldNames,
                                                      tVectorFieldNames,
                                                      tIntElementScalarFieldNames,
                                                      tSensitivity,
                                                      aOutputOptions.mInternalUseFlag);

        return tMeshData;
    }


    void package_elements_into_buckets(Output_Options<Integer> const & aOutputOptions,
                                       Cell<mesh::Bucket< Integer,Integer_Matrix >> & aXTKBuckets)
    {

        mesh::Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > const & tXTKMeshData = mXTKMesh.get_mesh_data();
        enum EntityRank tElementRank  = EntityRank::ELEMENT;
        Integer tElementId;
        Integer tElementIndex;
        Integer tChildMeshIndex;
        Cell<Integer>   tPartOrdinals;
        Cell<std::string> tPartNames;
        Cell<std::string> tAppendedPartNames;
        Cell<Cell<std::string>> tAppendedNonInterfacePartNames;
        Cell<Cell<std::string>> tAppendedInterfacePartNames;

        Integer tNumPhases  = mGeometryEngines.get_num_phases();
        Integer tNumBuckets = tXTKMeshData.get_num_buckets(tElementRank);

       moris::Matrix< Integer_Matrix > tPhaseIndex(1,1);
       moris::Matrix< Integer_Matrix > tNodeIndiceMat(1,1);
       moris::Matrix< Integer_Matrix > tElementsInBucket(1,1);

           /*
           * Iterate over element buckets and sort into XTK Buckets
           */
          for(Integer iBuck = 0; iBuck<tNumBuckets; iBuck++)
          {

              tElementsInBucket = tXTKMeshData.get_entities_in_bucket_loc_index(iBuck,tElementRank);
              Integer tNumElementsInBuckets = tElementsInBucket.n_cols();
              tXTKMeshData.get_entity_part_membership_ordinals(tElementsInBucket(0,0),tElementRank,tPartOrdinals);
              tXTKMeshData.get_part_name_from_part_ordinals(tPartOrdinals,tPartNames);

              // TODO: REMOVE THIS HACK (NOTE the unnamed_block... causes a truncation issue in STK)
              if(tPartNames(0).compare("unnamed_block_id:_1_type:_hex8") == 0)
              {

                  std::cout<<"Match"<<std::endl;
                  tPartNames(0) = "hex_bl_1";
              }

              /*
               * Add appendixes to the part names
               */
              append_interface_element_parts(tNumPhases,aOutputOptions,tPartNames,tAppendedInterfacePartNames);
              append_non_interface_parts(tNumPhases,aOutputOptions,tPartNames,tAppendedNonInterfacePartNames);



              /*
               * Setup a temporary bucket which takes the elements that do not have chidlren and sorts them into different phases
               */
              Cell<mesh::Bucket< Integer,Integer_Matrix >> tTemporaryBucketForNoChildrenElements(tNumPhases,mesh::Bucket< Integer,Integer_Matrix >(tNumElementsInBuckets,8,tAppendedNonInterfacePartNames.size()*tNumPhases,20));
              for(Integer iPhase = 0; iPhase<tNumPhases; iPhase++)
              {
                  if(aOutputOptions.output_phase(iPhase))
                  {
                      tTemporaryBucketForNoChildrenElements(iPhase).add_part_names(tAppendedNonInterfacePartNames(iPhase));
                  }
              }
              /*
               * Setup a temporary bucket for elements which do have children
               * Assuming all are tet 4s
               */

              Cell<mesh::Bucket< Integer,Integer_Matrix >> tTemporaryBucketForChildrenElements(tNumPhases,mesh::Bucket< Integer,Integer_Matrix >(tNumElementsInBuckets*24*6,8,tAppendedInterfacePartNames.size()*tNumPhases,20));
              /**
               * Setup a temporary bucket which takes the elements that do not have children and sorts them into different phases
               */
              for(Integer iPhase = 0; iPhase<tNumPhases; iPhase++)
              {
                  if(aOutputOptions.output_phase(iPhase))
                  {
                      tTemporaryBucketForChildrenElements(iPhase).add_part_names(tAppendedInterfacePartNames(iPhase));
                  }
              }



              for(Integer iElem = 0; iElem < tElementsInBucket.n_cols(); iElem++)
              {
                  tElementIndex = tElementsInBucket(0,iElem);

                  if(mXTKMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT))
                  {

                      tChildMeshIndex = mXTKMesh.child_mesh_index(tElementIndex,EntityRank::ELEMENT);

                      Cell<moris::Matrix< Integer_Matrix >> tElementCMInds;
                      Cell<moris::Matrix< Integer_Matrix >> tElementIds;
                      Cell<std::string> tPhaseNames;
                      mCutMesh.pack_cut_mesh_by_phase(tChildMeshIndex, tNumPhases, tElementCMInds,tElementIds);

                      moris::Matrix< Integer_Matrix > tElementToNode = mCutMesh.get_child_mesh(tChildMeshIndex).get_element_to_node_global();
                      for(Integer iPhase = 0; iPhase<tElementCMInds.size(); iPhase++ )
                      {
                          if(aOutputOptions.output_phase(iPhase))
                          {
                              tTemporaryBucketForChildrenElements(iPhase).add_entity_ids(tElementIds(iPhase));
                              for(Integer iE = 0; iE<tElementCMInds(iPhase).n_cols(); iE++)
                              {
                                  moris::Matrix< Integer_Matrix > tSingleElemToNode = tElementToNode.get_row(tElementCMInds(iPhase)(0,iE));
                                  tTemporaryBucketForChildrenElements(iPhase).add_entity(tSingleElemToNode);
                              }
                          }
                      }

                  }

                  else
                  {
                      moris::Matrix< Integer_Matrix > tElementNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tElementIndex,EntityRank::ELEMENT,EntityRank::NODE);

                      tNodeIndiceMat(0,0) = tElementNodes(0,0);

                      tElementId = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tElementIndex,EntityRank::ELEMENT);

                      mGeometryEngines.get_phase_index(tNodeIndiceMat,tPhaseIndex);

                      tElementNodes = tXTKMeshData.get_entity_connected_to_entity_glb_ids(tElementIndex,EntityRank::ELEMENT,EntityRank::NODE);

                      if(tElementNodes.n_cols() == 4 && mConvertedToTet10s)
                      {
                          tElementNodes.resize(1,10);
                          for(Integer iMN = 4; iMN<10; iMN++)
                          {
                              tElementNodes(0,iMN) = mMidsideElementToNode(tElementIndex,iMN-4);
                          }

                      }


                      if(aOutputOptions.output_phase(tPhaseIndex(0,0)))
                      {
                          tTemporaryBucketForNoChildrenElements(tPhaseIndex(0,0)).add_entity(tElementNodes);
                          tTemporaryBucketForNoChildrenElements(tPhaseIndex(0,0)).add_entity_ids(tElementId);
                      }
                  }

              }
              /*
               * Add temporary buckets to full bucket
               */

              aXTKBuckets.append(tTemporaryBucketForNoChildrenElements);
              aXTKBuckets.append(tTemporaryBucketForChildrenElements);
          }
    }

    void package_parent_side_sets_for_mesh_input(Cell<mesh::Side_Set_Input<Integer,Integer_Matrix>> & aMeshSideSets,
                                                 Output_Options<Integer> const & aOutputOptions)
    {
        mesh::Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > const & tXTKMeshData = mXTKMesh.get_mesh_data();

        Integer tNumPhases = mGeometryEngines.get_num_phases();
        bool tHasChildren;
        Integer tFaceIndex;
        Integer tElementIndex;
        Integer tElementId;
        Integer tChildMeshIndex;
        Integer tNumFacesInBucket;
        Integer tNumElementsAttachedToFace;
        Integer tPhaseIndex = 10;
        moris::Matrix< Integer_Matrix > tBucketFaces(1,1);
        moris::Matrix< Integer_Matrix > tElementsAttachedToFace(1,1);
        moris::Matrix< Integer_Matrix > tPhase(1,1);
        moris::Matrix< Integer_Matrix > tSingleNode(1,1);
        moris::Matrix< Integer_Matrix > tChildSideOrdinalofFace(1,1);
        moris::Matrix< Integer_Matrix > tChildElementsIndexConnectedToFace(1,1);
        moris::Matrix< Integer_Matrix > tChildElementsCMIndexConnectedToFace(1,1);

        Cell<std::string> tPartNames;
        Cell<Cell<std::string>> tAppendedPartNames;
        Cell<Integer> tPartOrdinals;

        Integer tNumSideBuckets = tXTKMeshData.get_num_buckets(EntityRank::FACE);

        aMeshSideSets.reserve(tNumSideBuckets);


        for(Integer iBuck = 0; iBuck<tNumSideBuckets; iBuck++)
        {
            Cell<Integer> tPhaseIndexInOutput(tNumPhases,INTEGER_MAX) ;
            tBucketFaces = tXTKMeshData.get_entities_in_bucket_loc_index(iBuck, EntityRank::FACE);

            if(tBucketFaces.n_cols() != 0)
            {
                tXTKMeshData.get_entity_part_membership_ordinals(tBucketFaces(0,0), EntityRank::FACE, tPartOrdinals);
                tXTKMeshData.get_part_name_from_part_ordinals(tPartOrdinals, tPartNames);

                if(tPartNames.size() != 0)
                {

                    // Append part names
                    this->append_non_interface_parts(tNumPhases, aOutputOptions, tPartNames, tAppendedPartNames);

                    for(Integer iPhase = 0; iPhase<tNumPhases; iPhase++)
                    {
                        if(aOutputOptions.output_phase(iPhase))
                        {
                            tPhaseIndexInOutput(iPhase) = aMeshSideSets.size();
                            aMeshSideSets.push_back(mesh::Side_Set_Input<Integer,Integer_Matrix>(tBucketFaces.n_cols(),30));

                            aMeshSideSets(tPhaseIndexInOutput(iPhase)).set_side_set_name(tAppendedPartNames(iPhase)(0));
                        }
                    }


                    /*
                     * Loop over faces and get elements attached to them
                     * then figure out which of these elements have children
                     */
                    tNumFacesInBucket = tBucketFaces.n_cols();
                    for(Integer iFace = 0; iFace<tNumFacesInBucket; iFace++)
                    {
                        tFaceIndex = tBucketFaces(0,iFace);
                        tElementsAttachedToFace = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndex,EntityRank::FACE,EntityRank::ELEMENT);

                        tNumElementsAttachedToFace = tElementsAttachedToFace.n_cols();
                        for( Integer iElem = 0; iElem<tNumElementsAttachedToFace; iElem++)
                        {

                            tElementIndex = tElementsAttachedToFace(0,iElem);
                            tHasChildren = mXTKMesh.entity_has_children(tElementIndex,EntityRank::ELEMENT);
                            // get the faces from the child mesh
                            if(tHasChildren)
                            {
                                tChildMeshIndex = mXTKMesh.child_mesh_index(tElementsAttachedToFace(0,iElem),EntityRank::ELEMENT);

                                Child_Mesh_Test< Real, Integer, Real_Matrix, Integer_Matrix > tChildMesh = mCutMesh.get_child_mesh(tChildMeshIndex);

                                tChildMesh.get_child_elements_connected_to_parent_face(tFaceIndex,
                                                                                       tChildElementsIndexConnectedToFace,
                                                                                       tChildElementsCMIndexConnectedToFace,
                                                                                       tChildSideOrdinalofFace);


                                // child mesh elemental phase vector
                                moris::Matrix< Integer_Matrix > const & tChildElementPhaseIndices = tChildMesh.get_element_phase_indices();
                                moris::Matrix< Integer_Matrix > const & tElementIds               = tChildMesh.get_element_ids();

                                // Get the child element phase
                                for(Integer iCElem  = 0; iCElem < tChildElementsCMIndexConnectedToFace.n_cols(); iCElem++)
                                {
                                    tPhaseIndex = tChildElementPhaseIndices(0,tChildElementsCMIndexConnectedToFace(0,iCElem));
                                    if(aOutputOptions.output_phase(tPhaseIndex))
                                    {

                                        // Child Element Id
                                        tElementId = tElementIds(0,tChildElementsCMIndexConnectedToFace(0,iCElem));
                                        aMeshSideSets(tPhaseIndexInOutput(tPhaseIndex)).add_element_id_and_side_ordinal(tElementId,
                                                                                                                        tChildSideOrdinalofFace(0,iCElem),
                                                                                                                        iCElem);
                                    }
                                }
                            }
                            // Get the face ordinal from the parent mesh
                            else
                            {
                                Integer tFaceOrdinal = tXTKMeshData.get_element_face_ordinal_loc_inds(tElementsAttachedToFace(0,iElem),tFaceIndex);
                                Integer tGlbElementId = tXTKMeshData.get_glb_entity_id_from_entity_loc_index(tElementsAttachedToFace(0,iElem), EntityRank::ELEMENT);
                                moris::Matrix< Integer_Matrix > tFaceNodes = tXTKMeshData.get_entity_connected_to_entity_loc_inds(tFaceIndex, EntityRank::FACE, EntityRank::NODE);

                                tSingleNode(0,0) = tFaceNodes(0,0);
                                mGeometryEngines.get_phase_index(tSingleNode,tPhase);

                                tPhaseIndex = tPhase(0,0);

                                if(aOutputOptions.output_phase(tPhaseIndex))
                                {
                                    aMeshSideSets(tPhaseIndexInOutput(tPhaseIndex)).add_element_id_and_side_ordinal(tGlbElementId,tFaceOrdinal,0);
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    void package_interface_side_sets_for_mesh_input(Output_Options<Integer> const & aOutputOptions,
                                                    Cell<mesh::Side_Set_Input<Integer,Integer_Matrix>> & aMeshSideSets,
                                                    Cell<std::string> const & aInterfaceSideNames)
    {

        mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> const & tXTKMeshData = mXTKMesh.get_mesh_data();
        Integer tNumChildMeshes = mCutMesh.get_num_simple_meshes();
        Integer tParentElementIndex = 0;
        Integer tElementOwner = 0 ;

        int tProcRank;
        MPI_Comm_rank(get_comm(), &tProcRank);


        moris::Matrix< Integer_Matrix > tElementIds(1,1);
        moris::Matrix< Integer_Matrix > tSideOrdinals(1,1);

        Integer tNumSides = tNumChildMeshes*35;

        // Pack into side sets
        mesh::Side_Set_Input<Integer,Integer_Matrix> tSideSets(tNumSides, 20);
        tSideSets.set_side_set_name(aInterfaceSideNames(0));

        for(Integer iElem = 0;  iElem<tNumChildMeshes; iElem++)
        {
            tParentElementIndex = mCutMesh.get_parent_element_index(iElem);
            tElementOwner = tXTKMeshData.get_entity_parallel_owner_rank(tParentElementIndex, EntityRank::ELEMENT);

            // See if I own the element
            if(tElementOwner == (Integer) tProcRank)
            {
                mCutMesh.pack_interface_sides(iElem,aOutputOptions,tElementIds,tSideOrdinals);
                tSideSets.add_element_id_and_side_ordinal(tElementIds, tSideOrdinals);
            }

        }

        aMeshSideSets.push_back(tSideSets);
    }

    void package_node_sets_for_mesh_input(Cell<mesh::Node_Set_Input<Real, Integer, Real_Matrix, Integer_Matrix>> & aMeshNodeSets)
    {
        mesh::Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > const & tXTKMeshData = mXTKMesh.get_mesh_data();
        Cell<Integer>     tPartOrdinals;
        Cell<std::string> tPartNames;
        enum EntityRank   tNodeRank  = EntityRank::NODE;


        Integer           tNumParts;
        Integer           tNumNodesInBucket;
        Integer           tMaxStringLength = 25;
        Integer           tNumBuckets = tXTKMeshData.get_num_buckets(tNodeRank);
        moris::Matrix< Integer_Matrix >  tNodesInBucketId;
        moris::Matrix< Integer_Matrix >  tNodesInBucketIndex(1,1);

        for(Integer iBuck = 0; iBuck<tNumBuckets; iBuck++)
        {
            // Get node indices and Ids
            tNodesInBucketIndex = tXTKMeshData.get_entities_in_bucket_loc_index(iBuck,tNodeRank);

            tNodesInBucketId = mesh::Mesh_Helper::get_glb_entity_id_from_entity_loc_index_range(tXTKMeshData, tNodesInBucketIndex, tNodeRank);

            // Get parts of bucket
            tXTKMeshData.get_entity_part_membership_ordinals(tNodesInBucketIndex(0,0),tNodeRank,tPartOrdinals);
            tXTKMeshData.get_part_name_from_part_ordinals(tPartOrdinals,tPartNames);

            // Initialize a temporary Node Set Input
            tNumParts = tPartNames.size();
            tNumNodesInBucket = tNodesInBucketIndex.n_cols();
            mesh::Node_Set_Input<Real, Integer, Real_Matrix, Integer_Matrix> tNodeSet( tNumNodesInBucket,tNumParts,tMaxStringLength );

            // Add the set names
            tNodeSet.add_node_set_names(tPartNames);

            // Add all nodes by ID
            for(Integer iNode = 0; iNode<tNumNodesInBucket; iNode++)
            {
                tNodeSet.add_node_id(tNodesInBucketId(0,iNode));
            }

            aMeshNodeSets.push_back(tNodeSet);
        }
    }

    void package_interface_node_sets_for_mesh_input(Cell<mesh::Node_Set_Input<Real, Integer, Real_Matrix, Integer_Matrix>> & aMeshNodeSets,
                                                    std::string const & aInterfaceAppendix,
                                                    std::string const & aInterfaceNodeBase,
                                                    Sensitivity<Real, Integer, Real_Matrix, Integer_Matrix> & aSensitivity,
                                                    Output_Options<Integer> const & aOutputOptions)
    {

        // Get interface node ids
        Cell<moris::Matrix< Integer_Matrix >> tInterfaceNodeIds = mXTKMesh.get_interface_nodes_glb_ids();

        // count interface nodes
        Integer tNumInterfaceNodes = 0;
        for(Integer iG = 0; iG<tInterfaceNodeIds.size(); iG++)
        {
            tNumInterfaceNodes = tNumInterfaceNodes + tInterfaceNodeIds(iG).n_cols();
        }

        Integer tNumDesignVars;
        if(mGeometryEngines.is_geometry_analytic())
        {
            tNumDesignVars = mGeometryEngines.get_num_design_vars_analytic();
        }

        else
        {
            tNumDesignVars = 2;
        }

        aSensitivity = Sensitivity<Real, Integer, Real_Matrix, Integer_Matrix> (tNumDesignVars,
                                                                                tNumInterfaceNodes,
                                                                                aOutputOptions.mDxDpName,
                                                                                aOutputOptions.mDxDpIndicesName,
                                                                                aOutputOptions.mDxDpNumIndicesName,
                                                                                aOutputOptions.mPackageDxDpSparsely);

        // Verify XTK mesh interface flag size matches the number of geometries in geometry engine
        XTK_ASSERT(tInterfaceNodeIds.size() == mGeometryEngines.get_num_geometries(),"Mismatch in number of geometries and expected geometries");

        for(Integer iG = 0; iG<mGeometryEngines.get_num_geometries(); iG++)
        {
            Integer tNumInterfaceNodes = tInterfaceNodeIds(iG).n_cols();

            mesh::Node_Set_Input<Real, Integer, Real_Matrix, Integer_Matrix>
            tNodeSet( tNumInterfaceNodes, 1, aInterfaceAppendix.length()+aInterfaceNodeBase.length()+ 3, mModelDimension, 0, 0, 0, 3);

            tNodeSet.add_node_set_name(aInterfaceNodeBase+aInterfaceAppendix + "_" + std::to_string(iG));


            // add interface node ids to node set
            tNodeSet.add_node_ids(tInterfaceNodeIds(iG));

            // append node set argument
            aMeshNodeSets.push_back(tNodeSet);
        }


        // FIXME:HARDCODED FOR 1 GEOMETRY
        if(mGeometryEngines.mComputeDxDp)
        {
            // iterate over interface nodes and get sensitivities
            moris::Matrix< Integer_Matrix > tInterfaceNodesInds = mXTKMesh.get_interface_nodes_loc_inds(0);
            for(Integer iS = 0; iS<tNumInterfaceNodes; iS++)
            {
                aSensitivity.add_node_sensitivity(tInterfaceNodeIds(0)(0,iS),
                                                  mGeometryEngines.get_node_adv_indices(tInterfaceNodesInds(0,iS)),
                                                  mGeometryEngines.get_node_dx_dp(tInterfaceNodesInds(0,iS)));
            }
        }

    }

    void append_face_parts(Output_Options<Integer> const & aOutputOptions,
                           Integer const & aNumPhases,
                           Cell<std::string> const & aPartNames,
                           Cell<std::string> & aAppendedPartNames)
    {
        for(Integer iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(Integer iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedPartNames.push_back(aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase));
                    }
                }
        }
    }

    void append_element_part_names_and_assign_topology(enum EntityTopology const & aParentElementTopology,
                                                       Output_Options<Integer> const & aOutputOptions,
                                                       Integer const & aNumPhases,
                                                       Cell<std::string> const & aPartNames,
                                                       Cell<std::string> & aAppendedPartNames,
                                                       Cell<enum EntityTopology> & aPartTopologies)
    {
        aAppendedPartNames.reserve(aPartNames.size()*aNumPhases*2*20);
        aPartTopologies.reserve(aPartNames.size()*aNumPhases*2);

        /*
         * Append the non interface part names
         */
        for(Integer iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(Integer iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedPartNames.push_back(aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase));
                        aPartTopologies.push_back(aParentElementTopology);
                    }
                }
        }

        enum EntityTopology tChildTopo = mCutMesh.get_child_element_topology();
        /*
         * Append the interface part names
         */
        for(Integer iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(Integer iPhase = 0; iPhase < aNumPhases; iPhase++)
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

    void append_non_interface_parts(Integer const & aNumPhases,
                                            Output_Options<Integer> const & aOutputOptions,
                                        Cell<std::string> const & aPartNames,
                                        Cell<Cell<std::string>> & aAppendedNonInterfacePartNames)
    {
        aAppendedNonInterfacePartNames = Cell<Cell<std::string>>(aNumPhases, xtk::Cell<std::string>(aPartNames.size()));
        for(Integer iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(Integer iPhase = 0; iPhase < aNumPhases; iPhase++)
                {
                    if(aOutputOptions.output_phase(iPhase))
                    {
                        aAppendedNonInterfacePartNames(iPhase)(iNames) = aPartNames(iNames)+aOutputOptions.mMaterialAppendix+std::to_string(iPhase);
                    }
                }
        }
    }

    void append_interface_element_parts(Integer const &           aNumPhases,
                                        Output_Options<Integer> const & aOutputOptions,
                                        Cell<std::string> const & aPartNames,
                                        Cell<Cell<std::string>> & aAppendedInterfacePartNames)
    {

        aAppendedInterfacePartNames.resize(aNumPhases, xtk::Cell<std::string>(aPartNames.size()));
        for(Integer iNames = 0; iNames<aPartNames.size(); iNames++  )
        {
                for(Integer iPhase = 0; iPhase < aNumPhases; iPhase++)
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
        if(get_rank(get_comm()) == 0)
        {
            for(Integer i = 0 ; i<aMethods.size(); i++)
            {
                std::cout<<"Decomposition Routine ("<<i<<"/"<<aMethods.size()-1<<"): "<<get_enum_str(aMethods(i))<<std::endl;

            }

            std::cout<<"Number of Geometries: " << mGeometryEngines.get_num_geometries()<<std::endl;
        }
    }


    Integer
    determine_element_phase_index(Integer aRowIndex,
                                  moris::Matrix< Integer_Matrix > const & aElementToNodeIndex)
    {
        Integer tNumGeom = mGeometryEngines.get_num_geometries();
        Integer tNumNodesPerElem = aElementToNodeIndex.n_cols();
        moris::Matrix< Integer_Matrix > tNodalPhaseVals(1,tNumGeom,INTEGER_MAX);

        for(Integer i = 0; i<tNumGeom; i++)
        {
            bool tFoundNonInterfaceNode = false;
            for( Integer j = 0; j<tNumNodesPerElem; j++)
            {
                if(!mXTKMesh.is_interface_node(aElementToNodeIndex(aRowIndex,j),i))
                {
                    tNodalPhaseVals(0,i) = mGeometryEngines.get_node_phase_index_wrt_a_geometry(aElementToNodeIndex(aRowIndex,j),i);
                    tFoundNonInterfaceNode = true;

                }
            }

            if(!tFoundNonInterfaceNode)
            {
                std::cout<<"Did not find a non-interface node for this element"<<std::endl;
                tNodalPhaseVals(0,i) = 1;
            }
        }


        Integer tElemPhaseVal = mGeometryEngines.get_elem_phase_index(tNodalPhaseVals);

        return tElemPhaseVal;
    }


};
}

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
