#include "cl_Model.hpp"

namespace moris
{
namespace xtk
{
// ----------------------------------------------------------------------------
Model::Model(mesh*                aMesh,
             GeometryEngine*      aGeometryEngine,
             uint                 aModelDimension)
{
    mParentMesh      = aMesh;
    mGeometryEngines = aGeometryEngine;
    mModelDimension  = aModelDimension;

}
// ----------------------------------------------------------------------------
Model::Model(uint    aModelDimension)
{
    mModelDimension = aModelDimension;
}
// ----------------------------------------------------------------------------
Model::Model(mesh*   aMesh,
             uint    aModelDimension)
{
    mModelDimension = aModelDimension;
    mParentMesh = aMesh;
}
// ----------------------------------------------------------------------------
Model::~Model()
{

}

// ----------------------------------------------------------------------------
// Decomposition Functions
Cell<MeshXTK>
Model::decompose(enum  Decomposition    aMethod)
{
    // Process for a decomposition
    // Note: Meat of decomposition code is inside subdivide function
    // This function should only initialize the XTK mesh, call the subdivision routines then return the XTK Mesh
    switch(aMethod)
    {
        case(Decomposition::REGULAR_HIER):
                {

            // Initialize XTK Mesh (these are modified within subdivide)
            // need to be initialized here for scope
            Cell<MeshXTK> tXTKMesh;

            // Subdivide using regular subdivision for a hex 8
            subdivide(RequestType::NC_REGULAR_SUBDIVISION_HEX8,tXTKMesh);

            // Perform conformal subdivision
            subdivide(RequestType::C_HIERARCHY_TET4,tXTKMesh);

            return tXTKMesh;

            break;
                }

        default:
        {
            MORIS_LOG_ERROR<<"Invalid Decomposition specified";
            Cell<MeshXTK> tCells(0,MeshXTK());

            return tCells;
            break;
        }
    }
}
// ----------------------------------------------------------------------------
void
Model::regular_subdivision_test_interface(Cell<MeshXTK>    & aXTKMeshList)
{
    //    regular_subdivision(aXTKMeshList);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void
Model::subdivide(enum RequestType   aReqType,
                 Cell<MeshXTK>    & aXTKMesh)
{
    // Please order cases in the same manner as the enum list for RequestType

    // 1. Get connectivity of entities of interests and ask geometry engine for intersection information
    // 2. Add to or initialize XTK Meshes for the geometry objects (an XTK Mesh starts with one element from the parent mesh
    // 3. Initialize entity requests with number of children allowed per parent entitity
    // 4. Loop over XTK meshes and set necessary requests and index pointers in XTK Meshs
    // 5. Tell the entity request object to handle requests once everything has been looped over (underhood MPI)
    // 6. Retrieve entity node Ids and Indices
    // 7. Generate/Modify templated mesh

    switch (aReqType)
    {
        case(RequestType::NC_REGULAR_SUBDIVISION_HEX8):
                    {

            uint tIntersectedCount = 0;

            // Package up node to element connectivity
            uint tNumElem = mParentMesh->get_num_elems();
            Mat<uint> tNodetoElemConnInd(tNumElem,8,UINT_MAX);
            Mat<uint> tNodetoElemConnId(tNumElem,8,UINT_MAX);

            for(uint i = 0; i<tNumElem; i++)
            {
                tNodetoElemConnInd.row(i) = mParentMesh->get_entity_local_ids_connected_to_entity(i,EntityRank::ELEMENT,EntityRank::NODE).row(0);
                tNodetoElemConnId.row(i) = mParentMesh->get_node_ids_from_local_map(tNodetoElemConnInd.row(i)).row(0);
            }

            // Intersected elements are flagged via the GeometryEngine
            Cell<GeometryEngine::Object> tGeoObjects = mGeometryEngines->is_intersected(mParentMesh->get_all_nodes_coords(),tNodetoElemConnInd,0);

            // Count number intersected
            tIntersectedCount = tGeoObjects.size();

            // Initialize simple mesh
            aXTKMesh.resize(tGeoObjects.size(),mModelDimension);       // Initialize simple xfem mesh


            // Give XTK Meshes node IDs by looping over geometry objects and selecting the
            // node indices from tNodetoElemConn
            uint      tempind = UINT_MAX;
            Mat<uint> temp(1,8,UINT_MAX);
            Cell<Mat<uint>> tEntities(mModelDimension+1); //
            for(uint j = 0; j<tIntersectedCount;j++)
            {
                tempind = tGeoObjects(j).get_parent_entity_index();
                temp    = tNodetoElemConnInd.row(tempind);
                tEntities(1) = mParentMesh->get_entity_local_ids_connected_to_entity(tempind, EntityRank::ELEMENT, EntityRank::EDGE);
                tEntities(2) = mParentMesh->get_entity_local_ids_connected_to_entity(tempind, EntityRank::ELEMENT, EntityRank::FACE);

                // Set parent element, nodes, and entity ancestry
                aXTKMesh(j).set_parent_element_index(tempind);   // Set the element index
                aXTKMesh(j).set_node_index(temp);                // Set the node indexes (only the node indexes from parent element)
                aXTKMesh(j).set_entity_ancestry(TemplateType::REGULAR_SUBDIVISION_HEX8,tEntities);
            }

            // Free up geometry objects
            tNodetoElemConnInd.resize(0,{0});
            tGeoObjects.resize(0,GeometryEngine::Object());

            // Create unique node IDs
            // Initialize request list for faces and elements
            Model::RequestHandler tFaceRequests(tIntersectedCount*6, 1, EntityRank::FACE,   EntityRank::NODE, mParentMesh); // 6 face requests per element
            Model::RequestHandler tElemRequests(tIntersectedCount,   1, EntityRank::ELEMENT,EntityRank::NODE, mParentMesh); // 1 face request per element

            // Loop over xtk meshes and place a node request on each face and at center
            for(uint i = 0; i< tIntersectedCount; i++)
            {
                // Initialize a cell of pointers to future node index
                Cell<uint*> tNodeInds(7);

                // Get element index
                uint tElemInd = aXTKMesh(i).get_parent_element_index();

                // Get index of faces connected to element
                Mat<uint> tFaceIndex = mParentMesh->get_entity_local_ids_connected_to_entity(tElemInd,EntityRank::ELEMENT,EntityRank::FACE);

                for(uint fi = 0; fi<6;fi++)
                {
                    tNodeInds(fi) = tFaceRequests.set_request_info(tFaceIndex(0,fi),mParentMesh->interpolate_to_location_on_entity(EntityRank::FACE,tFaceIndex(0,fi),{{0,0}}));
                }

                // Place node at center of element
                tNodeInds(6) = tElemRequests.set_request_info(tElemInd,mParentMesh->interpolate_to_location_on_entity(EntityRank::ELEMENT,tElemInd,{{0,0,0}}));

                aXTKMesh(i).set_pending_node_index_pointers(tNodeInds);
            }


            // Handle the requests (MPI communication here, happens in mesh communicate_mesh_info)
            tFaceRequests.handle_requests();
            tElemRequests.handle_requests();

            // Loop over xtk meshes again to commit the node indexes and generate a regular subdivided mesh
            Mat<uint> tConn;
            for(uint i = 0; i< tIntersectedCount; i++)
            {
                // Tell xtk meshes to save the local node indexes
                aXTKMesh(i).get_pending_node_inds();
                // Generate a regularly subdivided mesh
                aXTKMesh(i).generate_templated_mesh(TemplateType::REGULAR_SUBDIVISION_HEX8);
            }
            break;
                    }
        case(RequestType::C_HIERARCHY_TET4):
                {
            // Initialize
            uint                         tNumMesh   = aXTKMesh.size();
            uint                         tEdgeInd   = UINT_MAX;
            uint                         tEdgeId    = UINT_MAX;
            Mat<uint>                    tEdgetoElemConn(50,2,UINT_MAX);
            Mat<real>                    tNodeCoords(15,3,UINT_MAX);
            Mat<uint>                    tParentInfo(1,2,UINT_MAX);
            Mat<uint>                    tEdgeNodes(1,2,UINT_MAX);
            Mat<real>                    tLocalCoord(1,3,UINT_MAX);
            Mat<real>                    tGlobalCoord(1,3,UINT_MAX);
            Mat<real>                    tEdgeCoords(2,3,UINT_MAX);
            Cell<GeometryEngine::Object> tGeoObjects;



            // Initialize request list for faces and elements
            uint tEdgeChildren = 1;  // Number of children allowed on a parent mesh entity
            uint tFaceChildren = 4;
            uint tElemChildren = 24;

            uint tNumParentEdge = 12;
            uint tNumParentFace = 24;
            uint tNumParentElem = 14;
            Model::RequestHandler tEdgeRequests(tNumMesh*tNumParentEdge, tEdgeChildren, EntityRank::EDGE,   EntityRank::NODE, mParentMesh);
            Model::RequestHandler tFaceRequests(tNumMesh*tNumParentFace, tFaceChildren, EntityRank::FACE,   EntityRank::NODE, mParentMesh);
            Model::RequestHandler tElemRequests(tNumMesh*tNumParentElem, tElemChildren, EntityRank::ELEMENT,EntityRank::NODE, mParentMesh);

            // Specify I want intersection information on intersection
            uint tCheckType = 1;
            for(uint j = 0; j<aXTKMesh.size(); j++)
            {
                // Initialize auxiliary connectivity (needed to modify templated mesh)
                aXTKMesh(j).init_aux_connectivity(EntityRank::ELEMENT, EntityRank::NODE, EntityRank::EDGE);

                // Get full edge to element connectivity from XTK Mesh (probably slow)
                tEdgetoElemConn = aXTKMesh(j).get_full_connectivity(EntityRank::EDGE,EntityRank::NODE,0);

                // Get relevant node coordinates
                tNodeCoords     = mParentMesh->get_selected_nodes_coords_lcl_ind(aXTKMesh(j).get_all_node_inds());

                // Ask geometry engine which edges are intersected
                tGeoObjects     = mGeometryEngines->is_intersected(tNodeCoords,tEdgetoElemConn,tCheckType);
                Cell<uint*> tNodeInds(tGeoObjects.size());
                for(uint k = 0; k<tGeoObjects.size();k++)
                {
                    // Local index to XTK Mesh
                    tEdgeInd    = tGeoObjects(k).get_parent_entity_index();

                    // get a local coordinate [-1,1]
                    tLocalCoord = tGeoObjects(k).get_interface_lcl_coord();

                    // Get information about parent in stk mesh
                    tParentInfo = aXTKMesh(j).get_parent_entity(EntityRank::EDGE,tEdgeInd);

                    MORIS_ASSERT(tLocalCoord.n_cols()==1,"For C_HIERARCHY_TET4 subdivision, there should only be a local coordinate to a location on an edge corresponding to a levelset value =0");
                    MORIS_ASSERT((tLocalCoord(0,0)<1) || (tLocalCoord(0,0)<-1),"Local coordinate out of bounds (check geometry engine where local coordinate is set)");

                    // Add edge to the entity auxiliary connectivity
                    aXTKMesh(j).add_entity_to_aux_connectivity(k,tEdgeInd,0);

                    if(tParentInfo(0,0) == 1)
                    {
                        // Intersected edge is an existing stk edge
                        // Make request in edge requests
                        // This does not require a supplemental Id
                        tNodeInds(k) = tEdgeRequests.set_request_info(tParentInfo(0,1),mParentMesh->interpolate_to_location_on_entity(EntityRank::EDGE,tParentInfo(0,1),tLocalCoord));
                    }

                    else if (tParentInfo(0,0) == 2)
                    {
                        // Intersected edge was built on an stk face
                        // Make request in face requests
                        // This requires a supplemental Id

                        // tEdgeNodes is the local xtk index
                        tEdgeNodes   = aXTKMesh(j).get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tEdgeInd,0);

                        // Interpolate coordinate on edge using node coordinates
                        // Must happen using local XTK index
                        tEdgeCoords.row(0) = tNodeCoords.row(tEdgeNodes(0,0));
                        tEdgeCoords.row(1) = tNodeCoords.row(tEdgeNodes(0,1));
                        tGlobalCoord = Interpolation::linear_interpolation_location(tEdgeCoords,tLocalCoord);

                        // tEdge Nodes is now the processor local index
                        tEdgeNodes   = aXTKMesh(j).get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tEdgeInd,1);

                        // Convert to global id using mesh
                        tEdgeNodes   = mParentMesh->get_node_ids_from_local_map(tEdgeNodes);

                        // Create edge ID with nodes using cantor pairing function which also can be unpaired via an inverse function
                        // Currently, requires that the first node is a corner node and the second is an internal
                        tEdgeId = Pairing::cantor_pairing(tEdgeNodes);

                        tNodeInds(k) = tFaceRequests.set_request_info(tParentInfo(0,1), tEdgeId, tGlobalCoord);
                    }

                    else if (tParentInfo(0,0) == 3)
                    {
                        // Intersected edge was built in stk element
                        // Make request in element requests
                        // This requires a supplemental Id

                        // tEdgeNodes is the local xtk index
                        tEdgeNodes   = aXTKMesh(j).get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tEdgeInd,0);

                        // Interpolate coordinate on edge using node coordinates
                        // Must happen using local XTK index
                        tEdgeCoords.row(0) = tNodeCoords.row(tEdgeNodes(0,0));
                        tEdgeCoords.row(1) = tNodeCoords.row(tEdgeNodes(0,1));
                        tGlobalCoord = Interpolation::linear_interpolation_location(tEdgeCoords,tLocalCoord);

                        // tEdge Nodes is now the processor local index
                        tEdgeNodes   = aXTKMesh(j).get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tEdgeInd,1);

                        // Convert to global id using mesh
                        tEdgeNodes   = mParentMesh->get_node_ids_from_local_map(tEdgeNodes);

                        // Create edge ID with nodes using cantor pairing function which also can be unpaired via an inverse function
                        // Currently, requires that the first node is a corner node and the second is an internal
                        tEdgeId = Pairing::cantor_pairing(tEdgeNodes);

                        tNodeInds(k) = tElemRequests.set_request_info(tParentInfo(0,1),tEdgeId, tGlobalCoord);
                    }

                    else
                    {
                        MORIS_LOG_ERROR<<"Invalid ancestry returned from XTK Mesh";
                    }
                } // geometry object

                aXTKMesh(j).set_pending_node_index_pointers(tNodeInds);


            } // XTK Mesh loop

            //        MORIS_LOG_INFO<<"Edge Requests: " <<tEdgeRequests.get_num_requests();
            //        MORIS_LOG_INFO<<"Face Requests: " <<tFaceRequests.get_num_requests();
            //        MORIS_LOG_INFO<<"Elem Requests: " <<tElemRequests.get_num_requests();


            //        MORIS_LOG_INFO<<"Handling edge requests";
            tEdgeRequests.handle_requests();
            //        MORIS_LOG_INFO<<"Handling face requests";
            tFaceRequests.handle_requests();
            //        MORIS_LOG_INFO<<"Handling elem requests";
            tElemRequests.handle_requests();

            Mat<uint> tNodeIds;
            Mat<uint> tConn;
            for(uint j = 0; j<aXTKMesh.size(); j++)
            {
                aXTKMesh(j).get_pending_node_inds();
                tNodeIds = mParentMesh->get_node_ids_from_local_map(aXTKMesh(j).get_all_node_inds());
                aXTKMesh(j).set_node_ids(tNodeIds);
                aXTKMesh(j).modify_templated_mesh(TemplateType::HIERARCHY_TET4);
            }

            break;
                }
        default:
        {
            uint breaker = 0;
            MORIS_ASSERT(breaker!=0,"formulate_node_request should not enter the default case, check to see if your aCheckType is undefined.");
        }
    }
}
// ----------------------------------------------------------------------------
Cell<Mat<uint>>
Model::pack_XTK_mesh(Cell<MeshXTK>   & aXTKMesh,
                     uint              aPackType)
{
    uint tProcSize = par_size();

    if(tProcSize==1)
    {
        Cell<Mat<uint>> tConnPackage(2);
        uint tNumMesh = aXTKMesh.size();
        uint tNumElem = 0;
        for(uint i = 0; i<tNumMesh;i++)
        {
            tNumElem += aXTKMesh(i).get_num_entities(EntityRank::ELEMENT);
        }

        Mat<uint> tConn(tNumElem,4);
        tConnPackage(1) = Mat<uint>(1,tNumElem,0);
        uint k = 0;
        for(uint i = 0; i<tNumMesh;i++)
        {
            Mat<uint> temp = aXTKMesh(i).get_full_connectivity(EntityRank::ELEMENT,EntityRank::NODE,aPackType);

            for(uint r = 0; r<temp.n_rows(); r++)
            {
                tConn.row(k) = temp.row(r);
                k++;
            }
        }

        tConnPackage(0) = tConn;
        return tConnPackage;
    }

    else
    {
        Mat<uint> tOwnerVect;
        uint tNumMesh = aXTKMesh.size();
        uint tNumElem = 0;
        for(uint i = 0; i<tNumMesh;i++)
        {
                tNumElem += aXTKMesh(i).get_num_entities(EntityRank::ELEMENT);
        }

        Cell<Mat<uint>> tConnPackage(2,Mat<uint>(tNumElem,4));
        Mat<uint> tOwners(1,tNumElem);
        uint tOwner = 0;
        uint k = 0;
        for(uint i = 0; i<tNumMesh;i++)
        {
            tOwner = mParentMesh->parallel_owner_rank_by_entity_index(aXTKMesh(i).get_parent_element_index(),EntityRank::ELEMENT);

            Mat<uint> temp = aXTKMesh(i).get_full_connectivity(EntityRank::ELEMENT,EntityRank::NODE,aPackType);


            for(uint r = 0; r<temp.n_rows(); r++)
            {
                tConnPackage(0).row(k) = temp.row(r);
                tOwners(0,k) = tOwner;
                k++;
            }
        }
        tConnPackage(1) = tOwners;
        return tConnPackage;
    }
}

void
Model::handle_node_request(Cell<MeshXTK>    & aXTKMeshList)
{

}
//-----------------------------------------------------------------------------
void
Model::handle_node_request_test_interface(Cell<MeshXTK>    & aXTKMeshList)
{
    handle_node_request(aXTKMeshList);
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void Model::print_send_communication_Ids(Cell<Cell<Cell<uint>>>  & aSendIds,
                                         Cell<uint>                            & aActiveSendProc)
{
    MORIS_LOG_INFO << "========================================";
    MORIS_LOG_INFO<< "The following processors are active receivers:";
    for(uint ap = 0; ap<aActiveSendProc.size();ap++)
    {
        MORIS_LOG_INFO<< aActiveSendProc(ap);
    }


    // loop over active procs

    Mat<uint> temp(1,2,0);

    for(uint ap = 0; ap<aActiveSendProc.size();ap++)
    {
        uint tCurrentProc = aActiveSendProc(ap);
        MORIS_LOG_INFO << "The following information will be sent to processor:"<< tCurrentProc <<std::endl;

        // Unpack the information
        for(uint ui = 0;  ui<aSendIds(ap).size(); ui++)
        {
            Mat<uint> temp(1,aSendIds(ap)(ui).size(),0);
            // Loop through innermost cell
            for(uint in = 0;  in <aSendIds(ap)(ui).size();in++)
            {
                // Print the unsigned integers at the end
                temp(in) = aSendIds(ap)(ui)(in);
            }

            //MORIS_LOG_INFO << temp;
        }
    }
    MORIS_LOG_INFO << "========================================";
}
//-----------------------------------------------------------------------------

void Model::print_recv_communication_Ids(Cell<Cell<Cell<uint>>>  & aRecIds,
                                         Cell<uint>                            & aActiveRecProc)
{
    MORIS_LOG_INFO << "========================================";
    MORIS_LOG_INFO<< "The following processors are active senders:";

    for(uint ap = 0; ap<aActiveRecProc.size();ap++)
    {
        MORIS_LOG_INFO<< aActiveRecProc(ap);
    }

    // loop over active procs

    Mat<uint> temp(1,2,0);

    for(uint ap = 0; ap<aActiveRecProc.size();ap++)
    {
        uint tCurrentProc = aActiveRecProc(ap);
        MORIS_LOG_INFO << "The following information will be received from processor:"<< tCurrentProc <<std::endl;
        MORIS_LOG_INFO << "    Index       EID"<< tCurrentProc <<std::endl;
        // Unpack the information
        for(uint ui = 0;  ui<aRecIds(ap).size(); ui++)
        {
            Mat<uint> temp(1,aRecIds(ap)(ui).size(),0);
            // Loop through innermost cell
            for(uint in = 0;  in <aRecIds(ap)(ui).size();in++)
            {
                // Print the unsigned integers at the end
                temp(in) = aRecIds(ap)(ui)(in);
            }

            //MORIS_LOG_INFO << temp;
        }
    }
    MORIS_LOG_INFO << "========================================";

}


//-----------------------------------------------------------------------------

//Enrichment Functions---------------------------------------------------------

//-----------------------------------------------------------------------------



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Entity Tracker Object
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Request Handler Object
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
uint*
Model::RequestHandler::set_request_info(uint              aParentEntityIndex,
                                        Mat<real>  aChildCoords)
{
    // Ask entity tracker if this parent entity index has been used yet
    bool tUse = mEntityTracker.is_parent_entity_used(aParentEntityIndex);

    MORIS_ASSERT(aChildCoords.n_rows()==1,"Coordinates submitted to set_request_info need to be one row (x,y,z)");


    if(tUse!=true)
    {
        MORIS_ASSERT(mRequestCounter<mRequestInfo.n_rows(),"Not enough space allocated by constructor, check aNumExpectedRequests");
        mEntityTracker.mark_entity_as_used(aParentEntityIndex);
        mRequestInfo(mRequestCounter,0) = aParentEntityIndex;
        mRequestCounter++;
    }

    // Get pointer to where the index will be put
    return mEntityTracker.get_index_pointer(aParentEntityIndex);
}


uint*
Model::RequestHandler::set_request_info(uint              aParentEntityIndex,
                                        uint              aSecondaryEntityIndex,
                                        Mat<real>  aChildCoords)
{
    // Ask entity tracker if this parent entity index has been used yet
    // For this function the flag is embedded in pIndId as Cell 0
    Cell<uint*> pIdInd = mEntityTracker.is_parent_entity_used(aParentEntityIndex,aSecondaryEntityIndex);
    MORIS_ASSERT(mNumChildrenAllowed!=1,"Coordinates submitted to set_request_info need to be one row (x,y,z)");


    if(pIdInd(0)==NULL)
    {
        MORIS_ASSERT(mRequestCounter<mRequestInfo.n_rows(),"Not enough space allocated by constructor, check aNumExpectedRequests");
        mRequestInfo(mRequestCounter,0) = aParentEntityIndex;
        mRequestInfo(mRequestCounter,1) = aSecondaryEntityIndex;
        mRequestCounter++;
    }

    // Get pointer to where the index will be put
    return pIdInd(2);
}

//-----------------------------------------------------------------------------

void Model::RequestHandler::handle_requests()
{
    //     All requests have been placed in the request manager when this function is called
    //     resize must come before sort otherwise matrix will put zeros first
    mRequestInfo.resize(mRequestCounter,2);

    // tell mesh to communicate requests (mesh internally places ids and index in the entity tracker
//    mMeshPtr->communicate_entity_requests(mRequestInfo, &mEntityTracker, mParentEntityRank, mChildEntityRank ,0);
//
//    mMeshPtr->batch_create_new_entity(mPendingNodes, mChildEntityRank);
}

void
Model::create_regular_subdivided_mesh_output(Cell<MeshXTK>  & aXTKMesh)
{
    //    //        Look at subdivided mesh
    //            factory tMeshFactory;
    //            Mat<uint> tNodetoElemConnect(aXTKMesh.size()*24,4,UINT_MAX);
    //            Mat<real> tCoords = mParentMesh->get_all_nodes_coords();
    //            for(uint j = 0; j<aXTKMesh.size(); j++)
    //            {
    //                for(uint k = 0; k<24; k++)
    //                {
    //                    tNodetoElemConnect.row(j*24+k) = aXTKMesh(j).get_entities_connected_to_entity(EntityRank::ELEMENT,EntityRank::NODE,k,1).row(0);
    //                }
    //            }
    //
    //            Cell< uint >   aElemProcs = { 0, 1 };
    //            Cell< std::string >  aPartNames;
    //            aPartNames.push_back("block_1");
    //            aPartNames.push_back("block_2");
    //
    //            mesh* tMesh = tMeshFactory.create_mesh( MeshType::MTK, 3, tNodetoElemConnect, tCoords, aElemProcs, aPartNames );
    //            std::string OutputFileName1 = "test/src/mesh/TestOuputMesh.e";
    //            tMesh->create_output_mesh(OutputFileName1);

    //            delete tMesh;
}

}
}
