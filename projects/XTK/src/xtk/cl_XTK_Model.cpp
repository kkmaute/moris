/*
 * cl_XTK_Model.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Model.hpp"
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
        std::cout<<"deleting enrichment"<<std::endl;
        delete mEnrichment;
    }
}

// ----------------------------------------------------------------------------------
// Decomposition Source code
// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// Sensitivity Source code
// ----------------------------------------------------------------------------------

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

/**
 * Take the interface faces and create collapsed prisms
 */
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
        // initialize enrichment
        mEnrichment = new Enrichment(mGeometryEngine.get_num_phases(),
                                 this,
                                 &mCutMesh,
                                 &mBackgroundMesh);

        // Set verbose flag to match XTK.
        mEnrichment->mVerbose = mVerbose;

        // perform the enrichment
        mEnrichment->perform_enrichment();
    }


    // ----------------------------------------------------------------------------------
    // Export mesh Source code
    // ----------------------------------------------------------------------------------

}
