/*
 * cl_XTK_Model.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: doble
 */

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Background_Mesh.hpp"
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

    void
    Model::perform_multilevel_enrichment_internal()
    {
        mEnrichment->create_multilevel_enrichments();
    }
    // ----------------------------------------------------------------------------------
    // Export mesh Source code
    // ----------------------------------------------------------------------------------
    moris::mtk::Mesh*
    Model::get_xtk_as_mtk()
    {
        return new moris::mtk::XTK_Impl(this);
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
                     if(tSideSetData(i).numel() != 0)
                     {
                         tMtkMeshSets.add_side_set(&tBackgroundSideSets(i));
                     }
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
