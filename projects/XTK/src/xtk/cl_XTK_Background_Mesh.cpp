/*
 * cl_XTK_Background_Mesh.cpp
 *
 *  Created on: Mar 14, 2019
 *      Author: doble
 */

#include "cl_XTK_Background_Mesh.hpp"

namespace xtk
{
// ----------------------------------------------------------------------------------
// Constructor/Deconstructor Source code
// ----------------------------------------------------------------------------------
Background_Mesh::Background_Mesh(moris::mtk::Mesh* aMeshData):
    mMeshData(aMeshData),
    mChildMtkCells(0),
    mXtkMtkVertices(0),
    mXtkMtkVerticesInterpolation(0),
    mNodeIndexToChildMeshIndex(0,0)
{
    intialize_downward_inheritance();
    mExternalMeshData.set_up_external_entity_data(mMeshData);
    initialize_background_mesh_vertices();
}

// ----------------------------------------------------------------------------------

Background_Mesh::Background_Mesh(moris::mtk::Mesh* aMeshData,
                                 Geometry_Engine & aGeometryEngine):
    mMeshData(aMeshData),
    mChildMtkCells(0),
    mXtkMtkVertices(0),
    mXtkMtkVerticesInterpolation(0),
    mNodeIndexToChildMeshIndex(0,0)
{
    intialize_downward_inheritance();
    mExternalMeshData.set_up_external_entity_data(mMeshData);
    initialize_background_mesh_vertices();
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------
// Entity related functions
// ----------------------------------------------------------------------------------

moris::size_t
Background_Mesh::get_num_entities(enum EntityRank aEntityRank) const
{
    // Initialize
    moris::size_t tNumBackgroundEntities = mMeshData->get_num_entities((moris::EntityRank)aEntityRank);
    moris::size_t tExternalEntities = mExternalMeshData.get_num_entities_external_data(aEntityRank);

    return tNumBackgroundEntities + tExternalEntities;
}

// ----------------------------------------------------------------------------------

moris::size_t
Background_Mesh::get_num_entities_background(enum EntityRank aEntityRank) const
{
    return mMeshData->get_num_entities((moris::EntityRank)aEntityRank);
}

// ----------------------------------------------------------------------------------

moris::mtk::Vertex &
Background_Mesh::get_mtk_vertex(moris::moris_index aVertexIndex)
{
    MORIS_ASSERT(aVertexIndex <(moris::moris_index) mXtkMtkVertices.size(),"Vertex index is out of bounds");
    return mXtkMtkVertices(aVertexIndex);
}

// ----------------------------------------------------------------------------------
moris::mtk::Vertex_XTK &
Background_Mesh::get_mtk_vertex_xtk(moris::moris_index aVertexIndex)
{
    MORIS_ASSERT(aVertexIndex <(moris::moris_index) mXtkMtkVertices.size(),"Vertex index is out of bounds");
    return mXtkMtkVertices(aVertexIndex);
}
// ----------------------------------------------------------------------------------

moris::mtk::Vertex_Interpolation_XTK &
Background_Mesh::get_mtk_vertex_interpolation(moris::moris_index aVertexIndex)
{
    MORIS_ASSERT(aVertexIndex <(moris::moris_index) mXtkMtkVerticesInterpolation.size(),"Vertex index is out of bounds");
    return mXtkMtkVerticesInterpolation(aVertexIndex);
}

// ----------------------------------------------------------------------------------
moris::mtk::Cell &
Background_Mesh::get_mtk_cell(moris::moris_index aCellIndex)
{

    if(!this->is_background_cell(aCellIndex) || entity_has_children(aCellIndex,EntityRank::ELEMENT))
    {
        return get_child_element_mtk_cell(aCellIndex);
    }
    // otherwise the cell comes from the background mesh
    else
    {
        return mMeshData->get_mtk_cell(aCellIndex);
    }
}
// ----------------------------------------------------------------------------------

moris::mtk::Cell const &
Background_Mesh::get_mtk_cell(moris::moris_index aCellIndex) const
{

    if(!this->is_background_cell(aCellIndex) || entity_has_children(aCellIndex,EntityRank::ELEMENT))
    {
        return get_child_element_mtk_cell(aCellIndex);
    }
    // otherwise the cell comes from the background mesh
    else
    {
        return mMeshData->get_mtk_cell(aCellIndex);
    }
}

// ----------------------------------------------------------------------------------

moris::moris_id
Background_Mesh::allocate_entity_ids(moris::size_t         aNumReqs,
                    enum EntityRank aChildEntityRank)
{
    return mExternalMeshData.allocate_entity_ids_external_entity_data(aNumReqs,aChildEntityRank);
}

// ----------------------------------------------------------------------------------

moris::moris_index
Background_Mesh::get_first_available_index(enum EntityRank aEntityRank) const
{
    return mExternalMeshData.get_first_available_index_external_data(aEntityRank);
}


// ----------------------------------------------------------------------------------

void
Background_Mesh::update_first_available_index(moris::size_t         aNewFirstAvailableIndex,
                             enum EntityRank aEntityRank)
{
    mExternalMeshData.update_first_available_index_external_data(aNewFirstAvailableIndex, aEntityRank);
}

// ----------------------------------------------------------------------------------

//TODO: REMOVE THIS PENDING NODE STUFF WHEN FINISHED WITH DECOMP REFACTOR
void
Background_Mesh::batch_create_new_nodes(moris::Cell<xtk::Pending_Node> const & aPendingNodes)
{
    mExternalMeshData.batch_create_new_nodes_external_data(aPendingNodes);

    moris::uint tNumNewNodes = aPendingNodes.size();

    // Allocate space in the vertex interpolations
    mXtkMtkVerticesInterpolation.resize(mXtkMtkVerticesInterpolation.size() + tNumNewNodes);

    for(moris::uint i = 0; i <tNumNewNodes; i++)
    {
        mXtkMtkVertices.push_back(moris::mtk::Vertex_XTK( aPendingNodes(i).get_node_id(),
                                                          aPendingNodes(i).get_node_index(),
                                                          this));

        // add to map
        MORIS_ASSERT(mVertexGlbToLocalMap.find(aPendingNodes(i).get_node_id()) == mVertexGlbToLocalMap.end(),"Vertex already in map");
        mVertexGlbToLocalMap[aPendingNodes(i).get_node_id()] = aPendingNodes(i).get_node_index();
    }
}

void
Background_Mesh::batch_create_new_nodes(Cell<moris_index> const & aNewNodeIds,
                                        Cell<moris_index> const & aNewNodeIndices,
                                        Cell<moris::Matrix< moris::DDRMat >> const & aNewNodeCoordinates)
{
    // Batch create the new copied nodes in the mesh external data
    mExternalMeshData.batch_create_new_nodes_external_data(aNewNodeIds,aNewNodeIndices,aNewNodeCoordinates);

    moris::uint tNumNewNodes = aNewNodeIds.size();

    // Allocate space in the vertex interpolations
    mXtkMtkVerticesInterpolation.resize(mXtkMtkVerticesInterpolation.size() + tNumNewNodes);

    for(moris::uint i = 0; i <tNumNewNodes; i++)
    {
        mXtkMtkVertices.push_back(moris::mtk::Vertex_XTK( aNewNodeIds(i),
                                                          aNewNodeIndices(i),
                                                          this));

        // add to map
        MORIS_ASSERT(mVertexGlbToLocalMap.find(aNewNodeIds(i)) == mVertexGlbToLocalMap.end(),"Vertex already in map");
        mVertexGlbToLocalMap[aNewNodeIds(i)] = aNewNodeIndices(i);
    }
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::batch_create_new_nodes_as_copy_of_other_nodes(moris::Matrix< moris::IndexMat > const & aExistingNodeIndices,
                                                               moris::Matrix< moris::IndexMat > const & aNewNodeIds,
                                                               moris::Matrix< moris::IndexMat > const & aNewNodeIndices)
{
    // Collect node coordinates of nodes being copied
    moris::Matrix< moris::DDRMat > tNewNodeCoords = this->get_selected_node_coordinates_loc_inds(aExistingNodeIndices);

    // Batch create the new copied nodes in the mesh external data
    mExternalMeshData.batch_create_new_nodes_external_data(aNewNodeIds,aNewNodeIndices,tNewNodeCoords);

    moris::uint tNumNewNodes = aNewNodeIds.numel();

    // Allocate space in the vertex interpolations
    mXtkMtkVerticesInterpolation.resize(mXtkMtkVerticesInterpolation.size() + tNumNewNodes);

    for(moris::uint i = 0; i <tNumNewNodes; i++)
    {
        mXtkMtkVertices.push_back(moris::mtk::Vertex_XTK( aNewNodeIds(i),
                                                          aNewNodeIndices(i),
                                                          this));

        // add to map
        MORIS_ASSERT(mVertexGlbToLocalMap.find(aNewNodeIds(i)) == mVertexGlbToLocalMap.end(),"Vertex already in map");
        mVertexGlbToLocalMap[aNewNodeIds(i)] = aNewNodeIndices(i);
    }
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::allocate_external_node_to_child_mesh_associations()
{
    // Hard coded to 4 which is the case when a node is created on a parent edge
    // This should always be the max for XTK.
    const moris::size_t tNumCMPerNode = 4;

    // Number of external nodes
    moris::size_t tExtNumNodes = mExternalMeshData.get_num_entities_external_data(EntityRank::NODE);

    // Allocate matrix filled with moris_index max
    mNodeIndexToChildMeshIndex.resize(tExtNumNodes,tNumCMPerNode);
    mNodeIndexToChildMeshIndex.fill(std::numeric_limits<moris::moris_index>::max());
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::associate_external_nodes_to_child_mesh(moris::moris_index                       aChildMeshIndex,
                                                        moris::Matrix< moris::IndexMat > const & aNodeIndices)
{
    // Make sure the node to child mesh matrix has been allocated
    MORIS_ASSERT(mNodeIndexToChildMeshIndex.n_rows()==mExternalMeshData.get_num_entities_external_data(EntityRank::NODE),"mNodeIndexToChildMeshIndex has not been allocated");
    MORIS_ASSERT(moris::isvector(aNodeIndices), "Provided node indices need to be a vector");

    // Number of columns in mNodeIndexToChildMeshIndex
    size_t tNumCols = mNodeIndexToChildMeshIndex.n_cols();

    // Iterate over nodes and create associated between node index and child mesh index
    for(size_t i = 0 ; i < aNodeIndices.numel(); i++)
    {
        if(mExternalMeshData.is_external_entity(aNodeIndices(i),EntityRank::NODE))
        {
            moris::size_t  tExtIndex = mExternalMeshData.get_external_entity_index(aNodeIndices(i),EntityRank::NODE);
            for(size_t j = 0; j<tNumCols; j++)
            {
                if(mNodeIndexToChildMeshIndex(tExtIndex,j) == std::numeric_limits<moris::moris_index>::max())
                {
                    mNodeIndexToChildMeshIndex(tExtIndex , j ) = aChildMeshIndex;
                    break;
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------------

void
associate_external_nodes_to_child_mesh(moris::moris_index                       aChildMeshIndex,
                                       moris::Matrix< moris::IndexMat > const & aNodeIndices);

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat >
Background_Mesh::get_node_child_mesh_assocation( moris::moris_index aNodeIndex ) const
{
    MORIS_ASSERT(mExternalMeshData.is_external_entity(aNodeIndex,EntityRank::NODE),"Provided node index needs to be one created during the decomposition process");

    // External index
    moris::size_t  tExtIndex = mExternalMeshData.get_external_entity_index(aNodeIndex,EntityRank::NODE);


    moris::Matrix<moris::IndexMat> tChildMeshIndices = mNodeIndexToChildMeshIndex.get_row(tExtIndex);

    size_t tCount =0;
    for(size_t i = 0; i<tChildMeshIndices.n_cols(); i++)
    {
        if(tChildMeshIndices(i) == std::numeric_limits<moris::moris_index>::max())
        {
            break;
        }

        tCount++;
    }

    tChildMeshIndices.resize(1,tCount);
    return tChildMeshIndices;
}
// ----------------------------------------------------------------------------------

moris::size_t
Background_Mesh::get_glb_entity_id_from_entity_loc_index(moris::size_t         aEntityIndex,
                                                         enum EntityRank aEntityRank) const
{
    moris::size_t tGlbId = 0;
    if (mExternalMeshData.is_external_entity(aEntityIndex, aEntityRank))
    {
        tGlbId = mExternalMeshData.get_glb_entity_id_from_entity_loc_index_external_data(aEntityIndex, aEntityRank);
    }
    else
    {
        tGlbId = mMeshData->get_glb_entity_id_from_entity_loc_index(aEntityIndex,(moris::EntityRank)aEntityRank);
    }

    return tGlbId;
}


// ----------------------------------------------------------------------------------
moris::Matrix< moris::IdMat >
Background_Mesh::get_glb_entity_id_from_entity_loc_index_range(moris::Matrix< moris::IndexMat > const & tEntityIndices,
                                                               enum EntityRank aEntityRank) const
{
    MORIS_ASSERT(moris::isvector(tEntityIndices),"Entity indices are not provided in a vector");
    moris::size_t tNumEntities = tEntityIndices.numel();
    moris::Matrix< moris::IdMat > tEntityIds(1,tNumEntities);

    for(moris::size_t i =0; i<tNumEntities; i++)
    {
        tEntityIds(0,i) = this->get_glb_entity_id_from_entity_loc_index(tEntityIndices(i),aEntityRank);
    }
    return tEntityIds;
}

// ----------------------------------------------------------------------------------

void
Background_Mesh::convert_loc_entity_ind_to_glb_entity_ids(enum EntityRank aEntityRank,
                                             moris::Matrix< moris::IndexMat > & aEntityIndices) const
    {
        moris::uint tNumRows = aEntityIndices.n_rows();
        moris::uint tNumCols = aEntityIndices.n_cols();
        for(moris::uint j = 0; j<tNumCols; j++)
        {
            for(moris::uint i = 0; i < tNumRows; i++)
            {
                aEntityIndices(i,j) = this->get_glb_entity_id_from_entity_loc_index(aEntityIndices(i,j),aEntityRank);
            }
        }
    }
// ----------------------------------------------------------------------------------
bool
Background_Mesh::is_background_node(moris::moris_index aNodeIndex)
{
    return !mExternalMeshData.is_external_entity(aNodeIndex,EntityRank::NODE);
}

// ----------------------------------------------------------------------------------
bool
Background_Mesh::is_background_cell(moris::moris_index aNodeIndex) const
{
    return !mExternalMeshData.is_external_entity(aNodeIndex,EntityRank::ELEMENT);
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::DDRMat >
Background_Mesh::get_all_node_coordinates_loc_inds() const
{
    // Get counts inside of STK and in External Data (both happen in this function call)
    moris::size_t tNumNodes = this->get_num_entities(EntityRank::NODE);
    moris::size_t tNumBGNodes = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);

    moris::Matrix< moris::DDRMat > tAllNodeCoordinates = this->get_all_node_coordinates_loc_inds_background_mesh();
    tAllNodeCoordinates.resize(tNumNodes,3);
    // Get node coordinates from external entities
    mExternalMeshData.get_all_node_coordinates_loc_inds_external_data(tNumBGNodes,tAllNodeCoordinates);
    return tAllNodeCoordinates;
}

// ---------------------------------------------------------------------------------
moris::Matrix< moris::DDRMat >
Background_Mesh::get_selected_node_coordinates_loc_inds(moris::Matrix< moris::IndexMat > const & aNodeIndices) const
{

    MORIS_ERROR(moris::isvector(aNodeIndices),"Provided Node indices need to be a vector");
    // TODO: Add external entity check to see if xtk has the coordinate field or stk has it
    // Number of spatial dimensions
    moris::size_t tSpatialDimension = mMeshData->get_spatial_dim();

    // Get number of nodes provided
    moris::size_t tNumNodes = aNodeIndices.numel();

    // Initialize output matrix
    moris::Matrix< moris::DDRMat > tSelectedNodesCoords(tNumNodes, tSpatialDimension);

    // Loop over all nodes
    enum EntityRank tEntityRank = EntityRank::NODE;

    for (moris::size_t n = 0; n < tNumNodes; ++n)
    {
        if (mExternalMeshData.is_external_entity(aNodeIndices(n), tEntityRank))
        {
            moris::Matrix< moris::DDRMat > const & tNodeCoords = mExternalMeshData.get_selected_node_coordinates_loc_inds_external_data(aNodeIndices(n));
            tSelectedNodesCoords.set_row(n,tNodeCoords);
        }

        else
        {
            tSelectedNodesCoords.set_row(n,mMeshData->get_node_coordinate((moris_index)aNodeIndices(n)));
        }
    }

    return tSelectedNodesCoords;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IdMat >
Background_Mesh::get_local_to_global_map(enum EntityRank aEntityRank) const
{
    MORIS_ERROR(aEntityRank==EntityRank::NODE," This function is only implemented for node maps");
    size_t tNumNodes = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);

    moris::Matrix<moris::IndexMat> tLocalToGlobalBM(1,tNumNodes);

    for(size_t i = 0; i<tNumNodes; i++)
    {
        tLocalToGlobalBM(i) = (moris::size_t)mMeshData->get_glb_entity_id_from_entity_loc_index((moris_index)i,moris::EntityRank::NODE);
    }

    moris::Matrix<moris::IndexMat> const & tLocalToGlobalExt = mExternalMeshData.get_local_to_global_node_map();


    tNumNodes = this->get_num_entities(EntityRank::NODE);
    size_t tFirstExtNodeInd = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);
    MORIS_ERROR(tNumNodes = tLocalToGlobalBM.numel() + tLocalToGlobalExt.numel(),"Number of nodes returned does not match the number in the map");

    // combine the two maps
    moris::Matrix<moris::IdMat> tLocalToGlobal(1,tNumNodes);

    for(moris::size_t i = 0 ; i <tLocalToGlobalBM.numel(); i++)
    {
        tLocalToGlobal(i) = tLocalToGlobalBM(i);
    }

    for(moris::size_t i = 0 ; i <tLocalToGlobalExt.numel(); i++)
    {
        tLocalToGlobal(i+tFirstExtNodeInd) = tLocalToGlobalExt(i);
    }

    return tLocalToGlobal;
}

// ----------------------------------------------------------------------------------

moris::Matrix<moris::IdMat>
Background_Mesh::get_full_non_intersected_node_to_element_glob_ids() const
{
    moris::size_t tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);
    enum CellTopology tElemTopo = get_XTK_mesh_element_topology();

    moris::size_t tNumNodesPerElem = 0;
    if(tElemTopo == CellTopology::TET4)
    {
        tNumNodesPerElem = 4;
    }
    else if(tElemTopo == CellTopology::HEX8)
    {
        tNumNodesPerElem = 8;
    }
    else
    {
        tNumNodesPerElem = 8;
    }

    moris::Matrix<moris::IdMat> tElementToNode(tNumElementsBG,tNumNodesPerElem);
    moris::size_t tCount = 0;
    for(moris::size_t i = 0; i<tElementToNode.n_rows(); i++)
    {
        if(!this->entity_has_children(i,EntityRank::ELEMENT))
        {
            moris::Matrix<moris::IndexMat> tElementToNodeInd
            = mMeshData->get_entity_connected_to_entity_loc_inds((moris::moris_index)i,
                                                                  moris::EntityRank::ELEMENT,
                                                                  moris::EntityRank::NODE);

            moris::Matrix< moris::IdMat >tElementToNodeId
                = moris::mtk::convert_entity_indices_to_ids(tElementToNodeInd,moris::EntityRank::NODE,mMeshData);
            tElementToNode.set_row(tCount,tElementToNodeId);
            tCount++;
        }
    }
    tElementToNode.resize(tCount,tNumNodesPerElem);
    return tElementToNode;
}

// ----------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IdMat>>
Background_Mesh::get_non_intersected_element_to_node_by_phase(moris::uint aNumPhases)
{
    moris::size_t tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);
    enum CellTopology tElemTopo = get_XTK_mesh_element_topology();

    moris::size_t tNumNodesPerElem = 0;
    if(tElemTopo == CellTopology::TET4)
    {
        tNumNodesPerElem = 4;
    }
    else if(tElemTopo == CellTopology::HEX8)
    {
        tNumNodesPerElem = 8;
    }
    else
    {
        tNumNodesPerElem = 8;
    }

    Cell<moris::Matrix<moris::IdMat>> tElementToNodeByPhase(aNumPhases,moris::Matrix<moris::IdMat>(tNumElementsBG,tNumNodesPerElem));
    moris::Matrix<moris::DDUMat> tPhaseCount(1,aNumPhases,0);

    for(moris::size_t i = 0; i<tNumElementsBG; i++)
    {
        if(!this->entity_has_children(i,EntityRank::ELEMENT))
        {

            moris::moris_id tId = mMeshData->get_glb_entity_id_from_entity_loc_index(i,EntityRank::ELEMENT);
            moris::Matrix< moris::IdMat >tElementToNodeId = mMeshData->get_entity_connected_to_entity_glob_ids(tId,EntityRank::ELEMENT,EntityRank::NODE);

            moris::moris_index tPhaseIndex = this->get_element_phase_index(i);
            if(mMeshData->get_mesh_type() == MeshType::HMR)
            {
                tElementToNodeByPhase(tPhaseIndex).get_row(tPhaseCount(tPhaseIndex)) = moris::trans(tElementToNodeId);
            }

            else
            {
                tElementToNodeByPhase(tPhaseIndex).get_row(tPhaseCount(tPhaseIndex)) = tElementToNodeId.matrix_data();
            }
            tPhaseCount(tPhaseIndex)++;
        }
    }

    for(moris::uint i = 0; i <aNumPhases; i++)
    {
        tElementToNodeByPhase(i).resize(tPhaseCount(i),tNumNodesPerElem);
    }

    return tElementToNodeByPhase;
}

// ----------------------------------------------------------------------------------

moris::Matrix<moris::IdMat>
Background_Mesh::get_all_non_intersected_elements() const
{
    moris::size_t tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);

    moris::Matrix<moris::IdMat> tElementIds(tNumElementsBG,1);
    moris::size_t tCount = 0;
    for(moris::size_t i = 0; i<tElementIds.numel(); i++)
    {
        if(!this->entity_has_children(i,EntityRank::ELEMENT))
        {

            moris::moris_id tElementToNodeId = mMeshData->get_glb_entity_id_from_entity_loc_index((moris::moris_index)i,moris::EntityRank::ELEMENT);
            tElementIds(tCount) = tElementToNodeId;
            tCount++;
        }
    }
    tElementIds.resize(tCount,1);
    return tElementIds;
}

// ----------------------------------------------------------------------------------

Cell<moris::Matrix<moris::IdMat>>
Background_Mesh::get_all_non_intersected_elements_by_phase( uint aNumPhases ) const
{
    uint  tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);

    //Initialize output
    Cell<moris::Matrix<moris::IdMat>> tElementsByPhase(aNumPhases);
    moris::Matrix<moris::DDUMat> tPhaseCount(1,aNumPhases,0);
    for(uint i =0; i <aNumPhases; i++)
    {
        tElementsByPhase(i) = moris::Matrix<moris::IdMat>(1,tNumElementsBG);
    }

    for(uint iE = 0; iE<tNumElementsBG; iE++)
    {
        if(!this->entity_has_children(iE,EntityRank::ELEMENT))
        {
            moris::moris_index tPhaseInd = this->get_element_phase_index(iE);
            uint tCount = tPhaseCount(tPhaseInd);
            tElementsByPhase(tPhaseInd)(tCount) = mMeshData->get_glb_entity_id_from_entity_loc_index((moris::moris_index)iE,moris::EntityRank::ELEMENT);
            tPhaseCount(tPhaseInd)++;
        }
    }

    // Resize
    for(uint i =0; i <aNumPhases; i++)
    {
        tElementsByPhase(i).resize(1,tPhaseCount(i));
    }

    return  tElementsByPhase;
}
// ----------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Background_Mesh::get_all_non_intersected_elements_loc_inds() const
{
    moris::size_t tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);

    moris::Matrix<moris::IdMat> tElementInds(tNumElementsBG,1);
    moris::size_t tCount = 0;
    for(moris::size_t i = 0; i<tElementInds.numel(); i++)
    {
        if(!this->entity_has_children(i,EntityRank::ELEMENT))
        {
            tElementInds(tCount) = i;
            tCount++;
        }
    }
    tElementInds.resize(tCount,1);
    return tElementInds;
}

// ----------------------------------------------------------------------------------

void
Background_Mesh::register_new_downward_inheritance(Cell<std::pair<moris::moris_index,moris::moris_index>> const & aNewElementToChildMeshPairs)
{
    moris::size_t tNumNewPairs = aNewElementToChildMeshPairs.size();
    for(moris::size_t i = 0; i<tNumNewPairs; i++)
    {
        mElementDownwardInheritance.register_new_inheritance_pair(aNewElementToChildMeshPairs(i).first,aNewElementToChildMeshPairs(i).second);
    }
}

// ----------------------------------------------------------------------------------
bool
Background_Mesh::entity_has_children(moris::size_t aEntityIndex,
                                     enum EntityRank aEntityRank) const
{
    MORIS_ASSERT(aEntityRank==EntityRank::ELEMENT,"ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED");

    return mElementDownwardInheritance.has_inheritance(aEntityIndex);
}

// ----------------------------------------------------------------------------------
moris::moris_index const &
Background_Mesh::child_mesh_index(moris::size_t aEntityIndex,
                                  enum EntityRank aEntityRank)
{
    MORIS_ASSERT(aEntityRank==EntityRank::ELEMENT,"ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED");

    return mElementDownwardInheritance.get_inheritance(aEntityIndex);
}

// ---------------------------------------------------------------------------------

void
Background_Mesh::initialize_interface_node_flags(moris::size_t const & aNumNodes,
                                                 moris::size_t const & aNumGeometry)
{
    mInterfaceNodeFlag = moris::Matrix< moris::IndexMat >(aNumNodes,aNumGeometry,0);
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::allocate_space_in_interface_node_flags(moris::size_t aNumNodes,
                                       moris::size_t aNumGeometry)
{
    moris::size_t tCurrentSize = mInterfaceNodeFlag.n_rows();
    mInterfaceNodeFlag.resize(tCurrentSize + aNumNodes, aNumGeometry);

    for(moris::size_t i = tCurrentSize; i<tCurrentSize+aNumNodes; i++)
    {
        for(moris::size_t j = 0; j<aNumGeometry; j++)
        {
            mInterfaceNodeFlag(i,j) = 0;
        }
    }
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::mark_node_as_interface_node(moris::moris_index aNodeIndex,
                                             moris::size_t aGeomIndex)
{
    mInterfaceNodeFlag(aNodeIndex,aGeomIndex) = 1;
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::mark_nodes_as_interface_node_loc_inds(moris::Matrix<moris::IndexMat> aNodeIndices,
                                                       moris::size_t aGeomIndex)
{
    for(uint  i = 0; i <aNodeIndices.numel(); i++)
    {
        mInterfaceNodeFlag(aNodeIndices(i),aGeomIndex) = 1;
    }
}

// ----------------------------------------------------------------------------------

bool
Background_Mesh::is_interface_node(moris::moris_index aNodeIndex,
                                   moris::size_t aGeomIndex) const
{
    MORIS_ASSERT(aNodeIndex<(moris::moris_index)mInterfaceNodeFlag.n_rows(),"Attempting to access interface node flag for node index out of bounds. Have you called allocate_space_in_interface_node_flags?");

    if(mInterfaceNodeFlag(aNodeIndex,aGeomIndex) == 1)
    {
        return true;
    }

    else
    {
        return false;
    }
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat >
Background_Mesh::get_interface_nodes_loc_inds(moris::moris_index aGeometryIndex) const
{
    // initialize output
    moris::size_t tNumNodes = this->get_num_entities(EntityRank::NODE);
    moris::Matrix< moris::IndexMat > tInterfaceNodes(1,tNumNodes);

    // keep track of how many interface nodes
    moris::size_t tCount = 0;

    for(moris::size_t i = 0; i<tNumNodes; i++)
    {
        if(is_interface_node(i,aGeometryIndex))
        {
            tInterfaceNodes(0,tCount) = i;
            tCount++;
        }
    }

    tInterfaceNodes.resize(1,tCount);
    return tInterfaceNodes;
}

// ----------------------------------------------------------------------------------

Cell<moris::Matrix< moris::IdMat >>
Background_Mesh::get_interface_nodes_loc_inds() const
{
    // initialize output
    moris::size_t tNumGeoms = mInterfaceNodeFlag.n_cols();

    Cell<moris::Matrix< moris::IdMat >>tInterfaceNodes(tNumGeoms);

    for(moris::size_t i = 0 ; i <tNumGeoms; i++)
    {
        tInterfaceNodes(i) = get_interface_nodes_loc_inds(i);
    }

    return tInterfaceNodes;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IdMat >
Background_Mesh::get_interface_nodes_glob_ids(moris::moris_index aGeometryIndex) const
{
    // initialize output
    moris::size_t tNumNodes = this->get_num_entities(EntityRank::NODE);
    moris::Matrix< moris::IdMat > tInterfaceNodes(1,tNumNodes);

    // keep track of how many interface nodes
    moris::size_t tCount = 0;

    for(moris::size_t i = 0; i<tNumNodes; i++)
    {
        if(is_interface_node(i,aGeometryIndex))
        {
            tInterfaceNodes(0,tCount) = get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tCount++;
        }
    }

    tInterfaceNodes.resize(1,tCount);
    return tInterfaceNodes;
}

// ----------------------------------------------------------------------------------
Cell<moris::Matrix< moris::IdMat >>
Background_Mesh::get_interface_nodes_glob_ids()
{
    // initialize output
    moris::size_t tNumGeoms = mInterfaceNodeFlag.n_cols();

    Cell<moris::Matrix< moris::IdMat >>tInterfaceNodes(tNumGeoms);

    for(moris::size_t i = 0 ; i <tNumGeoms; i++)
    {
        tInterfaceNodes(i) = get_interface_nodes_glob_ids(i);
    }

    return tInterfaceNodes;
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::print_interface_node_flags()
{
    for(moris::size_t i = 0; i<mInterfaceNodeFlag.n_rows(); i++)
    {
        std::cout<<this->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE)<<" | ";
        for(moris::size_t j = 0; j<mInterfaceNodeFlag.n_cols(); j++)
        {
            std::cout<<mInterfaceNodeFlag(i,j)<<" ";
        }
        std::cout<<std::endl;
    }

}

// ----------------------------------------------------------------------------------

void
Background_Mesh::initialize_element_phase_indices(moris::size_t const & aNumElements)
{
    mElementPhaseIndex = moris::Matrix< moris::IndexMat >(aNumElements,1);
}

// ----------------------------------------------------------------------------------

void
Background_Mesh::set_element_phase_index(moris::size_t const & aElementIndex,
                                         moris::size_t const & aElementPhaseIndex)
{
    mElementPhaseIndex(aElementIndex,0) = aElementPhaseIndex;
}

// ----------------------------------------------------------------------------------

moris::moris_index const &
Background_Mesh::get_element_phase_index(moris::size_t const & aElementIndex) const
{
    return (mElementPhaseIndex)( aElementIndex, 0 );
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat >
Background_Mesh::get_element_phase_inds(moris::Matrix<moris::DDSTMat> const & aElementInds )
{
    moris::Matrix< moris::IndexMat > tElementPhasesInds(1,aElementInds.n_cols());

    for(moris::size_t i = 0; i < aElementInds.n_cols(); i++)
    {
        tElementPhasesInds(0,i) = this->get_element_phase_index(aElementInds(i));
    }

    return tElementPhasesInds;
}

// ----------------------------------------------------------------------------------

moris_index
Background_Mesh::get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                                       enum EntityRank aEntityRank) const
{

    MORIS_ERROR(aEntityRank == EntityRank::NODE,"Only a node map is implemented in XTK");

    auto tIter = mVertexGlbToLocalMap.find(aEntityId);

    MORIS_ERROR(tIter!=mVertexGlbToLocalMap.end(),
                "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u EntityRank = %u on process %u",aEntityId, (uint)aEntityRank, par_rank());

    return tIter->second;
}

// ----------------------------------------------------------------------------------

void
Background_Mesh::add_child_element_to_mtk_cells(moris::moris_index aElementIndex,
                                                moris::moris_index aElementId,
                                                moris::moris_index aCMElementIndex,
                                                Child_Mesh*        aChildMeshPtr)
{
    mChildMtkCells.push_back(moris::mtk::Cell_XTK(aElementId,
                                                  aElementIndex,
                                                  aCMElementIndex,
                                                  aChildMeshPtr,
                                                  this));

    MORIS_ASSERT(mChildMtkCellMap.find(aElementIndex) == mChildMtkCellMap.end(),"Element index already has an mtk cell associated with it");

    mChildMtkCellMap[aElementIndex] = mChildMtkCells.size()-1;

}

// ----------------------------------------------------------------------------------
const moris::mtk::Cell*
Background_Mesh::get_child_element_mtk_cell_ptr(moris::moris_index aElementIndex) const
{
    auto tIter = mChildMtkCellMap.find(aElementIndex);

    MORIS_ASSERT(mChildMtkCellMap.find(aElementIndex) != mChildMtkCellMap.end(),"Element index is not in the map");

    moris::moris_index tIndex = tIter->second;
    return &mChildMtkCells(tIndex);
}

// ----------------------------------------------------------------------------------

moris::mtk::Cell const &
Background_Mesh::get_child_element_mtk_cell(moris::moris_index aElementIndex) const
{
    auto tIter = mChildMtkCellMap.find(aElementIndex);

    MORIS_ASSERT(mChildMtkCellMap.find(aElementIndex) != mChildMtkCellMap.end(),"Element index is not in the map");

    moris::moris_index tIndex = tIter->second;
    return mChildMtkCells(tIndex);
}

// ----------------------------------------------------------------------------------

moris::mtk::Mesh &
Background_Mesh::get_mesh_data()
{
    return *mMeshData;
}

// ----------------------------------------------------------------------------------

moris::mtk::Mesh const &
Background_Mesh::get_mesh_data() const
{
    return *mMeshData;
}

// ----------------------------------------------------------------------------------
enum CellTopology
Background_Mesh::get_XTK_mesh_element_topology() const
{
    enum CellTopology tElementTopology = CellTopology::INVALID;
    moris::Matrix<  moris::IndexMat  > tElementNodes = mMeshData->get_entity_connected_to_entity_loc_inds(0,(moris::EntityRank)EntityRank::ELEMENT, (moris::EntityRank)EntityRank::NODE);
    if(tElementNodes.numel() == 8 && moris::isvector(tElementNodes))
    {
        tElementTopology = CellTopology::HEX8;

    }
    else if (tElementNodes.numel() == 4 && moris::isvector(tElementNodes))
    {
        tElementTopology = CellTopology::TET4;
    }
    else
    {
        std::cout<<"Topology not recognized in parent mesh";
    }

    return tElementTopology;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::DDRMat >
Background_Mesh::get_all_node_coordinates_loc_inds_background_mesh() const
{
    size_t tNumNodes = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);
    moris::Matrix<moris::DDRMat> tNodeCoords(tNumNodes,mMeshData->get_spatial_dim());
    for(size_t i = 0; i< tNumNodes; i++ )
    {
        tNodeCoords.set_row(i, mMeshData->get_node_coordinate(i));
    }
    return tNodeCoords;
}

// ----------------------------------------------------------------------------------

void
Background_Mesh::intialize_downward_inheritance()
{
    moris::size_t tNumElements = mMeshData->get_num_entities((moris::EntityRank)EntityRank::ELEMENT);
    MORIS_ASSERT(tNumElements!=0,"Empty Mesh Given to XTK Mesh");

    mElementDownwardInheritance = Downward_Inheritance<moris::moris_index,moris::moris_index>(tNumElements);
}

// ----------------------------------------------------------------------------------
void
Background_Mesh::initialize_background_mesh_vertices()
{
    moris::uint tNumNodes = mMeshData->get_num_entities(EntityRank::NODE);

    mXtkMtkVerticesInterpolation = moris::Cell<moris::mtk::Vertex_Interpolation_XTK>(tNumNodes);

    for(moris::uint  i = 0; i <tNumNodes; i++)
    {
        moris::mtk::Vertex * tVertex = &mMeshData->get_mtk_vertex(i);

        MORIS_ASSERT(mVertexGlbToLocalMap.find(tVertex->get_id()) == mVertexGlbToLocalMap.end(),"Vertex already in map");

        mVertexGlbToLocalMap[tVertex->get_id()] = tVertex->get_index();
        mXtkMtkVertices.push_back(tVertex);
    }

}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

}


