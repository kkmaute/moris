/*
 * cl_XTK_Mesh.hpp
 *
 *  Created on: Aug 10, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_MESH_HPP_
#define SRC_XTK_CL_XTK_MESH_HPP_

#include <unordered_map>
#include <utility>

// Mesh Includes:
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Enums.hpp"

// Assertion Includes:
#include "fn_assert.hpp"

// Linear Algebra Includes
#include "cl_Matrix.hpp"
#include "fn_print.hpp"
#include "fn_isvector.hpp"

//XTK Includes:
#include "cl_XTK_Node.hpp"
#include "cl_XTK_External_Mesh_Data.hpp"
#include "cl_XTK_Downward_Inheritance.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_MTK_Cell_XTK_Impl.hpp"

// Geometry Engine Includes
#include "cl_MGE_Geometry_Engine.hpp"

namespace xtk
{
class Background_Mesh
{
public:
    Background_Mesh(){};

    Background_Mesh(moris::mtk::Mesh* aMeshData):
        mMeshData(aMeshData),
        mChildCellPtrs(0),
        mNodeIndexToChildMeshIndex(0,0)
    {
        intialize_downward_inheritance();
        mExternalMeshData.set_up_external_entity_data(mMeshData);
    }

    Background_Mesh(moris::mtk::Mesh* aMeshData,
             Geometry_Engine & aGeometryEngine):
        mMeshData(aMeshData),
        mChildCellPtrs(0),
        mNodeIndexToChildMeshIndex(0,0)
    {
        intialize_downward_inheritance();
        mExternalMeshData.set_up_external_entity_data(mMeshData);
    }

    /*
     * Get number of entities in the background mesh and
     * the number of entities XTK has created
     *
     * NOTE: this function includes all nodes (background nodes and interface nodes)
     * but only includes background elements
     *
     */
    moris::size_t
    get_num_entities(enum EntityRank aEntityRank) const
    {
        // Initialize
        moris::size_t tNumBackgroundEntities = mMeshData->get_num_entities((moris::EntityRank)aEntityRank);
        moris::size_t tExternalEntities = mExternalMeshData.get_num_entities_external_data(aEntityRank);

        return tNumBackgroundEntities + tExternalEntities;
    }

    /*
     * Get number of entities in the background mesh
     *
     */
    moris::size_t
    get_num_entities_background(enum EntityRank aEntityRank) const
    {
        return mMeshData->get_num_entities((moris::EntityRank)aEntityRank);
    }

    moris::mtk::Cell &
    get_mtk_cell(moris::moris_index aCellIndex)
    {

        // if the cell index provided is not one which has children (we assume it is a child element)
        if(!entity_has_children(aCellIndex,EntityRank::ELEMENT))
        {
            return get_child_element_mtk_cell(aCellIndex);
        }
        // otherwise the cell comes from the background mesh
        else
        {
            return mMeshData->get_mtk_cell(aCellIndex);
        }
    }

    /*
     *  Get an offset which is the first id allocated (assumed ids are grouped)
     */
    moris::moris_id
    allocate_entity_ids(moris::size_t         aNumReqs,
                        enum EntityRank aChildEntityRank)
    {
        return mExternalMeshData.allocate_entity_ids_external_entity_data(aNumReqs,aChildEntityRank);
    }

    /*
     * Get the first available index
     */
    moris::moris_index
    get_first_available_index(enum EntityRank aEntityRank) const
    {
        return mExternalMeshData.get_first_available_index_external_data(aEntityRank);
    }


    /*
     * Increment the first available index by aNewFirstAvailableIndex
     */
    void
    update_first_available_index(moris::size_t         aNewFirstAvailableIndex,
                                 enum EntityRank aEntityRank)
    {
        mExternalMeshData.update_first_available_index_external_data(aNewFirstAvailableIndex, aEntityRank);
    }

    /*
     * Create a batch of new nodes
     */
    void
    batch_create_new_nodes(moris::Cell<xtk::Pending_Node> const & aPendingNodes)
    {
        mExternalMeshData.batch_create_new_nodes_external_data(aPendingNodes);
    }


    void
    batch_create_new_nodes_as_copy_of_other_nodes(moris::Matrix< moris::IndexMat > const & aExistingNodeIndices,
                                                moris::Matrix< moris::IndexMat > const & aNewNodeIds,
                                                moris::Matrix< moris::IndexMat > const & aNewNodeIndices)
    {
        // Collect node coordinates of nodes being copied
        moris::Matrix< moris::DDRMat > tNewNodeCoords = this->get_selected_node_coordinates_loc_inds(aExistingNodeIndices);

        // Batch create the new copied nodes in the mesh external data
        mExternalMeshData.batch_create_new_nodes_external_data(aNewNodeIds,aNewNodeIndices,tNewNodeCoords);
    }

    void
    allocate_external_node_to_child_mesh_associations()
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

    /*
     * Create association of external nodes to their child mesh index,
     * The node index vector does not necessarily need to be only external nodes
     * but only the ones which are external will be associated to a child mesh
     */
    void
    associate_external_nodes_to_child_mesh(moris::moris_index                       aChildMeshIndex,
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

    /*
     * Get the child mesh indices that a node belongs to.
     * Only implemented for nodes created during the decomposition process, called the
     * external nodes
     */
    moris::Matrix< moris::IndexMat >
    get_node_child_mesh_assocation( moris::moris_index aNodeIndex ) const
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


    /*
     * Get the global entity id from local index
     */
    moris::size_t
    get_glb_entity_id_from_entity_loc_index(moris::size_t         aEntityIndex,
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

    /*
     * From a vector of entity ids and ranks, return the global ids of these entities
     */
    moris::Matrix< moris::IdMat >
    get_glb_entity_id_from_entity_loc_index_range(moris::Matrix< moris::IndexMat > const & tEntityIndices,
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

    void
    convert_loc_entity_ind_to_glb_entity_ids(enum EntityRank aEntityRank,
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

    /*!
     * Returns whether a node was in the original background mesh
     */
    bool
    is_background_node(moris::moris_index aNodeIndex)
    {
        return !mExternalMeshData.is_external_entity(aNodeIndex,EntityRank::NODE);
    }

    /*
     * Return all node coordinates ordered by local indices
     */
    moris::Matrix< moris::DDRMat >
    get_all_node_coordinates_loc_inds() const
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

    /*
     * Return a coordinate matrix for the specified node indices
     */
    moris::Matrix< moris::DDRMat >
    get_selected_node_coordinates_loc_inds(
            moris::Matrix< moris::IndexMat > const & aNodeIndices) const
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

    moris::Matrix< moris::IdMat >
    get_local_to_global_map(enum EntityRank aEntityRank) const
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

    /*
     * Return a vector of all non-intersected elements'
     * element to node connectivity
     */
    moris::Matrix<moris::IdMat>
    get_full_non_intersected_node_to_element_glob_ids() const
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

    Cell<moris::Matrix<moris::IdMat>>
    get_non_intersected_element_to_node_by_phase(moris::uint aNumPhases)
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

    /*
     * Return all ids of non-intersected elements
     */
    moris::Matrix<moris::IdMat>
    get_all_non_intersected_elements() const
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

    /*
     * Return all ids of non-intersected elements
     */
    Cell<moris::Matrix<moris::IdMat>>
    get_all_non_intersected_elements_by_phase( uint aNumPhases ) const
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



    /*
     * Return all non-intersected element proc local indices
     */
    moris::Matrix<moris::IndexMat>
    get_all_non_intersected_elements_loc_inds() const
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

    // -------------------------------------------------------------------
    // Functions for setting up downard inheritance, where downward
    // downward inheritance is the relationship between a parent element
    // and its children elements
    // -------------------------------------------------------------------
    void
    register_new_downward_inheritance(Cell<std::pair<moris::moris_index,moris::moris_index>> const & aNewElementToChildMeshPairs)
    {
        moris::size_t tNumNewPairs = aNewElementToChildMeshPairs.size();
        for(moris::size_t i = 0; i<tNumNewPairs; i++)
        {
            mElementDownwardInheritance.register_new_inheritance_pair(aNewElementToChildMeshPairs(i).first,aNewElementToChildMeshPairs(i).second);
        }
    }

    /*
     * returns whether a given entity has any children entities
     */
    bool
    entity_has_children(moris::size_t aEntityIndex,
                        enum EntityRank aEntityRank) const
    {
        MORIS_ASSERT(aEntityRank==EntityRank::ELEMENT,"ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED");

        return mElementDownwardInheritance.has_inheritance(aEntityIndex);
    }

    /*
     * returns the child mesh index of entity with children
     */
    moris::moris_index const & child_mesh_index(moris::size_t aEntityIndex,
                                                enum EntityRank aEntityRank)
    {
        MORIS_ASSERT(aEntityRank==EntityRank::ELEMENT,"ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED");

        return mElementDownwardInheritance.get_inheritance(aEntityIndex);
    }

    // -------------------------------------------------------------------
    // Functions related to setting and accessing interface node information
    // -------------------------------------------------------------------

    void
    initialize_interface_node_flags(moris::size_t const & aNumNodes,
                                    moris::size_t const & aNumGeometry)
    {
        mInterfaceNodeFlag = moris::Matrix< moris::IndexMat >(aNumNodes,aNumGeometry,0);
    }

    /*
     * Allocates additional space in interface node flags
     */
    void
    allocate_space_in_interface_node_flags(moris::size_t aNumNodes,
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

    /*
     * Marks a node as an interface node for a given geometry index
     */
    void
    mark_node_as_interface_node(moris::moris_index aNodeIndex,
                                moris::size_t aGeomIndex)
    {
        mInterfaceNodeFlag(aNodeIndex,aGeomIndex) = 1;
    }

    /*
     * Marks a node as an interface node for a given geometry index
     */
    void
    mark_nodes_as_interface_node_loc_inds(moris::Matrix<moris::IndexMat> aNodeIndices,
                                          moris::size_t aGeomIndex)
    {
        for(uint  i = 0; i <aNodeIndices.numel(); i++)
        {
            mInterfaceNodeFlag(aNodeIndices(i),aGeomIndex) = 1;
        }
    }

    /*
     * Returns whether a node is an interface node for a given geometry index
     */
    bool
    is_interface_node(moris::moris_index aNodeIndex,
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


    /*
     * get the interface nodes with respect to a given geometry index
     */
    moris::Matrix< moris::IndexMat >
    get_interface_nodes_loc_inds(moris::moris_index aGeometryIndex) const
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

    /*
     * get the interface nodes with respect to a given geometry index
     */
    Cell<moris::Matrix< moris::IdMat >>
    get_interface_nodes_loc_inds() const
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

    /*
     * get the interface nodes with respect to a given geometry index
     */
    moris::Matrix< moris::IdMat >
    get_interface_nodes_glob_ids(moris::moris_index aGeometryIndex) const
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

    /*
     * get the interface nodes with respect to a given geometry index
     */
    Cell<moris::Matrix< moris::IdMat >>
    get_interface_nodes_glob_ids()
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


    void
    print_interface_node_flags()
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


    // -------------------------------------------------------------------
    // Functions related to setting and accessing element phase indices
    // -------------------------------------------------------------------
    /*
     * Allocate space for element phase indices
     */
    void
    initialize_element_phase_indices(moris::size_t const & aNumElements)
    {
        mElementPhaseIndex = moris::Matrix< moris::IndexMat >(aNumElements,1);
    }

    /*
     * Set the phase index value of element with element index. This is relative to each geometry.
     */
    void
    set_element_phase_index(moris::size_t const & aElementIndex,
                            moris::size_t const & aElementPhaseIndex)
    {
        mElementPhaseIndex(aElementIndex,0) = aElementPhaseIndex;
    }

    /*
     * Get the phase index value of element with element index
     */
    moris::moris_index const &
    get_element_phase_index(moris::size_t const & aElementIndex) const
    {
        return (mElementPhaseIndex)( aElementIndex, 0 );
    }

    moris::Matrix< moris::IndexMat >
    get_element_phase_inds(moris::Matrix<moris::DDSTMat> const & aElementInds )
    {
        moris::Matrix< moris::IndexMat > tElementPhasesInds(1,aElementInds.n_cols());

        for(moris::size_t i = 0; i < aElementInds.n_cols(); i++)
        {
            tElementPhasesInds(0,i) = this->get_element_phase_index(aElementInds(i));
        }

        return tElementPhasesInds;
    }


    void
    add_child_element_to_mtk_cells(moris::moris_index aElementIndex,
                                   moris::moris_index aElementId,
                                   moris::moris_index aCMElementIndex,
                                   Child_Mesh*        aChildMeshPtr)
    {
        mChildCellPtrs.push_back(moris::mtk::XTK_Cell(aElementId,
                                                      aElementIndex,
                                                      aCMElementIndex,
                                                      aChildMeshPtr,
                                                      this));

        MORIS_ASSERT(mChildCellPtrMap.find(aElementIndex) == mChildCellPtrMap.end(),"Element index already has an mtk cell associated with it");

        mChildCellPtrMap[aElementIndex] = mChildCellPtrs.size()-1;
    }

    const moris::mtk::Cell*
    get_child_element_mtk_cell_ptr(moris::moris_index aElementIndex) const
    {
        auto tIter = mChildCellPtrMap.find(aElementIndex);

        moris::moris_index tIndex = tIter->second;
        return &mChildCellPtrs(tIndex);
    }

    moris::mtk::Cell &
    get_child_element_mtk_cell(moris::moris_index aElementIndex)
    {
        auto tIter = mChildCellPtrMap.find(aElementIndex);

        moris::moris_index tIndex = tIter->second;
        return mChildCellPtrs(tIndex);
    }


    // -------------------------------------------------------------------
    // Access underlying mesh data functions
    moris::mtk::Mesh &
    get_mesh_data()
    {
        return *mMeshData;
    }

    moris::mtk::Mesh const &
    get_mesh_data() const
    {
        return *mMeshData;
    }
    // -------------------------------------------------------------------


    /*
     * Get the base topology of parent elements in the background mesh
     */
    enum CellTopology
    get_XTK_mesh_element_topology() const
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
private:
    // Background mesh data
    moris::mtk::Mesh* mMeshData;

    // External Entity information
    // The background mesh remains constant, and every new entity created is stored within XTK
    // child meshes
    Mesh_External_Entity_Data mExternalMeshData;

    // Downward inheritance pairs (links elements in XTK mesh to indices in Child Meshes)
    Downward_Inheritance<moris::moris_index, moris::moris_index> mElementDownwardInheritance;

    // Elements constructured by the decomposition process mtk Cells
    std::map < moris_id, moris_index > mChildCellPtrMap; /* To go from cell index to location in child cell ptrs*/
    moris::Cell<moris::mtk::XTK_Cell> mChildCellPtrs;

    // Associate external node indices to the child meshes they belong to
    // Row - External node index
    // Col - Child Mesh Index
    moris::Matrix< moris::IndexMat > mNodeIndexToChildMeshIndex;

    // Element Phase Index ordered by processor local indices
    moris::Matrix< moris::IndexMat > mElementPhaseIndex;

    // Nodal Phase Index
    // Note the exact phase value is located in the geometry index.
    // Columns - Geometry Index
    // Rows - Node Index
    // If Val = 0; This means the node is not an interface node for a given geometry
    // If Val = 1; This means the node is an interface node for a given geometry
    moris::Matrix< moris::IndexMat > mInterfaceNodeFlag;

    moris::Matrix< moris::DDRMat >
    get_all_node_coordinates_loc_inds_background_mesh() const
    {
        size_t tNumNodes = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);
        moris::Matrix<moris::DDRMat> tNodeCoords(tNumNodes,mMeshData->get_spatial_dim());
        for(size_t i = 0; i< tNumNodes; i++ )
        {
            tNodeCoords.set_row(i, mMeshData->get_node_coordinate(i));
        }
        return tNodeCoords;
    }


    /*
     * Allocate space in the downard inheritance map, one for each element in background mesh
     */
    void intialize_downward_inheritance()
    {
        moris::size_t tNumElements = mMeshData->get_num_entities((moris::EntityRank)EntityRank::ELEMENT);
        MORIS_ASSERT(tNumElements!=0,"Empty Mesh Given to XTK Mesh");

        mElementDownwardInheritance = Downward_Inheritance<moris::moris_index,moris::moris_index>(tNumElements);
    }


};
}

#endif /* SRC_XTK_CL_XTK_MESH_HPP_ */
