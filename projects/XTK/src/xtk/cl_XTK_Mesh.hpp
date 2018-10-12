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
#include "cl_XTK_External_Mesh_Data.hpp"

// Assertion Includes:
#include "fn_assert.hpp"

// Linear Algebra Includes
#include "cl_Matrix.hpp"
#include "fn_print.hpp"
#include "fn_isvector.hpp"

//XTK Includes:
#include "xtk/cl_XTK_Node.hpp"
#include "xtk/cl_XTK_External_Mesh_Data.hpp"
#include "xtk/cl_XTK_Downward_Inheritance.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"

// Geometry Engine Includes
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class XTK_Mesh
{
public:
    XTK_Mesh(){};

    XTK_Mesh(moris::mtk::Mesh* aMeshData):
        mMeshData(aMeshData),
        mNodeIndexToChildMeshIndex(0,0)
    {
        intialize_downward_inheritance();
        mExternalMeshData.set_up_external_entity_data(mMeshData);
    }

    XTK_Mesh(moris::mtk::Mesh* aMeshData,
             Geometry_Engine<Real, Integer, Real_Matrix, Integer_Matrix> & aGeometryEngine):
        mMeshData(aMeshData),
        mNodeIndexToChildMeshIndex(0,0)
    {
        intialize_downward_inheritance();
        mExternalMeshData.set_up_external_entity_data(mMeshData);
    }

    /*
     * Get number of entities in the background mesh and
     * the number of entities XTK has created
     */
    Integer
    get_num_entities(enum EntityRank aEntityRank) const
    {
        // Initialize
        Integer tNumBackgroundEntities = mMeshData->get_num_entities((moris::EntityRank)aEntityRank);
        Integer tExternalEntities = mExternalMeshData.get_num_entities_external_data(aEntityRank);

        return tNumBackgroundEntities + tExternalEntities;
    }

    /*
     *  Get an offset which is the first id allocate (assumed ids are grouped)
     */
    Integer
    allocate_entity_ids(Integer         aNumReqs,
                        enum EntityRank aChildEntityRank)
    {
        return mExternalMeshData.allocate_entity_ids_external_entity_data(aNumReqs,aChildEntityRank);
    }

    /*
     * Get the first available index
     */
    Integer
    get_first_available_index(enum EntityRank aEntityRank) const
    {
        return mExternalMeshData.get_first_available_index_external_data(aEntityRank);
    }

    /*
     * Increment the first available index by aNewFirstAvailableIndex
     */
    void
    update_first_available_index(Integer         aNewFirstAvailableIndex,
                                 enum EntityRank aEntityRank)
    {
        mExternalMeshData.update_first_available_index_external_data(aNewFirstAvailableIndex, aEntityRank);
    }

    /*
     * Create a batch of new nodes
     */
    void
    batch_create_new_nodes(xtk::Cell<xtk::Pending_Node<Real, Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes)
    {
        mExternalMeshData.batch_create_new_nodes_external_data(aPendingNodes);
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
        XTK_ASSERT(mExternalMeshData.is_external_entity(aNodeIndex,EntityRank::NODE),"Provided node index needs to be one created during the decomposition process");

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
    Integer
    get_glb_entity_id_from_entity_loc_index(Integer         aEntityIndex,
                                            enum EntityRank aEntityRank) const
    {
        Integer tGlbId = 0;
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
        Integer tNumEntities = tEntityIndices.numel();
        moris::Matrix< moris::IdMat > tEntityIds(1,tNumEntities);

        for(Integer i =0; i<tNumEntities; i++)
        {
            tEntityIds(0,i) = this->get_glb_entity_id_from_entity_loc_index(tEntityIndices(i),aEntityRank);
        }
        return tEntityIds;
    }


    /*
     * Return all node coordinates ordered by local indices
     */
    moris::Matrix< Real_Matrix >
    get_all_node_coordinates_loc_inds() const
    {
        // Get counts inside of STK and in External Data (both happen in this function call)
        Integer tNumNodes = this->get_num_entities(EntityRank::NODE);
        Integer tNumBGNodes = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);

        moris::Matrix< Real_Matrix > tAllNodeCoordinates = this->get_all_node_coordinates_loc_inds_background_mesh();
        tAllNodeCoordinates.resize(tNumNodes,3);
        // Get node coordinates from external entities
        mExternalMeshData.get_all_node_coordinates_loc_inds_external_data(tNumBGNodes,tAllNodeCoordinates);
        return tAllNodeCoordinates;
    }

    /*
     * Return a coordinate matrix for the specified node indices
     */
    moris::Matrix< Real_Matrix >
    get_selected_node_coordinates_loc_inds(
            moris::Matrix< moris::IndexMat > const & aNodeIndices) const
    {

        MORIS_ERROR(moris::isvector(aNodeIndices),"Provided Node indices need to be a vector");
        // TODO: Add external entity check to see if xtk has the coordinate field or stk has it
        // Number of spatial dimensions
        Integer tSpatialDimension = mMeshData->get_spatial_dim();

        // Get number of nodes provided
        Integer tNumNodes = aNodeIndices.numel();

        // Initialize output matrix
        moris::Matrix< Real_Matrix > tSelectedNodesCoords(tNumNodes, tSpatialDimension);

        // Loop over all nodes
        enum EntityRank tEntityRank = EntityRank::NODE;

        for (Integer n = 0; n < tNumNodes; ++n)
        {
            if (mExternalMeshData.is_external_entity(aNodeIndices(n), tEntityRank))
            {
                moris::Matrix< Real_Matrix > const & tNodeCoords = mExternalMeshData.get_selected_node_coordinates_loc_inds_external_data(aNodeIndices(n));
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

        moris::Matrix<Integer_Matrix> tLocalToGlobalBM(1,tNumNodes);

        for(size_t i = 0; i<tNumNodes; i++)
        {
            tLocalToGlobalBM(i) = (Integer)mMeshData->get_glb_entity_id_from_entity_loc_index((moris_index)i,moris::EntityRank::NODE);
        }

        moris::Matrix<Integer_Matrix>  tLocalToGlobalExt = mExternalMeshData.get_local_to_global_node_map();


        tNumNodes = this->get_num_entities(EntityRank::NODE);
        size_t tFirstExtNodeInd = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);
        MORIS_ERROR(tNumNodes = tLocalToGlobalBM.numel() + tLocalToGlobalExt.numel(),"Number of nodes returned does not match the number in the map");

        // combine the two maps
        moris::Matrix<moris::IdMat> tLocalToGlobal(1,tNumNodes);

        for(Integer i = 0 ; i <tLocalToGlobalBM.numel(); i++)
        {
            tLocalToGlobal(i) = tLocalToGlobalBM(i);
        }

        for(Integer i = 0 ; i <tLocalToGlobalExt.numel(); i++)
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
        Integer tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);
        enum EntityTopology tElemTopo = get_XTK_mesh_element_topology();

        Integer tNumNodesPerElem = 0;
        if(tElemTopo == EntityTopology::TET_4)
        {
            tNumNodesPerElem = 4;
        }
        else if(tElemTopo == EntityTopology::HEXA_8)
        {
            tNumNodesPerElem = 8;
        }
        else
        {
            tNumNodesPerElem = 8;
        }

        moris::Matrix<moris::IdMat> tElementToNode(tNumElementsBG,tNumNodesPerElem);
        Integer tCount = 0;
        for(Integer i = 0; i<tElementToNode.n_rows(); i++)
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

    /*
     * Return all ids of non-intersected elements
     */
    moris::Matrix<moris::IdMat>
    get_all_non_intersected_elements() const
    {
        Integer tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);

        moris::Matrix<moris::IdMat> tElementIds(tNumElementsBG,1);
        Integer tCount = 0;
        for(Integer i = 0; i<tElementIds.numel(); i++)
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
     * Return all non-intersected element proc local indices
     */
    moris::Matrix<moris::IndexMat>
    get_all_non_intersected_elements_loc_inds() const
    {
        Integer tNumElementsBG        = this->get_num_entities(EntityRank::ELEMENT);

        moris::Matrix<moris::IdMat> tElementInds(tNumElementsBG,1);
        Integer tCount = 0;
        for(Integer i = 0; i<tElementInds.numel(); i++)
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
    register_new_downward_inheritance(Cell<std::pair<Integer,Integer>> const & aNewElementToChildMeshPairs)
    {
        Integer tNumNewPairs = aNewElementToChildMeshPairs.size();
        for(Integer i = 0; i<tNumNewPairs; i++)
        {
            mElementDownwardInheritance.register_new_inheritance_pair(aNewElementToChildMeshPairs(i).first,aNewElementToChildMeshPairs(i).second);
        }
    }

    /*
     * returns whether a given entity has any children entities
     */
    bool
    entity_has_children(Integer aEntityIndex,
                        enum EntityRank aEntityRank) const
    {
        XTK_ASSERT(aEntityRank==EntityRank::ELEMENT,"ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED");

        return mElementDownwardInheritance.has_inheritance(aEntityIndex);
    }

    /*
     * returns the child mesh index of entity with children
     */
    Integer const & child_mesh_index(Integer aEntityIndex,
                             enum EntityRank aEntityRank)
    {
        XTK_ASSERT(aEntityRank==EntityRank::ELEMENT,"ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED");

        return mElementDownwardInheritance.get_inheritance(aEntityIndex);
    }

    // -------------------------------------------------------------------
    // Functions related to setting and accessing interface node information
    // -------------------------------------------------------------------

    void
    initialize_interface_node_flags(Integer const & aNumNodes,
                                    Integer const & aNumGeometry)
    {
        mInterfaceNodeFlag = moris::Matrix< Integer_Matrix >(aNumNodes,aNumGeometry,0);
    }

    /*
     * Allocates additional space in interface node flags
     */
    void
    allocate_space_in_interface_node_flags(Integer aNumNodes,
                                           Integer aNumGeometry)
    {
        Integer tCurrentSize = mInterfaceNodeFlag.n_rows();
        mInterfaceNodeFlag.resize(tCurrentSize + aNumNodes, aNumGeometry);

        for(Integer i = tCurrentSize; i<tCurrentSize+aNumNodes; i++)
        {
            for(Integer j = 0; j<aNumGeometry; j++)
            {
                mInterfaceNodeFlag(i,j) = 0;
            }
        }

    }

    /*
     * Marks a node as an interface node for a given geometry index
     */
    void
    mark_node_as_interface_node(Integer aNodeIndex,
                                Integer aGeomIndex)
    {
        mInterfaceNodeFlag(aNodeIndex,aGeomIndex) = 1;
    }

    /*
     * Returns whether a node is an interface node for a given geometry index
     */
    bool
    is_interface_node(Integer aNodeIndex,
                      Integer aGeomIndex)
    {
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
    get_interface_nodes_loc_inds(moris::moris_index aGeometryIndex)
    {
        // initialize output
        Integer tNumNodes = this->get_num_entities(EntityRank::NODE);
        moris::Matrix< moris::IndexMat > tInterfaceNodes(1,tNumNodes);

        // keep track of how many interface nodes
        Integer tCount = 0;

        for(Integer i = 0; i<tNumNodes; i++)
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


    void
    print_interface_node_flags()
    {
        for(Integer i = 0; i<mInterfaceNodeFlag.n_rows(); i++)
        {
            std::cout<<this->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE)<<" | ";
            for(Integer j = 0; j<mInterfaceNodeFlag.n_cols(); j++)
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
    initialize_element_phase_indices(Integer const & aNumElements)
    {
        mElementPhaseIndex = moris::Matrix< Integer_Matrix >(aNumElements,1);
    }

    /*
     * Set the phase index value of element with element index. This is relative to each geometry.
     */
    void
    set_element_phase_index(Integer const & aElementIndex,
                            Integer const & aElementPhaseIndex)
    {
        mElementPhaseIndex(aElementIndex,0) = aElementPhaseIndex;
    }

    /*
     * Get the phase index value of element with element index
     */
    Integer const &
    get_element_phase_index(Integer const & aElementIndex) const
    {
        return (mElementPhaseIndex)( aElementIndex, 0 );
    }

    moris::Matrix< Integer_Matrix >
    get_element_phase_inds(Cell<Integer> const & aElementInds )
    {
        moris::Matrix< Integer_Matrix > tElementPhasesInds(1,aElementInds.n_cols());

        for(Integer i = 0; i < aElementInds.n_cols(); i++)
        {
            tElementPhasesInds(0,i) = get_element_phase_index(aElementInds(0,i));
        }

        return tElementPhasesInds;
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
    enum EntityTopology
    get_XTK_mesh_element_topology() const
    {
        enum EntityTopology tElementTopology = EntityTopology::INVALID;
        moris::Matrix<  moris::IndexMat  > tElementNodes = mMeshData->get_entity_connected_to_entity_loc_inds(0,(moris::EntityRank)EntityRank::ELEMENT, (moris::EntityRank)EntityRank::NODE);
        if(tElementNodes.numel() == 8 && moris::isvector(tElementNodes))
        {
            tElementTopology = EntityTopology::HEXA_8;

        }
        else if (tElementNodes.numel() == 4 && moris::isvector(tElementNodes))
        {
            tElementTopology = EntityTopology::TET_4;
        }
        else
        {
            XTK_ERROR<<"Topology not recognized in parent mesh";
        }

        return tElementTopology;
    }
private:
    // Background mesh data
    moris::mtk::Mesh* mMeshData;

    // External Entity information
    // The background mesh remains constant, and every new entity created is stored within XTK
    // child meshes
    Mesh_External_Entity_Data<Real, Integer, Real_Matrix, Integer_Matrix> mExternalMeshData;

    // Downward inheritance pairs (links elements in XTK mesh to indices in Child Meshes)
    Downward_Inheritance<Integer,Integer> mElementDownwardInheritance;

    // Associate external node indices to the child meshes they belong to
    // Row - External node index
    // Col - Child Mesh Index
    moris::Matrix< moris::IndexMat > mNodeIndexToChildMeshIndex;

    // Element Phase Index ordered by processor local indices
    moris::Matrix< Integer_Matrix > mElementPhaseIndex;

    // Nodal Phase Index
    // Note the exact phase value is located in the geometry index.
    // Columns - Geometry Index
    // Rows - Node Index
    // If Val = 0; This means the node is not an interface node for a given geometry
    // If Val = 1; This means the node is an interface node for a given geometry
    moris::Matrix< Integer_Matrix > mInterfaceNodeFlag;

    moris::Matrix< Real_Matrix >
    get_all_node_coordinates_loc_inds_background_mesh() const
    {
        size_t tNumNodes = mMeshData->get_num_entities((moris::EntityRank)EntityRank::NODE);
        moris::Matrix<Real_Matrix> tNodeCoords(tNumNodes,mMeshData->get_spatial_dim());
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
        Integer tNumElements = mMeshData->get_num_entities((moris::EntityRank)EntityRank::ELEMENT);
        XTK_ASSERT(tNumElements!=0,"Empty Mesh Given to XTK Mesh");

        mElementDownwardInheritance = Downward_Inheritance<Integer,Integer>(tNumElements);
    }


};
}

#endif /* SRC_XTK_CL_XTK_MESH_HPP_ */
