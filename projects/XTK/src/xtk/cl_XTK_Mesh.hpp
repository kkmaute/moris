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
#include "mesh/cl_Mesh_Data.hpp"

// Assertion Includes:
#include "assert/fn_xtk_assert.hpp"

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

    XTK_Mesh(std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> aMeshData):
        mMeshData(aMeshData)
    {
        intialize_downward_inheritance();
    }

    XTK_Mesh(std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> aMeshData,
             Geometry_Engine<Real, Integer, Real_Matrix, Integer_Matrix> & aGeometryEngine):
        mMeshData(aMeshData)
    {
        intialize_downward_inheritance();
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
        mInterfaceNodeFlag = moris::Mat_New<Integer, Integer_Matrix>(aNumNodes,aNumGeometry,0);
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
     * Returns all interface node indices for a given geometry index
     *
     */
    moris::Mat_New<Integer, Integer_Matrix>
    get_interface_nodes_loc_inds(Integer aGeomIndex)
    {
        // initialize output
        Integer tNumNodes = mMeshData->get_num_entities(EntityRank::NODE);
        moris::Mat_New<Integer, Integer_Matrix> tInterfaceNodes(1,tNumNodes);

        // keep track of how many interface nodes
        Integer tCount = 0;
        for(Integer i = 0; i<tNumNodes; i++)
        {
            if(is_interface_node(i,aGeomIndex))
            {
                tInterfaceNodes(0,tCount) = i;
                tCount++;
            }
        }

        tInterfaceNodes.resize(1,tCount);
        return tInterfaceNodes;
    }


    /*
     * Get interface node for all geometries. returns the node ids rather than proc indices
     */
    Cell<moris::Mat_New<Integer, Integer_Matrix>>
    get_interface_nodes_glb_ids()
    {
        // initialize output
        Integer tNumNodes = mMeshData->get_num_entities(EntityRank::NODE);
        Integer tNumGeoms = mInterfaceNodeFlag.n_cols();
        Cell<moris::Mat_New<Integer, Integer_Matrix>> tInterfaceNodes(tNumGeoms);

        // iterate through geometries
        for(Integer iG = 0; iG<tNumGeoms; iG++)
        {
            // keep track of how many interface nodes
            Integer tCount = 0;

            tInterfaceNodes(iG) = moris::Mat_New<Integer, Integer_Matrix>(1,tNumNodes);

            for(Integer i = 0; i<mInterfaceNodeFlag.n_rows(); i++)
            {
                if(is_interface_node(i,iG))
                {
                    tInterfaceNodes(iG)(0,tCount) = mMeshData->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
                    tCount++;
                }
            }

            tInterfaceNodes(iG).resize(1,tCount);
        }
        return tInterfaceNodes;
    }

    Cell<moris::Mat_New<Integer, Integer_Matrix>>
    get_interface_nodes_loc_inds()
    {
        // initialize output
        Integer tNumNodes = mMeshData->get_num_entities(EntityRank::NODE);
        Integer tNumGeoms = mInterfaceNodeFlag.n_cols();
        Cell<moris::Mat_New<Integer,Integer_Matrix>> tInterfaceNodes(tNumGeoms);

        // iterate through geometries
        for(Integer iG = 0; iG<tNumGeoms; iG++)
        {
            // keep track of how many interface nodes
            Integer tCount = 0;

            tInterfaceNodes(iG) = moris::Mat_New<Integer, Integer_Matrix>(1,tNumNodes);

            for(Integer i = 0; i<mInterfaceNodeFlag.n_cols(); i++)
            {
                if(is_interface_node(i,iG))
                {
                    tInterfaceNodes(0,tCount) = i;
                    tCount++;
                }
            }

            tInterfaceNodes(iG).resize(1,tCount);
        }
        return tInterfaceNodes;
    }


    void
    print_interface_node_flags()
    {
        for(Integer i = 0; i<mInterfaceNodeFlag->n_rows(); i++)
        {
            std::cout<<mMeshData->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE)<<" | ";
            for(Integer j = 0; j<mInterfaceNodeFlag->n_cols(); j++)
            {
                std::cout<<(*mInterfaceNodeFlag)(i,j)<<" ";
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
        mElementPhaseIndex = moris::Mat_New<Integer, Integer_Matrix>(aNumElements,1);
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
    get_element_phase_index(Integer const & aElementIndex)
    {
        return (mElementPhaseIndex)( aElementIndex, 0 );
    }

    moris::Mat_New<Integer, Integer_Matrix>
    get_element_phase_inds(Cell<Integer> const & aElementInds )
    {
        moris::Mat_New<Integer, Integer_Matrix> tElementPhasesInds(1,aElementInds.n_cols());

        for(Integer i = 0; i < aElementInds.n_cols(); i++)
        {
            tElementPhasesInds(0,i) = get_element_phase_index(aElementInds(0,i));
        }

        return tElementPhasesInds;
    }


    // -------------------------------------------------------------------
    // Access underlying mesh data functions
    mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> &
    get_mesh_data()
    {
        return *mMeshData;
    }

    mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> const &
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
        moris::Mat_New<Integer, Integer_Matrix> tElementNodes = mMeshData->get_entity_connected_to_entity_loc_inds(0,EntityRank::ELEMENT, EntityRank::NODE);
        if(tElementNodes.n_cols() == 8)
        {
            tElementTopology = EntityTopology::HEXA_8;
        }
        else if (tElementNodes.n_cols() == 4)
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
    std::shared_ptr<mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix>> mMeshData;

    // Downward inheritance pairs (links elements in XTK mesh to indices in Child Meshes)
    Downward_Inheritance<Integer,Integer> mElementDownwardInheritance;

    // Element Phase Index ordered by processor local indices
    moris::Mat_New<Integer, Integer_Matrix> mElementPhaseIndex;

    // Nodal Phase Index
    // Note the exact phase value is located in the geometry index.
    // Columns - Geometry Index
    // Rows - Node Index
    // If Val = 0; This means the node is not an interface node for a given geometry
    // If Val = 1; This means the node is an interface node for a given geometry
    moris::Mat_New<Integer, Integer_Matrix> mInterfaceNodeFlag;


    /*
     * Allocate space in the downard inheritance map, one for each element in background mesh
     */
    void intialize_downward_inheritance()
    {
        Integer tNumElements = mMeshData->get_num_entities(EntityRank::ELEMENT);
        XTK_ASSERT(tNumElements!=0,"Empty Mesh Given to XTK Mesh");

        mElementDownwardInheritance = Downward_Inheritance<Integer,Integer>(tNumElements);
    }


};
}

#endif /* SRC_XTK_CL_XTK_MESH_HPP_ */
