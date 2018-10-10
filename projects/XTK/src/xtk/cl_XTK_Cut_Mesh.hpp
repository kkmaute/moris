/*
 * cl_XTK_Mesh.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_CUT_MESH_HPP_
#define SRC_XTK_CL_XTK_CUT_MESH_HPP_


// Standard Include
#include <memory>   // Shared ptrs
#include <mpi.h>

// XTKL: Linear Algebra includes
#include "linalg/cl_XTK_Matrix.hpp"

// XTKL: Logging and Assertion Includes
#include "assert/fn_xtk_assert.hpp"
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Mesh includes
#include "mesh/cl_Mesh_Data.hpp"  // For packaging only locally owned entities
#include "mesh/cl_Mesh_Enums.hpp" // For entity rank

// XTKL: Xtk includes
#include "xtk/cl_XTK_Child_Mesh.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "xtk/cl_XTK_Downward_Inheritance.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Cut_Mesh
{
public:
    Cut_Mesh(){};
    Cut_Mesh(Integer aModelDim) :
                 mNumberOfChildrenMesh(0),
                 mNumEntities(4,0),
                 mChildrenMeshes(1, Child_Mesh_Test<Real,Integer,Real_Matrix,Integer_Matrix>()),
                 mActiveSample(false),
                 mChildElementTopo(EntityTopology::TET_4)
    {

    }


    Cut_Mesh(Integer aNumSimpleMesh,
             Integer aModelDim) :
            mNumberOfChildrenMesh(aNumSimpleMesh),
            mChildrenMeshes(aNumSimpleMesh),
            mActiveSample(false),
            mChildElementTopo(EntityTopology::TET_4)
    {
        for(Integer i = 0; i <aNumSimpleMesh; i++)
        {
            mChildrenMeshes(i) = Child_Mesh_Test<Real,Integer,Real_Matrix,Integer_Matrix>();
        }
    }



    ~Cut_Mesh()
    {

    }

    /**
     * Allocate space for new children meshes
     * @param aNumSimpleMesh - number of simple meshes to allocate space for
     * @param aModelDim - model dimension (1D, 2D or 3D)
     */
    void inititalize_new_child_meshes(Integer aNumNewChildMesh, Integer aModelDim)
    {

        mNumberOfChildrenMesh = aNumNewChildMesh + mNumberOfChildrenMesh;
        mChildrenMeshes.resize(mNumberOfChildrenMesh, Child_Mesh_Test<Real,Integer,Real_Matrix,Integer_Matrix>());
    }

    /**
     * generate_templated_mesh genereates a templated mesh on a given child mesh with the specified edge intersection permutation (if required by a template)
     *
     * @param[in] aChildMeshIndex  - specifies which simple mesh to generate the template on
     * @param[in] aTemplate - Specifies the template used
     * @param[in] aPermutation - If required, specifies the edge intersection pattern
     */
    void generate_templated_mesh(Integer           aChildMeshIndex,
                                 enum TemplateType aTemplate)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds. Consider allocating more space");
        mChildrenMeshes(aChildMeshIndex).modify_child_mesh(aTemplate);


        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }

    /**
     * generate_templated_mesh tells all simple meshes to generate a specified templated mesh
     *
     * @param[in] aTemplate    - Specifies the template to use
     */
    void generate_templated_mesh(enum TemplateType aTemplate)
    {
        for (Integer i = 0; i < mNumberOfChildrenMesh; i++)
        {
            mChildrenMeshes(i).modify_child_mesh(aTemplate);
        }

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }

    /*
     * Generate a tempalted mesh for a subset of the children meshes.
     */
    void generate_templated_mesh(moris::Matrix< Integer_Matrix > const &aChildMeshIndices,
                                 enum TemplateType aTemplate)
    {
        for (Integer i = 0; i < aChildMeshIndices.n_cols(); i++)
        {
            mChildrenMeshes(aChildMeshIndices(0,i)).modify_child_mesh(aTemplate);
        }

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }


    /*
     * Converts existing tet4 child mesh to tet10s
     *
     */
    void convert_cut_mesh_to_tet10s()
    {
        for(Integer i = 0; i<get_num_simple_meshes(); i++)
        {
            mChildrenMeshes(i).convert_tet4_to_tet10_child();
        }

        mChildElementTopo = EntityTopology::TET_10;
    }

    /**
     * @ brief Sets up template ancestry with parametric information
     * @param[in] aChildMeshIndex        - simple mesh index
     * @param[in] aTemplate       - specifies the template ancestry to use
     * @param[in] aParentEntities - cell of row vectors of parent entity indices
     */
    void initialize_new_mesh_from_parent_element(Integer                                  aChildMeshIndex,
                                                 enum TemplateType                        aTemplate,
                                                 moris::Matrix< moris::IndexMat >       & aNodeIndices,
                                                 Cell<moris::Matrix< moris::IndexMat >> & aParentEntities)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");

        // Construct a template and initialize this new mesh with the template
        moris::Matrix< Integer_Matrix > tParentEdgeRanks(1,aParentEntities(1).numel(),1);

        moris::Matrix< Integer_Matrix > tParentFaceRanks(1,aParentEntities(2).numel(),2);

        moris::Matrix< Integer_Matrix > tInterfaceSides(1,1,std::numeric_limits<Integer>::max());

        mChildrenMeshes(aChildMeshIndex) = Child_Mesh_Test<Real,Integer,Real_Matrix,Integer_Matrix>(aParentEntities(3)(0,0),
                                                                                                    aNodeIndices,
                                                                                                    aNodeIndices,
                                                                                                    aParentEntities(1),
                                                                                                    tParentEdgeRanks,
                                                                                                    aParentEntities(2),
                                                                                                    tParentFaceRanks,
                                                                                                    tInterfaceSides);

        // set parent element parametric coordinate
        switch(aTemplate)
        {
            case TemplateType::HEX_8:
            {
                const moris::Matrix< moris::DDRMat > tParamCoords(
                {{-1.0, -1.0, -1.0},
                 { 1.0, -1.0, -1.0},
                 { 1.0,  1.0, -1.0},
                 {-1.0,  1.0, -1.0},
                 {-1.0, -1.0,  1.0},
                 { 1.0, -1.0,  1.0},
                 { 1.0,  1.0,  1.0},
                 {-1.0,  1.0,  1.0}});

                // Add hex 8 parametric coordinates
                mChildrenMeshes(aChildMeshIndex).allocate_parametric_coordinates(8);
                mChildrenMeshes(aChildMeshIndex).add_node_parametric_coordinate(aNodeIndices,tParamCoords);


                break;
            }
            case TemplateType::TET_4:
            {
                MORIS_ERROR(0,"TET_4 parametric coordinates not implemented");
                break;
            }

            default:
            {
                MORIS_ERROR(0,"Parent element template type not recognized");
            }

         }



    }



    //Modify Template mesh is only useful for unit tests without node vars
    void modify_templated_mesh(Integer aChildMeshIndex,
                               enum TemplateType aTemplate)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).modify_child_mesh(aTemplate);

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    void
    init_intersect_connectivity(moris::Matrix< Integer_Matrix > const & aChildMeshIndices)
    {
        Integer tSizeActive =aChildMeshIndices.n_cols();
        for(Integer i = 0; i <tSizeActive; i++)
        {
            mChildrenMeshes(aChildMeshIndices(0,i)).init_intersect_connectivity();
        }

    }

    void
    init_intersect_connectivity()
    {
        for(Integer i = 0; i <mNumberOfChildrenMesh; i++)
        {
            mChildrenMeshes(i).init_intersect_connectivity();
        }

    }


    /**
     * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
     *       - 1 means the provided aDPrime1Ind is an XTK index
     *
     * aDPrime2Ind must be XTK local index
     */
    void
    add_entity_to_intersect_connectivity(Integer aChildMeshIndex,
                                   Integer aDPrime1Ind,
                                   Integer aDPrime2Ind,
                                   Integer aReturnType)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).add_entity_to_intersect_connectivity(aDPrime1Ind, aDPrime2Ind, aReturnType);
    }

    /*
     * Set node indices in a child mesh
     * @param[in] aChildMeshIndex - Child mesh index
     * @param[in] aNodeInd - Node indices
     */
    void
    set_node_index(Integer const &                aChildMeshIndex,
                   moris::Matrix< Integer_Matrix > & aNodeInd)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).add_node_indices(aNodeInd);

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }

    /*
     * Set node ids in a child mesh
     * @param[in] aChildMeshIndex - Child mesh index
     * @param[in] aNodeInd - Node ids
     */
    void
    set_node_ids(Integer const & aChildMeshIndex,
                 moris::Matrix< moris::IdMat > & aNodeIds)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_node_ids(aNodeIds);

        // Meshes have changed so counts are wrong and need to be updated before use
        mConsistentCounts = false;
    }



    /*
     * Set node ids in a child mesh
     * @param[in] aChildMeshIndex - Child mesh index
     * @param[in] aNodeInd - Node ids
     */
    void
    add_node_ids(Integer const & aChildMeshIndex,
                 moris::Matrix< Integer_Matrix > & aNodeIds)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).add_node_ids(aNodeIds);

        // Meshes have changed so counts are wrong and need to be updated before use
        mConsistentCounts = false;
    }

    /*
     * Set element ids in a child mesh
     * @param[in] aChildMeshIndex - Child mesh index
     * @param[in] aElementIdOffset - Element Id offset
     */
    void
    set_child_element_ids(Integer const &   aChildMeshIndex,
                          moris::moris_id & aElementIdOffset )
    {
        mChildrenMeshes(aChildMeshIndex).set_child_element_ids(aElementIdOffset);
    }

    /*
     * Set element indicess in a child mesh
     * @param[in] aChildMeshIndex - Child mesh index
     * @param[in] aElementIdOffset - Element Ind offset
     */
    void
    set_child_element_inds(Integer const &      aChildMeshIndex,
                           moris::moris_index & aElementIndOffset )
    {
        mChildrenMeshes(aChildMeshIndex).set_child_element_inds(aElementIndOffset);
    }


    /*
     * Get element Ids in a child mesh
     */
    moris::Matrix< moris::IdMat > const &
    get_element_ids(Integer const & aChildMeshIndex)
    {
        return  mChildrenMeshes(aChildMeshIndex).get_element_ids();
    }

    /*
     * Get element Ids in a child mesh
     */
    moris::Matrix< moris::IdMat >
    get_all_element_ids()
    {
        Integer tNumElems = this->get_num_entities(EntityRank::ELEMENT);
        moris::Matrix< moris::IdMat > tElementIds(1,tNumElems);

        Integer tCount = 0;
        for(Integer iCM = 0; iCM<this->get_num_simple_meshes(); iCM++)
        {
            moris::Matrix< moris::IdMat > const & tCMIds = this->get_element_ids(iCM);
            for(Integer iE = 0; iE<tCMIds.numel(); iE++)
            {
                tElementIds(tCount) = tCMIds(iE);
                tCount++;
            }

        }
        return  tElementIds;
    }

    /*
     * Get element Inds from a child mesh
     */
    Cell<Integer> const &
    get_element_inds(Integer const & aChildMeshIndex) const
    {
        return  mChildrenMeshes(aChildMeshIndex).get_element_inds();
    }



    void set_pending_node_index_pointers(Integer aChildMeshIndex, Integer* aNodeIndPtr)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_pending_node_index_pointers(aNodeIndPtr);
    }


    void set_pending_node_index_pointers(Integer                              aChildMeshIndex,
                                         Cell<Integer*> &                     aNodeIndPtr,
                                         moris::Matrix< Real_Matrix > const & aParametricCoordinates)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_pending_node_index_pointers(aNodeIndPtr,aParametricCoordinates);
    }


    void set_pending_node_index_pointers_with_dx_dp(Integer aChildMeshIndex,
                                                    Integer* & aNodeIndPtr,
                                                    moris::Matrix< Integer_Matrix >* & aDxDpNode,
                                                    moris::Matrix< Integer_Matrix >* & aADVIndices)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_pending_node_index_pointers_with_dx_dp(aNodeIndPtr,aDxDpNode,aADVIndices);
    }

    void set_pending_node_index_pointers_with_dx_dp(Integer aChildMeshIndex,
                                                    Cell<Integer*> & aNodeIndPtr,
                                                    Cell<moris::Matrix< Integer_Matrix >*> & aDxDpNode,
                                                    Cell<moris::Matrix< Integer_Matrix >*> & aADVIndices)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_pending_node_index_pointers_with_dx_dp(aNodeIndPtr,aDxDpNode,aADVIndices);
    }

   void retrieve_pending_node_inds()
    {
        for(Integer i = 0; i < mNumberOfChildrenMesh; i++)
        {
            mChildrenMeshes(i).retrieve_pending_node_inds();
        }
        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }


    /*
     * Add node indices to a child mesh
     */
    void
    add_node_index(Integer aChildMeshIndex,
                   moris::Matrix< Integer_Matrix > & aNodeInd)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).add_node_ind(aNodeInd);
    }


    /*
     * Get node indices of a child mesh
     */
    Integer
    get_node_index(Integer aChildMeshIndex,
                   Integer aIndex) const
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_node_ind(aIndex);
    }

    /*
     * Get node ids from a child mesh
     */
    Integer get_node_id(Integer aChildMeshIndex,
                        Integer aIndex) const
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_node_id(aIndex);
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_indices(Integer aChildMeshIndex)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_node_indices();
    }

    EntityTopology
    get_child_element_topology()
    {
        return mChildElementTopo;
    }

    Integer get_num_simple_meshes() const
    {
        return mNumberOfChildrenMesh;
    }

    Integer get_num_entities(Integer aChildMeshIndex, enum EntityRank aEntityRank) const
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_num_entities(aEntityRank);
    }

    Integer get_num_entities(enum EntityRank aEntityRank) const
    {
        // Make sure counts are up to date
        get_entity_counts();

        // Make sure the counts are up to date
        XTK_ASSERT(mConsistentCounts, "Make sure to call get_entity_counts otherwise the mNumEntities variable is out dated and garbage is returned");

        return mNumEntities((Integer) aEntityRank);
    }

    Integer get_num_interface_nodes() const
    {
        Integer tNumInterfaceNodes = 0;
        for(Integer i = 0; i<get_num_simple_meshes();i++)
        {
            tNumInterfaceNodes+=mChildrenMeshes(i).get_num_interface_nodes();
        }

        return tNumInterfaceNodes;
    }

    Integer get_num_entities_connected_to_entity(Integer aChildMeshIndex,
                                                 enum EntityRank aEntity,
                                                 enum EntityRank aEntityPrime,
                                                 Integer aEntityInd) const
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_num_entities_connected_to_entity(aEntity, aEntityPrime, aEntityInd);
    }

    moris::Matrix< Integer_Matrix >
    get_entities_connected_to_entity(Integer aChildMeshIndex,
                                     enum EntityRank aEntity,
                                     enum EntityRank aEntityPrime,
                                     Integer aEntityIndex,
                                     Integer aReturnType) const
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_entities_connected_to_entity(aEntity, aEntityPrime, aEntityIndex, aReturnType);
    }


    Integer get_parent_element_index(Integer aChildMeshIndex)
    {
        XTK_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_parent_element_index();
    }

    /*
     * Returns all the children elements connected to a provided face index in a single child mesh
     *
     * @param[in]  aChildMeshIndex         - Child Mesh Index
     * @param[in]  aParentFaceIndex        - Parent Face Index
     * @param[out] aChildrenElementId      - Children Element Global Id
     * @param[out] aChildrenElementCMInd   - Child element index local to child mesh
     * @param[out] aFaceOrdinal            - Face Ordinal relative to element
     */
    void get_child_elements_connected_to_parent_face(Integer const & aChildMeshIndex,
                                                     Integer const & aParentFaceIndex,
                                                     moris::Matrix< Integer_Matrix > & aChildrenElementId,
                                                     moris::Matrix< Integer_Matrix > & aChildrenElementCMInd,
                                                     moris::Matrix< Integer_Matrix > & aFaceOrdinal) const
    {
        mChildrenMeshes(aChildMeshIndex).get_child_elements_connected_to_parent_face(aParentFaceIndex,aChildrenElementId,aChildrenElementCMInd,aFaceOrdinal);
    }

    /*
     * Get element processor local index from a Child mesh local index
     */
    moris::Matrix< Integer_Matrix > const &
    get_child_element_inds(Integer const & aChildMeshIndex) const
    {
        return mChildrenMeshes(aChildMeshIndex).get_element_inds();
    }


    moris::Matrix< Integer_Matrix >
    get_child_element_inds(Integer const & aChildMeshIndex,
                         moris::Matrix< Integer_Matrix > const &  aElementIndex) const
    {

        Integer tNumElements = aElementIndex.n_cols();
        moris::Matrix< Integer_Matrix > tElementInds = aElementIndex.create(1,tNumElements);

        for(Integer i= 0;  i<tNumElements; i++)
        {
            (*tElementInds)(0,i) = mChildrenMeshes(aChildMeshIndex).get_child_element_ind(aElementIndex(0,i));
        }

        return tElementInds;
    }


    Integer get_entity_phase_index( Integer    const & aMeshIndex,
                                    EntityRank const & aEntityRank,
                                    Integer    const & aEntityIndex) const
    {
        return mChildrenMeshes(aMeshIndex).get_entity_phase_index(aEntityRank,aEntityIndex);
    }

    // Outputting Mesh information


    void pack_cut_mesh_by_phase(Integer const & aMeshIndex,
                                Integer const & aNumPhases,
                                Cell<moris::Matrix< Integer_Matrix >> & aElementIds,
                                Cell<moris::Matrix< Integer_Matrix >> & aElementCMInds) const
    {

        mChildrenMeshes(aMeshIndex).pack_child_mesh_by_phase(aNumPhases,aElementCMInds,aElementIds);
    }

    void
    pack_interface_sides(Integer const & aMeshIndex,
                         Output_Options<Integer> const & aOutputOptions,
                         moris::Matrix< Integer_Matrix > & aElementIds,
                         moris::Matrix< Integer_Matrix > & aSideOrdinals) const
    {
        mChildrenMeshes(aMeshIndex).pack_interface_sides(aOutputOptions,aElementIds,aSideOrdinals);
    }

    /*
     * Get full element to node connectivity (ids). Full here means for all children meshes
     */
    moris::Matrix<moris::IdMat>
    get_full_element_to_node_glob_ids()
    {
        EntityTopology tChildElementTopo = this->get_child_element_topology();
        Integer        tNumElements = this->get_num_entities(EntityRank::ELEMENT);
        Integer tNumNodesPerElem = 0;
        if(tChildElementTopo == EntityTopology::TET_4)
        {
            tNumNodesPerElem = 4;
        }
        else
        {
            MORIS_ERROR(0,"Not implemented");
        }

        moris::Matrix<moris::IdMat> tElementToNodeIds(tNumElements,tNumNodesPerElem);
        Integer tCount = 0;
        for(Integer i = 0; i<this->get_num_simple_meshes(); i++)
        {
           moris::Matrix<moris::IdMat> tElementToNodeIdsCM = mChildrenMeshes(i).get_element_to_node_global();
           for(Integer j = 0; j<tElementToNodeIdsCM.n_rows(); j++)
           {
               tElementToNodeIds.set_row(tCount, tElementToNodeIdsCM.get_row(j));
               tCount++;
           }

        }

        return tElementToNodeIds;
    }


    void
    set_aux_connectivity(Integer aChildMeshIndex,
    moris::Matrix< Integer_Matrix > const & aAuxConn)
    {
        XTK_ASSERT(aChildMeshIndex<mNumberOfChildrenMesh,"The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_aux_connectivity(aAuxConn);
    }

    moris::Matrix< Integer_Matrix > const &
    get_aux_connectivity(Integer aChildMeshIndex)
    {
        XTK_ASSERT(aChildMeshIndex<mNumberOfChildrenMesh,"The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_aux_connectivity();
    }



    void print_node_to_entity_connectivity(Integer const & aChildMeshIndex) const
    {
        mChildrenMeshes(aChildMeshIndex).print_node_to_entity_connectivity();
    }

    void print_node_to_entity_connectivity_with_ancestry(Integer const & aChildMeshIndex) const
    {
        mChildrenMeshes(aChildMeshIndex).print_node_to_entity_connectivity_with_ancestry();
    }

    Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> const & get_child_mesh(Integer const & aChildMeshIndex) const
        {
            return mChildrenMeshes(aChildMeshIndex);
        }

    Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> & get_child_mesh(Integer const & aChildMeshIndex)
        {
            return mChildrenMeshes(aChildMeshIndex);
        }

private:
    Integer mNumberOfChildrenMesh;

    mutable Cell<Integer> mNumEntities;
    Cell<Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix>> mChildrenMeshes;

    // Number of entities total and if the current tally is accurate
    mutable bool mConsistentCounts;

    // Sampling member variables
    bool mActiveSample;
    Integer mActiveSampleRank;

    // TOPOLOGY TYPE
    EntityTopology mChildElementTopo;

private:

    void
    get_entity_counts() const
    {
        // If the counts aren't up to date then recount them.
        if(!mConsistentCounts)
        {
            enum EntityRank tRank = EntityRank::NODE;
            for(Integer r = 0; r<(Integer)EntityRank::END_ENUM; r++)
            {
                Integer tCount = 0;
                // Cast r to enum
                tRank = static_cast<EntityRank>(r);
                for(Integer m = 0; m<mNumberOfChildrenMesh; m++)
                {
                    tCount += mChildrenMeshes(m).get_num_entities(tRank);
                }
                mNumEntities(r) = tCount;
            }

            // This should be the only spot where this is changed to true
            mConsistentCounts = true;
        }
    }
};
}



#endif /* SRC_XTK_CL_XTK_CUT_MESH_HPP_ */
