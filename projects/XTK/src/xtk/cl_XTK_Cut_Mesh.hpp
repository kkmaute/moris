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
#include "cl_Matrix.hpp"

// XTKL: Logging and Assertion Includes

#include "cl_Cell.hpp"

// XTKL: Mesh includes
#include "cl_Mesh_Enums.hpp" // For entity rank

// XTKL: Xtk includes
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Output_Options.hpp"
#include "cl_XTK_Downward_Inheritance.hpp"
#include "cl_XTK_Interface_Element.hpp"

namespace xtk
{
class Cut_Mesh
{
public:
    Cut_Mesh(){};
    Cut_Mesh(moris::uint aModelDim) :
                 mNumberOfChildrenMesh(0),
                 mChildrenMeshes(1, Child_Mesh()),
                 mNumEntities(4,0),
                 mChildElementTopo(CellTopology::TET4)
    {

    }


    Cut_Mesh(moris::size_t aNumSimpleMesh,
             moris::size_t aModelDim) :
            mNumberOfChildrenMesh(aNumSimpleMesh),
            mChildrenMeshes(aNumSimpleMesh),
            mChildElementTopo(CellTopology::TET4)
    {
        for(moris::size_t i = 0; i <aNumSimpleMesh; i++)
        {
            mChildrenMeshes(i) = Child_Mesh();
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
    void inititalize_new_child_meshes(moris::size_t aNumNewChildMesh, moris::size_t aModelDim)
    {

        mNumberOfChildrenMesh = aNumNewChildMesh + mNumberOfChildrenMesh;
        mChildrenMeshes.resize(mNumberOfChildrenMesh, Child_Mesh());
    }

    /**
     * generate_templated_mesh genereates a templated mesh on a given child mesh with the specified edge intersection permutation (if required by a template)
     *
     * @param[in] aChildMeshIndex  - specifies which simple mesh to generate the template on
     * @param[in] aTemplate - Specifies the template used
     * @param[in] aPermutation - If required, specifies the edge intersection pattern
     */
    void generate_templated_mesh(moris::size_t           aChildMeshIndex,
                                 enum TemplateType aTemplate)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds. Consider allocating more space");
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
        for (moris::size_t i = 0; i < mNumberOfChildrenMesh; i++)
        {
            mChildrenMeshes(i).modify_child_mesh(aTemplate);
        }

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }

    /*
     * Generate a tempalted mesh for a subset of the children meshes.
     */
    void generate_templated_mesh(moris::Matrix< moris::DDSTMat > const &aChildMeshIndices,
                                 enum TemplateType aTemplate)
    {
        for (moris::size_t i = 0; i < aChildMeshIndices.n_cols(); i++)
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
        for(moris::size_t i = 0; i<get_num_child_meshes(); i++)
        {
            mChildrenMeshes(i).convert_tet4_to_tet10_child();
        }

        mChildElementTopo = CellTopology::TET10;
    }

    /**
     * @ brief Sets up template ancestry with parametric information
     * @param[in] aChildMeshIndex        - simple mesh index
     * @param[in] aTemplate       - specifies the template ancestry to use
     * @param[in] aParentEntities - cell of row vectors of parent entity indices
     */
    void initialize_new_mesh_from_parent_element(moris::size_t                            aChildMeshIndex,
                                                 enum TemplateType                        aTemplate,
                                                 moris::Matrix< moris::IndexMat >       & aNodeIndices,
                                                 Cell<moris::Matrix< moris::IndexMat >> & aParentEntities)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");

        // Construct a template and initialize this new mesh with the template

        // Set all node parent ranks to EntityRank::NODE which = 0;
        moris::Matrix< moris::DDSTMat > tElementNodeParentRanks(1,8);
        tElementNodeParentRanks.fill(0);

        // Set all edge parent ranks to EntityRank::EDGE which = 1;
        moris::Matrix< moris::DDSTMat > tParentEdgeRanks(1,aParentEntities(1).numel());
        tParentEdgeRanks.fill(1);

        // Set all node parent ranks to EntityRank::FACE which = 2;
        moris::Matrix< moris::DDSTMat > tParentFaceRanks(1,aParentEntities(2).numel());
        tParentFaceRanks.fill(2);

        // No interface sides
        moris::Matrix< moris::DDSTMat > tInterfaceSides(1,1);
        tInterfaceSides.fill(std::numeric_limits<moris::size_t>::max());

        // Note for this: child mesh node indices, parent node indices, and element to node connectivity are
        // the same thing. This is why aNodeIndices appears 3 times in this call.
        mChildrenMeshes(aChildMeshIndex) = Child_Mesh(aParentEntities(3)(0,0),
                                                           aNodeIndices,
                                                           aNodeIndices,
                                                           tElementNodeParentRanks,
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
                mChildrenMeshes(aChildMeshIndex).allocate_parametric_coordinates(8,3);
                mChildrenMeshes(aChildMeshIndex).add_node_parametric_coordinate(aNodeIndices,tParamCoords);


                break;
            }
            case TemplateType::TET_4:
            {
                const moris::Matrix< moris::DDRMat > tParamCoords(
                {{ 1.0, 0.0, 0.0, 0.0},
                 { 0.0, 1.0, 0.0, 0.0},
                 { 0.0, 0.0, 1.0, 0.0},
                 { 0.0, 0.0, 0.0, 1.0}});

                // Add tetra parametric coordinates
                mChildrenMeshes(aChildMeshIndex).allocate_parametric_coordinates(4,4);
                mChildrenMeshes(aChildMeshIndex).add_node_parametric_coordinate(aNodeIndices,tParamCoords);
                break;
            }

            default:
            {
                MORIS_ERROR(0,"Parent element template type not recognized");
            }

         }



    }



    //Modify Template mesh is only useful for unit tests without node vars
    void modify_templated_mesh(moris::size_t aChildMeshIndex,
                               enum TemplateType aTemplate)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).modify_child_mesh(aTemplate);

        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }
    void
    init_intersect_connectivity(moris::Matrix< moris::IndexMat > const & aChildMeshIndices)
    {
        moris::size_t tSizeActive =aChildMeshIndices.n_cols();
        for(moris::size_t i = 0; i <tSizeActive; i++)
        {
            mChildrenMeshes(aChildMeshIndices(0,i)).init_intersect_connectivity();
        }

    }

    void
    init_intersect_connectivity()
    {
        for(moris::size_t i = 0; i <mNumberOfChildrenMesh; i++)
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
    add_entity_to_intersect_connectivity(moris::size_t aChildMeshIndex,
                                         moris::size_t aDPrime1Ind,
                                         moris::size_t aDPrime2Ind,
                                         moris::size_t aReturnType)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).add_entity_to_intersect_connectivity(aDPrime1Ind, aDPrime2Ind, aReturnType);
    }


    /**
     * add an interface element
     */
    void
    add_interface_element(Interface_Element const & aInterfaceElement)
    {
        mInterfaceElements.push_back(aInterfaceElement);
    }

    moris::Cell<Interface_Element> &
    get_interface_elements()
    {
        return mInterfaceElements;
    }

    moris::Matrix<moris::IdMat>
    get_interface_element_ids()
    {
        moris::Matrix<moris::IdMat> tInterfaceElementIds(1,mInterfaceElements.size());

        for(moris::uint i = 0; i <mInterfaceElements.size(); i++)
        {
            tInterfaceElementIds(i) = mInterfaceElements(i).get_element_id();
        }

        return tInterfaceElementIds;
    }



    moris::Matrix<moris::IndexMat>
    get_extracted_interface_elements_loc_inds()
    {
        MORIS_ASSERT(mChildElementTopo == CellTopology::TET4,"Interface unzipping only tested on child meshes with tet4s");

        // hardcoded for tet4s
        moris::uint tNumNodesPerElement = 6;

        moris::uint tNumInterfaceElements = mInterfaceElements.size();

        moris::Matrix<moris::IndexMat> tInterfaceElemToNode(tNumInterfaceElements,tNumNodesPerElement);

        for(moris::uint iInt =0; iInt < tNumInterfaceElements; iInt++)
        {
            tInterfaceElemToNode.get_row(iInt) = mInterfaceElements(iInt).extract_as_standard_element_loc_inds().get_row(0);
        }


        return tInterfaceElemToNode;

    }

    /*
     * Set node indices in a child mesh
     * @param[in] aChildMeshIndex - Child mesh index
     * @param[in] aNodeInd - Node indices
     */
    void
    set_node_index(moris::size_t const &                aChildMeshIndex,
                   moris::Matrix< moris::IndexMat > & aNodeInd)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
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
    set_node_ids(moris::size_t const & aChildMeshIndex,
                 moris::Matrix< moris::IdMat > & aNodeIds)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
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
    add_node_ids(moris::size_t const & aChildMeshIndex,
                 moris::Matrix< moris::IdMat > & aNodeIds)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
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
    set_child_element_ids(moris::size_t const &   aChildMeshIndex,
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
    set_child_element_inds(moris::size_t const &      aChildMeshIndex,
                           moris::moris_index & aElementIndOffset )
    {
        mChildrenMeshes(aChildMeshIndex).set_child_element_inds(aElementIndOffset);
    }


    /*
     * Get element Ids in a child mesh
     */
    moris::Matrix< moris::IdMat > const &
    get_element_ids(moris::size_t const & aChildMeshIndex)
    {
        return  mChildrenMeshes(aChildMeshIndex).get_element_ids();
    }

    /*
     * Get element Ids in a child mesh
     */
    moris::Matrix< moris::IdMat >
    get_all_element_ids()
    {
        moris::size_t tNumElems = this->get_num_entities(EntityRank::ELEMENT);
        moris::Matrix< moris::IdMat > tElementIds(1,tNumElems);

        moris::size_t tCount = 0;
        for(moris::size_t iCM = 0; iCM<this->get_num_child_meshes(); iCM++)
        {
            moris::Matrix< moris::IdMat > const & tCMIds = this->get_element_ids(iCM);
            for(moris::size_t iE = 0; iE<tCMIds.numel(); iE++)
            {
                tElementIds(tCount) = tCMIds(iE);
                tCount++;
            }

        }
        return  tElementIds;
    }


    /*
     * Get element Ids in the cut mesh of a given id
     */
    Cell<moris::Matrix< moris::IdMat >>
    get_child_elements_by_phase(uint aNumPhases)
    {
        moris::size_t tNumElems = this->get_num_entities(EntityRank::ELEMENT);

        //Initialize output
        Cell<moris::Matrix<moris::IdMat>> tElementsByPhase(aNumPhases);
        moris::Matrix<moris::DDUMat> tPhaseCount(1,aNumPhases,0);
        for(uint i =0; i <aNumPhases; i++)
        {
            tElementsByPhase(i) = moris::Matrix<moris::IdMat>(1,tNumElems);
        }

        for(moris::size_t iCM = 0; iCM<this->get_num_child_meshes(); iCM++)
        {
            Child_Mesh const & tCM = get_child_mesh(iCM);
            moris::Matrix< moris::IdMat >    const & tCMIds     = tCM.get_element_ids();
            moris::Matrix< moris::IndexMat > const & tElemPhase = tCM.get_element_phase_indices();
            for(moris::size_t iE = 0; iE<tCMIds.numel(); iE++)
            {
                moris::moris_index tPhaseInd = tElemPhase(iE);
                uint tCount = tPhaseCount(tElemPhase(iE));
                tElementsByPhase(tPhaseInd)(tCount) = tCMIds(iE);
                tPhaseCount(tElemPhase(iE))++;
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
     * Get element Inds from a child mesh
     */
    moris::Matrix<moris::IndexMat> const &
    get_element_inds(moris::size_t const & aChildMeshIndex) const
    {
        return  mChildrenMeshes(aChildMeshIndex).get_element_inds();
    }


    void set_pending_node_index_pointers(moris::moris_index                   aChildMeshIndex,
                                         Cell<moris::moris_index*> &          aNodeIndPtr,
                                         moris::Matrix< moris::DDRMat > const & aParametricCoordinates)
    {
        MORIS_ASSERT(aChildMeshIndex < (moris::moris_index)mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        mChildrenMeshes(aChildMeshIndex).set_pending_node_index_pointers(aNodeIndPtr,aParametricCoordinates);
    }


   void retrieve_pending_node_inds()
    {
        for(moris::size_t i = 0; i < mNumberOfChildrenMesh; i++)
        {
            mChildrenMeshes(i).retrieve_pending_node_inds();
        }
        // Meshes have changed so counts are wrong
        mConsistentCounts = false;
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_indices(moris::size_t aChildMeshIndex)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_node_indices();
    }

    enum CellTopology
    get_child_element_topology()
    {
        return mChildElementTopo;
    }

    moris::size_t get_num_child_meshes() const
    {
        return mNumberOfChildrenMesh;
    }

    moris::size_t get_num_entities(moris::size_t aChildMeshIndex, enum EntityRank aEntityRank) const
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
        return mChildrenMeshes(aChildMeshIndex).get_num_entities(aEntityRank);
    }

    moris::size_t get_num_entities(enum EntityRank aEntityRank) const
    {
        // Make sure counts are up to date
        get_entity_counts();

        // Make sure the counts are up to date
        MORIS_ASSERT(mConsistentCounts, "Make sure to call get_entity_counts otherwise the mNumEntities variable is out dated and garbage is returned");

        return mNumEntities((moris::size_t) aEntityRank);
    }


    moris::moris_index get_parent_element_index(moris::size_t aChildMeshIndex)
    {
        MORIS_ASSERT(aChildMeshIndex < mNumberOfChildrenMesh, "The requested mesh index is out of bounds.");
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
    void get_child_elements_connected_to_parent_face(moris::moris_index const & aChildMeshIndex,
                                                     moris::moris_index const & aParentFaceIndex,
                                                     moris::Matrix< moris::IdMat > & aChildrenElementId,
                                                     moris::Matrix< moris::IndexMat > & aChildrenElementCMInd,
                                                     moris::Matrix< moris::IndexMat > & aFaceOrdinal) const
    {
        mChildrenMeshes(aChildMeshIndex).get_child_elements_connected_to_parent_face(aParentFaceIndex,aChildrenElementId,aChildrenElementCMInd,aFaceOrdinal);
    }

    // Outputting Mesh information


    void pack_cut_mesh_by_phase(moris::size_t const & aMeshIndex,
                                moris::size_t const & aNumPhases,
                                Cell<moris::Matrix< moris::DDSTMat >> & aElementIds,
                                Cell<moris::Matrix< moris::DDSTMat >> & aElementCMInds) const
    {

        mChildrenMeshes(aMeshIndex).pack_child_mesh_by_phase(aNumPhases,aElementCMInds,aElementIds);
    }

    moris::Matrix< moris::IdMat >
    pack_interface_sides() const
    {
        uint tNumElem = this->get_num_entities(EntityRank::ELEMENT);

        moris::Matrix< moris::IdMat > tFullElementIdAndSideOrdinals(tNumElem,2);

        uint tCount = 0;
        uint tNumElemsFromCM = 0;
        for(uint i = 0; i <this->get_num_child_meshes(); i++)
        {
           moris::Matrix< moris::IdMat > tSingleCMElementIdAndSideOrds = mChildrenMeshes(i).pack_interface_sides();

           tNumElemsFromCM = tSingleCMElementIdAndSideOrds.n_rows();
           if(tNumElemsFromCM!=0)
           {
               tFullElementIdAndSideOrdinals({tCount,tCount+tNumElemsFromCM-1},{0,1}) = tSingleCMElementIdAndSideOrds({0,tNumElemsFromCM-1},{0,1});
           }

           tCount = tCount + tNumElemsFromCM;
        }

        tFullElementIdAndSideOrdinals.resize(tCount,2);
        return tFullElementIdAndSideOrdinals;
    }

    /*
     * Get full element to node connectivity (ids). Full here means for all children meshes
     */
    moris::Matrix<moris::IdMat>
    get_full_element_to_node_glob_ids()
    {
        enum CellTopology tChildElementTopo = this->get_child_element_topology();
        moris::size_t     tNumElements = this->get_num_entities(EntityRank::ELEMENT);
        moris::size_t tNumNodesPerElem = 0;
        if(tChildElementTopo == CellTopology::TET4)
        {
            tNumNodesPerElem = 4;
        }
        else
        {
            MORIS_ERROR(0,"Not implemented");
        }

        moris::Matrix<moris::IdMat> tElementToNodeIds(tNumElements,tNumNodesPerElem);
        moris::size_t tCount = 0;
        for(moris::size_t i = 0; i<this->get_num_child_meshes(); i++)
        {
           moris::Matrix<moris::IdMat> tElementToNodeIdsCM = mChildrenMeshes(i).get_element_to_node_global();
           for(moris::size_t j = 0; j<tElementToNodeIdsCM.n_rows(); j++)
           {
               tElementToNodeIds.set_row(tCount, tElementToNodeIdsCM.get_row(j));
               tCount++;
           }

        }

        return tElementToNodeIds;
    }

    Child_Mesh const & get_child_mesh(moris::size_t const & aChildMeshIndex) const
        {
            return mChildrenMeshes(aChildMeshIndex);
        }

    Child_Mesh & get_child_mesh(moris::size_t const & aChildMeshIndex)
        {
            return mChildrenMeshes(aChildMeshIndex);
        }

private:
    // TODO: REMOVE THIS MEMBER VARIABLE
    moris::size_t mNumberOfChildrenMesh;

    Cell<Child_Mesh> mChildrenMeshes;

    // Interface elements
    moris::Cell<Interface_Element> mInterfaceElements;

    // Number of entities total in child meshes and if current count is accurate
    // mutable because some const function need this information, and if the counts
    // are not consistent we need to be able to update these vars
    mutable bool mConsistentCounts;
    mutable Cell<moris::size_t> mNumEntities;

    // topology of child elements (i.e. TET4)
    enum CellTopology mChildElementTopo;

    // Interface

private:

    void
    get_entity_counts() const
    {
        // If the counts aren't up to date then recount them.
        if(!mConsistentCounts)
        {
            enum EntityRank tRank = EntityRank::NODE;
            for(moris::size_t r = 0; r<(moris::size_t)EntityRank::ELEMENT+1; r++)
            {
                moris::size_t tCount = 0;
                // Cast r to enum
                tRank = static_cast<EntityRank>(r);
                for(moris::size_t m = 0; m<mNumberOfChildrenMesh; m++)
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
