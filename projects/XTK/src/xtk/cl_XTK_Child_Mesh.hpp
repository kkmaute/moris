/*
 * cl_XTK_Child_Mesh_Test.hpp
 *
 *  Created on: Jun 21, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_CHILD_MESH_HPP_
#define SRC_XTK_CL_XTK_CHILD_MESH_HPP_
#include <unordered_map>

#include "cl_Matrix.hpp"
#include "fn_print.hpp"
#include "fn_isvector.hpp"
#include "fn_iscol.hpp"
#include "fn_trans.hpp"

#include "containers/cl_XTK_Cell.hpp"

#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/fn_create_faces_from_element_to_node.hpp"
#include "xtk/fn_create_edges_from_element_to_node.hpp"
#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"

#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/fn_verify_tet_topology.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{
//TODO:CHANGE NAME TO CHILD_MESH from CHILD_MESH_TEST
class Child_Mesh_Test
{
public:

    Child_Mesh_Test():
            mElementTopology(EntityTopology::TET_4),
            mChildElementIds(0,0),
            mChildElementInds(0,0),
            mNodeIds(0,0),
            mNodeInds(0,0),
            mNodeParametricCoord(0,3),
            mHasFaceConn(false),
            mFaceToNode(0,0),
            mNodeToFace(0,0),
            mFaceToElement(0,0),
            mElementToFace(0,0),
            mFaceParentInds(0,0),
            mFaceParentRanks(0,0),
            mHasEdgeConn(false),
            mEdgeToNode(0,0),
            mNodeToEdge(0,0),
            mEdgeToElement(0,0),
            mElementToEdge(0,0),
            mEdgeParentInds(0,0),
            mEdgeParentRanks(0,0),
            mHasElemToElem(false),
            mElementToElement(0,0),
            mHasCoincidentEdges(false),
            mHasPhaseInfo(false),
            mElementPhaseIndices(0,0),
            mElementBinIndex(0,0),
            mBinBulkPhase(0),
            mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
    {};

    Child_Mesh_Test(moris::moris_index                 aParentElementIndex,
                    moris::Matrix< moris::IndexMat > & aNodeInds,
                    moris::Matrix< moris::IndexMat > & aElementToNode,
                    moris::Matrix< moris::IndexMat > & aElementEdgeParentInds,
                    moris::Matrix< moris::DDSTMat >  & aElementEdgeParentRanks,
                    moris::Matrix< moris::IndexMat > & aElementFaceParentInds,
                    moris::Matrix< moris::DDSTMat >  & aElementFaceParentRanks,
                    moris::Matrix< moris::DDSTMat >  & aElementInferfaceSides):
                        mElementTopology(EntityTopology::TET_4),
                        mChildElementIds(0,0),
                        mChildElementInds(0,0),
                        mNodeIds(0,0),
                        mNodeParametricCoord(0,3),
                        mHasFaceConn(false),
                        mFaceToNode(0,0),
                        mNodeToFace(0,0),
                        mFaceToElement(0,0),
                        mElementToFace(0,0),
                        mFaceParentInds(0,0),
                        mFaceParentRanks(0,0),
                        mHasEdgeConn(false),
                        mEdgeToNode(0,0),
                        mNodeToEdge(0,0),
                        mEdgeToElement(0,0),
                        mElementToEdge(0,0),
                        mEdgeParentInds(0,0),
                        mEdgeParentRanks(0,0),
                        mHasElemToElem(false),
                        mElementToElement(0,0),
                        mHasCoincidentEdges(false),
                        mHasPhaseInfo(false),
                        mElementPhaseIndices(0,0),
                        mElementBinIndex(0,0),
                        mBinBulkPhase(0),
                        mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
{
        // Check for row vector connectivity (if not it is transposed)
        row_vector_connectivity_check( aNodeInds );
        row_vector_connectivity_check( aElementEdgeParentInds );
        row_vector_connectivity_check( aElementEdgeParentRanks);
        row_vector_connectivity_check( aElementFaceParentInds );
        row_vector_connectivity_check( aElementFaceParentRanks );


        mParentElementIndex     = aParentElementIndex;
        mNumElem                = aElementToNode.n_rows();
        mNodeInds               = aNodeInds.copy();
        mElementToNode          = aElementToNode.copy();
        mElementEdgeParentInds  = aElementEdgeParentInds.copy();
        mElementEdgeParentRanks = aElementEdgeParentRanks.copy();
        mElementFaceParentInds  = aElementFaceParentInds.copy();
        mElementFaceParentRanks = aElementFaceParentRanks.copy();
        mElementInferfaceSides  = aElementInferfaceSides.copy();

        set_up_node_map();
}


    Child_Mesh_Test(Mesh_Modification_Template & aMeshModTemplate):
                        mElementTopology(EntityTopology::TET_4),
                        mChildElementIds(0,0),
                        mChildElementInds(0,0),
                        mNodeIds(0,0),
                        mNodeParametricCoord(0,3),
                        mHasFaceConn(false),
                        mFaceToNode(0,0),
                        mNodeToFace(0,0),
                        mFaceToElement(0,0),
                        mElementToFace(0,0),
                        mFaceParentInds(0,0),
                        mFaceParentRanks(0,0),
                        mHasEdgeConn(false),
                        mEdgeToNode(0,0),
                        mNodeToEdge(0,0),
                        mEdgeToElement(0,0),
                        mElementToEdge(0,0),
                        mEdgeParentInds(0,0),
                        mEdgeParentRanks(0,0),
                        mHasElemToElem(false),
                        mElementToElement(0,0),
                        mIntersectConnectivity(0,0),
                        mHasCoincidentEdges(false),
                        mHasPhaseInfo(false),
                        mElementPhaseIndices(0,0),
                        mElementBinIndex(0,0),
                        mBinBulkPhase(0),
                        mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
{

        reindex_template(aMeshModTemplate);

        //TODO: Try using Noahs move operator here.
        mParentElementIndex     = aMeshModTemplate.mParentElemInd;
        mNumElem                = aMeshModTemplate.mNumNewElem;
        mNodeInds               = aMeshModTemplate.mNodeInds.copy();
        mElementToNode          = aMeshModTemplate.mNewElementToNode.copy();
        mElementEdgeParentInds  = aMeshModTemplate.mNewParentEdgeOrdinals.copy();
        mElementEdgeParentRanks = aMeshModTemplate.mNewParentEdgeRanks.copy();
        mElementFaceParentInds  = aMeshModTemplate.mNewParentFaceOrdinals.copy();
        mElementFaceParentRanks = aMeshModTemplate.mNewParentFaceRanks.copy();
        mElementInferfaceSides  = aMeshModTemplate.mNewElementInterfaceSides.copy();

        set_up_node_map();

        generate_connectivities(true,true,true);

}

    ~Child_Mesh_Test(){}


    // Declare iterator type
    typename std::unordered_map<moris::size_t,moris::size_t>::iterator Map_Iterator;


    // --------------------------------------------------------------
    // Functions to access connectivity
    // --------------------------------------------------------------

    /*
     * Get number of a entity of a given rank
     */
    moris::size_t
    get_num_entities(enum EntityRank aEntityRank) const
    {
        moris::size_t tNumEntities = 0;
        if(aEntityRank == EntityRank::NODE)
        {
            tNumEntities = mNodeInds.n_cols();
        }
        else if(aEntityRank == EntityRank::EDGE)
        {
            XTK_ASSERT(mHasEdgeConn,"Without Edge connectivity it is unclear how many edges there are. Please call generate_connectivities w/ the edge flag on");
            tNumEntities = mEdgeToNode.n_rows();
        }

        else if(aEntityRank == EntityRank::FACE)
        {
            XTK_ASSERT(mHasFaceConn,"Without Edge connectivity it is unclear how many edges there are. Please call generate_connectivities w/ the edge flag on");
            tNumEntities = mFaceToNode.n_rows();
        }
        else if(aEntityRank == EntityRank::ELEMENT)
        {
            tNumEntities = mElementToNode.n_rows();
        }

        else
        {
            XTK_ASSERT(tNumEntities!=0,"Invalid entity rank specified");
        }

        return tNumEntities;
    }

    /*
     * return element to node connectivity matrix (processor local indices)
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_node() const
    {
        return mElementToNode;
    }

    /*
     * Converts the existing element to node connectivity (which contains proc local indices)
     * to child mesh local indices
     */
    moris::Matrix< moris::IndexMat >
    get_element_to_node_local() const
    {
        return convert_to_local_indices(mElementToNode);
    }


    /*
     * Converts the existing element to node connectivity (which contains proc local indices)
     * to child mesh local indices
     */
    moris::Matrix< moris::IdMat >
    get_element_to_node_global() const
    {
        moris::Matrix< moris::IdMat > tElementToNodeCMLoc = convert_to_local_indices(mElementToNode);
        return convert_to_glob_ids(tElementToNodeCMLoc);
    }

    /*
     * Return edge to node
     */
    moris::Matrix< moris::IndexMat > &
    get_edge_to_node()
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return mEdgeToNode;
    }



    /*
     * Return edge to node
     */
    moris::Matrix< moris::IndexMat >
    get_edge_to_node_local() const
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return convert_to_local_indices(mEdgeToNode);
    }

    /*
     * Return element to edge
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_edge() const
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return mElementToEdge;
    }

    /*
     * Return edges connected to elements
     */
    moris::Matrix< moris::IndexMat > const &
    get_edge_to_element() const
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return mEdgeToElement;
    }

    /*
     * return faces to node
     */
    moris::Matrix< moris::IndexMat > const &
    get_face_to_node() const
    {
        XTK_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return mFaceToNode;
    }

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    moris::Matrix< moris::IndexMat >
    get_face_to_node_local() const
    {
        XTK_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return convert_to_local_indices(mFaceToNode);
    }

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_face() const
    {
        XTK_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return mElementToFace;
    }

    /*
     * Return element to element connectivity (cm local element indices)
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_element() const
    {
        XTK_ASSERT(mHasElemToElem,"Element to element connectivity has not been generated with call to generate_element_to_element_connectivity");
        return mElementToElement;
    }

    /*
     * Return the parametric coordinates of the nodes (ordered by cm local index)
     */
    moris::Matrix< moris::DDRMat > const &
    get_parametric_coordinates() const
    {
        return mNodeParametricCoord;
    }

    /*
     * Return the parametric coordinates of a node using a node index
     */
    moris::Matrix< moris::DDRMat >
    get_parametric_coordinates(moris::moris_index aNodeIndex) const
    {

        // get child mesh local index
        auto tIter = mNodeIndsToCMInd.find(aNodeIndex);

        MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

        return mNodeParametricCoord.get_row(tIter->second);
    }



    // Functions to access ancestry

    moris::moris_index
    get_parent_element_index() const
    {
        return mParentElementIndex;
    }

    moris::Matrix< moris::IndexMat > const &
    get_face_parent_inds() const
    {
        return mFaceParentInds;
    }

    moris::Matrix< moris::DDSTMat > const &
    get_face_parent_ranks() const
    {
        return mFaceParentRanks;
    }



    moris::Matrix< moris::IndexMat > const &
    get_edge_parent_inds() const
    {
        return mEdgeParentInds;
    }

    moris::Matrix< moris::DDSTMat > const &
    get_edge_parent_ranks() const
    {
        return mEdgeParentRanks;
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_indices() const
    {
        return mNodeInds;
    }


    moris::Matrix< moris::IdMat > const &
    get_node_ids() const
    {
        return mNodeIds;
    }

    moris::Matrix< moris::IdMat > const &
    get_element_ids() const
    {
        return mChildElementIds;
    }


    moris::Matrix< moris::IndexMat > const &
    get_element_inds() const
    {
        return mChildElementInds;
    }


    // Function to access ordinals
    moris::Matrix< moris::IndexMat >
    get_edge_ordinal_from_element_and_edge_indices(moris::moris_index const & aElementIndex,
                                                   moris::Matrix< moris::IndexMat > const & aEdgeIndices) const
    {

        // get the elemnt to edge connectivity
        moris::Matrix< moris::IndexMat > const & tElemToEdge = get_element_to_edge();

        // Initialize bounds on loops
        moris::size_t tNumEdges = tElemToEdge.n_cols();
        moris::size_t tNumOrdinalstoFind = aEdgeIndices.n_cols();

        // Initialize output
        moris::Matrix< moris::IndexMat > tEdgeOrdinals(1,tNumOrdinalstoFind);

        // find the edge ordinals
        moris::size_t tCount = 0;
        for(moris::size_t iOrd = 0; iOrd<tNumOrdinalstoFind; iOrd++)
        {

            for(moris::size_t iEdge = 0; iEdge < tNumEdges; iEdge++)
            {
                if(aEdgeIndices(0,iOrd) == tElemToEdge(aElementIndex,iEdge))
                {
                    tEdgeOrdinals(0,tCount) = iEdge;
                    tCount++;
                    break;
                }
            }
        }

        XTK_ASSERT(tNumOrdinalstoFind == tCount,"All edge ordinals not found");

        return tEdgeOrdinals;

    }

    // Function to access ordinals
    moris::moris_index
    get_edge_ordinal_from_element_and_edge_indices(moris::moris_index const & aElementIndex,
                                                   moris::moris_index const & aEdgeIndex) const
    {

        // get the elemnt to edge connectivity
        moris::Matrix< moris::IndexMat > const & tElemToEdge = get_element_to_edge();

        // Initialize bounds on loops
        moris::size_t tNumEdges = tElemToEdge.n_cols();

        // Initialize output
        moris::moris_index tEdgeOrdinals = 1000;

        // find the edge ordinals
        moris::size_t tCount = 0;

        for(moris::size_t iEdge = 0; iEdge < tNumEdges; iEdge++)
        {
            if(aEdgeIndex == tElemToEdge(aElementIndex,iEdge))
            {
                tEdgeOrdinals = iEdge;
                tCount++;
                break;
            }
        }

        XTK_ASSERT(1 == tCount,"All edge ordinals not found");

        return tEdgeOrdinals;

    }


    // Function to access ordinals
    moris::moris_index
    get_face_ordinal_from_element_and_face_index(moris::moris_index const & aElementIndex,
                                                 moris::moris_index         aFaceIndex) const
    {
        // get the elemnt to edge connectivity
        moris::Matrix< moris::IndexMat > const & tElemToFace = get_element_to_face();

        moris::size_t tNumEdges = tElemToFace.n_cols();

        // Initialize output
        moris::moris_index tFaceOrdinal = 1000;

        bool tSuccess = false;

        for(moris::size_t iEdge = 0; iEdge < tNumEdges; iEdge++)
        {

            if(aFaceIndex == tElemToFace(aElementIndex,iEdge))
            {
                tFaceOrdinal = iEdge;
                tSuccess = true;
                continue;
            }
        }

        XTK_ASSERT(tSuccess
                   ,"Face ordinal not found");

        return tFaceOrdinal;

    }

    /*
     * Returns the child element and face ordinal connected to a provided parent face
     */
    void
    get_child_elements_connected_to_parent_face(moris::moris_index         const & aParentFaceIndex,
                                                moris::Matrix< moris::IdMat >    & aChildElemsIdsOnFace,
                                                moris::Matrix< moris::IndexMat > & aChildElemsCMIndOnFace,
                                                moris::Matrix< moris::IndexMat > & aChildElemOnFaceOrdinal) const
    {
        // First determine which of this meshes face live on the parent face
        moris::size_t tNumFaces = get_num_entities(EntityRank::FACE);
        moris::size_t tNumElemsOnFace = 0;
        moris::Matrix< moris::IndexMat > tLocFacesOnParentFace(1,tNumFaces);


        // Iterate through face ancestry ranks
        for(moris::size_t i = 0; i < tNumFaces; i++)
        {
            if(mFaceParentInds(0,i) == aParentFaceIndex && mFaceParentRanks(0,i) == 2)
            {
                tLocFacesOnParentFace(0,tNumElemsOnFace) = i;
                tNumElemsOnFace++;
            }
        }

        // Iterate through face to element connectivities of faces connected to the parent face
        aChildElemsIdsOnFace    = moris::Matrix< moris::IdMat >(1,tNumElemsOnFace*2);
        aChildElemsCMIndOnFace  = moris::Matrix< moris::IndexMat >(1,tNumElemsOnFace*2);
        aChildElemOnFaceOrdinal = moris::Matrix< moris::IndexMat >(1,tNumElemsOnFace*2);
        moris::size_t tCount = 0;

        for(moris::size_t i = 0; i<tNumElemsOnFace; i++)
        {
            moris::size_t tFaceInd = tLocFacesOnParentFace(0,i);
            for(moris::size_t j = 0; j<mFaceToElement.n_cols(); j++)
            {
                if(mFaceToElement(tFaceInd,j) != std::numeric_limits<moris::moris_index>::max())
                {
                    aChildElemsCMIndOnFace(0,tCount)  = mFaceToElement(tFaceInd,j);
                    aChildElemsIdsOnFace(0,tCount)    = mChildElementIds(0,mFaceToElement(tFaceInd,j));
                    aChildElemOnFaceOrdinal(0,tCount) = get_face_ordinal_from_element_and_face_index(mFaceToElement(tFaceInd,j),tFaceInd);
                    tCount++;
                }
            }
        }
        aChildElemsCMIndOnFace.resize(1,tCount);
        aChildElemsIdsOnFace.resize(1,tCount);
        aChildElemOnFaceOrdinal.resize(1,tCount);
    }

    // --------------------------------------------------------------
    // Functions to modify the mesh
    // --------------------------------------------------------------

    /*
     * Modify child mesh by selecting template using Intersection connectivity
     * to determine correct template and insert the template
     */
    void
    modify_child_mesh(enum TemplateType aTemplate)
    {

        modify_child_mesh_internal(aTemplate);

        // Verify the topology before modifying
        MORIS_ASSERT(verify_tet4_topology(this->get_element_to_node(),
                                          this->get_element_to_edge(),
                                          this->get_element_to_face(),
                                          this->get_edge_to_node(),
                                          this->get_face_to_node()),
                     "The generated mesh has an invalid topology");

    }

    void
    convert_tet4_to_tet10_child()
    {
        std::cerr<<"Warning conversion tet4 to tet10 not implemented in new child mesh structure"<<std::endl;
    }

    /*
     * Add more nodes to the existing node indices list
     */
    void
    add_node_indices(moris::Matrix< moris::IndexMat > const & aNewNodeInds)
    {
        // Allocate space
        moris::size_t tNumNewNodes = aNewNodeInds.n_cols();
        moris::size_t tNumCurrNodes = mNodeInds.n_cols();

        // add nodes to map
        add_nodes_to_map(aNewNodeInds);

        // resize the node index list
        mNodeInds.resize(1,tNumCurrNodes+tNumNewNodes);

        // add nodes to matrix
        for(moris::size_t i = 0; i<tNumNewNodes; i++)
        {
            mNodeInds(0,i+tNumCurrNodes) = aNewNodeInds(0,i);
        }

    }

    /*
     * Add more nodes to the existing node indices list
     */
    void
    add_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds)
    {
        // Allocate space
        moris::size_t tNumNewNodes = aNewNodeIds.n_cols();
        moris::size_t tNumCurrNodes = mNodeIds.n_cols();

        mNodeIds.resize(1,tNumCurrNodes+tNumNewNodes);

        // add nodes to matrix
        for(moris::size_t i = 0; i<tNumNewNodes; i++)
        {
            mNodeIds(0,i+tNumCurrNodes) = aNewNodeIds(0,i);
        }
    }


    /*
     * sets the node ids. note this overwrites existing mNodeIds data
     */

    void
    set_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds)
    {
        XTK_ASSERT(aNewNodeIds.n_cols() == get_num_entities(EntityRank::NODE),"Number of node ids does not match the number of nodes in this child mesh");
        mNodeIds = aNewNodeIds.copy();
    }

    /*
     * Sets the globally unique element Ids for the child mesh. This is important for mesh from data creation
     * @param[in] aElementId - First element Id (note: this is incremented as the id is assigned)
     */
    void set_child_element_ids(moris::moris_id & aElementId)
    {
        XTK_ASSERT(mChildElementIds.n_cols() == 0, "Element Ids already set");
        moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
        mChildElementIds = moris::Matrix< moris::IdMat >(1,tNumElements);

        for(moris::size_t iElem = 0; iElem<tNumElements; iElem++)
        {
            mChildElementIds(0,iElem) = aElementId;
            aElementId++;
        }
    }

    /*
     * Sets the processor unique element indices for the child mesh.
     * @param[in] aElementInd - First element Ind (note: this is incremented as the ind is assigned)
     */
    void set_child_element_inds(moris::moris_index & aElementInd)
    {
        XTK_ASSERT(mChildElementInds.n_cols() == 0, "Element Inds already set");


        moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
        mChildElementInds = moris::Matrix< moris::IndexMat >(1,tNumElements);

        for(moris::size_t iElem = 0; iElem<tNumElements; iElem++)
        {
            mChildElementInds(0,iElem) = aElementInd;
            aElementInd++;
        }
    }

    /*!
     * Add node parametric coordinate for a single node
     */
    void
    add_node_parametric_coordinate( moris::size_t aNodeIndex,
                                    moris::Matrix< moris::DDRMat > const & aParamCoord )
    {
        MORIS_ASSERT(aParamCoord(0) >= -1.0 && aParamCoord(0) <= -1.0, "Parametric coordinate is out of bounds - zeta");
        MORIS_ASSERT(aParamCoord(1) >= -1.0 && aParamCoord(1) <= -1.0, "Parametric coordinate is out of bounds - eta");
        MORIS_ASSERT(aParamCoord(2) >= -1.0 && aParamCoord(2) <= -1.0, "Parametric coordinate is out of bounds - xsi");

        // get child mesh local index
        auto tIter = mNodeIndsToCMInd.find(aNodeIndex);

        MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

        moris::size_t tNodeCMIndex = tIter->second;
        mNodeParametricCoord.set_row(tNodeCMIndex,aParamCoord);

    }

    /*!
     * Add node parametric coordinate for multiple nodes
     * @param[in] aNodeIndices - List of node indices with parametric coordinates
     * @param[in] aParamCoord  - Parametric coordinates for nodes (note: row 1 of aParamCoord is the Parametric coordinate of aNodeIndices(1))
     */
    void
    add_node_parametric_coordinate( moris::Matrix< moris::IndexMat> const & aNodeIndices,
                                    moris::Matrix< moris::DDRMat >  const & aParamCoord )
    {

        moris::size_t tNumNodes = aNodeIndices.numel();
        // Only do the checks if there are nodes to be added (this is only an issue for a child mesh
        // intersected by multipl geometries.
        if(tNumNodes !=0 )
        {
            MORIS_ASSERT(moris::isvector(aNodeIndices),"Node indices need to be a vector");
            MORIS_ASSERT(aParamCoord.n_rows() == tNumNodes ,"Number of nodes and parametric coordinates size");
            MORIS_ASSERT(aParamCoord.max() <= 1.0, "At least one of the parametric coordinates provided is out of bound");
            MORIS_ASSERT(aParamCoord.min() >= -1.0, "At least one of the parametric coordinates provided is out of bound");

            // Iterate over nodes and add the parametric coordinate
            for(moris::size_t i = 0; i <tNumNodes; i++)
            {
                // get child mesh local index
                auto tIter = mNodeIndsToCMInd.find(aNodeIndices(i));

                MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

                moris::size_t tNodeCMIndex = tIter->second;

                mNodeParametricCoord.set_row(tNodeCMIndex,aParamCoord.get_row(i));
            }
        }

    }

    void
    allocate_parametric_coordinates( moris::size_t aNumNewNodes )
    {
        moris::size_t tCurrentRow = mNodeParametricCoord.n_rows();
        moris::size_t tCurrentCol = mNodeParametricCoord.n_cols();
        mNodeParametricCoord.resize(tCurrentRow+aNumNewNodes,tCurrentCol);
    }

    /*
     * Resizes element related matrices
     */
    void
    allocate_more_elements(moris::size_t const & aNumMoreElements)
    {
        moris::size_t tCurrNumElem = mNumElem;
        moris::size_t tNewNumElem = tCurrNumElem + aNumMoreElements;

        mElementToNode.resize(tNewNumElem,mElementToNode.n_cols());
        mElementEdgeParentInds.resize(tNewNumElem,mElementEdgeParentInds.n_cols());
        mElementEdgeParentRanks.resize(tNewNumElem,mElementEdgeParentRanks.n_cols());
        mElementFaceParentInds.resize(tNewNumElem,mElementFaceParentInds.n_cols());
        mElementFaceParentRanks.resize(tNewNumElem,mElementFaceParentRanks.n_cols());
        mElementInferfaceSides.resize(tNewNumElem,mElementInferfaceSides.n_cols());
    }

    void
    insert_child_mesh_template(Mesh_Modification_Template & aMeshModTemplate)
    {
        reindex_template(aMeshModTemplate);


        // Replace the element which is marked as replaceable
        replace_element(aMeshModTemplate.mElemIndToReplace,
                        0,
                        aMeshModTemplate.mNewElementToNode,
                        aMeshModTemplate.mNewParentEdgeOrdinals,
                        aMeshModTemplate.mNewParentEdgeRanks,
                        aMeshModTemplate.mNewParentFaceOrdinals,
                        aMeshModTemplate.mNewParentFaceRanks,
                        aMeshModTemplate.mNewElementInterfaceSides);

        // Number of elements that are new
        moris::size_t tNumElemToAdd = aMeshModTemplate.mNumNewElem-aMeshModTemplate.mNumElemToReplace;
        for(moris::size_t i = 0; i<tNumElemToAdd; i++)
        {
            add_element(i+aMeshModTemplate.mNumElemToReplace,
                        aMeshModTemplate.mNewElementToNode,
                        aMeshModTemplate.mNewParentEdgeOrdinals,
                        aMeshModTemplate.mNewParentEdgeRanks,
                        aMeshModTemplate.mNewParentFaceOrdinals,
                        aMeshModTemplate.mNewParentFaceRanks,
                        aMeshModTemplate.mNewElementInterfaceSides);
        }


    }

    /*
     * Convert the element to node, parent information to  the local indexing scheme rather than ordinal
     */
    void reindex_template(Mesh_Modification_Template & aMeshModTemplate)
    {
        aMeshModTemplate.mNewElementToNode = reindex_matrix(aMeshModTemplate.mNewElementToNode, 0, aMeshModTemplate.mNodeInds);

        // Reindex template edge and face ordinals
        reindex_template_parent_information(aMeshModTemplate);
    }


    /*
     * Replace the information associated with an element
     */
    void
    replace_element(moris::size_t                     const & aElementIndexToReplace,
                    moris::size_t                     const & aRowIndex,
                    moris::Matrix< moris::IndexMat > const & aElementToNode,
                    moris::Matrix< moris::IndexMat > const & aElementEdgeParentInds,
                    moris::Matrix< moris::DDSTMat > const & aElementEdgeParentRanks,
                    moris::Matrix< moris::IndexMat > const & aElementFaceParentInds,
                    moris::Matrix< moris::DDSTMat > const & aElementFaceParentRanks,
                    moris::Matrix< moris::DDSTMat > const & aElementInterfaceFaces)
    {
        mHasFaceConn   = false;
        mHasEdgeConn   = false;
        mHasElemToElem = false;
        replace_row(aRowIndex, aElementToNode         , aElementIndexToReplace, mElementToNode);
        replace_row(aRowIndex, aElementEdgeParentInds , aElementIndexToReplace, mElementEdgeParentInds);
        replace_row(aRowIndex, aElementEdgeParentRanks, aElementIndexToReplace, mElementEdgeParentRanks);
        replace_row(aRowIndex, aElementFaceParentInds , aElementIndexToReplace, mElementFaceParentInds);
        replace_row(aRowIndex, aElementFaceParentRanks, aElementIndexToReplace, mElementFaceParentRanks);
        replace_row(aRowIndex, aElementInterfaceFaces , aElementIndexToReplace, mElementInferfaceSides);
    }


    /*
     * Replace the information associated with an element
     */
    void
    add_element(moris::size_t                     const & aRowIndex,
                moris::Matrix< moris::IndexMat > const & aElementToNode,
                moris::Matrix< moris::IndexMat > const & aElementEdgeParentInds,
                moris::Matrix< moris::DDSTMat > const & aElementEdgeParentRanks,
                moris::Matrix< moris::IndexMat > const & aElementFaceParentInds,
                moris::Matrix< moris::DDSTMat > const & aElementFaceParentRanks,
                moris::Matrix< moris::DDSTMat> const & aElementInterfaceFaces)
    {
        XTK_ASSERT(mNumElem<mElementToNode.n_rows(),"Not enough space allocated in call to allocate_more_elements");
        mHasFaceConn = false;
        mHasEdgeConn = false;
        mHasElemToElem = false;

        replace_row(aRowIndex, aElementToNode         , mNumElem, mElementToNode);
        replace_row(aRowIndex, aElementEdgeParentInds , mNumElem, mElementEdgeParentInds);
        replace_row(aRowIndex, aElementEdgeParentRanks, mNumElem, mElementEdgeParentRanks);
        replace_row(aRowIndex, aElementFaceParentInds , mNumElem, mElementFaceParentInds);
        replace_row(aRowIndex, aElementFaceParentRanks, mNumElem, mElementFaceParentRanks);
        replace_row(aRowIndex, aElementInterfaceFaces , mNumElem, mElementInferfaceSides);

        mNumElem++;
    }


    /*
     * Generates face connectivities, edge connectivities, element to element graph
     */

    void
    generate_connectivities(bool aGenerateFaceConn,
                            bool aGenerateEdgeConn,
                            bool aGenerateElemToElem)
    {
        moris::Matrix< moris::IndexMat > tElementToNodeLoc = get_element_to_node_local();

        if( aGenerateFaceConn )
        {
            generate_face_connectivity_and_ancestry(tElementToNodeLoc);
        }

        if( aGenerateEdgeConn )
        {
            generate_edge_connectivity_and_ancestry(tElementToNodeLoc);
        }

        if(aGenerateElemToElem)
        {
            generate_element_to_element_connectivity();
        }
    }

    /**
      * Tracks the intersection connectivty prior to adding templates to the mesh. Each intersected edge of an element tracks is flagged to have an
      * intersection nod
     */
    void init_intersect_connectivity()
    {
        moris::size_t tNumD = get_num_entities(EntityRank::ELEMENT);
        moris::size_t tNumEdgesToElem = get_element_to_edge().n_cols();
        mIntersectConnectivity   = moris::Matrix< moris::IndexMat >(tNumD, tNumEdgesToElem * 2 + 1, 0); // Needs to be zero for use column
        mEdgeOnInterface = moris::Matrix< moris::IndexMat >(this->get_num_entities(EntityRank::EDGE), 1, 0);
    }

    /**
      * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
      *       - 1 means the provided aDPrime1Ind is an XTK index
      *
      * aDPrime2Ind must be XTK local index
      */
     void
     add_entity_to_intersect_connectivity(moris::moris_index aCMNodeInd,
                                          moris::moris_index aCMEdgeInd,
                                          moris::size_t aFlag)
     {
         if(aFlag == 0)
         {
             aCMNodeInd = aCMNodeInd + get_num_entities(EntityRank::NODE);
         }

         // Get entities of dimension edge connected to elements
         moris::Matrix< moris::IndexMat > const tEdgeToElem = get_edge_to_element();


         moris::size_t tNumElemsConnected = 0;
         for(moris::size_t i = 0; i<tEdgeToElem.n_cols(); i++)
         {
             if(tEdgeToElem(aCMEdgeInd, i)== std::numeric_limits<moris::moris_index>::max())
             {
                 break;
             }
             tNumElemsConnected++;
         }

         for(moris::size_t i = 0; i < tNumElemsConnected; i++)
         {
             moris::size_t tElemCMInd = tEdgeToElem(aCMEdgeInd, i);

             XTK_ASSERT(mIntersectConnectivity(tElemCMInd, 0) < (moris::moris_index)this->get_element_to_edge().n_cols(),
                        "Entity corresponding to provided aDInd has exceeded allocated space");
             XTK_ASSERT(tElemCMInd < mIntersectConnectivity.n_rows(),
                        "aDInd is outside of bounds. Has auxiliary connectivity been initialized?");

             mIntersectConnectivity(tElemCMInd, mIntersectConnectivity(tElemCMInd, 0) + 1) = aCMNodeInd;
             mIntersectConnectivity(tElemCMInd, mIntersectConnectivity(tElemCMInd, 0) + mElementToEdge.n_cols() + 1) = aCMEdgeInd;
             mIntersectConnectivity(tElemCMInd, 0)++;
         }

     }

     void
     mark_edge_as_on_interface(moris::moris_index aEdgeIndex)
     {
         XTK_ASSERT(aEdgeIndex<(moris::moris_index)get_num_entities(EntityRank::EDGE),"Edge index out of bounds");
         mEdgeOnInterface(aEdgeIndex) = 1;
         mHasCoincidentEdges = true;
     }

    /*
     * Tracks intersection connectivity prior to the child mesh being modified
     */
      void
      set_intersect_connectivity(moris::Matrix< moris::IndexMat > const & aIntConnectivity)
      {
          // Copy so the data is completely owned by the child mesh
          mIntersectConnectivity = aIntConnectivity.copy();
      }

      /**
       * Tells the child mesh where a node index will be placed once it has been communicated
       */
      void set_pending_node_index_pointers(Cell<moris::moris_index*>            aNodeIndPtr,
                                           moris::Matrix< moris::DDRMat > const & aNodeParamCoordinate)
      {
          mPtrPendingNodeIndex     = aNodeIndPtr;
          mPendingParamCoordinates = aNodeParamCoordinate.copy();
      }

      /**
       *  retrieve_pending_node_inds
       * XTK mesh retrieves pending node indices. Prior to this call the unique node assignments
       * must be complete via the request structure in XTK model
       */
      void retrieve_pending_node_inds()
      {
          // Number of nodes before adding the new nodes
          moris::size_t tNumNodesPre = this->get_num_entities(EntityRank::NODE);

          // number of new nodes
          moris::size_t tNumNodesToAdd = mPtrPendingNodeIndex.size();

          // Collect node indices from their address
          moris::Matrix< moris::IndexMat > tNodeMat(1, mPtrPendingNodeIndex.size());
          for(moris::size_t i = 0; i < mPtrPendingNodeIndex.size(); i++)
          {
              tNodeMat(0, i) = *mPtrPendingNodeIndex(i);
          }

          // Add node indices to child mesh
          add_node_indices(tNodeMat);

          // Resize the parametric coordinate information
          mNodeParametricCoord.resize(tNumNodesPre + tNumNodesToAdd ,3);

          // Add parametric information
          add_node_parametric_coordinate(tNodeMat,mPendingParamCoordinates);

          // size out since the nodes have been retrieved
          mPtrPendingNodeIndex.resize(0, NULL);
          mPendingParamCoordinates.resize(0,0);
      }


      /*
       * Take the information of edges on the interface and
       * figure out which faces are on the interface, then
       * mark element edges as on interface.
       */
      void
      mark_interface_faces_from_interface_coincident_faces()
      {
          if(mHasCoincidentEdges)
          {
              uint tNumElements         = this->get_num_entities(EntityRank::ELEMENT);
              uint tNumEdgesPerElem     = mElementToEdge.n_cols();
              uint tNumFacesPerElem     = mElementToFace.n_cols();
              uint tNumEdgesPerFace     = 3;
              moris::moris_index tDummy = std::numeric_limits<moris::moris_index>::max();
              moris::Matrix<moris::IndexMat> tElementEdgeOrdsOnInterface(tNumElements,tNumEdgesPerElem);
              tElementEdgeOrdsOnInterface.fill(tDummy);
              moris::Matrix<moris::IndexMat> tElemCounter(tNumElements,1,0);

              // Get reference to edge to element connectivity
              moris::Matrix< moris::IndexMat > const & tEdgeToElement = get_edge_to_element();

              for(uint iEdge = 0; iEdge<mEdgeOnInterface.numel(); iEdge++)
              {
                  if(mEdgeOnInterface(iEdge) == 1)
                  {
                      for(uint i = 0; i <tEdgeToElement.n_cols(); i++)
                      {
                          if(tEdgeToElement(iEdge,i) ==  tDummy)
                          {
                              break;
                          }
                          moris::moris_index tElemInd = tEdgeToElement(iEdge,i);
                          moris::moris_index tEdgeOrd = this->get_edge_ordinal_from_element_and_edge_indices(tElemInd,iEdge);
                          tElementEdgeOrdsOnInterface(tElemInd,tElemCounter(tElemInd)) = tEdgeOrd;
                          tElemCounter(tElemInd)++;
                      }
                  }
              }

              // Now that we know all the element edges that are intersected
              // figure out which faces are on the interface
              uint tNumFacesPerEdge = 2;
              moris::Matrix<moris::IndexMat> const tEdgeOrdinalToFaceOrdMap = {{0,3},{1,3},{2,3},{0,2},{0,1},{1,2}};
              moris::Matrix<moris::DDUMat>       tFaceCounter(1,tNumFacesPerElem,0);

              for(uint iElem = 0; iElem < tNumElements; iElem++)
              {
                  for(uint iEdge = 0; iEdge<tNumEdgesPerElem; iEdge++)
                  {

                      moris::moris_index tEdgeOrd = tElementEdgeOrdsOnInterface(iElem,iEdge);
                      if(tEdgeOrd==tDummy)
                      {
                          break;
                      }

                      for(uint iEO = 0; iEO<tNumFacesPerEdge; iEO++)
                      {
                          moris::moris_index tFaceOrd = tEdgeOrdinalToFaceOrdMap(tEdgeOrd,iEO);
                          tFaceCounter(tFaceOrd)++;

                          if(tFaceCounter(tFaceOrd) == tNumEdgesPerFace)
                          {
                              mElementInferfaceSides(iElem) = tFaceOrd;
                          }
                      }



                  }
                  tFaceCounter.fill(0);
              }
          }

          mHasCoincidentEdges = false;
          mEdgeOnInterface.resize(0,0);
      }


      // --------------------------------------------------------------
      // Functions for elemental phase information and subphase bins
      // --------------------------------------------------------------

      void
      initialize_element_phase_mat()
      {
          mElementPhaseIndices = moris::Matrix< moris::IndexMat >(1,get_num_entities(EntityRank::ELEMENT),std::numeric_limits<moris::moris_index>::max());
      }

      void set_element_phase_index(moris::size_t aEntityIndex,
                                   moris::size_t aEntityPhaseIndex)
      {
          XTK_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
          mHasPhaseInfo = true;
          mElementPhaseIndices(0,aEntityIndex) = aEntityPhaseIndex;
      }

      moris::size_t get_element_phase_index( moris::size_t const & aEntityIndex) const
      {
          XTK_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
          XTK_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
          return mElementPhaseIndices(0,aEntityIndex);
      }

      moris::Matrix< moris::IndexMat > const &
      get_element_phase_indices() const
      {
          XTK_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
          return mElementPhaseIndices;
      }

      moris::size_t
      get_num_subphase_bins()
      {
          return mSubPhaseBins.size();
      }

      Cell<moris::moris_index> const &
      get_subphase_bin_bulk_phase()
      {
          return mBinBulkPhase;
      }



      /*
       * Sets the elemental subphase value in this child mesh
       */
      void
      set_elemental_subphase(moris::Matrix< moris::IndexMat > const & aElementSubPhase)
      {
          construct_subphase_bins(aElementSubPhase);
      }

      /*
       * Returns the elemental subphase values
       */
      moris::Matrix< moris::IndexMat > const &
      get_elemental_subphase_bin_membership()
      {
          return mElementBinIndex;
      }

      // --------------------------------------------------------------
      // Functions IO
      // --------------------------------------------------------------

      void pack_child_mesh_by_phase(moris::size_t const & aNumPhases,
                                    Cell<moris::Matrix< moris::DDSTMat >> & aElementIds,
                                    Cell<moris::Matrix< moris::DDSTMat >> & aElementCMInds) const
      {
          moris::size_t tNumElems = get_num_entities(EntityRank::ELEMENT);

          aElementIds = Cell<moris::Matrix< moris::DDSTMat >>(aNumPhases);
          aElementCMInds = Cell<moris::Matrix< moris::DDSTMat >>(aNumPhases);

          for(moris::size_t i = 0; i<aNumPhases; i++)
          {
              aElementIds(i) = moris::Matrix< moris::DDSTMat >(1,tNumElems);
              aElementCMInds(i) = moris::Matrix< moris::DDSTMat >(1,tNumElems);
          }

          Cell<moris::size_t> tPhaseCounter(aNumPhases,0);

          moris::Matrix< moris::IndexMat > const & tElementPhaseIndices = get_element_phase_indices();
          moris::Matrix< moris::IdMat > const & tElementIds  = get_element_ids();

          for(moris::size_t i = 0; i < tNumElems; i++)
          {
              moris::size_t tPhaseIndex = tElementPhaseIndices(0,i);
              moris::size_t tPhaseCount = tPhaseCounter(tPhaseIndex);
              aElementIds(tPhaseIndex)(0,tPhaseCount) = tElementIds(0,i);
              aElementCMInds(tPhaseIndex)(0,tPhaseCount) = i;
              tPhaseCounter(tPhaseIndex)++;
          }

          for(moris::size_t i = 0; i<aNumPhases; i++)
          {
              aElementIds(i).resize(1,tPhaseCounter(i));
              aElementCMInds(i).resize(1,tPhaseCounter(i));
          }
      }


      moris::Matrix< moris::IdMat >
      pack_interface_sides() const
      {
          // Loop bound and sizing
          moris::size_t tNumElem = get_num_entities(EntityRank::ELEMENT);

          moris::Matrix< moris::IdMat > tInterfaceSideSetInfo(tNumElem,2);

          // Keep track of the number of interface sides
          moris::size_t tCount = 0;

         // Iterate over each element and if the element has an interface side it will be in mElementInferfaceSides vector
         for(moris::size_t iEl =0 ; iEl<tNumElem; iEl++)
         {
             //TODO: NOTE THIS WILL NOT WORK WITH MULTI-MATERIAL YET
             if(mElementInferfaceSides(iEl) != std::numeric_limits<moris::size_t>::max())
             {
                 tInterfaceSideSetInfo(tCount,0) = mChildElementIds(iEl);
                 tInterfaceSideSetInfo(tCount,1) = mElementInferfaceSides(iEl);
                 tCount++;
             }
         }

         // Size out space
         tInterfaceSideSetInfo.resize(tCount,2);

         return tInterfaceSideSetInfo;
      }

private:
    // Parent element index
    moris::moris_index                     mParentElementIndex;

    // Element To Node and Ancestry Information (This is the only data that is set with templates.
    // all other is generated with an algorithm)
    // All node connectivity is indexed by proc local indexs
    enum EntityTopology              mElementTopology;
    moris::size_t                          mNumElem;
    moris::Matrix< moris::IndexMat > mElementToNode; // node indices correspond to the local child mesh index
    moris::Matrix< moris::IndexMat > mElementEdgeParentInds;
    moris::Matrix< moris::DDSTMat  > mElementEdgeParentRanks;
    moris::Matrix< moris::IndexMat > mElementFaceParentInds;
    moris::Matrix< moris::DDSTMat  > mElementFaceParentRanks;
    moris::Matrix< moris::DDSTMat  > mElementInferfaceSides;

    // Child element information ---------------------------
    moris::Matrix< moris::IdMat >    mChildElementIds;
    moris::Matrix< moris::IndexMat > mChildElementInds;

    // Node information ------------------------------------
    moris::Matrix< moris::IdMat >    mNodeIds;
    moris::Matrix< moris::IndexMat > mNodeInds;

    // Map where  key - proc local ind, val - local child mesh index
    std::unordered_map<moris::size_t, moris::size_t> mNodeIndsToCMInd;

    // Parametric coordinate relative to parent element
    moris::Matrix< moris::DDRMat >   mNodeParametricCoord;

    // Face Connectivity -----------------------------------
    bool mHasFaceConn;
    moris::Matrix< moris::IndexMat > mFaceToNode;
    moris::Matrix< moris::IndexMat > mNodeToFace;
    moris::Matrix< moris::IndexMat > mFaceToElement;
    moris::Matrix< moris::IndexMat > mElementToFace;
    moris::Matrix< moris::IndexMat > mFaceParentInds;
    moris::Matrix< moris::DDSTMat >  mFaceParentRanks;

    // Edge connectivity -----------------------------------
    bool mHasEdgeConn;
    moris::Matrix< moris::IndexMat > mEdgeToNode;
    moris::Matrix< moris::IndexMat > mNodeToEdge;
    moris::Matrix< moris::IndexMat > mEdgeToElement;
    moris::Matrix< moris::IndexMat > mElementToEdge;
    moris::Matrix< moris::IndexMat > mEdgeParentInds;
    moris::Matrix< moris::DDSTMat >  mEdgeParentRanks;

    // Element to Element graph ----------------------------
    bool mHasElemToElem;
    moris::Matrix< moris::IndexMat > mElementToElement;

    // Auxiliary connectivity data and pending nodes (mesh modification data)
    moris::Matrix< moris::IndexMat > mIntersectConnectivity;
    bool                             mHasCoincidentEdges;
    moris::Matrix< moris::IndexMat > mEdgeOnInterface;
    Cell<moris::moris_index*>        mPtrPendingNodeIndex;
    moris::Matrix< moris::DDRMat >   mPendingParamCoordinates;

    // Phase member variables (structured) -----------------
    bool                                   mHasPhaseInfo;
    moris::Matrix< moris::IndexMat >       mElementPhaseIndices;
    moris::Matrix< moris::IndexMat >       mElementBinIndex;
    Cell<moris::moris_index>               mBinBulkPhase;
    Cell<moris::Matrix< moris::IndexMat >> mSubPhaseBins;

private:
    void
    generate_face_connectivity_and_ancestry(moris::Matrix< moris::IndexMat > const & aElementToNodeLocal)
    {
        create_faces_from_element_to_node(mElementTopology,
                                          get_num_entities(EntityRank::NODE),
                                          aElementToNodeLocal,
                                          mElementToFace,
                                          mFaceToNode,
                                          mNodeToFace,
                                          mFaceToElement);

        // Convert face to node from cm indices to proc indices
        mFaceToNode = convert_to_proc_indices(mFaceToNode);

        mHasFaceConn = true;

        setup_face_ancestry();
    }

    void
    generate_edge_connectivity_and_ancestry(moris::Matrix< moris::IndexMat > const & aElementToNodeLocal)
    {
        create_edges_from_element_to_node(mElementTopology,
                                          get_num_entities(EntityRank::NODE),
                                          aElementToNodeLocal,
                                          mElementToEdge,
                                          mEdgeToNode,
                                          mNodeToEdge,
                                          mEdgeToElement);


        // convert back from local to proc indices
        mEdgeToNode = convert_to_proc_indices(mEdgeToNode);

        mHasEdgeConn = true;

        setup_edge_ancestry();

    }

    void
    generate_element_to_element_connectivity()
    {
        XTK_ASSERT(mHasFaceConn,"Face conn needs to be consistent prior to generating element to element");
        mHasElemToElem          = true;
        moris::size_t tNumFacePerElem = 4;
        mElementToElement      = generate_element_to_element(mFaceToElement,
                                                              mNumElem,
                                                              tNumFacePerElem,
                                                              std::numeric_limits<moris::moris_index>::max());
    }

    void
    set_up_node_map()
    {
        moris::size_t tNumNodes = get_num_entities(EntityRank::NODE);

        for( moris::size_t i = 0; i<tNumNodes; i++)
        {
            if(mNodeIndsToCMInd.find(mNodeInds(i)) == mNodeIndsToCMInd.end())
            {
                mNodeIndsToCMInd[mNodeInds(i)] = i;
            }

            else
            {
                moris::size_t tBreaker = 0;
                XTK_ASSERT(tBreaker!=0," Attempted to add duplicate node in constructor, are your node indices correct?");
            }
        }
    }


    /*
     * Note: this function only adds nodes which do not already exist to the map
     * and ignores nodes which are already in the map.
     */
    void
    add_nodes_to_map(moris::Matrix< moris::IndexMat > const & aNodesToAdd)
    {
        moris::size_t tNumNodes = get_num_entities(EntityRank::NODE);
        moris::size_t tNumNewNodes = aNodesToAdd.n_cols();

        for( moris::size_t i = 0; i<tNumNewNodes; i++)
        {

            if(mNodeIndsToCMInd.find(aNodesToAdd(0,i)) == mNodeIndsToCMInd.end())
            {
                mNodeIndsToCMInd[aNodesToAdd(0,i)] = i+tNumNodes;
            }
        }

    }

    /*
     * Takes the element face ordinal ancestry and creates a face ancestry
     */
    void
    setup_face_ancestry()
    {
        moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
        moris::size_t tNumFacePerElement = mElementToFace.n_cols();
        moris::size_t tNumFaces = get_num_entities(EntityRank::FACE);

        mFaceParentInds  = moris::Matrix< moris::IndexMat >(1,tNumFaces);
        mFaceParentRanks = moris::Matrix< moris::DDSTMat >(1,tNumFaces);

        for( moris::size_t i = 0; i<tNumElements; i++)
        {
            for(moris::size_t j = 0; j<tNumFacePerElement; j++)
            {
                moris::size_t tFaceIndex = mElementToFace(i,j);

                mFaceParentInds(0,tFaceIndex) = mElementFaceParentInds(i,j);
                mFaceParentRanks(0,tFaceIndex) = mElementFaceParentRanks(i,j);
            }
        }

    }

    void
    setup_edge_ancestry()
    {
        moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
        moris::size_t tNumEdgePerElement = mElementToEdge.n_cols();
        moris::size_t tNumFaces = get_num_entities(EntityRank::EDGE);

        mEdgeParentInds  = moris::Matrix< moris::IndexMat >(1,tNumFaces);
        mEdgeParentRanks = moris::Matrix< moris::DDSTMat >(1,tNumFaces);


        for( moris::size_t i = 0; i<tNumElements; i++)
        {
            for(moris::size_t j = 0; j<tNumEdgePerElement; j++)
            {
                moris::size_t tFaceIndex = mElementToEdge(i,j);

                mEdgeParentInds(0,tFaceIndex) = mElementEdgeParentInds(i,j);
                mEdgeParentRanks(0,tFaceIndex) = mElementEdgeParentRanks(i,j);
            }
        }

    }



    moris::Matrix< moris::IndexMat >
    convert_to_local_indices(moris::Matrix< moris::IndexMat > const & aEntityRankToNode) const
    {
        moris::size_t tNumRows = aEntityRankToNode.n_rows();
        moris::size_t tNumCols = aEntityRankToNode.n_cols();

        moris::Matrix< moris::IndexMat > tLocalEntityRankToNode(tNumRows,tNumCols);


        for(moris::size_t i = 0; i < tNumRows; i++)
        {
            for(moris::size_t j = 0; j < tNumCols; j++)
            {

                auto tIter = mNodeIndsToCMInd.find(aEntityRankToNode(i,j));

                XTK_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

                tLocalEntityRankToNode(i,j) = tIter->second;
            }
        }

        return tLocalEntityRankToNode;
    }


    moris::Matrix< moris::IndexMat >
    convert_to_proc_indices(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
    {
        moris::size_t tNumRows = aEntityRankToNodeLocal.n_rows();
        moris::size_t tNumCols = aEntityRankToNodeLocal.n_cols();

        moris::Matrix< moris::IndexMat > tProcEntityRankToNode(tNumRows,tNumCols);

        for(moris::size_t i = 0; i < tNumRows; i++)
        {
            for(moris::size_t j = 0; j < tNumCols; j++)
            {
                tProcEntityRankToNode(i,j) = mNodeInds(0,aEntityRankToNodeLocal(i,j));
            }
        }

        return tProcEntityRankToNode;
    }

    moris::Matrix< moris::IdMat >
    convert_to_glob_ids(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
    {
        moris::size_t tNumRows = aEntityRankToNodeLocal.n_rows();
        moris::size_t tNumCols = aEntityRankToNodeLocal.n_cols();

        moris::Matrix<  moris::IdMat  > tProcEntityRankToNode(tNumRows,tNumCols);

        for(moris::size_t i = 0; i < tNumRows; i++)
        {
            for(moris::size_t j = 0; j < tNumCols; j++)
            {
                tProcEntityRankToNode(i,j) = mNodeIds(0,aEntityRankToNodeLocal(i,j));
            }
        }

        return tProcEntityRankToNode;
    }

    void
    modify_child_mesh_internal(enum TemplateType aTemplate)
    {
        if(aTemplate == TemplateType::HIERARCHY_TET4)
        {
            // Verify the topology before modifying
            MORIS_ASSERT(verify_tet4_topology(this->get_element_to_node(),
                                              this->get_element_to_edge(),
                                              this->get_element_to_face(),
                                              this->get_edge_to_node(),
                                              this->get_face_to_node()),
                                              "The generated mesh has an invalid topology");
            XTK_ASSERT(mIntersectConnectivity.n_rows() == mNumElem, "There needs to be a row for each element in the child mesh in the intersect connectivity");
            XTK_ASSERT(aTemplate == TemplateType::HIERARCHY_TET4,"This function needs to be abstracted and has only been tested with HIERARCHY_TET4.");

            // Iterate over number of elements (only ones that existed at the beginning of the modification)
            moris::size_t tNumExistElem = mNumElem;

            // Container for all the templates to add
            Cell<Mesh_Modification_Template> tTemplatesToAdd(tNumExistElem);

            // Keep track of the number of intersected elements
            moris::size_t tNumIntersected = 0;
            moris::size_t tNumNewElem     = 0;

            moris::Matrix< moris::IndexMat > tEdgeToNodeCMLoc = get_edge_to_node_local();

            moris::Matrix< moris::IndexMat > tElemToNodeCMLoc = get_element_to_node_local();
            moris::Matrix< moris::IndexMat > tElemToEdgeCMLoc = get_element_to_edge();

            for(moris::size_t iE = 0; iE<tNumExistElem; iE++)
            {
                if(mIntersectConnectivity(iE,0) == 3 || mIntersectConnectivity(iE,0) == 4)
                {
                    moris::size_t tPermutationId = std::numeric_limits<moris::size_t>::max();

                    moris::Matrix< moris::IndexMat > tIntersectConnRow = mIntersectConnectivity.get_row(iE);

                    moris::Matrix< moris::IndexMat > tSortedNodes
                    = sort_nodes(aTemplate, tIntersectConnRow, tEdgeToNodeCMLoc, iE, tPermutationId);

                    // Select more specific tet4 hierarchy template
                    // TODO: Clean these if statements up
                    enum TemplateType tNarrowTemplate = TemplateType::INVALID_TEMPLATE_TYPE;
                    if(mIntersectConnectivity(iE,0) == 4)
                    {
                        if(tSortedNodes(0, 8) == 0)
                        {
                            tNarrowTemplate = TemplateType::HIERARCHY_TET4_4Na;
                        }
                        else if(tSortedNodes(0, 8) == 1)
                        {
                            tNarrowTemplate = TemplateType::HIERARCHY_TET4_4Nb;
                        }
                        else if(tSortedNodes(0, 8) == 2)
                        {
                            tNarrowTemplate = TemplateType::HIERARCHY_TET4_4Nc;
                        }

                        tSortedNodes.resize(1,8);
                    }
                    else if(mIntersectConnectivity(iE,0) == 3)
                    {
                        tSortedNodes.resize(1,7);
                        tNarrowTemplate = TemplateType::HIERARCHY_TET4_3N;
                    }

                    // Get parent element information
                    moris::Matrix< moris::IndexMat >  tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    moris::Matrix< moris::IndexMat > tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentFaceRanks = mElementFaceParentRanks.get_row(iE);



                    tSortedNodes = convert_to_proc_indices(tSortedNodes);


                    // Setup template with this information

                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template
                    (tElementsAncestry(0,0),
                     iE,
                     tSortedNodes,
                     tParentEdgeInds,
                     tParentEdgeRanks,
                     tParentFaceInds,
                     tParentFaceRanks,
                     tNarrowTemplate,
                     tPermutationId);

                    // Increment the count of number of intersected elements and number of new elements
                    tNumNewElem = tNumNewElem + tTemplatesToAdd(tNumIntersected).mNumNewElem - tTemplatesToAdd(tNumIntersected).mNumElemToReplace;
                    tNumIntersected++;

                }

                else if (mIntersectConnectivity(iE,0) == 1)
                {

                    // Figure out which nodes are which in the template

                    moris::moris_index tEdgeOrd = this->get_edge_ordinal_from_element_and_edge_indices(iE,mIntersectConnectivity(iE,7));
                    moris::Matrix< moris::IndexMat > tSortedNodes = {{mElementToNode(iE,0),
                                                                     mElementToNode(iE,1),
                                                                     mElementToNode(iE,2),
                                                                     mElementToNode(iE,3),
                                                                     mNodeInds(mIntersectConnectivity(iE,1))}};

                    // Get parent element information
                    moris::Matrix< moris::IndexMat >  tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    moris::Matrix< moris::IndexMat > tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentFaceRanks = mElementFaceParentRanks.get_row(iE);

                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template
                    (tElementsAncestry(0,0),
                     iE,
                     tSortedNodes,
                     tParentEdgeInds,
                     tParentEdgeRanks,
                     tParentFaceInds,
                     tParentFaceRanks,
                     TemplateType::BISECTED_TET4,
                     tEdgeOrd);

                    // Increment the count of number of intersected elements and number of new elements
                    tNumNewElem = tNumNewElem + tTemplatesToAdd(tNumIntersected).mNumNewElem - tTemplatesToAdd(tNumIntersected).mNumElemToReplace;
                    tNumIntersected++;
                    continue;
                }
                else if(mIntersectConnectivity(iE,0) == 2 )
                {
                    moris::moris_index tNodeH = std::numeric_limits<moris::moris_index>::max();
                    moris::moris_index tEdgeH = std::numeric_limits<moris::moris_index>::max();
                    moris::moris_index tNodeL = std::numeric_limits<moris::moris_index>::max();
                    moris::moris_index tEdgeL = std::numeric_limits<moris::moris_index>::max();
                    moris::Matrix< moris::IdMat > const & tNodeIds = this->get_node_ids();
                    if(tNodeIds(mIntersectConnectivity(iE,1)) > tNodeIds(mIntersectConnectivity(iE,2)))
                    {
                        tNodeH = mIntersectConnectivity(iE,1);
                        tEdgeH = mIntersectConnectivity(iE,7);
                        tNodeL = mIntersectConnectivity(iE,2);
                        tEdgeL = mIntersectConnectivity(iE,8);
                    }
                    else
                    {
                        tNodeL = mIntersectConnectivity(iE,1);
                        tEdgeL = mIntersectConnectivity(iE,7);
                        tNodeH = mIntersectConnectivity(iE,2);
                        tEdgeH = mIntersectConnectivity(iE,8);
                    }


                    moris::moris_index tEdgeOrdL = this->get_edge_ordinal_from_element_and_edge_indices(iE,tEdgeL);
                    moris::moris_index tEdgeOrdH = this->get_edge_ordinal_from_element_and_edge_indices(iE,tEdgeH);

                    moris::moris_index tPermutation = 10*tEdgeOrdL + tEdgeOrdH;

                    moris::Matrix< moris::IndexMat > tSortedNodes = {{mElementToNode(iE,0),
                                                                      mElementToNode(iE,1),
                                                                      mElementToNode(iE,2),
                                                                      mElementToNode(iE,3),
                                                                      mNodeInds(tNodeL),
                                                                      mNodeInds(tNodeH)}};
                    // Get parent element information
                    moris::Matrix< moris::IndexMat > tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    moris::Matrix< moris::IndexMat > tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentFaceRanks = mElementFaceParentRanks.get_row(iE);

                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template
                    (tElementsAncestry(0,0),
                     iE,
                     tSortedNodes,
                     tParentEdgeInds,
                     tParentEdgeRanks,
                     tParentFaceInds,
                     tParentFaceRanks,
                     TemplateType::HIERARCHY_TET4_2,
                     tPermutation);

                    // Increment the count of number of intersected elements and number of new elements
                    tNumNewElem = tNumNewElem + tTemplatesToAdd(tNumIntersected).mNumNewElem - tTemplatesToAdd(tNumIntersected).mNumElemToReplace;
                    tNumIntersected++;
                    continue;

                }

                else if (mIntersectConnectivity(iE,0) == 0)
                {
                    continue;
                }
                else
                {
                    XTK_ERROR << "Invalid connectivity for nodal hierarchy template, (should be 3 or 4 nodes)\n";
                }
            }

            // Allocate space for new elements
            allocate_more_elements(tNumNewElem);

            // Add templates to the connectivity
            for(moris::size_t i = 0; i<tNumIntersected; i++)
            {
                insert_child_mesh_template(tTemplatesToAdd(i));
            }

            // generate face, edge and element to element connectivity
            generate_connectivities(true,true,true);

            // Clear the intersection connectivity
            mIntersectConnectivity.resize(1,1);
        }

        else if(aTemplate == TemplateType::REGULAR_SUBDIVISION_HEX8)
        {

            Mesh_Modification_Template tRegSubTemplate(mParentElementIndex,
                                                                                                0,
                                                                                                mNodeInds,
                                                                                                mElementEdgeParentInds,
                                                                                                mElementEdgeParentRanks,
                                                                                                mElementFaceParentInds,
                                                                                                mElementFaceParentRanks,
                                                                                                TemplateType::REGULAR_SUBDIVISION_HEX8);

            mElementToNode.resize(1,4);
            mElementEdgeParentInds.resize(1,6);
            mElementEdgeParentRanks.resize(1,6);
            mElementFaceParentInds.resize(1,4);
            mElementFaceParentRanks.resize(1,4);


            // Allocate space for new elements
            allocate_more_elements(tRegSubTemplate.mNumNewElem - tRegSubTemplate.mNumElemToReplace);

            // Insert the template in the mesh
            insert_child_mesh_template(tRegSubTemplate);

            // generate face, edge and element to element connectivity
            generate_connectivities(true,true,true);

        }
    }

    /*
     * Edges needed to get permutation
     */
    moris::Matrix< moris::IndexMat >
    sort_nodes(enum TemplateType                    aTemplate,
               moris::Matrix< moris::IndexMat > const & aIntConnectivity,
               moris::Matrix< moris::IndexMat > const & aEdgeToNodeCMLoc,
               moris::size_t const &                      aElementIndex,
               moris::size_t &                            aPermutation)
    {
        //Locate highest node in intersection connectivity
        switch(aTemplate)
        {
        case (TemplateType::HIERARCHY_TET4):
                         {

            if(aIntConnectivity(0, 0) == 3)
            {
                moris::moris_index tIntersectionCase = std::numeric_limits<moris::moris_index>::max();
                moris::Matrix< moris::IndexMat > tHigh({{ aIntConnectivity(0, 1), aIntConnectivity(0, 7)}});
                moris::Matrix< moris::IndexMat > tMid({{ aIntConnectivity(0, 1), aIntConnectivity(0, 7)}});
                moris::Matrix< moris::IndexMat > tLow({{ aIntConnectivity(0, 1), aIntConnectivity(0, 7)}});

                for(moris::size_t i = 0; i < 2; i++)
                {
                    if(mNodeIds(0,aIntConnectivity(0, i + 2)) > mNodeIds(0,tHigh(0, 0)))
                    {
                        tMid.set_row(0, tHigh);
                        tHigh(0, 0) = aIntConnectivity(0, i + 2);
                        tHigh(0, 1) = aIntConnectivity(0, i + 8);
                    }
                    else if(mNodeIds(0,aIntConnectivity(0, i + 2)) < mNodeIds(0,tLow(0, 0)))
                    {

                        tMid.set_row(0, tLow);
                        tLow(0, 0) = aIntConnectivity(0, i + 2);
                        tLow(0, 1) = aIntConnectivity(0, i + 8);
                    }
                    else
                    {
                        tMid(0, 0) = aIntConnectivity(0, i + 2);
                        tMid(0, 1) = aIntConnectivity(0, i + 8);
                    }
                }

                moris::Matrix< moris::IndexMat > tNodes13 = aEdgeToNodeCMLoc.get_row(tMid(0, 1));
                moris::Matrix< moris::IndexMat > tNodes12 = aEdgeToNodeCMLoc.get_row(tLow(0, 1));
                moris::Matrix< moris::IndexMat > tNodes14 = aEdgeToNodeCMLoc.get_row(tHigh(0, 1));

                // Initialize Edge matrix for the 3 node intersection
                moris::Matrix< moris::IndexMat > tEdgeIndices({{tLow(0, 1),tMid(0, 1),tHigh(0, 1)}});

                // Get edge ordinals of the edge indices
                moris::Matrix< moris::IndexMat > tEdgeOrdinals = get_edge_ordinal_from_element_and_edge_indices(aElementIndex,tEdgeIndices);
                get_intersection_permutation(tEdgeOrdinals,aPermutation);

                // Find the shared node in all the above lists
                // Decide which node is which (using intersections)
                moris::moris_index tN1 = std::numeric_limits<moris::moris_index>::max();
                moris::moris_index tN2 = std::numeric_limits<moris::moris_index>::max();
                moris::moris_index tN3 = std::numeric_limits<moris::moris_index>::max();
                moris::moris_index tN4 = std::numeric_limits<moris::moris_index>::max();

                // Find nodes 1,4,3
                if(tNodes14(0, 0) == tNodes13(0, 0))
                {
                    tN1 = tNodes14(0, 0);
                    tN4 = tNodes14(0, 1);
                    tN3 = tNodes13(0, 1);
                }
                else if(tNodes14(0, 0) == tNodes13(0, 1))
                {
                    tN1 = tNodes14(0, 0);
                    tN4 = tNodes14(0, 1);
                    tN3 = tNodes13(0, 0);
                }
                else if(tNodes14(0, 1) == tNodes13(0, 0))
                {
                    tN1 = tNodes14(0, 1);
                    tN4 = tNodes14(0, 0);
                    tN3 = tNodes13(0, 1);
                }
                else if(tNodes14(0, 1) == tNodes13(0, 1))
                {
                    tN1 = tNodes14(0, 1);
                    tN4 = tNodes14(0, 0);
                    tN3 = tNodes13(0, 0);
                }
                else
                    XTK_ERROR << "Duplicate node not found, invalid edge intersection configuration";

                // Find node 2
                if(tN1 == tNodes12(0, 0))
                {
                    tN2 = tNodes12(0, 1);
                }
                else if(tN1 == tNodes12(0, 1))
                {
                    tN2 = tNodes12(0, 0);
                }
                else
                    XTK_ERROR << "Node 2 not found, invalid edge intersection configuration";

                moris::Matrix< moris::IndexMat > tSortedNodes(
                        {{tN1, tN2, tN3, tN4, tLow(0, 0), tMid(0, 0), tHigh(0, 0), tIntersectionCase}});

                return tSortedNodes;
                         }

            else if(aIntConnectivity(0, 0) == 4)
            {
                // Sort Auxiliary nodes from highest to lowest using bubble sort (but based on first column then swap row

                moris::Matrix< moris::IndexMat > tSortedNodes(1, 9);
                moris::Matrix< moris::IndexMat > tNodes(
                {
                    {aIntConnectivity(0, 1), aIntConnectivity(0, 7)},
                    {aIntConnectivity(0, 2), aIntConnectivity(0, 8)},
                    {aIntConnectivity(0, 3), aIntConnectivity(0, 9)},
                    {aIntConnectivity(0, 4), aIntConnectivity(0, 10)}});

                moris::size_t j = 0;
                moris::size_t n = 4; // 4 numbers to sort
                bool swapped = true;

                // Temporary row storage
                moris::Matrix< moris::IndexMat > tRowStorage(1, 2);
                moris::Matrix< moris::IndexMat > tRowSwapper(1, 2);

                while(swapped)
                {
                    swapped = false;
                    j++;
                    for(moris::size_t i = 0; i < n - j; i++)
                    {
                        if(mNodeIds(0,tNodes(i, 0)) > mNodeIds(0,tNodes(i + 1, 0)))
                        {
                            tRowStorage = tNodes.get_row(i);
                            tRowSwapper = tNodes.get_row(i + 1);
                            tNodes.set_row(i, tRowSwapper);
                            tNodes.set_row(i + 1, tRowStorage);
                            swapped = true;
                        }
                    }
                }
                // Determine the relationship between high and low

                moris::Matrix< moris::IndexMat > tEdgeToNode = get_edge_to_node_local();

                moris::Matrix< moris::IndexMat > tNodesL  = tEdgeToNode.get_row(tNodes(0, 1));
                moris::Matrix< moris::IndexMat > tNodesML = tEdgeToNode.get_row(tNodes(1, 1));
                moris::Matrix< moris::IndexMat > tNodesMH = tEdgeToNode.get_row(tNodes(2, 1));
                moris::Matrix< moris::IndexMat > tNodesH  = tEdgeToNode.get_row(tNodes(3, 1));

                moris::size_t tHLOppFlag  = 1;
                moris::size_t tHMHOppFlag = 1;
                moris::size_t tHMLOppFlag = 1;

                // Initialize Edge matrix for the 3 node intersection
                moris::Matrix< moris::IndexMat > tEdgeIndices({{tNodes(0, 1),tNodes(1, 1),tNodes(2, 1),tNodes(3, 1)}});

                moris::Matrix< moris::IndexMat > tEdgeOrdinals = get_edge_ordinal_from_element_and_edge_indices(aElementIndex,tEdgeIndices);

                get_intersection_permutation(tEdgeOrdinals,aPermutation);


                // If L and H share a node then it is not case a
                for(moris::size_t i = 0; i < 2; i++)
                {
                    if(tNodesL(0, i) == tNodesH(0, 0))
                    {
                        tHLOppFlag = 0;
                    }

                    else if(tNodesL(0, i) == tNodesH(0, 1))
                    {
                        tHLOppFlag = 0;
                    }
                }

                // If MH and H share a node then its not case b
                for(moris::size_t i = 0; i < 2; i++)
                {
                    if(tNodesMH(0, i) == tNodesH(0, 0))
                    {
                        tHMHOppFlag = 0;
                    }

                    else if(tNodesMH(0, i) == tNodesH(0, 1))
                    {
                        tHMHOppFlag = 0;
                    }
                }

                // If ML and H share a node then its not case c
                for(moris::size_t i = 0; i < 2; i++)
                {
                    if(tNodesML(0, i) == tNodesH(0, 0))
                    {
                        tHMLOppFlag = 0;
                    }

                    else if(tNodesML(0, i) == tNodesH(0, 1))
                    {
                        tHMLOppFlag = 0;
                    }
                }
                if(tHLOppFlag)
                {
                    tSortedNodes(0, 8) = 0; // Indicating that you need to use template node_hier_4_node_a.inc
                    tSortedNodes(0, 4) = tNodes(0, 0);
                    tSortedNodes(0, 5) = tNodes(1, 0);
                    tSortedNodes(0, 6) = tNodes(2, 0);
                    tSortedNodes(0, 7) = tNodes(3, 0);

                    // Get node shared by MH and H
                    j = 1;
                    moris::size_t tSuccess = 0;
                    for(moris::size_t i = 0; i < 2; i++)
                    {
                        if(tNodesMH(0, i) == tNodesH(0, 0))
                        {
                            tSortedNodes(0, 1) = tNodesH(0, 1);  // High nodes independent node
                            tSortedNodes(0, 0) = tNodesMH(0, i);// Shared Node
                            tSortedNodes(0, 2) = tNodesMH(0, j);// Mid highs ind node
                            tSuccess = 1;
                        }

                        else if(tNodesMH(0, i) == tNodesH(0, 1))
                        {
                            tSortedNodes(0, 1) = tNodesH(0, 0);  // High nodes independent node
                            tSortedNodes(0, 0) = tNodesMH(0, i);// Shared Node
                            tSortedNodes(0, 2) = tNodesMH(0, j);// Mid highs ind node
                            tSuccess = 1;
                        }

                        else
                        {
                            j = i;
                        }
                    }

                    XTK_ASSERT(tSuccess = 1, "Sorting to find nodes 1,2,3 unsuccessful");

                    // Use Midlow and Low to find last node
                    j = 1;
                    tSuccess = 0;
                    for(moris::size_t i = 0; i < 2; i++)
                    {
                        if(tNodesML(0, i) == tNodesL(0, 0))
                        {
                            tSortedNodes(0, 3) = tNodesML(0, i);
                            tSuccess = 1;
                        }

                        else if(tNodesML(0, i) == tNodesL(0, 1))
                        {
                            tSortedNodes(0, 3) = tNodesML(0, i);
                            tSuccess = 1;
                        }

                        else
                        {
                            j = i;
                        }
                    }

                    XTK_ASSERT(tSuccess = 1, "Sorting to find node 4 unsuccessful");

                }
                else if(tHMHOppFlag)
                {
                    tSortedNodes(0, 8) = 1; // Indicating that you need to use template node_hier_4_node_b.inc
                    tSortedNodes(0, 4) = tNodes(0, 0);// Low
                    tSortedNodes(0, 5) = tNodes(1, 0);// Mid low
                    tSortedNodes(0, 6) = tNodes(2, 0);// Mid high
                    tSortedNodes(0, 7) = tNodes(3, 0);// High

                    // Get node shared by H and L
                    j = 1;
                    moris::size_t tSuccess = 0;
                    for(moris::size_t i = 0; i < 2; i++)
                    {
                        if(tNodesL(0, i) == tNodesH(0, 0))
                        {
                            tSortedNodes(0, 0) = tNodesL(0, i); // Shared Node
                            tSortedNodes(0, 1) = tNodesH(0, 1);// High nodes independent node
                            tSortedNodes(0, 2) = tNodesL(0, j);// Mid highs ind node
                            tSuccess = 1;
                        }

                        else if(tNodesL(0, i) == tNodesH(0, 1))
                        {
                            tSortedNodes(0, 0) = tNodesL(0, i); // Shared Node
                            tSortedNodes(0, 1) = tNodesH(0, 0);// High nodes independent node
                            tSortedNodes(0, 2) = tNodesL(0, j);// Mid highs ind node
                            tSuccess = 1;
                        }

                        else
                        {
                            j = i;
                        }
                    }
                    XTK_ASSERT(tSuccess = 1, "Sorting to find nodes 1,2,3 unsuccessful");

                    // Midlow and MidHigh to find last node
                    j = 1;
                    tSuccess = 0;
                    for(moris::size_t i = 0; i < 2; i++)
                    {
                        if(tNodesML(0, i) == tNodesMH(0, 0))
                        {
                            tSortedNodes(0, 3) = tNodesML(0, i);
                            tSuccess = 1;
                        }

                        else if(tNodesML(0, i) == tNodesMH(0, 1))
                        {
                            tSortedNodes(0, 3) = tNodesML(0, i);
                            tSuccess = 1;
                        }

                        else
                            j = i;
                    }
                    XTK_ASSERT(tSuccess = 1, "Sorting to find node 4 unsuccessful");

                }
                else if(tHMLOppFlag)
                {
                    tSortedNodes(0, 8) = 2; // Indicating that you need to use template node_hier_4_node_b.inc
                    tSortedNodes(0, 4) = tNodes(0, 0);// Low
                    tSortedNodes(0, 5) = tNodes(1, 0);// Mid low
                    tSortedNodes(0, 6) = tNodes(2, 0);// Mid high
                    tSortedNodes(0, 7) = tNodes(3, 0);// High

                    // Get node shared by H and L
                    j = 1;
                    moris::size_t tSuccess = 0;
                    for(moris::size_t i = 0; i < 2; i++)
                    {
                        if(tNodesL(0, i) == tNodesH(0, 0))
                        {
                            tSortedNodes(0, 0) = tNodesL(0, i); // Shared Node
                            tSortedNodes(0, 1) = tNodesH(0, 1);// High nodes independent node
                            tSortedNodes(0, 2) = tNodesL(0, j);// Mid highs ind node
                            tSuccess = 1;
                        }

                        else if(tNodesL(0, i) == tNodesH(0, 1))
                        {
                            tSortedNodes(0, 0) = tNodesL(0, i); // Shared Node
                            tSortedNodes(0, 1) = tNodesH(0, 0);// High nodes independent node
                            tSortedNodes(0, 2) = tNodesL(0, j);// Mid highs ind node
                            tSuccess = 1;
                        }

                        else
                        {
                            j = i;
                        }
                    }
                    XTK_ASSERT(tSuccess = 1, "Sorting to find nodes 1,2,3 unsuccessful");

                    // Midlow and MidHigh to find last node
                    j = 1;
                    tSuccess = 0;
                    for(moris::size_t i = 0; i < 2; i++)
                    {
                        if(tNodesML(0, i) == tNodesMH(0, 0))
                        {
                            tSortedNodes(0, 3) = tNodesML(0, i);
                            tSuccess = 1;
                        }

                        else if(tNodesML(0, i) == tNodesMH(0, 1))
                        {
                            tSortedNodes(0, 3) = tNodesML(0, i);
                            tSuccess = 1;
                        }

                        else
                        {
                            j = i;
                        }
                    }
                    XTK_ASSERT(tSuccess = 1, "Sorting to find node 4 unsuccessful");
                }
                else
                    XTK_ERROR << "Sorting Failed (invalid flagging). Did a node appear twice?";
                return tSortedNodes;
            }

            else
            {
                XTK_ERROR << "SORTING NOT COMPLETED! Check to see if this function is called for a non-intersected element";
                moris::Matrix< moris::IndexMat > dummy(1, 1);
                return dummy;
            }

            break;
                         }

        default:
        {
            XTK_ERROR << "Sorting for specified template type not implemented";

            moris::Matrix< moris::IndexMat > dummy(1, 1);
            return dummy;

            break;
        }

        }

    }


    /*
    * This assumes the edge ordinals provided are ordered in the same order as is found in the auxiliary connectivity
    */
    void
    get_intersection_permutation(moris::Matrix< moris::IndexMat > const & aOrderedEdgeOrdinals,
                                 moris::size_t & aPermutation)
    {
        if(aOrderedEdgeOrdinals.n_cols() == 3)
        {
            // Determine permutation
            // Rule:  1   * edge index containing the lowest node ID
            //      + 10  * edge index containing the middle node ID
            //      + 100 * edge index containing the highest node ID

            /*
             * Determine the permutation using the provided edge ordinals
             */
            aPermutation = aOrderedEdgeOrdinals(0, 0) + 10 * aOrderedEdgeOrdinals(0, 1)  + 100 *aOrderedEdgeOrdinals(0, 2) ;
        }

        else if (aOrderedEdgeOrdinals.n_cols() == 4)
        {
            // Determine permutation
            // Rule:  1    * edge ordinal containing the lowest node ID
            //      + 10   * edge ordinal containing the middle lowest node ID
            //      + 100  * edge ordinal containing the middle highest node ID
            //      + 1000 * edge ordinal containing the highest node ID

            /*
             * Determine the permutation using the provided edge ordinals
             */
            aPermutation = aOrderedEdgeOrdinals(0, 0) + 10 * aOrderedEdgeOrdinals(0, 1)  + 100 * aOrderedEdgeOrdinals(0, 2) + 1000 * aOrderedEdgeOrdinals(0, 3);
        }
        else
        {
            XTK_ERROR<<"Permutation rule not implemented";
        }


    }


    void
    reindex_template_parent_information(Mesh_Modification_Template & aMeshModTemplate)
    {
        // Reindex parent ordinals for edges
        moris::size_t tNumEdgePerElem = aMeshModTemplate.mNewParentEdgeOrdinals.n_cols();
        moris::size_t tNumFacePerElem = aMeshModTemplate.mNewParentFaceOrdinals.n_cols();
        moris::size_t tNumElems       = aMeshModTemplate.mNumNewElem;
        moris::Matrix< moris::IndexMat > tParentElem({{aMeshModTemplate.mParentElemInd}});
        moris::Matrix< moris::DDSTMat > tParentElemRank({{3}});

        // Place ptrs in cell to avoid if statement checking rank in for loop
        Cell<moris::Matrix< moris::IndexMat >*> tParentEntitiesInds({&aMeshModTemplate.mParentEdgeInds,
                                                            &aMeshModTemplate.mParentFaceInds,
                                                            &tParentElem});

        Cell<moris::Matrix< moris::DDSTMat >*> tParentEntitiesRanks({& aMeshModTemplate.mParentEdgeRanks,
                                                                 & aMeshModTemplate.mParentFaceRanks,
                                                                 & tParentElemRank});


        for(moris::size_t i = 0; i<tNumElems; i++)
        {
            // Reindex edges
            for(moris::size_t j = 0; j<tNumEdgePerElem; j++)
            {
                moris::size_t tEdgeParentRank = aMeshModTemplate.mNewParentEdgeRanks(i,j);
                aMeshModTemplate.mNewParentEdgeRanks(i,j) = (*tParentEntitiesRanks(tEdgeParentRank-1))(0,aMeshModTemplate.mNewParentEdgeOrdinals(i,j));
                aMeshModTemplate.mNewParentEdgeOrdinals(i,j) = (*tParentEntitiesInds(tEdgeParentRank-1))(0,aMeshModTemplate.mNewParentEdgeOrdinals(i,j));

            }

            for(moris::size_t j = 0; j<tNumFacePerElem; j++)
            {
                moris::size_t tFaceParentRank = aMeshModTemplate.mNewParentFaceRanks(i,j);
                aMeshModTemplate.mNewParentFaceRanks(i,j) = (*tParentEntitiesRanks(tFaceParentRank-1))(0,aMeshModTemplate.mNewParentFaceOrdinals(i,j));
                aMeshModTemplate.mNewParentFaceOrdinals(i,j) = (*tParentEntitiesInds(tFaceParentRank-1))(0,aMeshModTemplate.mNewParentFaceOrdinals(i,j));
            }
        }
    }

    void
    construct_subphase_bins(moris::Matrix< moris::IndexMat > const & aElementSubPhase)
    {
        // Number of bins corresponds to the maxmimum value in the element sub phase vector
        moris::size_t tNumBins = aElementSubPhase.max() + 1;
        moris::size_t tNumElements = this->get_num_entities(EntityRank::ELEMENT);

        // Initialize member variables
        // Element Sub-phase bins
        mBinBulkPhase = Cell<moris::moris_index>(tNumBins);
        mElementBinIndex = aElementSubPhase.copy();

        moris::Matrix< moris::DDSTMat > tBinSizeCounter(1,tNumBins,0);
        mSubPhaseBins = Cell<moris::Matrix< moris::IndexMat >>(tNumBins);
        for(moris::size_t i = 0; i<tNumBins; i++)
        {
            mSubPhaseBins(i) = moris::Matrix< moris::IndexMat >(1,tNumElements);
        }

        // Place the elements into bins;
        for(moris::size_t i = 0; i<tNumElements; i++)
        {
            moris::size_t tBinIndex = aElementSubPhase(0,i);
            moris::size_t tBinCount = tBinSizeCounter(0,tBinIndex);

            tBinSizeCounter(0,tBinIndex)++;

            mSubPhaseBins(tBinIndex)(0,tBinCount) = i;
        }

        // Size out extra space and set bin element bulk phase
        for(moris::size_t i = 0; i<tNumBins; i++)
        {
            moris::size_t tBinCount = tBinSizeCounter(0,i);
            mBinBulkPhase(i) = get_element_phase_index(mSubPhaseBins(i)(0,0));
            mSubPhaseBins(i).resize(1,tBinCount);
        }

    }

    template<typename Mat_Type>
    void
    row_vector_connectivity_check(moris::Matrix<Mat_Type> & aConnectivity)
    {
        MORIS_ASSERT(moris::isvector(aConnectivity),"Provided connectivity is not a vector");

        if(moris::iscol(aConnectivity))
        {
            aConnectivity = moris::trans(aConnectivity);
        }
    }
};

}
#endif /* SRC_XTK_CL_XTK_CHILD_MESH_HPP_ */
