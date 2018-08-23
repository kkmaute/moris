/*
 * cl_XTK_Child_Mesh_Test.hpp
 *
 *  Created on: Jun 21, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_CHILD_MESH_HPP_
#define SRC_XTK_CL_XTK_CHILD_MESH_HPP_
#include <unordered_map>

#include "linalg/cl_XTK_Matrix.hpp"

#include "containers/cl_XTK_Cell.hpp"

#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/fn_create_faces_from_element_to_node.hpp"
#include "xtk/fn_create_edges_from_element_to_node.hpp"
#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"

#include "mesh/cl_Mesh_Enums.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Child_Mesh_Test
{
public:

    Child_Mesh_Test():
            mElementTopology(EntityTopology::TET_4),
            mChildElementIds(0,0),
            mChildElementInds(0,0),
            mNodeIds(0,0),
            mNodeInds(0,0),
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
            mHasPhaseInfo(false),
            mElementPhaseIndices(0,0),
            mElementBinIndex(0,0),
            mBinBulkPhase(0),
            mSubPhaseBins(0,Mat<Integer,Integer_Matrix>(0,0))
    {};

    Child_Mesh_Test(Integer                     const & aParentElementIndex,
                    Mat<Integer,Integer_Matrix> const & aNodeInds,
                    Mat<Integer,Integer_Matrix> const & aElementToNode,
                    Mat<Integer,Integer_Matrix> const & aElementEdgeParentInds,
                    Mat<Integer,Integer_Matrix> const & aElementEdgeParentRanks,
                    Mat<Integer,Integer_Matrix> const & aElementFaceParentInds,
                    Mat<Integer,Integer_Matrix> const & aElementFaceParentRanks,
                    Mat<Integer,Integer_Matrix> const & aElementInferfaceSides):
                        mElementTopology(EntityTopology::TET_4),
                        mChildElementIds(0,0),
                        mChildElementInds(0,0),
                        mNodeIds(0,0),
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
                        mHasPhaseInfo(false),
                        mElementPhaseIndices(0,0),
                        mElementBinIndex(0,0),
                        mBinBulkPhase(0),
                        mSubPhaseBins(0,Mat<Integer,Integer_Matrix>(0,0))
{
        mParentElementIndex     = aParentElementIndex;
        mNumElem                = aElementToNode.get_num_rows();
        mNodeInds               = aNodeInds.copy();
        mElementToNode          = aElementToNode.copy();
        mElementEdgeParentInds  = aElementEdgeParentInds.copy();
        mElementEdgeParentRanks = aElementEdgeParentRanks.copy();
        mElementFaceParentInds  = aElementFaceParentInds.copy();
        mElementFaceParentRanks = aElementFaceParentRanks.copy();
        mElementInferfaceSides  = aElementInferfaceSides.copy();

        set_up_node_map();
}


    Child_Mesh_Test(Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix> & aMeshModTemplate):
                        mElementTopology(EntityTopology::TET_4),
                        mChildElementIds(0,0),
                        mChildElementInds(0,0),
                        mNodeIds(0,0),
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
                        mHasPhaseInfo(false),
                        mElementPhaseIndices(0,0),
                        mElementBinIndex(0,0),
                        mBinBulkPhase(0),
                        mSubPhaseBins(0,Mat<Integer,Integer_Matrix>(0,0))
{

        reindex_template(aMeshModTemplate);
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
    typename std::unordered_map<Integer,Integer>::iterator Map_Iterator;


    // --------------------------------------------------------------
    // Functions to access connectivity
    // --------------------------------------------------------------

    /*
     * Get number of a entity of a given rank
     */
    Integer
    get_num_entities(enum EntityRank aEntityRank) const
    {
        Integer tNumEntities = 0;
        if(aEntityRank == EntityRank::NODE)
        {
            tNumEntities = mNodeInds.get_num_columns();
        }
        else if(aEntityRank == EntityRank::EDGE)
        {
            XTK_ASSERT(mHasEdgeConn,"Without Edge connectivity it is unclear how many edges there are. Please call generate_connectivities w/ the edge flag on");
            tNumEntities = mEdgeToNode.get_num_rows();
        }

        else if(aEntityRank == EntityRank::FACE)
        {
            XTK_ASSERT(mHasFaceConn,"Without Edge connectivity it is unclear how many edges there are. Please call generate_connectivities w/ the edge flag on");
            tNumEntities = mFaceToNode.get_num_rows();
        }
        else if(aEntityRank == EntityRank::ELEMENT)
        {
            tNumEntities = mElementToNode.get_num_rows();
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
    Mat<Integer,Integer_Matrix> const &
    get_element_to_node() const
    {
        return mElementToNode;
    }

    /*
     * Converts the existing element to node connectivity (which contains proc local indices)
     * to child mesh local indices
     */
    Mat<Integer,Integer_Matrix>
    get_element_to_node_local() const
    {
        return convert_to_local_indices(mElementToNode);
    }


    /*
     * Converts the existing element to node connectivity (which contains proc local indices)
     * to child mesh local indices
     */
    Mat<Integer,Integer_Matrix>
    get_element_to_node_global() const
    {
        Mat<Integer,Integer_Matrix> tElementToNodeCMLoc = convert_to_local_indices(mElementToNode);
        return convert_to_glob_ids(tElementToNodeCMLoc);
    }

    /*
     * Return edge to node
     */
    Mat<Integer,Integer_Matrix> &
    get_edge_to_node()
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return mEdgeToNode;
    }



    /*
     * Return edge to node
     */
    Mat<Integer,Integer_Matrix>
    get_edge_to_node_local() const
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return convert_to_local_indices(mEdgeToNode);
    }

    /*
     * Return element to edge
     */
    Mat<Integer,Integer_Matrix> const &
    get_element_to_edge() const
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return mElementToEdge;
    }

    /*
     * Return element to edge
     */
    Mat<Integer,Integer_Matrix> const &
    get_edge_to_element() const
    {
        XTK_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
        return mEdgeToElement;
    }

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    Mat<Integer,Integer_Matrix> const &
    get_face_to_node() const
    {
        XTK_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return mFaceToNode;
    }

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    Mat<Integer,Integer_Matrix>
    get_face_to_node_local() const
    {
        XTK_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return convert_to_local_indices(mFaceToNode);
    }

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    Mat<Integer,Integer_Matrix> const &
    get_element_to_face() const
    {
        XTK_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return mElementToFace;
    }

    /*
     * Return element to element connectivity (cm local element indices)
     */
    Mat<Integer,Integer_Matrix> const &
    get_element_to_element() const
    {
        XTK_ASSERT(mHasElemToElem,"Element to element connectivity has not been generated with call to generate_element_to_element_connectivity");
        return mElementToElement;
    }


    // Functions to access ancestry

    Integer
    get_parent_element_index() const
    {
        return mParentElementIndex;
    }

    Mat<Integer,Integer_Matrix> const &
    get_face_parent_inds() const
    {
        return mFaceParentInds;
    }

    Mat<Integer,Integer_Matrix> const &
    get_face_parent_ranks() const
    {
        return mFaceParentRanks;
    }



    Mat<Integer,Integer_Matrix> const &
    get_edge_parent_inds() const
    {
        return mEdgeParentInds;
    }

    Mat<Integer,Integer_Matrix> const &
    get_edge_parent_ranks() const
    {
        return mEdgeParentRanks;
    }

    Mat<Integer,Integer_Matrix> const &
    get_node_indices() const
    {
        return mNodeInds;
    }


    Mat<Integer,Integer_Matrix> const &
    get_node_ids() const
    {
        return mNodeIds;
    }

    Mat<Integer,Integer_Matrix> const &
    get_element_ids() const
    {
        return mChildElementIds;
    }


    Mat<Integer,Integer_Matrix> const &
    get_element_inds() const
    {
        return mChildElementInds;
    }


    // Function to access ordinals
    Mat<Integer,Integer_Matrix>
    get_edge_ordinal_from_element_and_edge_indices(Integer const & aElementIndex,
                                                   Mat<Integer,Integer_Matrix> const & aEdgeIndices) const
    {

        // get the elemnt to edge connectivity
        Mat<Integer,Integer_Matrix> const & tElemToEdge = get_element_to_edge();

        // Initialize bounds on loops
        Integer tNumEdges = tElemToEdge.get_num_columns();
        Integer tNumOrdinalstoFind = aEdgeIndices.get_num_columns();

        // Initialize output
        Mat<Integer,Integer_Matrix> tEdgeOrdinals(1,tNumOrdinalstoFind);

        // find the edge ordinals
        Integer tCount = 0;
        for(Integer iOrd = 0; iOrd<tNumOrdinalstoFind; iOrd++)
        {

            for(Integer iEdge = 0; iEdge < tNumEdges; iEdge++)
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
    Integer
    get_face_ordinal_from_element_and_face_index(Integer const & aElementIndex,
                                                 Integer         aFaceIndex) const
    {
        // get the elemnt to edge connectivity
        Mat<Integer,Integer_Matrix> const & tElemToFace = get_element_to_face();

        Integer tNumEdges = tElemToFace.get_num_columns();

        // Initialize output
        Integer tFaceOrdinal = 1000;

        bool tSuccess = false;

        for(Integer iEdge = 0; iEdge < tNumEdges; iEdge++)
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
    get_child_elements_connected_to_parent_face(Integer const & aParentFace,
                                                Mat<Integer,Integer_Matrix> & aChildElemsIdsOnFace,
                                                Mat<Integer,Integer_Matrix> & aChildElemsCMIndOnFace,
                                                Mat<Integer,Integer_Matrix> & aChildElemOnFaceOrdinal) const
    {
        // First determine which of this meshes face live on the parent face
        Integer tNumFaces = get_num_entities(EntityRank::FACE);
        Integer tNumElemsOnFace = 0;
        Mat<Integer,Integer_Matrix> tLocFacesOnParentFace(1,tNumFaces);


        // Iterate through face ancestry ranks
        for(Integer i = 0; i < tNumFaces; i++)
        {
            if(mFaceParentInds(0,i) == aParentFace && mFaceParentRanks(0,i) == 2)
            {
                tLocFacesOnParentFace(0,tNumElemsOnFace) = i;
                tNumElemsOnFace++;
            }
        }

        // Iterate through face to element connectivities of faces connected to the parent face
        aChildElemsIdsOnFace    = Mat<Integer,Integer_Matrix>(1,tNumElemsOnFace*2);
        aChildElemsCMIndOnFace  = Mat<Integer,Integer_Matrix>(1,tNumElemsOnFace*2);
        aChildElemOnFaceOrdinal = Mat<Integer,Integer_Matrix>(1,tNumElemsOnFace*2);
        Integer tCount = 0;

        for(Integer i = 0; i<tNumElemsOnFace; i++)
        {
            Integer tFaceInd = tLocFacesOnParentFace(0,i);
            for(Integer j = 0; j<mFaceToElement.get_num_columns(); j++)
            {
                if(mFaceToElement(tFaceInd,j) != std::numeric_limits<Integer>::max())
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
    add_node_indices(Mat<Integer,Integer_Matrix> const & aNewNodeInds)
    {
        // Allocate space
        Integer tNumNewNodes = aNewNodeInds.get_num_columns();
        Integer tNumCurrNodes = mNodeInds.get_num_columns();

        // add nodes to map
        add_nodes_to_map(aNewNodeInds);

        mNodeInds.resize(1,tNumCurrNodes+tNumNewNodes);

        // add nodes to matrix
        for(Integer i = 0; i<tNumNewNodes; i++)
        {
            mNodeInds(0,i+tNumCurrNodes) = aNewNodeInds(0,i);
        }

    }

    /*
     * Add more nodes to the existing node indices list
     */
    void
    add_node_ids(Mat<Integer,Integer_Matrix> const & aNewNodeIds)
    {
        // Allocate space
        Integer tNumNewNodes = aNewNodeIds.get_num_columns();
        Integer tNumCurrNodes = mNodeIds.get_num_columns();

        mNodeIds.resize(1,tNumCurrNodes+tNumNewNodes);

        // add nodes to matrix
        for(Integer i = 0; i<tNumNewNodes; i++)
        {
            mNodeIds(0,i+tNumCurrNodes) = aNewNodeIds(0,i);
        }
    }


    /*
     * sets the node ids. note this overwrites existing mNodeIds data
     */

    void
    set_node_ids(Mat<Integer,Integer_Matrix> const & aNewNodeIds)
    {
        XTK_ASSERT(aNewNodeIds.get_num_columns() == get_num_entities(EntityRank::NODE),"Number of node ids does not match the number of nodes in this child mesh");
        mNodeIds = aNewNodeIds.copy();
    }

    /*
     * Sets the globally unique element Ids for the child mesh. This is important for mesh from data creation
     * @param[in] aElementId - First element Id (note: this is incremented as the id is assigned)
     */
    void set_child_element_ids(Integer & aElementId)
    {
        XTK_ASSERT(mChildElementIds.get_num_columns() == 0, "Element Ids already set");
        Integer tNumElements = get_num_entities(EntityRank::ELEMENT);
        mChildElementIds = Mat<Integer,Integer_Matrix>(1,tNumElements);

        for(Integer iElem = 0; iElem<tNumElements; iElem++)
        {
            mChildElementIds(0,iElem) = aElementId;
            aElementId++;
        }
    }

    /*
     * Sets the processor unique element indices for the child mesh.
     * @param[in] aElementInd - First element Ind (note: this is incremented as the ind is assigned)
     */
    void set_child_element_inds(Integer & aElementInd)
    {
        XTK_ASSERT(mChildElementInds.get_num_columns() == 0, "Element Inds already set");


        Integer tNumElements = get_num_entities(EntityRank::ELEMENT);
        mChildElementInds = Mat<Integer,Integer_Matrix>(1,tNumElements);

        for(Integer iElem = 0; iElem<tNumElements; iElem++)
        {
            mChildElementInds(0,iElem) = aElementInd;
            aElementInd++;
        }
    }

    /*
     * Resizes element related matrices
     */
    void
    allocate_more_elements(Integer const & aNumMoreElements)
    {
        Integer tCurrNumElem = mNumElem;
        Integer tNewNumElem = tCurrNumElem + aNumMoreElements;

        mElementToNode.resize(tNewNumElem,mElementToNode.get_num_columns());
        mElementEdgeParentInds.resize(tNewNumElem,mElementEdgeParentInds.get_num_columns());
        mElementEdgeParentRanks.resize(tNewNumElem,mElementEdgeParentRanks.get_num_columns());
        mElementFaceParentInds.resize(tNewNumElem,mElementFaceParentInds.get_num_columns());
        mElementFaceParentRanks.resize(tNewNumElem,mElementFaceParentRanks.get_num_columns());
        mElementInferfaceSides.resize(tNewNumElem,mElementInferfaceSides.get_num_columns());
    }

    void
    insert_child_mesh_template(Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix> & aMeshModTemplate)
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
        Integer tNumElemToAdd = aMeshModTemplate.mNumNewElem-aMeshModTemplate.mNumElemToReplace;
        for(Integer i = 0; i<tNumElemToAdd; i++)
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
    void reindex_template(Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix> & aMeshModTemplate)
    {
        aMeshModTemplate.mNewElementToNode = reindex_matrix(aMeshModTemplate.mNewElementToNode, 0, aMeshModTemplate.mNodeInds);

        // Reindex template edge and face ordinals
        reindex_template_parent_information(aMeshModTemplate);
    }


    /*
     * Replace the information associated with an element
     */
    void
    replace_element(Integer                     const & aElementIndexToReplace,
                    Integer                     const & aRowIndex,
                    Mat<Integer,Integer_Matrix> const & aElementToNode,
                    Mat<Integer,Integer_Matrix> const & aElementEdgeParentInds,
                    Mat<Integer,Integer_Matrix> const & aElementEdgeParentRanks,
                    Mat<Integer,Integer_Matrix> const & aElementFaceParentInds,
                    Mat<Integer,Integer_Matrix> const & aElementFaceParentRanks,
                    Mat<Integer,Integer_Matrix> const & aElementInterfaceFaces)
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
    add_element(Integer                     const & aRowIndex,
                Mat<Integer,Integer_Matrix> const & aElementToNode,
                Mat<Integer,Integer_Matrix> const & aElementEdgeParentInds,
                Mat<Integer,Integer_Matrix> const & aElementEdgeParentRanks,
                Mat<Integer,Integer_Matrix> const & aElementFaceParentInds,
                Mat<Integer,Integer_Matrix> const & aElementFaceParentRanks,
                Mat<Integer,Integer_Matrix> const & aElementInterfaceFaces)
    {
        XTK_ASSERT(mNumElem<mElementToNode.get_num_rows(),"Not enough space allocated in call to allocate_more_elements");
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
    generate_connectivities(bool mGenerateFaceConn,
                            bool mGenerateEdgeConn,
                            bool mGenerateElemToElem)
    {
        Mat<Integer,Integer_Matrix> tElementToNodeLoc = get_element_to_node_local();

        if( mGenerateFaceConn )
        {
            generate_face_connectivity_and_ancestry(tElementToNodeLoc);
        }

        if( mGenerateEdgeConn )
        {
            generate_edge_connectivity_and_ancestry(tElementToNodeLoc);
        }

        if(mGenerateElemToElem)
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
        Integer tNumD = get_num_entities(EntityRank::ELEMENT);
        Integer tNumEdgesToElem = get_element_to_edge().get_num_columns();
        mIntersectConnectivity = Mat<Integer,Integer_Matrix>(tNumD, tNumEdgesToElem * 2 + 1, 0); // Needs to be zero for use column
    }

    /**
      * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
      *       - 1 means the provided aDPrime1Ind is an XTK index
      *
      * aDPrime2Ind must be XTK local index
      */
     void add_entity_to_intersect_connectivity(Integer aCMNodeInd,
                                               Integer aCMEdgeInd,
                                               Integer aFlag)
     {
         if(aFlag == 0)
         {
             aCMNodeInd = aCMNodeInd + get_num_entities(EntityRank::NODE);
         }

         // Get entities of dimension mAuxDPrime2 connected to entities of mAuxD (i.e. elements connected to given edge aDPrime2Ind)
         Mat<Integer,Integer_Matrix> const tEdgeToElem = get_edge_to_element();


         Integer tNumElemsConnected = 0;
         for(Integer i = 0; i<tEdgeToElem.get_num_columns(); i++)
         {
             if(tEdgeToElem(aCMEdgeInd, i)== std::numeric_limits<Integer>::max())
             {
                 break;
             }
             tNumElemsConnected++;
         }

         for(Integer i = 0; i < tNumElemsConnected; i++)
         {
             Integer tElemCMInd = tEdgeToElem(aCMEdgeInd, i);

             XTK_ASSERT(mIntersectConnectivity(tElemCMInd, 0) < this->get_element_to_edge().get_num_columns(),
                        "Entity corresponding to provided aDInd has exceeded allocated space");
             XTK_ASSERT(tElemCMInd < mIntersectConnectivity.get_num_rows(),
                        "aDInd is outside of bounds. Has auxiliary connectivity been initialized?");

             mIntersectConnectivity(tElemCMInd, mIntersectConnectivity(tElemCMInd, 0) + 1) = aCMNodeInd;
             mIntersectConnectivity(tElemCMInd, mIntersectConnectivity(tElemCMInd, 0) + mElementToEdge.get_num_columns() + 1) = aCMEdgeInd;
             mIntersectConnectivity(tElemCMInd, 0)++;
         }

     }

    /*
     * Tracks intersection connectivity prior to the child mesh being modified
     */
      void
      set_intersect_connectivity(Mat<Integer,Integer_Matrix> const & aIntConnectivity)
      {
          // Copy so the data is completely owned by the child mesh
          mIntersectConnectivity = aIntConnectivity.copy();
      }

      /**
       * Tells the child mesh where a node index will be placed once it has been communicated
       */
      void set_pending_node_index_pointers(Cell<Integer*> aNodeIndPtr)
      {
          mPtrPendingNodeIndex = aNodeIndPtr;
      }

      /**
       *  retrieve_pending_node_inds
       * @param aInterface - indicates the node is on the interface (set to false by default
       * XTK mesh retrieves pending node indices. Prior to this call the unique node assignments
       * must be complete via the request structure in XTK model
       */
      void retrieve_pending_node_inds()
      {
          Mat<Integer,Integer_Matrix> tNodeMat(1, mPtrPendingNodeIndex.size(), 0);

          for(Integer i = 0; i < mPtrPendingNodeIndex.size(); i++)
          {
              tNodeMat(0, i) = *mPtrPendingNodeIndex(i);
          }

          add_node_indices(tNodeMat);

          // size out since the nodes have been retrieved
          mPtrPendingNodeIndex.resize(0, NULL);
      }


      // --------------------------------------------------------------
      // Functions for elemental phase information and subphase bins
      // --------------------------------------------------------------

      void
      initialize_element_phase_mat()
      {
          mElementPhaseIndices = Mat<Integer,Integer_Matrix>(1,get_num_entities(EntityRank::ELEMENT),std::numeric_limits<Integer>::max());
      }

      void set_element_phase_index(Integer aEntityIndex,
                                   Integer aEntityPhaseIndex)
      {
          XTK_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
          mHasPhaseInfo = true;
          mElementPhaseIndices(0,aEntityIndex) = aEntityPhaseIndex;
      }

      Integer get_element_phase_index( Integer const & aEntityIndex) const
      {
          XTK_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
          XTK_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
          return mElementPhaseIndices(0,aEntityIndex);
      }

      Mat<Integer,Integer_Matrix> const &
      get_element_phase_indices() const
      {
          XTK_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
          return mElementPhaseIndices;
      }

      Integer
      get_num_subphase_bins()
      {
          return mSubPhaseBins.size();
      }

      Cell<Integer> const &
      get_subphase_bin_bulk_phase()
      {
          return mBinBulkPhase;
      }



      /*
       * Sets the elemental subphase value in this child mesh
       */
      void
      set_elemental_subphase(Mat<Integer,Integer_Matrix> const & aElementSubPhase)
      {
          construct_subphase_bins(aElementSubPhase);
      }

      /*
       * Returns the elemental subphase values
       */
      Mat<Integer,Integer_Matrix> const &
      get_elemental_subphase_bin_membership()
      {
          return mElementBinIndex;
      }

      // --------------------------------------------------------------
      // Functions IO
      // --------------------------------------------------------------

      void pack_child_mesh_by_phase(Integer const & aNumPhases,
                                    Cell<Mat<Integer,Integer_Matrix>> & aElementIds,
                                    Cell<Mat<Integer,Integer_Matrix>> & aElementCMInds) const
      {
          Integer tNumElems = get_num_entities(EntityRank::ELEMENT);

          aElementIds = Cell<Mat<Integer,Integer_Matrix>>(aNumPhases);
          aElementCMInds = Cell<Mat<Integer,Integer_Matrix>>(aNumPhases);

          for(Integer i = 0; i<aNumPhases; i++)
          {
              aElementIds(i) = Mat<Integer,Integer_Matrix>(1,tNumElems);
              aElementCMInds(i) = Mat<Integer,Integer_Matrix>(1,tNumElems);
          }

          Cell<Integer> tPhaseCounter(aNumPhases,0);

          Mat<Integer,Integer_Matrix> const & tElementPhaseIndices = get_element_phase_indices();
          Mat<Integer,Integer_Matrix> const & tElementIds  = get_element_ids();

          for(Integer i = 0; i < tNumElems; i++)
          {
              Integer tPhaseIndex = tElementPhaseIndices(0,i);
              Integer tPhaseCount = tPhaseCounter(tPhaseIndex);
              aElementIds(tPhaseIndex)(0,tPhaseCount) = tElementIds(0,i);
              aElementCMInds(tPhaseIndex)(0,tPhaseCount) = i;
              tPhaseCounter(tPhaseIndex)++;
          }

          for(Integer i = 0; i<aNumPhases; i++)
          {
              aElementIds(i).resize(1,tPhaseCounter(i));
              aElementCMInds(i).resize(1,tPhaseCounter(i));
          }
      }


      void
      pack_interface_sides(Output_Options<Integer> const & aOutputOptions,
                           Mat<Integer, Integer_Matrix> & aElementIds,
                           Mat<Integer, Integer_Matrix> & aSideOrdinals) const
      {
          // Loop bound and sizing
          Integer tNumElem = get_num_entities(EntityRank::ELEMENT);
          aElementIds   = Mat<Integer,Integer_Matrix>(1,tNumElem);
          aSideOrdinals = Mat<Integer,Integer_Matrix>(1,tNumElem);

          // Keep track of the number of interface sides
          Integer tCount = 0;

          // Reference to elemental phase vector
          Mat<Integer,Integer_Matrix> const & tElemPhases = get_element_phase_indices();

         // Iterate over each element and if the element has an interface side it will be in mElementInferfaceSides vector
         for(Integer iEl =0 ; iEl<tNumElem; iEl++)
         {
             //TODO: NOTE THIS WILL NOT WORK WITH MULTI-MATERIAL YET
             if(aOutputOptions.output_phase(tElemPhases(0,iEl)) && mElementInferfaceSides(iEl,0) != std::numeric_limits<Integer>::max())
             {
                 aElementIds(0,tCount)   = mChildElementIds(0,iEl);
                 aSideOrdinals(0,tCount) = mElementInferfaceSides(iEl,0);
                 tCount++;
             }
         }

         // Size out space
         aElementIds.resize(1,tCount);
         aSideOrdinals.resize(1,tCount);
      }

private:
    // Parent element index
    Integer                     mParentElementIndex;

    // Element To Node and Ancestry Information (This is the only data that is set with templates. All other is generated with
    // an algorithm)
    // All node connectivity is indexed by proc local indexs
    enum EntityTopology         mElementTopology;
    Integer                     mNumElem;
    Mat<Integer,Integer_Matrix> mElementToNode; // node indices correspond to the local child mesh index
    Mat<Integer,Integer_Matrix> mElementEdgeParentInds;
    Mat<Integer,Integer_Matrix> mElementEdgeParentRanks;
    Mat<Integer,Integer_Matrix> mElementFaceParentInds;
    Mat<Integer,Integer_Matrix> mElementFaceParentRanks;
    Mat<Integer,Integer_Matrix> mElementInferfaceSides;

    // Child element information
    Mat<Integer,Integer_Matrix> mChildElementIds;
    Mat<Integer,Integer_Matrix> mChildElementInds;

    // Nodes related to mesh index
    Mat<Integer,Integer_Matrix>          mNodeIds;
    Mat<Integer,Integer_Matrix>          mNodeInds;
    std::unordered_map<Integer, Integer> mNodeIndsToCMInd; // key - proc local ind, val - local child mesh index

    // Face Connectivity
    bool mHasFaceConn;
    Mat<Integer,Integer_Matrix> mFaceToNode;
    Mat<Integer,Integer_Matrix> mNodeToFace;
    Mat<Integer,Integer_Matrix> mFaceToElement;
    Mat<Integer,Integer_Matrix> mElementToFace;
    Mat<Integer,Integer_Matrix> mFaceParentInds;
    Mat<Integer,Integer_Matrix> mFaceParentRanks;

    // Edge connectivity
    bool mHasEdgeConn;
    Mat<Integer,Integer_Matrix> mEdgeToNode;
    Mat<Integer,Integer_Matrix> mNodeToEdge;
    Mat<Integer,Integer_Matrix> mEdgeToElement;
    Mat<Integer,Integer_Matrix> mElementToEdge;
    Mat<Integer,Integer_Matrix> mEdgeParentInds;
    Mat<Integer,Integer_Matrix> mEdgeParentRanks;

    // Element to Element graph
    bool mHasElemToElem;
    Mat<Integer,Integer_Matrix> mElementToElement;

    // Auxiliary connectivity data and pending nodes (mesh modification data)
    Mat<Integer,Integer_Matrix> mIntersectConnectivity;
    Cell<Integer*>              mPtrPendingNodeIndex;

    // Phase member variables (structured)
    bool                              mHasPhaseInfo;
    Mat<Integer,Integer_Matrix>       mElementPhaseIndices;
    Mat<Integer,Integer_Matrix>       mElementBinIndex;
    Cell<Integer>                     mBinBulkPhase;
    Cell<Mat<Integer,Integer_Matrix>> mSubPhaseBins;

private:
    void
    generate_face_connectivity_and_ancestry(Mat<Integer,Integer_Matrix> const & aElementToNodeLocal)
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
    generate_edge_connectivity_and_ancestry(Mat<Integer,Integer_Matrix> const & aElementToNodeLocal)
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
        Integer tNumFacePerElem = 4;
        mElementToElement      = generate_element_to_element(mFaceToElement,
                                                              mNumElem,
                                                              tNumFacePerElem,
                                                              std::numeric_limits<Integer>::max());
    }

    void
    set_up_node_map()
    {
        Integer tNumNodes = get_num_entities(EntityRank::NODE);

        for( Integer i = 0; i<tNumNodes; i++)
        {
            if(mNodeIndsToCMInd.find(mNodeInds(0,i)) == mNodeIndsToCMInd.end())
            {
                mNodeIndsToCMInd[mNodeInds(0,i)] = i;
            }

            else
            {
                Integer tBreaker = 0;
                XTK_ASSERT(tBreaker!=0," Attempted to add duplicate node in constructor, are your node indices correct?");
            }
        }
    }


    /*
     * Note: this function only adds nodes which do not already exist to the map
     * and ignores nodes which are already in the map.
     */
    void
    add_nodes_to_map(Mat<Integer,Integer_Matrix> const & aNodesToAdd)
    {
        Integer tNumNodes = get_num_entities(EntityRank::NODE);
        Integer tNumNewNodes = aNodesToAdd.get_num_columns();

        for( Integer i = 0; i<tNumNewNodes; i++)
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
        Integer tNumElements = get_num_entities(EntityRank::ELEMENT);
        Integer tNumFacePerElement = mElementToFace.get_num_columns();
        Integer tNumFaces = get_num_entities(EntityRank::FACE);

        mFaceParentInds  = Mat<Integer,Integer_Matrix>(1,tNumFaces);
        mFaceParentRanks = Mat<Integer,Integer_Matrix>(1,tNumFaces);

        for( Integer i = 0; i<tNumElements; i++)
        {
            for(Integer j = 0; j<tNumFacePerElement; j++)
            {
                Integer tFaceIndex = mElementToFace(i,j);

                mFaceParentInds(0,tFaceIndex) = mElementFaceParentInds(i,j);
                mFaceParentRanks(0,tFaceIndex) = mElementFaceParentRanks(i,j);
            }
        }

    }

    void
    setup_edge_ancestry()
    {
        Integer tNumElements = get_num_entities(EntityRank::ELEMENT);
        Integer tNumEdgePerElement = mElementToEdge.get_num_columns();
        Integer tNumFaces = get_num_entities(EntityRank::EDGE);

        mEdgeParentInds  = Mat<Integer,Integer_Matrix>(1,tNumFaces);
        mEdgeParentRanks = Mat<Integer,Integer_Matrix>(1,tNumFaces);


        for( Integer i = 0; i<tNumElements; i++)
        {
            for(Integer j = 0; j<tNumEdgePerElement; j++)
            {
                Integer tFaceIndex = mElementToEdge(i,j);

                mEdgeParentInds(0,tFaceIndex) = mElementEdgeParentInds(i,j);
                mEdgeParentRanks(0,tFaceIndex) = mElementEdgeParentRanks(i,j);
            }
        }

    }



    Mat<Integer,Integer_Matrix>
    convert_to_local_indices(Mat<Integer,Integer_Matrix> const & aEntityRankToNode) const
    {
        Integer tNumRows = aEntityRankToNode.get_num_rows();
        Integer tNumCols = aEntityRankToNode.get_num_columns();

        Mat<Integer,Integer_Matrix> tLocalEntityRankToNode(tNumRows,tNumCols);


        for(Integer i = 0; i < tNumRows; i++)
        {
            for(Integer j = 0; j < tNumCols; j++)
            {

                auto tIter = mNodeIndsToCMInd.find(aEntityRankToNode(i,j));

                XTK_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

                tLocalEntityRankToNode(i,j) = tIter->second;
            }
        }
        return tLocalEntityRankToNode;
    }


    Mat<Integer,Integer_Matrix>
    convert_to_proc_indices(Mat<Integer,Integer_Matrix> const & aEntityRankToNodeLocal) const
    {
        Integer tNumRows = aEntityRankToNodeLocal.get_num_rows();
        Integer tNumCols = aEntityRankToNodeLocal.get_num_columns();

        Mat<Integer,Integer_Matrix> tProcEntityRankToNode(tNumRows,tNumCols);

        for(Integer i = 0; i < tNumRows; i++)
        {
            for(Integer j = 0; j < tNumCols; j++)
            {
                tProcEntityRankToNode(i,j) = mNodeInds(0,aEntityRankToNodeLocal(i,j));
            }
        }

        return tProcEntityRankToNode;
    }

    Mat<Integer,Integer_Matrix>
    convert_to_glob_ids(Mat<Integer,Integer_Matrix> const & aEntityRankToNodeLocal) const
    {
        Integer tNumRows = aEntityRankToNodeLocal.get_num_rows();
        Integer tNumCols = aEntityRankToNodeLocal.get_num_columns();

        Mat<Integer,Integer_Matrix> tProcEntityRankToNode(tNumRows,tNumCols);

        for(Integer i = 0; i < tNumRows; i++)
        {
            for(Integer j = 0; j < tNumCols; j++)
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

            XTK_ASSERT(mIntersectConnectivity.get_num_rows() == mNumElem, "There needs to be a row for each element in the child mesh in the intersect connectivity");
            XTK_ASSERT(aTemplate == TemplateType::HIERARCHY_TET4,"This function needs to be abstracted and has only been tested with HIERARCHY_TET4.");

            // Iterate over number of elements (only ones that existed at the beginning of the modification)
            Integer tNumExistElem = mNumElem;

            // Container for all the templates to add
            Cell<Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix>> tTemplatesToAdd(tNumExistElem);

            // Keep track of the number of intersected elements
            Integer tNumIntersected = 0;
            Integer tNumNewElem     = 0;

            Mat<Integer, Integer_Matrix> tEdgeToNodeCMLoc = get_edge_to_node_local();

            for(Integer iE = 0; iE<tNumExistElem; iE++)
            {
                if(mIntersectConnectivity(iE,0) == 3 || mIntersectConnectivity(iE,0) == 4)
                {
                    Integer tPermutationId = std::numeric_limits<Integer>::max();

                    Mat<Integer,Integer_Matrix> tIntersectConnRow = mIntersectConnectivity.get_row(iE);

                    Mat<Integer,Integer_Matrix> tSortedNodes = sort_nodes(aTemplate, tIntersectConnRow, tEdgeToNodeCMLoc, iE, tPermutationId);



                    // Select more specific tet4 hierarchy template
                    // TODO: Clean these if statements up
                    enum TemplateType tNarrowTemplate;
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
                    Mat<Integer,Integer_Matrix> tElementsAncestry({{mParentElementIndex}}); // Not used
                    Mat<Integer,Integer_Matrix> tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    Mat<Integer,Integer_Matrix> tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    Mat<Integer,Integer_Matrix> tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    Mat<Integer,Integer_Matrix> tParentFaceRanks = mElementFaceParentRanks.get_row(iE);



                    tSortedNodes = convert_to_proc_indices(tSortedNodes);


                    // Setup template with this information
                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix>
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

                else if (mIntersectConnectivity(iE,0) == 0)
                {
                    continue;
                }
                else
                {
//                    XTK_ERROR << "Invalid connectivity for nodal hierarchy template, (should be 3 or 4 nodes)\n";
                }
            }

            // Allocate space for new elements
            allocate_more_elements(tNumNewElem);

            // Add templates to the connectivity
            for(Integer i = 0; i<tNumIntersected; i++)
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

            Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix> tRegSubTemplate(mParentElementIndex,
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
    Mat<Integer,Integer_Matrix>
    sort_nodes(enum TemplateType                    aTemplate,
               Mat<Integer, Integer_Matrix> const & aIntConnectivity,
               Mat<Integer, Integer_Matrix> const & aEdgeToNodeCMLoc,
               Integer const &                      aElementIndex,
               Integer &                            aPermutation)
    {
        //Locate highest node in intersection connectivity
        switch(aTemplate)
        {
        case (TemplateType::HIERARCHY_TET4):
                         {

            if(aIntConnectivity(0, 0) == 3)
            {
                Integer tIntersectionCase = std::numeric_limits<Integer>::max();
                Mat<Integer, Integer_Matrix> tHigh({{ aIntConnectivity(0, 1), aIntConnectivity(0, 7)}});
                Mat<Integer, Integer_Matrix> tMid({{ aIntConnectivity(0, 1), aIntConnectivity(0, 7)}});
                Mat<Integer, Integer_Matrix> tLow({{ aIntConnectivity(0, 1), aIntConnectivity(0, 7)}});

                for(Integer i = 0; i < 2; i++)
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


                Mat<Integer, Integer_Matrix> tNodes13 = aEdgeToNodeCMLoc.get_row(tMid(0, 1));
                Mat<Integer, Integer_Matrix> tNodes12 = aEdgeToNodeCMLoc.get_row(tLow(0, 1));
                Mat<Integer, Integer_Matrix> tNodes14 = aEdgeToNodeCMLoc.get_row(tHigh(0, 1));

                // Initialize Edge matrix for the 3 node intersection
                Mat<Integer, Integer_Matrix> tEdgeIndices({{tLow(0, 1),tMid(0, 1),tHigh(0, 1)}});

                // Get edge ordinals of the edge indices
                Mat<Integer, Integer_Matrix> tEdgeOrdinals = get_edge_ordinal_from_element_and_edge_indices(aElementIndex,tEdgeIndices);
                get_intersection_permutation(tEdgeOrdinals,aPermutation);

                // Find the shared node in all the above lists
                // Decide which node is which (using intersections)
                Integer tN1 = std::numeric_limits<Integer>::max();
                Integer tN2 = std::numeric_limits<Integer>::max();
                Integer tN3 = std::numeric_limits<Integer>::max();
                Integer tN4 = std::numeric_limits<Integer>::max();

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

                Mat<Integer, Integer_Matrix> tSortedNodes(
                        {{tN1, tN2, tN3, tN4, tLow(0, 0), tMid(0, 0), tHigh(0, 0), tIntersectionCase}});

                return tSortedNodes;
                         }

            else if(aIntConnectivity(0, 0) == 4)
            {
                // Sort Auxiliary nodes from highest to lowest using bubble sort (but based on first column then swap row

                Mat<Integer, Integer_Matrix> tSortedNodes(1, 9);
                Mat<Integer, Integer_Matrix> tNodes(
                {
                    {aIntConnectivity(0, 1), aIntConnectivity(0, 7)},
                    {aIntConnectivity(0, 2), aIntConnectivity(0, 8)},
                    {aIntConnectivity(0, 3), aIntConnectivity(0, 9)},
                    {aIntConnectivity(0, 4), aIntConnectivity(0, 10)}});

                Integer j = 0;
                Integer n = 4; // 4 numbers to sort
                bool swapped = true;

                // Temporary row storage
                Mat<Integer, Integer_Matrix> tRowStorage(1, 2);
                Mat<Integer, Integer_Matrix> tRowSwapper(1, 2);

                while(swapped)
                {
                    swapped = false;
                    j++;
                    for(Integer i = 0; i < n - j; i++)
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

                Mat<Integer, Integer_Matrix> tEdgeToNode = get_edge_to_node_local();

                Mat<Integer,Integer_Matrix> tNodesL  = tEdgeToNode.get_row(tNodes(0, 1));
                Mat<Integer,Integer_Matrix> tNodesML = tEdgeToNode.get_row(tNodes(1, 1));
                Mat<Integer,Integer_Matrix> tNodesMH = tEdgeToNode.get_row(tNodes(2, 1));
                Mat<Integer,Integer_Matrix> tNodesH  = tEdgeToNode.get_row(tNodes(3, 1));

                Integer tHLOppFlag  = 1;
                Integer tHMHOppFlag = 1;
                Integer tHMLOppFlag = 1;

                // Initialize Edge matrix for the 3 node intersection
                Mat<Integer,Integer_Matrix> tEdgeIndices({{tNodes(0, 1),tNodes(1, 1),tNodes(2, 1),tNodes(3, 1)}});

                Mat<Integer, Integer_Matrix> tEdgeOrdinals = get_edge_ordinal_from_element_and_edge_indices(aElementIndex,tEdgeIndices);

                get_intersection_permutation(tEdgeOrdinals,aPermutation);


                // If L and H share a node then it is not case a
                for(Integer i = 0; i < 2; i++)
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
                for(Integer i = 0; i < 2; i++)
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
                for(Integer i = 0; i < 2; i++)
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
                    Integer tSuccess = 0;
                    for(Integer i = 0; i < 2; i++)
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
                    for(Integer i = 0; i < 2; i++)
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
                    Integer tSuccess = 0;
                    for(Integer i = 0; i < 2; i++)
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
                    for(Integer i = 0; i < 2; i++)
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
                    Integer tSuccess = 0;
                    for(Integer i = 0; i < 2; i++)
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
                    for(Integer i = 0; i < 2; i++)
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
                Mat<Integer, Integer_Matrix> dummy(1, 1);
                return dummy;
            }

            break;
                         }

        default:
        {
            XTK_ERROR << "Sorting for specified template type not implemented";

            Mat<Integer, Integer_Matrix> dummy(1, 1);
            return dummy;

            break;
        }

        }

    }


    /*
    * This assumes the edge ordinals provided are ordered in the same order as is found in the auxiliary connectivity
    */
    void
    get_intersection_permutation(Mat<Integer,Integer_Matrix> const & aOrderedEdgeOrdinals,
                                 Integer & aPermutation)
    {
        if(aOrderedEdgeOrdinals.get_num_columns() == 3)
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

        else if (aOrderedEdgeOrdinals.get_num_columns() == 4)
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
    reindex_template_parent_information(Mesh_Modification_Template<Real,Integer,Real_Matrix,Integer_Matrix> & aMeshModTemplate)
    {
        // Reindex parent ordinals for edges
        Integer tNumEdgePerElem = aMeshModTemplate.mNewParentEdgeOrdinals.get_num_columns();
        Integer tNumFacePerElem = aMeshModTemplate.mNewParentFaceOrdinals.get_num_columns();
        Integer tNumElems       = aMeshModTemplate.mNumNewElem;
        Mat<Integer,Integer_Matrix> tParentElem({{aMeshModTemplate.mParentElemInd}});
        Mat<Integer,Integer_Matrix> tParentElemRank({{3}});

        // Place ptrs in cell to avoid if statement checking rank in for loop
        Cell<Mat<Integer,Integer_Matrix>*> tParentEntitiesInds({&aMeshModTemplate.mParentEdgeInds,
                                                            &aMeshModTemplate.mParentFaceInds,
                                                            &tParentElem});

        Cell<Mat<Integer,Integer_Matrix>*> tParentEntitiesRanks({& aMeshModTemplate.mParentEdgeRanks,
                                                                 & aMeshModTemplate.mParentFaceRanks,
                                                                 & tParentElemRank});


        for(Integer i = 0; i<tNumElems; i++)
        {
            // Reindex edges
            for(Integer j = 0; j<tNumEdgePerElem; j++)
            {
                Integer tEdgeParentRank = aMeshModTemplate.mNewParentEdgeRanks(i,j);
                aMeshModTemplate.mNewParentEdgeRanks(i,j) = (*tParentEntitiesRanks(tEdgeParentRank-1))(0,aMeshModTemplate.mNewParentEdgeOrdinals(i,j));
                aMeshModTemplate.mNewParentEdgeOrdinals(i,j) = (*tParentEntitiesInds(tEdgeParentRank-1))(0,aMeshModTemplate.mNewParentEdgeOrdinals(i,j));

            }

            for(Integer j = 0; j<tNumFacePerElem; j++)
            {
                Integer tFaceParentRank = aMeshModTemplate.mNewParentFaceRanks(i,j);
                aMeshModTemplate.mNewParentFaceRanks(i,j) = (*tParentEntitiesRanks(tFaceParentRank-1))(0,aMeshModTemplate.mNewParentFaceOrdinals(i,j));
                aMeshModTemplate.mNewParentFaceOrdinals(i,j) = (*tParentEntitiesInds(tFaceParentRank-1))(0,aMeshModTemplate.mNewParentFaceOrdinals(i,j));
            }
        }
    }

    void
    construct_subphase_bins(Mat<Integer,Integer_Matrix> const & aElementSubPhase)
    {
        // Number of bins corresponds to the maxmimum value in the element sub phase vector
        Integer tNumBins = aElementSubPhase.get_max_value() + 1;
        Integer tNumElements = this->get_num_entities(EntityRank::ELEMENT);

        // Initialize member variables
        // Element Sub-phase bins
        mBinBulkPhase = Cell<Integer>(tNumBins);
        mElementBinIndex = aElementSubPhase.copy();

        Mat<Integer,Integer_Matrix> tBinSizeCounter(1,tNumBins,0);
        mSubPhaseBins = Cell<Mat<Integer,Integer_Matrix>>(tNumBins);
        for(Integer i = 0; i<tNumBins; i++)
        {
            mSubPhaseBins(i) = Mat<Integer,Integer_Matrix>(1,tNumElements);
        }

        // Place the elements into bins;
        for(Integer i = 0; i<tNumElements; i++)
        {
            Integer tBinIndex = aElementSubPhase(0,i);
            Integer tBinCount = tBinSizeCounter(0,tBinIndex);

            tBinSizeCounter(0,tBinIndex)++;

            mSubPhaseBins(tBinIndex)(0,tBinCount) = i;
        }

        // Size out extra space and set bin element bulk phase
        for(Integer i = 0; i<tNumBins; i++)
        {
            Integer tBinCount = tBinSizeCounter(0,i);
            mBinBulkPhase(i) = get_element_phase_index(mSubPhaseBins(i)(0,0));
            mSubPhaseBins(i).resize(1,tBinCount);
        }

    }
};

}
#endif /* SRC_XTK_CL_XTK_CHILD_MESH_HPP_ */
