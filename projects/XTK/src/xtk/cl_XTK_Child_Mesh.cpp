/*
 * cl_XTK_Child_Mesh.cpp
 *
 *  Created on: Mar 14, 2019
 *      Author: doble
 */

#include "cl_XTK_Child_Mesh.hpp"
#include "assert.hpp"
#include "cl_MTK_Vertex.hpp"

using namespace moris;
namespace xtk
{
Child_Mesh::Child_Mesh():
                    mElementTopology(CellTopology::TET4),
                    mConnectivity(nullptr),
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
                    mSubPhaseBinIndices(0),
                    mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
            {}


// ----------------------------------------------------------------------------------


Child_Mesh::Child_Mesh(moris::uint                        aSpatialDimension,
                       moris::moris_index                 aParentElementIndex,
                       moris::Matrix< moris::IndexMat > & aNodeInds,
                       moris::Matrix< moris::IndexMat > & aElementNodeParentInds,
                       moris::Matrix< moris::DDSTMat >  & aElementNodeParentRanks,
                       moris::Matrix< moris::IndexMat > & aElementToNode,
                       moris::Matrix< moris::IndexMat > & aElementEdgeParentInds,
                       moris::Matrix< moris::DDSTMat >  & aElementEdgeParentRanks,
                       moris::Matrix< moris::IndexMat > & aElementFaceParentInds,
                       moris::Matrix< moris::DDSTMat >  & aElementFaceParentRanks,
                       moris::Matrix< moris::DDSTMat >  & aElementInferfaceSides ):
                        mElementTopology(CellTopology::TET4),
                        mConnectivity(nullptr),
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
                        mSubPhaseBinIndices(0),
                        mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
{
        // Check for row vector connectivity (if not it is transposed)
        row_vector_connectivity_check( aNodeInds );
        row_vector_connectivity_check( aElementNodeParentInds  );
        row_vector_connectivity_check( aElementNodeParentRanks );
        row_vector_connectivity_check( aElementEdgeParentInds  );
        row_vector_connectivity_check( aElementEdgeParentRanks );
        row_vector_connectivity_check( aElementFaceParentInds  );
        row_vector_connectivity_check( aElementFaceParentRanks );

        // copy all the information
        mSpatialDimension       = aSpatialDimension;
        mParentElementIndex     = aParentElementIndex;
        mNumElem                = aElementToNode.n_rows();
        mNodeInds               = aNodeInds.copy();
        mNodeParentInds         = aElementNodeParentInds.copy();
        mNodeParentRank         = aElementNodeParentRanks.copy();
        mElementToNode          = aElementToNode.copy();
        mElementEdgeParentInds  = aElementEdgeParentInds.copy();
        mElementEdgeParentRanks = aElementEdgeParentRanks.copy();
        mElementFaceParentInds  = aElementFaceParentInds.copy();
        mElementFaceParentRanks = aElementFaceParentRanks.copy();
        mElementInterfaceSides  = aElementInferfaceSides.copy();

        set_up_proc_local_to_cm_local_node_map();

        // construct connectivity
        moris::mtk::Cell_Info_Factory tConnFact;
        mConnectivity = tConnFact.create_cell_info_sp(mElementTopology);

}

Child_Mesh::Child_Mesh(Mesh_Modification_Template & aMeshModTemplate):
                    mElementTopology(CellTopology::TET4),
                    mConnectivity(nullptr),
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
                    mSubPhaseBinIndices(0),
                    mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
{

    reindex_template(aMeshModTemplate);

    MORIS_ASSERT(aMeshModTemplate.mIsReindexed,"Prior to inserting the template, reindexing must occur, in Child_Mesh from Child_Mesh_Template constructor");

    mSpatialDimension       = aMeshModTemplate.mSpatialDimension;
    mElementTopology        = aMeshModTemplate.mElementTopology;
    mParentElementIndex     = aMeshModTemplate.mParentElemInd;
    mNumElem                = aMeshModTemplate.mNumNewElem;
    mNodeInds               = aMeshModTemplate.mNodeInds.copy();
    mElementToNode          = aMeshModTemplate.mNewElementToNode.copy();
    mElementEdgeParentInds  = aMeshModTemplate.mNewParentEdgeOrdinals.copy();
    mElementEdgeParentRanks = aMeshModTemplate.mNewParentEdgeRanks.copy();
    mElementFaceParentInds  = aMeshModTemplate.mNewParentFaceOrdinals.copy();
    mElementFaceParentRanks = aMeshModTemplate.mNewParentFaceRanks.copy();
    mElementInterfaceSides  = aMeshModTemplate.mNewElementInterfaceSides.copy();

    // construct connectivity
    moris::mtk::Cell_Info_Factory tConnFact;
    mConnectivity = tConnFact.create_cell_info_sp(mElementTopology);

    set_up_proc_local_to_cm_local_node_map();

    if(mElementTopology == CellTopology::TET4)
    {
    	generate_connectivities(true,true,true);
    }
    else if(mElementTopology == CellTopology::TRI3)
    {
    	// 2D elements do not have faces
    	generate_connectivities(false,true,true);
    }
    else
    {
    	MORIS_ASSERT(0, "Only hex8 reg sub, tet4, quad4 reg sub, and tri3 are accepted presently.");
    }

}

Child_Mesh::~Child_Mesh()
{
}

moris::size_t
Child_Mesh::get_num_entities(enum EntityRank aEntityRank) const
{
    moris::size_t tNumEntities = 0;
    if(aEntityRank == EntityRank::NODE)
    {
        tNumEntities = mNodeInds.n_cols();
    }
    else if(aEntityRank == EntityRank::EDGE)
    {
        MORIS_ASSERT(mHasEdgeConn,"Without Edge connectivity it is unclear how many edges there are. Please call generate_connectivities w/ the edge flag on");
        tNumEntities = mEdgeToNode.n_rows();
    }

    else if(aEntityRank == EntityRank::FACE)
    {
        MORIS_ASSERT(mHasFaceConn,"Without face connectivity it is unclear how many faces there are. Please call generate_connectivities w/ the edge flag on");
        tNumEntities = mFaceToNode.n_rows();
    }
    else if(aEntityRank == EntityRank::ELEMENT)
    {
        tNumEntities = mElementToNode.n_rows();
    }

    else
    {
        MORIS_ASSERT(tNumEntities!=0,"Invalid entity rank specified");
    }

    return tNumEntities;
}
// ----------------------------------------------------------------------------------
enum EntityRank
Child_Mesh::get_facet_rank() const
{
    if(mSpatialDimension == 2)
    {
        return EntityRank::EDGE;
    }
    else if(mSpatialDimension == 3)
    {
        return EntityRank::FACE;
    }
    else
    {
        MORIS_ASSERT(0,"Facet rank only implemented in 2 and 3 spatial dimensions");
        return EntityRank::END_ENUM;
    }
}

// ----------------------------------------------------------------------------------
moris::mtk::Cell_Info const *
Child_Mesh::get_cell_info() const
{
    return mConnectivity.get();
}
// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_element_to_node() const
{
    return mElementToNode;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IdMat >
Child_Mesh::get_element_to_node_glob_ids(moris::moris_index aCMElemIndex) const
{
    moris::Matrix< moris::IdMat > tElementToNode = this->get_element_to_node().get_row(aCMElemIndex);

    tElementToNode = this->convert_proc_to_glob_ids(tElementToNode);

    return tElementToNode;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat >
Child_Mesh::get_node_to_element_local() const
{
    // Get the element to node locally indexed connectivity
    moris::Matrix< moris::IndexMat > tElementToNodeLocal = this->get_element_to_node_local();

    // Transpose the element to node local connectivity to get node to element
    // Allocate max size matrix (each node belongs to each element)
    moris::Matrix< moris::IndexMat > tNodeToElementLocal(this->get_num_entities(EntityRank::NODE),this->get_num_entities(EntityRank::ELEMENT));
    tNodeToElementLocal.fill(MORIS_INDEX_MAX);

    // keep track of number of nodes a node is connected to.
    moris::Matrix< moris::IndexMat > tNodetoElementCounter(this->get_num_entities(EntityRank::NODE),1,0);

    // Iterate through elements
    for(moris::uint iEl = 0; iEl<this->get_num_entities(EntityRank::ELEMENT); iEl++)
    {
        // iterate through nodes on the element
        for(moris::uint iN = 0; iN<tElementToNodeLocal.n_cols(); iN++)
        {
            // Get the node index
            moris::moris_index tNodeIndex = tElementToNodeLocal(iEl,iN);

            // add node to local node to element connectivity
            tNodeToElementLocal(tNodeIndex,tNodetoElementCounter(tNodeIndex)) = iEl;

            // Add to count
            tNodetoElementCounter(tNodeIndex)++;
        }
    }

    // Resize
    tNodeToElementLocal.resize(this->get_num_entities(EntityRank::NODE),tNodetoElementCounter.max());

    return tNodeToElementLocal;
}

// ----------------------------------------------------------------------------------

moris::moris_index
Child_Mesh::get_element_node_ordinal(moris::moris_index aCMLocElemIndex,
                         moris::moris_index aProcLocNodeIndex)
{
    for(moris::uint i = 0; i <mElementToNode.n_cols(); i++)
    {
        if(mElementToNode(aCMLocElemIndex,i) == aProcLocNodeIndex)
        {
            return i;
        }
    }

    return MORIS_INDEX_MAX;
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat >
Child_Mesh::get_element_to_node_local() const
{
    return convert_to_cm_local_indices(mElementToNode);
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IdMat >
Child_Mesh::get_element_to_node_global() const
{
    return convert_proc_to_glob_ids(mElementToNode);
}


// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > &
Child_Mesh::get_edge_to_node()
{
   MORIS_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
   return mEdgeToNode;
}

// ----------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_edge_to_node() const
{
   MORIS_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
   return mEdgeToNode;
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat >
Child_Mesh::get_edge_to_node_local() const
{
    MORIS_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
    return convert_to_cm_local_indices(mEdgeToNode);
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_element_to_edge() const
{
    MORIS_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
    return mElementToEdge;
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_edge_to_element() const
{
    MORIS_ASSERT(mHasEdgeConn,"Edge connectivity has not been generated with call to generate_edge_connectivity");
    return mEdgeToElement;
}

// ----------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_face_to_node() const
{
    MORIS_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
    return mFaceToNode;
}

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_facet_to_node() const
{
    if(mSpatialDimension == 3)
    {
        return this->get_face_to_node();
    }
    else if(mSpatialDimension == 2)
    {
        return this->get_edge_to_node();
    }
    else
    {
        MORIS_ERROR(0,"Not implemented");
        return mFaceToNode;
    }
}

// ---------------------------------------------------------------------------------

moris::Matrix<moris::IndexMat> const &
Child_Mesh::get_facet_to_element() const
{
    if(mSpatialDimension == 3)
    {
        return this->get_face_to_element();
    }
    else if(mSpatialDimension == 2)
    {
        return this->get_edge_to_element();
    }
    else
    {
        MORIS_ERROR(0,"Not implemented");
        return this->get_face_to_element();
    }
}


// ---------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_face_to_element() const
{
    MORIS_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
    return mFaceToElement;
}

// ---------------------------------------------------------------------------------


moris::Matrix< moris::IndexMat >
Child_Mesh::get_face_to_node_local() const
{
    MORIS_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
    return convert_to_cm_local_indices(mFaceToNode);
}

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_element_to_face() const
{
    MORIS_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
    return mElementToFace;
}

// ---------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat> const &
Child_Mesh::get_element_to_facet() const
{
    if(mSpatialDimension == 2)
    {
        return mElementToEdge;
    }
    else if(mSpatialDimension == 3)
    {
        return mElementToFace;
    }
    else
    {
        MORIS_ASSERT(0,"get_element_to_facet only implemented in 2 and 3 spatial dimensions");
        return mElementToFace;
    }
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_element_to_element() const
{
    MORIS_ASSERT(mHasElemToElem,"Element to element connectivity has not been generated with call to generate_element_to_element_connectivity");
    return mElementToElement;
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::DDRMat > const &
Child_Mesh::get_parametric_coordinates() const
{
    return mNodeParametricCoord;
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::DDRMat >
Child_Mesh::get_parametric_coordinates(moris::moris_index aNodeIndex) const
{

    // get child mesh local index
    auto tIter = mNodeIndsToCMInd.find(aNodeIndex);

    MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

    return mNodeParametricCoord.get_row(tIter->second);
}

// ---------------------------------------------------------------------------------

moris::mtk::Geometry_Type
Child_Mesh::get_child_geometry_type() const
{
    return mConnectivity->get_cell_geometry();
}
// ---------------------------------------------------------------------------------

moris::mtk::Interpolation_Order
Child_Mesh::get_child_interpolation_order() const
{
    return mConnectivity->get_cell_interpolation_order();
}
// ---------------------------------------------------------------------------------
moris::moris_index
Child_Mesh::get_entity_parent_entity_rank(enum EntityRank aEntityRank,
                                          moris::moris_index aCMLocIndex)
 {
     switch(aEntityRank)
     {
         case( EntityRank::NODE ):
         {
             return mNodeParentRank(aCMLocIndex);
             break;
         }
         case( EntityRank::EDGE ):
         {
             MORIS_ASSERT(mHasEdgeConn, "Attempting to access edge ancestry in a child mesh that did not generate edges (in get_entity_parent_entity_rank)");
             return mEdgeParentRanks(aCMLocIndex);
             break;
         }
         case( EntityRank::FACE ):
         {
             MORIS_ASSERT(mHasFaceConn, "Attempting to access face ancestry in a child mesh that did not generate face (in get_entity_parent_entity_rank)");
             return mFaceParentRanks(aCMLocIndex);
             break;
         }

         case( EntityRank::ELEMENT ):
         {
             return (moris::moris_index) EntityRank::ELEMENT;
             break;
         }

         default:
         {
             MORIS_ERROR(0,"Invalid entity rank specified in get_entity_parent_entity_rank");
             return MORIS_INDEX_MAX;
             break;
         }

     }
 }

// ---------------------------------------------------------------------------------
moris::moris_index
Child_Mesh::get_entity_parent_entity_proc_ind(enum EntityRank  aEntityRank,
                                              moris::moris_index aCMLocIndex)
{
    switch(aEntityRank)
    {
        case( EntityRank::NODE ):
        {
            return mNodeParentInds(aCMLocIndex);
            break;
        }
        case( EntityRank::EDGE ):
        {
            MORIS_ASSERT(mHasEdgeConn, "Attempting to access edge ancestry in a child mesh that did not generate edges (in get_entity_parent_entity_rank)");
            return mEdgeParentInds(aCMLocIndex);
            break;
        }
        case( EntityRank::FACE ):
        {
            MORIS_ASSERT(mHasFaceConn, "Attempting to access face ancestry in a child mesh that did not generate face (in get_entity_parent_entity_rank)");
            return mFaceParentInds(aCMLocIndex);
            break;
        }

        case( EntityRank::ELEMENT ):
        {
            return mParentElementIndex;
            break;
        }

        default:
        {
            MORIS_ERROR(0,"Invalid entity rank specified in get_entity_parent_entity_rank");
            return MORIS_INDEX_MAX;
            break;
        }

    }
}

// ---------------------------------------------------------------------------------

moris::moris_index
Child_Mesh::get_parent_element_index() const
{
    return mParentElementIndex;
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_face_parent_inds() const
{
    return mFaceParentInds;
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::DDSTMat > const &
Child_Mesh::get_face_parent_ranks() const
{
    return mFaceParentRanks;
}


moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_facet_parent_inds() const
{
    MORIS_ASSERT(mSpatialDimension == 2 || mSpatialDimension == 3,"get facet parent ranks only supported in 2 or 3 dimension");
    if(mSpatialDimension == 2)
    {
        return mEdgeParentInds;
    }
    else
    {
        return mFaceParentInds;
    }
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::DDSTMat > const &
Child_Mesh::get_facet_parent_ranks() const
{
    MORIS_ASSERT(mSpatialDimension == 2 || mSpatialDimension == 3,"get facet parent ranks only supported in 2 or 3 dimension");
    if(mSpatialDimension == 2)
    {
        return mEdgeParentRanks;
    }
    else
    {
        return mFaceParentRanks;
    }

}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_edge_parent_inds() const
{
    return mEdgeParentInds;
}


// ---------------------------------------------------------------------------------

moris::Matrix< moris::DDSTMat > const &
Child_Mesh::get_edge_parent_ranks() const
{
    return mEdgeParentRanks;
}

// ---------------------------------------------------------------------------------
moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_node_indices() const
{
    return mNodeInds;
}
// ---------------------------------------------------------------------------------

moris::Matrix< moris::IdMat > const &
Child_Mesh::get_node_ids() const
{
    return mNodeIds;
}
moris::Cell<moris::mtk::Vertex const *> const &
Child_Mesh::get_vertices() const
{
    return mVertices;
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::IdMat > const &
Child_Mesh::get_element_ids() const
{
    return mChildElementIds;
}


// ---------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat > const &
Child_Mesh::get_element_inds() const
{
    return mChildElementInds;
}

// ---------------------------------------------------------------------------------
moris::Matrix<moris::IndexMat>
Child_Mesh::convert_to_proc_local_elem_inds(moris::Matrix< moris::IndexMat > aConnOfElementCMIndices)
{
    moris::Matrix<moris::IndexMat> tConnWithElementsProcIndices(aConnOfElementCMIndices.n_rows(),aConnOfElementCMIndices.n_cols());

    for(moris::uint i = 0; i <aConnOfElementCMIndices.n_rows(); i++)
    {
        for(moris::uint j =0; j <aConnOfElementCMIndices.n_cols(); j++)
        {
            tConnWithElementsProcIndices(i,j) = mChildElementInds(aConnOfElementCMIndices(i,j));
        }
    }

    return tConnWithElementsProcIndices;
}

// ---------------------------------------------------------------------------------

moris::Matrix< moris::IndexMat >
Child_Mesh::get_edge_ordinal_from_element_and_edge_indices(moris::moris_index const & aElementIndex,
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

    MORIS_ASSERT(tNumOrdinalstoFind == tCount,"All edge ordinals not found");

    return tEdgeOrdinals;

}

// ---------------------------------------------------------------------------------


moris::moris_index
Child_Mesh::get_edge_ordinal_from_element_and_edge_indices(moris::moris_index const & aElementIndex,
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

    MORIS_ASSERT(1 == tCount,"All edge ordinals not found");

    return tEdgeOrdinals;

}

// ---------------------------------------------------------------------------------

moris::moris_index
Child_Mesh::get_face_ordinal_from_element_and_face_index(moris::moris_index const & aElementIndex,
                                                         moris::moris_index         aFacet) const
{
    // get the elemnt to edge connectivity
    moris::Matrix< moris::IndexMat > const & tElemToFacet = get_element_to_facet();

    moris::size_t tNumFacets = tElemToFacet.n_cols();

    // Initialize output
    moris::moris_index tFaceOrdinal = 1000;

    bool tSuccess = false;

    for(moris::size_t iFacet = 0; iFacet < tNumFacets; iFacet++)
    {

        if(aFacet == tElemToFacet(aElementIndex,iFacet))
        {
            tFaceOrdinal = iFacet;
            tSuccess = true;
            continue;
        }
    }
    MORIS_ERROR(tSuccess,"Face ordinal not found");
    return tFaceOrdinal;
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::get_child_elements_connected_to_parent_facet(moris::moris_index         const & aParentFaceIndex,
                                                         moris::Matrix< moris::IdMat >    & aChildElemsIdsOnFacet,
                                                         moris::Matrix< moris::IndexMat > & aChildElemsCMIndOnFacet,
                                                         moris::Matrix< moris::IndexMat > & aChildElemOnFacetOrdinal) const
{
    // First determine which of this meshes face live on the parent face
    enum EntityRank tFacetRank = this->get_facet_rank();

    moris::size_t tNumFacets = get_num_entities(tFacetRank);
    moris::size_t tNumElemsOnFace = 0;
    moris::Matrix< moris::IndexMat > tLocFacesOnParentFace(1,tNumFacets);

    moris::Matrix< moris::IndexMat > const & tFacetParentInds = this->get_facet_parent_inds();
    moris::Matrix< moris::DDSTMat  > const & tFacetParentRanks = this->get_facet_parent_ranks();
    moris::Matrix< moris::IndexMat > const & tFacetToElement = this->get_facet_to_element();

    // Iterate through face ancestry ranks
    for(moris::size_t i = 0; i < tNumFacets; i++)
    {
        if(tFacetParentInds(0,i) == aParentFaceIndex && tFacetParentRanks(0,i) == (uint)(tFacetRank))
        {
            tLocFacesOnParentFace(0,tNumElemsOnFace) = i;
            tNumElemsOnFace++;
        }
    }

    // Iterate through face to element connectivities of faces connected to the parent face
    aChildElemsIdsOnFacet    = moris::Matrix< moris::IdMat >(1,tNumElemsOnFace*2);
    aChildElemsCMIndOnFacet  = moris::Matrix< moris::IndexMat >(1,tNumElemsOnFace*2);
    aChildElemOnFacetOrdinal = moris::Matrix< moris::IndexMat >(1,tNumElemsOnFace*2);
    moris::size_t tCount = 0;

    for(moris::size_t i = 0; i<tNumElemsOnFace; i++)
    {
        moris::size_t tFacetInd = tLocFacesOnParentFace(0,i);
        for(moris::size_t j = 0; j<tFacetToElement.n_cols(); j++)
        {
            if(tFacetToElement(tFacetInd,j) != std::numeric_limits<moris::moris_index>::max())
            {
                aChildElemsCMIndOnFacet(0,tCount)  = tFacetToElement(tFacetInd,j);
                aChildElemsIdsOnFacet(0,tCount)    = mChildElementIds(0,tFacetToElement(tFacetInd,j));
                aChildElemOnFacetOrdinal(0,tCount) = this->get_face_ordinal_from_element_and_face_index(tFacetToElement(tFacetInd,j),tFacetInd);
                tCount++;
            }
        }
    }
    aChildElemsCMIndOnFacet.resize(1,tCount);
    aChildElemsIdsOnFacet.resize(1,tCount);
    aChildElemOnFacetOrdinal.resize(1,tCount);
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::modify_child_mesh(enum TemplateType aTemplate)
{
    modify_child_mesh_internal(aTemplate);
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::initialize_unzipping()
{
    MORIS_ASSERT(mUnzippingFlag == false,"Error in child mesh, unzipping has already been initialized");
    mUnzippingFlag = true;
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::finalize_unzipping()
{
    MORIS_ASSERT(mUnzippingFlag == true,"Error in child mesh, finalize_unzipping called on a child mesh which is not unzipping");
    mUnzippingFlag = false;
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::unzip_child_mesh_interface_get_interface_element_pairs(moris::uint aGeometryIndex,
                                                       bool & aNoPairFoundFlag,
                                                       moris::Matrix<moris::IndexMat> & aInterfaceElementPairs,
                                                       moris::Matrix<moris::IndexMat> & aInterfacePairSideOrds) const
{
    aNoPairFoundFlag = false;
    moris::uint tNumElements = this->get_num_entities(EntityRank::ELEMENT);

    // Figure out which faces are interface faces
    moris::Cell<moris::moris_index> tInterfaceFaceIndices;

    moris::Matrix<moris::IndexMat> const & tElementToFacet = this->get_element_to_facet();

    for(moris::uint iEl = 0; iEl<tNumElements; iEl++)
    {
        // Get the side ordinal from our element interface side matrix
        moris::size_t tInterfaceSideOrdinal = mElementInterfaceSides(iEl,aGeometryIndex);
        if(tInterfaceSideOrdinal!=std::numeric_limits<moris::size_t>::max())
        {
            // Add the child mesh local face index to the interface indice cell
            tInterfaceFaceIndices.push_back(tElementToFacet(iEl,tInterfaceSideOrdinal));
        }
    }

    // Remove duplicates
    unique(tInterfaceFaceIndices);

    // Number of interface nodes
    moris::uint tNumInterfaceFaces = tInterfaceFaceIndices.size();

    // Construct interface pairs
    moris::uint tNumElementsPerFace = 2;
    aInterfaceElementPairs.resize(tNumElementsPerFace,tInterfaceFaceIndices.size());
    aInterfacePairSideOrds.resize(tNumElementsPerFace,tInterfaceFaceIndices.size());

    // Face to element matrix
    moris::Matrix<moris::IndexMat> const & tFacetToElem = this->get_facet_to_element();

    for( moris::uint iF = 0; iF< tNumInterfaceFaces; iF++)
    {
        moris::uint tInterfaceFaceIndex = tInterfaceFaceIndices(iF);

        // iterate through elements attached to face
        for(moris::uint iEl = 0; iEl<tNumElementsPerFace; iEl++)
        {
            moris::uint tCMElementIndex = tFacetToElem(tInterfaceFaceIndex,iEl);

            if( tCMElementIndex != MORIS_INDEX_MAX )
            {
                aInterfaceElementPairs(iEl,iF)  = tCMElementIndex;
                aInterfacePairSideOrds(iEl,iF)  = mElementInterfaceSides(tCMElementIndex,aGeometryIndex);
            }
            else
            {
                aNoPairFoundFlag = true;
                aInterfaceElementPairs(iEl,iF) = MORIS_INDEX_MAX;
                aInterfacePairSideOrds(iEl,iF) = MORIS_INDEX_MAX;
            }
        }
    }

}

// ---------------------------------------------------------------------------------

void
Child_Mesh::unzip_child_mesh_interface(moris::moris_index                       aGeometryIndex,
                           moris::Matrix< moris::IndexMat > const & aInterfaceElementPairsCMIndex,
                           moris::Matrix< moris::IndexMat > const & aElementUsingZippedNodes,
                           moris::Matrix< moris::IndexMat > const & aInterfaceNodeIndices,
                           moris::Matrix< moris::IndexMat > const & aUnzippedInterfaceNodeIndices,
                           moris::Matrix< moris::IndexMat > const & aUnzippedInterfaceNodeIds)
{
    MORIS_ASSERT(mUnzippingFlag == true,"Error in child mesh, unzip_child_mesh_interface called on a child mesh which is not unzipping. Please call initialize_unzipping.");
    MORIS_ASSERT(aGeometryIndex == 0 ,"unzip_child_mesh_interface_get_interface_element_pairs not implemented for multiple geometries because mElementInterfaceFaces needs to accomodate this");

    this->add_node_indices(aUnzippedInterfaceNodeIndices);
    this->add_node_ids(aUnzippedInterfaceNodeIds);

    // Mark connectivity as being inconsistent since the nodes have changed
    mHasFaceConn   = false;
    mHasEdgeConn   = false;
    mHasElemToElem = false;

    // Element face to node map for tet4
    moris::Matrix< moris::IndexMat > tNodeToFaceMap = mConnectivity->get_node_to_face_map();

    // Copy starting element to node connectivity
    moris::Matrix< moris::IndexMat > tOldElementToNode = mElementToNode.copy();

    // Get the node to element connectivity
    moris::Matrix< moris::IndexMat > tNodeToElement = get_node_to_element_local();

    // Number of interface nodes
    moris::uint tNumInterfaceNodes = aInterfaceNodeIndices.numel();

    // Construct map between interface node index and index in the aInterfaceNodeIndices
    // Prefer not to use map but cannot see how without repeatedly finding values in aInterfaceNodeIndices
    std::unordered_map<moris::size_t, moris::size_t> tInterfaceNodesMap;
    for(moris::uint iN = 0; iN<tNumInterfaceNodes; iN++)
    {
        MORIS_ASSERT(tInterfaceNodesMap.find(aInterfaceNodeIndices(iN))==tInterfaceNodesMap.end()," interface node index already in the map");

        tInterfaceNodesMap[aInterfaceNodeIndices(iN)] = iN;
    }

    moris::Matrix<moris::IndexMat> tElementsConnectedToNodeTraversed(tNumInterfaceNodes,1,0);

    // Verify number of interface pairs and get the number we have
    moris::uint tNumInterfacePairs = aInterfaceElementPairsCMIndex.n_cols();
    MORIS_ASSERT(aElementUsingZippedNodes.numel() == aInterfaceElementPairsCMIndex.n_cols()," Dimension mismatch between interface element pair matrix and element using unzipped nodes" );


    for(moris::uint iP =0; iP<tNumInterfacePairs; iP++)
    {
        // Cm mesh element index for element which needs to change nodes
        moris::moris_index tElementUsingUnzippedNodes = aInterfaceElementPairsCMIndex(aElementUsingZippedNodes(iP),iP);

        // Side ordinal on interface
        moris::size_t tElementInterfaceSideOrdinal = mElementInterfaceSides(tElementUsingUnzippedNodes,aGeometryIndex);

        MORIS_ASSERT(tElementInterfaceSideOrdinal!=std::numeric_limits<moris::size_t>::max(),"Invalid interface side ordinal found");

        // change the element to node connectivity
        for(moris::moris_index iN =0; iN<(moris::moris_index)tNodeToFaceMap.n_cols(); iN++)
        {
            // Node ordinal relative to element's node connectivity
            moris::moris_index tNodeOrd   = tNodeToFaceMap(tElementInterfaceSideOrdinal,iN);

            // Get the node index using the above ordinal
            moris::moris_index tNodeIndex = tOldElementToNode(tElementUsingUnzippedNodes,tNodeOrd);

            // Find the index in the interface nods
            auto tInterfaceOrd = tInterfaceNodesMap.find(tNodeIndex);

            // Check if its actually there
            MORIS_ASSERT(tInterfaceOrd != tInterfaceNodesMap.end(),"Node not in interface node map, conversion to local interface node ordinal failed");

            // Use the index to get the new node index
            moris::moris_index tUnzippedNodeIndex = aUnzippedInterfaceNodeIndices(tInterfaceOrd->second);

            // Set new node index in the element to node connectivity
            mElementToNode(tElementUsingUnzippedNodes,tNodeOrd) = tUnzippedNodeIndex;

            // traverse over all elements connected to this interface node if it has not been done yet
            if(tElementsConnectedToNodeTraversed(tInterfaceOrd->second) == 0)
            {
                // My phase index
                moris::size_t tUnzippedPhase = get_element_phase_index(tElementUsingUnzippedNodes);

                // the child mesh local node index
                moris::moris_index tNodeIndexCM = get_cm_local_node_index(tNodeIndex);

                // iterate through node to element (NTE)
                for(moris::uint iNTE = 0; iNTE<tNodeToElement.n_cols(); iNTE++)
                    {
                    moris::moris_index tElementInd = tNodeToElement(tNodeIndexCM,iNTE);

                    if(tElementInd != MORIS_INDEX_MAX)
                    {
                        if(get_element_phase_index(tElementInd) == tUnzippedPhase)
                        {
                            moris::moris_index tNodeOrd = get_element_node_ordinal(tElementInd,tNodeIndex);

                            if(tNodeOrd!=MORIS_INDEX_MAX)
                            {
                                mElementToNode(tElementInd,tNodeOrd) = tUnzippedNodeIndex;
                            }
                        }

                    }

                    else
                    {
                        break;
                    }

                }

                tElementsConnectedToNodeTraversed(tInterfaceOrd->second) = 1;
            }

        }

    }

    // regenerate the connectivities with the changes
    generate_connectivities(true,true,true);


    // allocate space
    allocate_parametric_coordinates( tNumInterfaceNodes, 3 );

    // iterate over unzipped nodes and copy the parametric coordinates of the original node
    for(moris::uint iInt =0; iInt<tNumInterfaceNodes; iInt++)
    {
        add_node_parametric_coordinate(aUnzippedInterfaceNodeIndices(iInt),
                                       this->get_parametric_coordinates(aInterfaceNodeIndices(iInt)));
    }


}

void
Child_Mesh::convert_tet4_to_tet10_child()
{
    std::cerr<<"Warning conversion tet4 to tet10 not implemented in new child mesh structure"<<std::endl;
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::add_node_indices(moris::Matrix< moris::IndexMat > const & aNewNodeInds)
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

    // Add node inheritance space for new nodes
    mNodeParentInds.resize(1,tNumNewNodes + tNumCurrNodes);
    mNodeParentRank.resize(1,tNumNewNodes + tNumCurrNodes);

}

// ---------------------------------------------------------------------------------

void
Child_Mesh::add_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds)
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

void
Child_Mesh::add_vertices(moris::Cell<moris::mtk::Vertex const *> const & aVertices)
{
    mVertices.append(aVertices);
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::set_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds)
{
    MORIS_ASSERT(aNewNodeIds.n_cols() == get_num_entities(EntityRank::NODE),"Number of node ids does not match the number of nodes in this child mesh");
    mNodeIds = aNewNodeIds.copy();
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::set_child_element_ids(moris::moris_id & aElementId)
{
    MORIS_ASSERT(mChildElementIds.n_cols() == 0, "Element Ids already set");
    moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
    mChildElementIds = moris::Matrix< moris::IdMat >(1,tNumElements);

    for(moris::size_t iElem = 0; iElem<tNumElements; iElem++)
    {
        mChildElementIds(0,iElem) = aElementId;
        aElementId++;
    }
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::set_child_element_inds(moris::moris_index & aElementInd)
{
    MORIS_ASSERT(mChildElementInds.n_cols() == 0, "Element Inds already set");


    moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
    mChildElementInds = moris::Matrix< moris::IndexMat >(1,tNumElements);

    for(moris::size_t iElem = 0; iElem<tNumElements; iElem++)
    {
            mChildElementInds(0,iElem) = aElementInd;
            aElementInd++;
    }
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::add_node_parametric_coordinate( moris::size_t aNodeIndex,
                                            moris::Matrix< moris::DDRMat > const & aParamCoord )
{
    // get child mesh local index
    auto tIter = mNodeIndsToCMInd.find(aNodeIndex);

    MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

    moris::size_t tNodeCMIndex = tIter->second;
    mNodeParametricCoord.set_row(tNodeCMIndex,aParamCoord);

}

// ---------------------------------------------------------------------------------

void
Child_Mesh::add_node_parametric_coordinate( moris::Matrix< moris::IndexMat> const & aNodeIndices,
                                            moris::Matrix< moris::DDRMat >  const & aParamCoord )
{

    moris::size_t tNumNodes = aNodeIndices.numel();
    // Only do the checks if there are nodes to be added (this is only an issue for a child mesh
    // intersected by multiple geometries.
    if(tNumNodes !=0 )
    {
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

// ---------------------------------------------------------------------------------

void
Child_Mesh::allocate_parametric_coordinates( moris::size_t aNumNewNodes,
                                             moris::size_t aDimOfParmCoord)
{
    moris::size_t tCurrentRow = mNodeParametricCoord.n_rows();
    mNodeParametricCoord.resize(tCurrentRow+aNumNewNodes,aDimOfParmCoord);
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::allocate_more_elements(moris::size_t const & aNumMoreElements)
{
    moris::size_t tCurrNumElem = mNumElem;
    moris::size_t tNewNumElem  = tCurrNumElem + aNumMoreElements;

    mElementToNode.resize(tNewNumElem,mElementToNode.n_cols());
    mElementEdgeParentInds.resize(tNewNumElem,mElementEdgeParentInds.n_cols());
    mElementEdgeParentRanks.resize(tNewNumElem,mElementEdgeParentRanks.n_cols());
    mElementFaceParentInds.resize(tNewNumElem,mElementFaceParentInds.n_cols());
    mElementFaceParentRanks.resize(tNewNumElem,mElementFaceParentRanks.n_cols());
    mElementInterfaceSides.resize(tNewNumElem,mElementInterfaceSides.n_cols());
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::insert_child_mesh_template(Mesh_Modification_Template & aMeshModTemplate)
  {
      // Reindex the child mehs template
      reindex_template(aMeshModTemplate);

      mElementTopology = aMeshModTemplate.mElementTopology;

      // construct connectivity
      moris::mtk::Cell_Info_Factory tConnFact;
      mConnectivity = tConnFact.create_cell_info_sp(mElementTopology);


      // Assert that reindex template actually happened
      MORIS_ASSERT(aMeshModTemplate.mIsReindexed,"Prior to inserting the template, reindexing must occur, in insert_child_mesh_template.");

      // index of current geometry
      moris_index tGeomIndex = *(mGeometryIndex.end() - 1);


      // Replace the element which is marked as replaceable (if in the future multiple elements are replaced at one time then this could be a for loop)
      // But for now only one element is replaced

      replace_element(aMeshModTemplate.mElemIndToReplace,
                      0,
                      tGeomIndex,
                      aMeshModTemplate.mNewElementToNode,
                      aMeshModTemplate.mNewParentEdgeOrdinals,
                      aMeshModTemplate.mNewParentEdgeRanks,
                      aMeshModTemplate.mNewParentFaceOrdinals,
                      aMeshModTemplate.mNewParentFaceRanks,
                      aMeshModTemplate.mNewElementInterfaceSides);

      // Number of elements that are new
      moris::size_t tNumElemToAdd = aMeshModTemplate.mNumNewElem-aMeshModTemplate.mNumElemToReplace;

      // Iterate over new elements to add and add them
      for(moris::size_t i = 0; i<tNumElemToAdd; i++)
      {
          add_element(i+aMeshModTemplate.mNumElemToReplace,
                      tGeomIndex,
                      aMeshModTemplate.mNewElementToNode,
                      aMeshModTemplate.mNewParentEdgeOrdinals,
                      aMeshModTemplate.mNewParentEdgeRanks,
                      aMeshModTemplate.mNewParentFaceOrdinals,
                      aMeshModTemplate.mNewParentFaceRanks,
                      aMeshModTemplate.mNewElementInterfaceSides);
      }

      // Add node inheritance from the child mesh template
      if(aMeshModTemplate.mHasNodeInheritance)
      {
          for(moris::size_t i = 0; i < aMeshModTemplate.mNodeInds.numel(); i++)
          {
              moris::size_t tLocalProcNodeInd = aMeshModTemplate.mNodeInds(i);
              moris::size_t tNodeParentInds   = aMeshModTemplate.mNewNodeParentOrdinals(i);
              moris::size_t tNodeParentRank   = aMeshModTemplate.mNewNodeParentRanks(i);
              moris::size_t tLocalCMNodeInd   = this->get_cm_local_node_index(tLocalProcNodeInd);

              set_node_inheritance(tLocalCMNodeInd,tNodeParentInds,tNodeParentRank);
          }
      }
  }

// ---------------------------------------------------------------------------------

void
Child_Mesh::reindex_template(Mesh_Modification_Template & aMeshModTemplate)
{
    MORIS_ASSERT(!aMeshModTemplate.mIsReindexed,"reindex_template() called on a child mesh modification template which has already been reindexed");

    aMeshModTemplate.mNewElementToNode = reindex_matrix(aMeshModTemplate.mNewElementToNode, 0, aMeshModTemplate.mNodeInds);

    // Reindex template edge and face ordinals
    reindex_template_parent_information(aMeshModTemplate);
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::replace_element(moris::size_t                    const & aElementIndexToReplace,
                            moris::size_t                    const & aRowIndex,
                            moris::size_t                    const & aGeomIndex,
                            moris::Matrix< moris::IndexMat > const & aElementToNode,
                            moris::Matrix< moris::IndexMat > const & aElementEdgeParentInds,
                            moris::Matrix< moris::DDSTMat >  const & aElementEdgeParentRanks,
                            moris::Matrix< moris::IndexMat > const & aElementFaceParentInds,
                            moris::Matrix< moris::DDSTMat >  const & aElementFaceParentRanks,
                            moris::Matrix< moris::DDSTMat >  const & aElementInterfaceFaces)
{
    mHasFaceConn   = false;
    mHasEdgeConn   = false;
    mHasElemToElem = false;
    replace_row(aRowIndex, aElementToNode         , aElementIndexToReplace, mElementToNode);
    replace_row(aRowIndex, aElementEdgeParentInds , aElementIndexToReplace, mElementEdgeParentInds);
    replace_row(aRowIndex, aElementEdgeParentRanks, aElementIndexToReplace, mElementEdgeParentRanks);
    replace_row(aRowIndex, aElementFaceParentInds , aElementIndexToReplace, mElementFaceParentInds);
    replace_row(aRowIndex, aElementFaceParentRanks, aElementIndexToReplace, mElementFaceParentRanks);
    mElementInterfaceSides(aElementIndexToReplace,aGeomIndex) = aElementInterfaceFaces(aRowIndex);
}


// ---------------------------------------------------------------------------------

void
Child_Mesh::add_element(moris::size_t                    const & aRowIndex,
                        moris::size_t                    const & aGeomIndex,
                        moris::Matrix< moris::IndexMat > const & aElementToNode,
                        moris::Matrix< moris::IndexMat > const & aElementEdgeParentInds,
                        moris::Matrix< moris::DDSTMat >  const & aElementEdgeParentRanks,
                        moris::Matrix< moris::IndexMat > const & aElementFaceParentInds,
                        moris::Matrix< moris::DDSTMat >  const & aElementFaceParentRanks,
                        moris::Matrix< moris::DDSTMat>   const & aElementInterfaceFaces)
{
    MORIS_ASSERT(mNumElem<mElementToNode.n_rows(),"Not enough space allocated in call to allocate_more_elements");
    mHasFaceConn = false;
    mHasEdgeConn = false;
    mHasElemToElem = false;

    replace_row(aRowIndex, aElementToNode         , mNumElem, mElementToNode);
    replace_row(aRowIndex, aElementEdgeParentInds , mNumElem, mElementEdgeParentInds);
    replace_row(aRowIndex, aElementEdgeParentRanks, mNumElem, mElementEdgeParentRanks);
    replace_row(aRowIndex, aElementFaceParentInds , mNumElem, mElementFaceParentInds);
    replace_row(aRowIndex, aElementFaceParentRanks, mNumElem, mElementFaceParentRanks);

    mElementInterfaceSides(mNumElem,aGeomIndex) = aElementInterfaceFaces(aRowIndex);

    mNumElem++;
}
// ---------------------------------------------------------------------------------

void
Child_Mesh::set_node_inheritance(moris::moris_index aCMLocNodeIndex,
                                 moris::moris_index aProcLocParentEntityInd,
                                 moris::moris_index aParentEntityRank)
{
    MORIS_ASSERT(aCMLocNodeIndex<(moris::moris_index)this->get_num_entities(EntityRank::NODE),
                 "Trying to set node inheritance on a node that doesn't exist in child mesh");

    mNodeParentRank(aCMLocNodeIndex) = aParentEntityRank;
    mNodeParentInds(aCMLocNodeIndex) = aProcLocParentEntityInd;
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::generate_connectivities(bool aGenerateFaceConn,
                                    bool aGenerateEdgeConn,
                                    bool aGenerateElemToElem)
{

    moris::mtk::Cell_Info_Factory tConnFact;
    mConnectivity = tConnFact.create_cell_info_sp(mElementTopology);

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

// ---------------------------------------------------------------------------------

void
Child_Mesh::add_new_geometry_interface(moris_index aGeomIndex)
{

    mGeometryIndex.push_back(aGeomIndex);
    mElementInterfaceSides.resize(this->get_num_entities(EntityRank::ELEMENT),mGeometryIndex.size());
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::add_entity_to_intersect_connectivity(moris::moris_index aCMNodeInd,
                                                 moris::moris_index aCMEdgeInd,
                                                 moris::size_t aFlag)
{
    if(aFlag == 0)
    {
        aCMNodeInd = aCMNodeInd + get_num_entities(EntityRank::NODE);
    }

    // add node information
    mIntersectedCMNodeIndex.push_back(aCMNodeInd);
    mIntersectedEdges.push_back(aCMEdgeInd);

}

// ---------------------------------------------------------------------------------

void
Child_Mesh::mark_edge_as_on_interface(moris::moris_index aEdgeIndex)
{
    MORIS_ASSERT(aEdgeIndex<(moris::moris_index)get_num_entities(EntityRank::EDGE),"Edge index out of bounds");
    MORIS_ASSERT(mEdgeOnInterface.numel() == get_num_entities(EntityRank::EDGE),"mEdgeOnInterface has not been properly initialized");
    mEdgeOnInterface(aEdgeIndex) = 1;
    mHasCoincidentEdges = true;
}

// ---------------------------------------------------------------------------------


void
Child_Mesh::mark_interface_faces_from_interface_coincident_faces()
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
        moris::Matrix<moris::IndexMat> const tEdgeOrdinalToFaceOrdMap = mConnectivity->get_edge_to_face_map();


        moris::Matrix<moris::DDUMat>       tFaceCounter(1,tNumFacesPerElem,0);
        moris_index tGeomIndex = mGeometryIndex.size()-1;


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
                        mElementInterfaceSides(iElem,tGeomIndex) = tFaceOrd;
                    }
                }



            }
            tFaceCounter.fill(0);
        }
    }

    mHasCoincidentEdges = false;
    mEdgeOnInterface.resize(0,0);
}

// ---------------------------------------------------------------------------------

void
Child_Mesh::initialize_element_phase_mat()
{
    mElementPhaseIndices = moris::Matrix< moris::IndexMat >(1,get_num_entities(EntityRank::ELEMENT),std::numeric_limits<moris::moris_index>::max());
}


// ---------------------------------------------------------------------------------

void
Child_Mesh::set_element_phase_index(moris::size_t aEntityIndex,
                                    moris::size_t aEntityPhaseIndex)
{
    MORIS_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
    mHasPhaseInfo = true;
    mElementPhaseIndices(0,aEntityIndex) = aEntityPhaseIndex;
}

// ---------------------------------------------------------------------------------

moris::size_t
Child_Mesh::get_element_phase_index( moris::size_t const & aEntityIndex) const
{
    MORIS_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
    MORIS_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
    return mElementPhaseIndices(0,aEntityIndex);
}

// ---------------------------------------------------------------------------------
moris_index
Child_Mesh::get_element_subphase_index(moris::size_t const & aEntityIndex) const
{
    MORIS_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
    MORIS_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
    return mElementBinIndex(aEntityIndex);
}
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------


}

