/*
 * cl_XTK_Child_Mesh.cpp
 *
 *  Created on: Mar 14, 2019
 *      Author: doble
 */

#include "cl_XTK_Child_Mesh.hpp"
#include "assert.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_TOL_Memory_Map.hpp"
#include "cl_MTK_Mesh_Core.hpp"

using namespace moris;
namespace xtk
{
    Child_Mesh::Child_Mesh()
    : mElementTopology(CellTopology::TET4),
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
      mHasInterChildMeshInterface(false),
      mHasPhaseInfo(false),
      mElementPhaseIndices(0,0),
      mElementBinIndex(0,0),
      mBinBulkPhase(0),
      mSubPhaseBinIndices(0),
      mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
    {}

    // ----------------------------------------------------------------------------------

    Child_Mesh::Child_Mesh(
            moris::uint                        aSpatialDimension,
            moris::moris_index                 aParentElementIndex,
            moris::Matrix< moris::IndexMat > & aNodeInds,
            moris::Matrix< moris::IndexMat > & aElementNodeParentInds,
            moris::Matrix< moris::DDSTMat >  & aElementNodeParentRanks,
            moris::Matrix< moris::IndexMat > & aElementToNode,
            moris::Matrix< moris::IndexMat > & aElementEdgeParentInds,
            moris::Matrix< moris::DDSTMat >  & aElementEdgeParentRanks,
            moris::Matrix< moris::IndexMat > & aElementFaceParentInds,
            moris::Matrix< moris::DDSTMat >  & aElementFaceParentRanks,
            moris::Matrix< moris::DDSTMat >  & aElementInferfaceSides )
    : mElementTopology(CellTopology::TET4),
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
      mHasInterChildMeshInterface(false),
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

    Child_Mesh::Child_Mesh(
            moris_index        aChildMeshIndex,
            moris::uint        aSpatialDimension,
            moris::mtk::Cell*  aParentCell,
            moris::mtk::Mesh*  aBackgroundMesh,
            Cell<moris_index> const & aGeometryIndices)
    : mElementTopology(CellTopology::TET4),
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
      mHasInterChildMeshInterface(false),
      mHasPhaseInfo(false),
      mElementPhaseIndices(0,0),
      mElementBinIndex(0,0),
      mBinBulkPhase(0),
      mSubPhaseBinIndices(0),
      mSubPhaseBins(0,moris::Matrix< moris::IndexMat >(0,0))
    {
        mChildMeshIndex = aChildMeshIndex;
        // Check for row vector connectivity (if not it is transposed)
        Matrix<IndexMat> tVertexIndices = aParentCell->get_vertex_inds();
        row_vector_connectivity_check( tVertexIndices );

        Matrix<IndexMat> tVertexIds = aParentCell->get_vertex_ids();
        row_vector_connectivity_check( tVertexIndices );
        
        // get the number of cells
        mtk::Cell_Info  const * tParentCellInfo = aParentCell->get_cell_info(); 


        Matrix< IndexMat > tElementToEdge = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(aParentCell->get_index(), moris::EntityRank::ELEMENT, moris::EntityRank::EDGE);
        row_vector_connectivity_check( tElementToEdge  );

        Matrix< IndexMat > tFacetoElemConnInd = aBackgroundMesh->get_entity_connected_to_entity_loc_inds(
                        aParentCell->get_index(),
                        moris::EntityRank::ELEMENT,
                        moris::EntityRank::FACE);
        row_vector_connectivity_check( tFacetoElemConnInd );

        // copy all the information
        mSpatialDimension       = aSpatialDimension;
        mParentElementIndex     = aParentCell->get_index();
        mNumElem                = 1;
        mNodeInds               = tVertexIndices;
        mNodeIds                = tVertexIds;
        mNodeParentInds         = tVertexIndices;
        mNodeParentRank         = moris::Matrix< moris::DDSTMat >(1,tVertexIndices.numel());
        mElementToNode          = tVertexIndices;
        mElementEdgeParentInds  = tElementToEdge;
        mElementEdgeParentRanks = moris::Matrix< moris::DDSTMat >(1,tFacetoElemConnInd.numel(),1);
        mElementFaceParentInds  = tFacetoElemConnInd;
        mElementFaceParentRanks = moris::Matrix< moris::DDSTMat >(1,tFacetoElemConnInd.numel(),2);
        mElementInterfaceSides  = moris::Matrix< moris::DDSTMat > (1,aGeometryIndices.size(),std::numeric_limits<moris::size_t>::max() );
        
        tParentCellInfo->get_loc_coords_of_cell(mNodeParametricCoord);
        
        mGeometryIndex = aGeometryIndices;

        // add nodes to the
        set_up_proc_local_to_cm_local_node_map();

        // construct connectivity
        moris::mtk::Cell_Info_Factory tConnFact;
        mConnectivity = tConnFact.create_cell_info_sp(mElementTopology);
    }

    // ----------------------------------------------------------------------------------

    Child_Mesh::Child_Mesh(Mesh_Modification_Template & aMeshModTemplate)
    : mElementTopology(CellTopology::TET4),
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
      mHasInterChildMeshInterface(false),
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

    // ----------------------------------------------------------------------------------

    Child_Mesh::~Child_Mesh()
    {
    }

    // ----------------------------------------------------------------------------------

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
        if(mSpatialDimension == 2 && !mHMR)
        {
            return EntityRank::EDGE;
        }
        else if(mSpatialDimension == 2 && mHMR)
        {
            return EntityRank::FACE;
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

    enum EntityRank
    Child_Mesh::get_facet_rank_internal() const
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
    std::shared_ptr<moris::mtk::Cell_Info>
    Child_Mesh::get_cell_info_sp() const
    {
        return mConnectivity;
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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

    moris_index
    Child_Mesh::get_cell_facet_ordinal(moris_index aCellIndex,
            moris_index aFacetIndex) const
    {
        moris::Matrix<moris::IndexMat> const & tCellToFacet = this->get_element_to_facet();
        for(moris::uint i = 0; i < tCellToFacet.n_cols(); i++)
        {
            if(tCellToFacet(aCellIndex,i) == aFacetIndex)
            {
                return (moris_index)i;
            }
        }

        return MORIS_INDEX_MAX;
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
            case EntityRank::NODE:
            {
                return mNodeParentRank(aCMLocIndex);
                break;
            }
            case EntityRank::EDGE:
            {
                MORIS_ASSERT(mHasEdgeConn, "Attempting to access edge ancestry in a child mesh that did not generate edges (in get_entity_parent_entity_rank)");
                return mEdgeParentRanks(aCMLocIndex);
                break;
            }
            case EntityRank::FACE:
            {
                MORIS_ASSERT(mHasFaceConn, "Attempting to access face ancestry in a child mesh that did not generate face (in get_entity_parent_entity_rank)");
                return mFaceParentRanks(aCMLocIndex);
                break;
            }
            case EntityRank::ELEMENT:
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
    Child_Mesh::get_entity_parent_entity_proc_ind(
            enum EntityRank    aEntityRank,
            moris::moris_index aCMLocIndex)
    {
        switch(aEntityRank)
        {
            case EntityRank::NODE:
            {
                return mNodeParentInds(aCMLocIndex);
                break;
            }
            case EntityRank::EDGE:
            {
                MORIS_ASSERT(mHasEdgeConn, "Attempting to access edge ancestry in a child mesh that did not generate edges (in get_entity_parent_entity_rank)");
                return mEdgeParentInds(aCMLocIndex);
                break;
            }
            case EntityRank::FACE:
            {
                MORIS_ASSERT(mHasFaceConn, "Attempting to access face ancestry in a child mesh that did not generate face (in get_entity_parent_entity_rank)");
                return mFaceParentInds(aCMLocIndex);
                break;
            }
            case EntityRank::ELEMENT:
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

    moris::moris_index
    Child_Mesh::get_child_mesh_index() const
    {
        return mChildMeshIndex;
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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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
    Child_Mesh::get_edge_ordinal_from_element_and_edge_indices(
            moris::moris_index               const & aElementIndex,
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
    Child_Mesh::get_edge_ordinal_from_element_and_edge_indices(
            moris::moris_index const & aElementIndex,
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
    Child_Mesh::get_face_ordinal_from_element_and_face_index(
            moris::moris_index const & aElementIndex,
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
    Child_Mesh::get_child_elements_connected_to_parent_facet(
            moris::moris_index         const & aParentFaceIndex,
            moris::uint                      & aNumElemsOnFace,
            moris::Matrix< moris::IdMat >    & aChildElemsIdsOnFacet,
            moris::Matrix< moris::IndexMat > & aChildElemsCMIndOnFacet,
            moris::Matrix< moris::IndexMat > & aChildElemOnFacetOrdinal) const
    {
        // First determine which of this meshes face live on the parent face
        enum EntityRank tFacetRank = this->get_facet_rank();

        moris::size_t tNumFacets = get_num_entities(this->get_facet_rank_internal());
        moris::size_t tNumElemsOnFace = 0;
        moris::Matrix< moris::IndexMat > tLocFacesOnParentFace(1,tNumFacets);

        moris::Matrix< moris::IndexMat > const & tFacetParentInds  = this->get_facet_parent_inds();
        moris::Matrix< moris::DDSTMat  > const & tFacetParentRanks = this->get_facet_parent_ranks();
        moris::Matrix< moris::IndexMat > const & tFacetToElement   = this->get_facet_to_element();

        // Iterate through face ancestry ranks
        for(moris::size_t i = 0; i < tNumFacets; i++)
        {
            //FIXME: FIgure out how to handle the 2d case in HMR since it calls edges faces ><
            if(tFacetParentInds(0,i) == aParentFaceIndex && tFacetParentRanks(0,i) ==(size_t) tFacetRank)
            {
                tLocFacesOnParentFace(0,tNumElemsOnFace) = i;
                tNumElemsOnFace++;
            }
        }

        // Check that arrays are sufficiently large
        MORIS_ERROR( aChildElemsCMIndOnFacet.numel() > 2*tNumElemsOnFace,
                "Child_Mesh::get_child_elements_connected_to_parent_facet - output vectors are not sufficiently large.\n" );

        // initialize counter
        aNumElemsOnFace = 0;

        for(moris::size_t i = 0; i<tNumElemsOnFace; i++)
        {
            moris::size_t tFacetInd = tLocFacesOnParentFace(0,i);

            for(moris::size_t j = 0; j<tFacetToElement.n_cols(); j++)
            {
                if(tFacetToElement(tFacetInd,j) != std::numeric_limits<moris::moris_index>::max())
                {
                    aChildElemsCMIndOnFacet(0,aNumElemsOnFace)  = tFacetToElement(tFacetInd,j);
                    aChildElemsIdsOnFacet(0,aNumElemsOnFace)    = mChildElementIds(0,tFacetToElement(tFacetInd,j));

                    aChildElemOnFacetOrdinal(0,aNumElemsOnFace) =
                            this->get_face_ordinal_from_element_and_face_index(tFacetToElement(tFacetInd,j),tFacetInd);

                    aNumElemsOnFace++;
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------

    moris::moris_index
    Child_Mesh::get_cm_local_node_index(moris::moris_index aNodeProcIndex) const
    {
        auto tIter = mNodeIndsToCMInd.find(aNodeProcIndex);

        MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(), "Node not in map, conversion to local indices failed");

        return tIter->second;
    }

    // ---------------------------------------------------------------------------------

    bool Child_Mesh::has_interface_along_geometry(moris_index aGeomIndex) const
    {
        uint tNumGeoms = mGeometryIndex.size();
        for (uint iG = 0; iG < tNumGeoms; iG++)
        {
            if (mGeometryIndex(iG) == aGeomIndex)
            {
                return true;
            }
        }

        return false;
    }

    // ---------------------------------------------------------------------------------

    moris_index
    Child_Mesh::get_local_geom_index(moris_index aGeomIndex) const
    {
        uint tNumGeoms = mGeometryIndex.size();
        for (uint iG = 0; iG < tNumGeoms; iG++)
        {
            if (mGeometryIndex(iG) == aGeomIndex)
            {
                return (moris_index)iG;
            }
        }

        return MORIS_INDEX_MAX;
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
    Child_Mesh::unzip_child_mesh_interface_get_interface_element_pairs(
            moris::uint                      aGeometryIndex,
            bool                           & aNoPairFoundFlag,
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
    Child_Mesh::unzip_child_mesh_interface(
            moris::moris_index                       aGeometryIndex,
            moris::Matrix< moris::IndexMat > const & aInterfaceElementPairsCMIndex,
            moris::Matrix< moris::IndexMat > const & aElementUsingZippedNodes,
            moris::Matrix< moris::IndexMat > const & aInterfaceNodeIndices,
            moris::Matrix< moris::IndexMat > const & aUnzippedInterfaceNodeIndices,
            moris::Matrix< moris::IndexMat > const & aUnzippedInterfaceNodeIds)
    {
        MORIS_ASSERT(mUnzippingFlag == true,
                "Error in child mesh, unzip_child_mesh_interface called on a child mesh which is not unzipping. Please call initialize_unzipping.");
        MORIS_ASSERT(aGeometryIndex == 0 ,
                "unzip_child_mesh_interface_get_interface_element_pairs not implemented for multiple geometries because mElementInterfaceFaces needs to accommodate this");

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

        MORIS_ASSERT(aElementUsingZippedNodes.numel() == aInterfaceElementPairsCMIndex.n_cols(),
                " Dimension mismatch between interface element pair matrix and element using unzipped nodes" );

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

    void
    Child_Mesh::add_vertices(moris::Cell<moris::mtk::Vertex const *> const & aVertices)
    {
        mVertices.append(aVertices);
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::set_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds)
    {
        MORIS_ASSERT(aNewNodeIds.n_cols() == get_num_entities(EntityRank::NODE),
                "Number of node ids does not match the number of nodes in this child mesh");

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
    Child_Mesh::add_node_parametric_coordinate(
            moris::Matrix< moris::IndexMat> const & aNodeIndices,
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

        // set the new entries in element interface sides to max
        for(moris::uint i = tCurrNumElem; i < tNewNumElem; i++)
        {
            for(moris::uint j = 0; j < mElementInterfaceSides.n_cols(); j++)
            {
                mElementInterfaceSides(i,j) = std::numeric_limits<moris::size_t>::max();
            }
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::insert_child_mesh_template(Mesh_Modification_Template & aMeshModTemplate)
    {
        // if we have more than 1 geometry we need to figure out if the intersected cell also had sides on the interface
        Matrix<moris::DDSTMat> tParentCellInterfaceSideOrdinals;
        this->setup_template_interface_facets( aMeshModTemplate, tParentCellInterfaceSideOrdinals );

        // Reindex the child mehs template
        // this takes it from ordinals and such to indices
        reindex_template(aMeshModTemplate);

        mElementTopology = aMeshModTemplate.mElementTopology;

        // construct connectivity
        moris::mtk::Cell_Info_Factory tConnFact;
        mConnectivity = tConnFact.create_cell_info_sp(mElementTopology);

        // Assert that reindex template actually happened
        MORIS_ASSERT(aMeshModTemplate.mIsReindexed,"Prior to inserting the template, reindexing must occur, in insert_child_mesh_template.");

        // index of current geometry
        moris_index tGeomIndex = mGeometryIndex.size()-1;

        // Replace the element which is marked as replaceable (if in the future multiple elements are replaced at one time then this could be a for loop)
        // But for now only one element is replaced

        replace_element(
                aMeshModTemplate.mElemIndToReplace,
                0,
                tGeomIndex,
                aMeshModTemplate.mNewElementToNode,
                aMeshModTemplate.mNewParentEdgeOrdinals,
                aMeshModTemplate.mNewParentEdgeRanks,
                aMeshModTemplate.mNewParentFaceOrdinals,
                aMeshModTemplate.mNewParentFaceRanks,
                tParentCellInterfaceSideOrdinals);

        // Number of elements that are new
        moris::size_t tNumElemToAdd = aMeshModTemplate.mNumNewElem-aMeshModTemplate.mNumElemToReplace;

        // Iterate over new elements to add and add them
        for(moris::size_t i = 0; i<tNumElemToAdd; i++)
        {
            add_element(
                    i+aMeshModTemplate.mNumElemToReplace,
                    tGeomIndex,
                    aMeshModTemplate.mNewElementToNode,
                    aMeshModTemplate.mNewParentEdgeOrdinals,
                    aMeshModTemplate.mNewParentEdgeRanks,
                    aMeshModTemplate.mNewParentFaceOrdinals,
                    aMeshModTemplate.mNewParentFaceRanks,
                    tParentCellInterfaceSideOrdinals);
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
    Child_Mesh::setup_template_interface_facets(
            Mesh_Modification_Template & aMeshModTemplate,
            Matrix<DDSTMat> &            aInterfaceSideOrdinals)
    {
        if(mGeometryIndex.size() > 1)
        {
            // resize and fill the output
            aInterfaceSideOrdinals.resize(aMeshModTemplate.mNumNewElem,mGeometryIndex.size());
            aInterfaceSideOrdinals.fill(std::numeric_limits<moris::size_t>::max());

            // last column always just the new interface sides
            aInterfaceSideOrdinals.get_column(aInterfaceSideOrdinals.n_cols()-1) = aMeshModTemplate.mNewElementInterfaceSides.get_column(0);

            // get the facet rank
            enum EntityRank tFacetRank = EntityRank::INVALID;

            // current index
            moris::Matrix<moris::IndexMat> * tParentFacetOrds = nullptr;
            moris::Matrix<moris::DDSTMat> * tParentFacetRanks = nullptr;
            if(mSpatialDimension == 3)
            {
                tParentFacetOrds =  & aMeshModTemplate.mNewParentFaceOrdinals;
                tParentFacetRanks =  & aMeshModTemplate.mNewParentFaceRanks;
                tFacetRank = EntityRank::FACE;
            }
            else if(mSpatialDimension == 2)
            {
                tParentFacetOrds =  & aMeshModTemplate.mNewParentEdgeOrdinals;
                tParentFacetRanks =  & aMeshModTemplate.mNewParentEdgeRanks;
                tFacetRank = EntityRank::EDGE;
            }
            else
            {
                MORIS_ASSERT(0,"Only 2d and 3d meshes supported");
            }

            // iterate through previous geometry indices
            for(moris::uint  iG = 0 ; iG < mGeometryIndex.size() - 1; iG++ )
            {
                // side ordinal of base cell on interface with respect to a previous geometry iG
                size_t tInterfaceSideOrd = mElementInterfaceSides(aMeshModTemplate.mElemIndToReplace,iG);
                if(tInterfaceSideOrd != std::numeric_limits<size_t>::max())
                {
                    // iterate through the facet ords ods 
                    for(moris::uint iCell = 0; iCell < tParentFacetOrds->n_rows(); iCell++ )
                    {
                        for(moris::uint iF = 0; iF < tParentFacetOrds->n_cols(); iF++ )
                        {
                            if((*tParentFacetRanks)(iCell,iF) == (size_t) tFacetRank)
                            {
                                if((*tParentFacetOrds)(iCell,iF) == (moris_index) tInterfaceSideOrd)
                                {
                                    aInterfaceSideOrdinals(iCell,iG) = iF;
                                }
                            }
                        }
                    }
                }
            }

        }
        else
        {
            aInterfaceSideOrdinals = aMeshModTemplate.mNewElementInterfaceSides;
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::reindex_template(Mesh_Modification_Template & aMeshModTemplate)
    {
        MORIS_ASSERT(!aMeshModTemplate.mIsReindexed,
                "reindex_template() called on a child mesh modification template which has already been reindexed");

        aMeshModTemplate.mNewElementToNode = reindex_matrix(aMeshModTemplate.mNewElementToNode, 0, aMeshModTemplate.mNodeInds);

        // Reindex template edge and face ordinals
        reindex_template_parent_information(aMeshModTemplate);
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::replace_element(
            moris::size_t                    const & aElementIndexToReplace,
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
        replace_row(aRowIndex, aElementInterfaceFaces , aElementIndexToReplace, mElementInterfaceSides);
    }


    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::add_element(
            moris::size_t                    const & aRowIndex,
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

        replace_row(aRowIndex, aElementInterfaceFaces , mNumElem, mElementInterfaceSides);

        mNumElem++;
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::set_node_inheritance(
            moris::moris_index aCMLocNodeIndex,
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
    Child_Mesh::generate_connectivities(
            bool aGenerateFaceConn,
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

        // last col index
        moris_index tLastCol = mGeometryIndex.size()-1;
        for(moris::uint i = 0; i < this->get_num_entities(EntityRank::ELEMENT); i++)
        {
            mElementInterfaceSides(i,tLastCol) = std::numeric_limits<moris::size_t>::max();
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::add_entity_to_intersect_connectivity(
            moris::moris_index aCMNodeInd,
            moris::moris_index aCMEdgeInd,
            moris::size_t      aFlag)
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
    Child_Mesh::mark_facet_as_on_interface(moris_index aFacetCMIndex,
                               moris_index aGeometryIndex)
    {
        // facet to element
        moris::Matrix<moris::IndexMat> const & tFacetToElement = this->get_facet_to_element();

        // local geometry index
        moris_index tLocalGeomIndex = this->get_local_geom_index(aGeometryIndex);

        // iterate through elements attached to this facet
        for(moris::uint iCell = 0; iCell < tFacetToElement.n_cols(); iCell++)
        {
            // cell index
            moris_index tCellIndex = tFacetToElement(aFacetCMIndex,iCell);
            moris_index tSideOrdinal = MORIS_INDEX_MAX;

            // only do this if we have a valid cell attached to this facet
            if(tCellIndex != MORIS_INDEX_MAX)
            {

                if(mSpatialDimension == 2)
                {
                    tSideOrdinal = this->get_edge_ordinal_from_element_and_edge_indices(tCellIndex,aFacetCMIndex);
                }
                else if (mSpatialDimension == 3)
                {
                    tSideOrdinal = this->get_face_ordinal_from_element_and_face_index(tCellIndex,aFacetCMIndex);
                }
                else
                {
                    MORIS_ERROR(0,"two and three dimensions only");
                }

                mElementInterfaceSides(tCellIndex,tLocalGeomIndex) = (size_t)tSideOrdinal;
            }
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::mark_interface_faces_from_interface_coincident_faces()
    {
        if(mHasCoincidentEdges && mSpatialDimension == 3)
        {
            // do it differently in 3d because it is a more involved process
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

        else if(mHasCoincidentEdges && mSpatialDimension == 2)
        {
            moris::Matrix< moris::IndexMat > const tEdgeToElem = this->get_edge_to_element();
            moris_index tGeomIndex = mGeometryIndex.size()-1;
            moris::moris_index tDummy = std::numeric_limits<moris::moris_index>::max();

            // iterate through edges on interface
            for(moris::uint iEdge = 0; iEdge < mEdgeOnInterface.numel(); iEdge++ )
            {
                if(mEdgeOnInterface(iEdge) == 1)
                {
                    for(uint i = 0; i <tEdgeToElem.n_cols(); i++)
                    {
                        if(tEdgeToElem(iEdge,i) ==  tDummy)
                        {
                            break;
                        }
                        moris::moris_index tElemInd = tEdgeToElem(iEdge,i);
                        moris::moris_index tEdgeOrd = this->get_edge_ordinal_from_element_and_edge_indices(tElemInd,iEdge);
                        mElementInterfaceSides(tElemInd,tGeomIndex) = tEdgeOrd;
                    }
                }
            }
        }

        mHasCoincidentEdges = false;
        mEdgeOnInterface.set_size(0,0);
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
    moris::Matrix< moris::IndexMat > const &
    Child_Mesh::get_element_phase_indices() const
    {
        MORIS_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
        return mElementPhaseIndices;
    }
    // ---------------------------------------------------------------------------------
    moris_index
    Child_Mesh::get_element_subphase_id(moris::size_t const & aEntityIndex) const
    {
        MORIS_ASSERT(aEntityIndex<get_num_entities(EntityRank::ELEMENT),"EntityIndex out of bounds, aEntityIndex should be a child mesh local index");
        MORIS_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
        return mSubPhaseBinId(mElementBinIndex(aEntityIndex));
    }

    // ---------------------------------------------------------------------------------

    moris::size_t
    Child_Mesh::get_num_subphase_bins() const
    {
        return mSubPhaseBins.size();
    }

    // ---------------------------------------------------------------------------------

    Cell<moris::moris_index> const &
    Child_Mesh::get_subphase_bin_bulk_phase() const
    {
        return mBinBulkPhase;
    }

    // ---------------------------------------------------------------------------------

    Cell<moris::Matrix< moris::IndexMat >> const &
    Child_Mesh::get_subphase_groups() const
    {
        return mSubPhaseBins;
    }

    // ---------------------------------------------------------------------------------

    Cell<moris_index> const &
    Child_Mesh::get_subphase_indices( ) const
    {
        return mSubPhaseBinIndices;
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat> const &
    Child_Mesh::get_subphase_ids( ) const
    {
        return mSubPhaseBinId;
    }

    // ---------------------------------------------------------------------------------

    moris_index
    Child_Mesh::get_subphase_loc_index(moris_index aSubPhaseIndex) const
    {
        moris_index tLocSubIndex = MORIS_INDEX_MAX;
        for(moris::moris_index i = 0; i < (moris_index)mSubPhaseBinIndices.size(); i++)
        {   
            if(mSubPhaseBinIndices(i) == aSubPhaseIndex)
            {
                tLocSubIndex = i;
                break;
            }
        }

        MORIS_ASSERT(tLocSubIndex != MORIS_INDEX_MAX,"Subphase index does not show up in child mesh");
        return tLocSubIndex;
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::get_subphases_attached_to_facet(
            moris_index         aFacetIndex,
            Cell<moris_index> & aSubPhaseCMIndex,
            Cell<moris_index> & aRepresentativeChildCellInd,
            Cell<moris_index> & aRepresentativeChildCellSideOrdinal) const
    {

        aRepresentativeChildCellInd.clear();
        aRepresentativeChildCellSideOrdinal.clear();
        // estimate maximum number of elements on face
         const uint tMaxElemOnFace = 100;

        // get child cells on facet
        Matrix<IndexMat> tChildElemsIdsOnFace(1,tMaxElemOnFace);
        Matrix<IndexMat> tChildElemsCMIndOnFace(1,tMaxElemOnFace);
        Matrix<IndexMat> tChildElemOnFaceOrdinal(1,tMaxElemOnFace);

        // define varialbe for actual number of child elements on face
        uint tNumberOfChildElemsOnFace;

        this->get_child_elements_connected_to_parent_facet(
                aFacetIndex,
                tNumberOfChildElemsOnFace,
                tChildElemsIdsOnFace,
                tChildElemsCMIndOnFace,
                tChildElemOnFaceOrdinal);

        // get the subphases of cells
        moris::Matrix< moris::IndexMat > const & tElementSubphase = this->get_elemental_subphase_bin_membership();
        moris::Matrix<moris::IndexMat> const & tElementInds = this->get_element_inds();

        std::unordered_map<moris_index,moris_index> tUniqueMap;

        for(moris::uint i = 0; i < tNumberOfChildElemsOnFace; i++)
        {
            moris_index tSubphase = tElementSubphase(tChildElemsCMIndOnFace(i));
            if(tUniqueMap.find(tSubphase) == tUniqueMap.end())
            {
                // get the subphase membership of cell on facet
                aSubPhaseCMIndex.push_back(tSubphase);
                aRepresentativeChildCellInd.push_back(tElementInds(tChildElemsCMIndOnFace(i)));
                aRepresentativeChildCellSideOrdinal.push_back(tChildElemOnFaceOrdinal(i));
                tUniqueMap[tSubphase] = 1; // value not important
            }
        }
    }

    // ---------------------------------------------------------------------------------

    bool
    Child_Mesh::construct_internal_double_sides_between_subphases()
    {
        Matrix<IndexMat> const & tFacetToCell = this->get_facet_to_element();

        Matrix<IndexMat> const & tCellFacets  = this->get_element_to_facet();

        // information about subphases
        moris::Cell<moris::Matrix< moris::IndexMat >> const &  tSubphaseClusters = this->get_subphase_groups();
        Cell<moris::moris_index> const & tSubphaseBulkPhases  = this->get_subphase_bin_bulk_phase();

        // facet rank
        enum EntityRank tFacetRank = this->get_facet_rank();

        // get the interface side ords for the cells in this child mesh (rows local cell index in cm, col geometry)
        moris::Matrix< moris::DDSTMat  > const & tCellInterfaceSideOrds = this->get_cell_interface_side_ords();

        // facet parent inds and ranks
        // moris::Matrix< moris::IndexMat > const & tFacetParentInds  = tChildMesh->get_facet_parent_inds();
        moris::Matrix< moris::DDSTMat >  const & tFacetParentRanks = this->get_facet_parent_ranks();
        moris::Matrix<moris::IndexMat>   const & tCellFacetInds    = this->get_element_to_facet();

        // determine if there is an interchild mesh interface
        for (moris::uint iLC = 0; iLC < tCellInterfaceSideOrds.n_rows(); iLC++)
        {
            for (moris::uint iG = 0; iG < tCellInterfaceSideOrds.n_cols(); iG++)
            {
                moris::size_t tInterfaceOrdinal = tCellInterfaceSideOrds(iLC, iG);

                // if this is an interface facet
                if (tInterfaceOrdinal != std::numeric_limits<moris::size_t>::max())
                {
                    moris_index tFacetIndex = tCellFacetInds(iLC, tInterfaceOrdinal);
                    if (tFacetParentRanks(tFacetIndex) == (size_t)tFacetRank)
                    {
                        mHasInterChildMeshInterface = true;
                    }
                }
            }
        }

        //iterate through subphases and create neighborhoods between subphases
        // from low bulk phase to high bulk phase always
        for(moris::uint iSP0 = 0; iSP0 < tSubphaseClusters.size(); iSP0++)
        {
            for(moris::uint iSP1 = 0; iSP1 < tSubphaseClusters.size(); iSP1++)
            {
                // skip the same index and if the first subphase bulk phase is higher
                if(iSP1 == iSP0 || tSubphaseBulkPhases(iSP0) > tSubphaseBulkPhases(iSP1))
                {
                    continue;
                }

                // get cells in subphase isp0
                Matrix<IndexMat> const & tSubphaseCells0 = mSubPhaseBins(iSP0);

                moris_index tDblSideIndex = MORIS_INDEX_MAX;

                // count how many pairs are between this subphase combo
                moris_index tCount = 0;

                for(moris::uint i = 0; i < tSubphaseCells0.numel(); i++)
                {
                    moris_index tCMLocCellInd = tSubphaseCells0(i);

                    // iterate through geometries
                    for(moris::uint iG =0; iG< this->mGeometryIndex.size(); iG++)
                    {
                        if(mElementInterfaceSides(tCMLocCellInd,iG) != std::numeric_limits<size_t>::max())
                        {
                            moris_index tSideOrd = mElementInterfaceSides(tCMLocCellInd,iG);
                            moris_index tFacetCmIndex = tCellFacets(tCMLocCellInd,tSideOrd);

                            // figure out the cell index of the neighbor
                            moris_index tInterfaceNeighbor        = MORIS_INDEX_MAX;
                            moris_index tInterfaceNeighborSideOrd = MORIS_INDEX_MAX;
                            for(moris::uint iS = 0; iS < tFacetToCell.n_cols(); iS++)
                            {
                                if(tFacetToCell(tFacetCmIndex,iS) != tCMLocCellInd )
                                {
                                    tInterfaceNeighbor = tFacetToCell(tFacetCmIndex,iS);
                                    if(tInterfaceNeighbor != std::numeric_limits<moris_index>::max())
                                    {
                                        tInterfaceNeighborSideOrd = this->get_cell_facet_ordinal(tInterfaceNeighbor,tFacetCmIndex);
                                        break;
                                    }
                                    else
                                    {
                                        mHasInterChildMeshInterface = true;
                                    }
                                }
                            }

                            // if this cell belongs to the other subphase add the pair and shared side ordinal

                            if(tInterfaceNeighbor != std::numeric_limits<moris_index>::max())
                            {
                                if(this->get_element_subphase_index(tInterfaceNeighbor) == (moris_index) iSP1)
                                {
                                    // add this subphase pair to the double side sets if we havent done so
                                    if(tCount == 0)
                                    {
                                        tDblSideIndex = mDoubleSideSetSubphaseInds.size();
                                        mDoubleSideSetSubphaseInds.push_back({(moris_index)iSP0,(moris_index)iSP1});
                                        mDoubleSideSetCellPairs.push_back(Cell<Cell< moris_index >>(0));
                                        mDoubleSideSetFacetPairs.push_back(Cell<Cell< moris_index >>(0));
                                    }

                                    // add the pair information to the set
                                    mDoubleSideSetCellPairs(tDblSideIndex).push_back({tCMLocCellInd,tInterfaceNeighbor});
                                    mDoubleSideSetFacetPairs(tDblSideIndex).push_back({tSideOrd,tInterfaceNeighborSideOrd});

                                    tCount++;
                                }
                            }

                            else
                            {
                              mHasInterChildMeshInterface = true;
                            }
                        }
                    }
                }
            }
        }

        return mHasInterChildMeshInterface;
    }

    bool
    Child_Mesh::has_inter_child_mesh_interfaces()
    {
        return mHasInterChildMeshInterface;
    }

    // ---------------------------------------------------------------------------------

    uint
    Child_Mesh::get_num_double_side_interfaces() const
    {
        return mDoubleSideSetSubphaseInds.size();
    }

    // ---------------------------------------------------------------------------------

    Cell< moris_index > const &
    Child_Mesh::get_double_side_interface_subphase_indices(moris_index aDblSideCMIndex) const
    {
        return mDoubleSideSetSubphaseInds(aDblSideCMIndex);
    }

    // ---------------------------------------------------------------------------------

    Cell<Cell< moris_index >> const &
    Child_Mesh::get_double_side_interface_cell_pairs(moris_index aDblSideCMIndex) const
    {
        return mDoubleSideSetCellPairs(aDblSideCMIndex);
    }

    // ---------------------------------------------------------------------------------

    Cell<Cell< moris_index >> const &
    Child_Mesh::get_double_side_interface_cell_pairs_facet_ords(moris_index aDblSideCMIndex) const
    {
        return mDoubleSideSetFacetPairs(aDblSideCMIndex);
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::delete_double_sides_interface_sets()
    {
        mDoubleSideSetSubphaseInds.resize(0);
        mDoubleSideSetCellPairs.resize(0);
        mDoubleSideSetFacetPairs.resize(0);
    }
    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::print_double_sides_between_subphases(moris_index aVerboseLevel)
    {
        moris::uint tNumDoubleSideSets = mDoubleSideSetSubphaseInds.size();

        std::cout<<"\n--------------------------------------------------------"<<std::endl;
        std::cout<<"Child Mesh Double Sides:"<<std::endl;
        std::cout<<" Num Double Side Sets: "<<std::setw(9)<<tNumDoubleSideSets<<std::endl;

        for(moris::uint i = 0 ; i <tNumDoubleSideSets; i++)
        {
            std::cout<<" Subphase Pairs: "<<std::setw(9)<<mSubPhaseBinIndices(mDoubleSideSetSubphaseInds(i)(0));
            std::cout<<std::setw(9)<<mSubPhaseBinIndices(mDoubleSideSetSubphaseInds(i)(1));
            std::cout<<" | Num Cell Pairs: "<<mDoubleSideSetCellPairs(i).size();

            if(aVerboseLevel > 0)
            {
                std::cout<<"\n Cell Pairs: "<<std::endl;
                for(moris::uint j = 0 ; j < mDoubleSideSetCellPairs(i).size(); j++)
                {
                    std::cout<<"  Cell 1 ID/Ord: "<<std::setw(9)<<mChildElementIds(mDoubleSideSetCellPairs(i)(j)(0))<<std::setw(9)<<mDoubleSideSetFacetPairs(i)(j)(0);
                    std::cout<<"  Cell 2 ID/Ord: "<<std::setw(9)<<mChildElementIds(mDoubleSideSetCellPairs(i)(j)(1))<<std::setw(9)<<mDoubleSideSetFacetPairs(i)(j)(1)<<std::endl;
                }
            }
            std::cout<<std::endl;
        }
        std::cout<<"--------------------------------------------------------"<<std::endl;
    }    \

    // ---------------------------------------------------------------------------------

    Cell<moris_index> const &
    Child_Mesh::get_subphase_basis_enrichment_levels(moris_index aSubphaseBin) const
    {
        MORIS_ASSERT(aSubphaseBin < (moris_index)mSubphaseBasisEnrichmentLevel.size(), "Subphase group index out of bounds");

        return mSubphaseBasisEnrichmentLevel(aSubphaseBin);
    }
    
    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::reindex_cells(Cell<moris_index> & aOldIndexToNewCellIndex)
    {
        for(moris::uint iC = 0; iC < mChildElementInds.numel(); iC++)
        {
            moris_index tOldIndex = mChildElementInds(iC);
            moris_index tNewIndex = aOldIndexToNewCellIndex(tOldIndex);
            MORIS_ASSERT(tNewIndex != MORIS_INDEX_MAX,"Trying to reindex with a max value. Was this cell deleted?");
            mChildElementInds(iC) = tNewIndex;
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::set_elemental_subphase(
            moris::moris_index                     & aSubPhaseIndex,
            moris::Matrix< moris::IndexMat > const & aElementSubPhase)
    {
        this->construct_subphase_bins(aSubPhaseIndex,aElementSubPhase);
    }

    // ----------------------------------------------------------------------------------

    void
    Child_Mesh::set_subphase_id(moris_id const & aSubphaseIndex,
            moris_id & aSubphaseId)
    {
        mSubPhaseBinId(aSubphaseIndex) = aSubphaseId;
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat > const &
    Child_Mesh::get_elemental_subphase_bin_membership() const
    {
        return mElementBinIndex;
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::add_basis_and_enrichment_to_subphase_group(
            moris_index aSubphaseBinIndex,
            moris_index aBasisIndex,
            moris_index aBasisEnrLev)
    {
        moris_index tLocSubIndex = this->get_subphase_loc_index(aSubphaseBinIndex);
        mSubphaseBasisIndices(tLocSubIndex).push_back(aBasisIndex);
        mSubphaseBasisEnrichmentLevel(tLocSubIndex).push_back(aBasisEnrLev);
    }

    // ---------------------------------------------------------------------------------

    Cell<moris_index> const &
    Child_Mesh::get_subphase_basis_indices(moris_index aSubphaseBin) const
    {
        MORIS_ASSERT(aSubphaseBin < (moris_index)mSubphaseBasisIndices.size(),"Subphase group index out of bounds");

        return mSubphaseBasisIndices(aSubphaseBin);
    }

    // ---------------------------------------------------------------------------------
    moris::Matrix< moris::DDSTMat  > const & 
    Child_Mesh::get_cell_interface_side_ords()
    {
        return mElementInterfaceSides;
    }

    // ---------------------------------------------------------------------------------
    void
    Child_Mesh::pack_child_mesh_by_phase(
            moris::size_t                 const & aNumPhases,
            Cell<moris::Matrix< moris::IdMat >> & aElementIds,
            Cell<moris::Matrix< moris::IdMat >> & aElementCMInds) const
    {
        moris::size_t tNumElems = get_num_entities(EntityRank::ELEMENT);

        aElementIds = Cell<moris::Matrix< moris::IdMat >>(aNumPhases);
        aElementCMInds = Cell<moris::Matrix< moris::IdMat >>(aNumPhases);

        for(moris::size_t i = 0; i<aNumPhases; i++)
        {
            aElementIds(i) = moris::Matrix< moris::IdMat >(1,tNumElems);
            aElementCMInds(i) = moris::Matrix< moris::IdMat >(1,tNumElems);
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

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Child_Mesh::pack_interface_sides(
            moris_index aGeometryIndex,
            moris_index aPhaseIndex0,
            moris_index aPhaseIndex1,
            moris_index aIndexFlag) const
    {
        // if this geometry index does not show up in this child mesh
        if(!this->has_interface_along_geometry(aGeometryIndex))
        {
            return Matrix<moris::IdMat>(0,0);
        }
        // Loop bound and sizing
        moris::size_t tNumElem = get_num_entities(EntityRank::ELEMENT);

        // local geometry index
        moris_index tLocGeomIndex = this->get_local_geom_index(aGeometryIndex);

        moris::Matrix< moris::IdMat > tInterfaceCMInfo(tNumElem,2,MORIS_ID_MAX);

        // Keep track of the number of interface sides
        moris::size_t tCount = 0;

        // construct a membership map in interface side set
        Matrix<IndexMat> tCellsIndexInSet(this->get_num_entities(EntityRank::ELEMENT),1,MORIS_INDEX_MAX);

        // Iterate over each element and if the element has an interface side it will be in mElementInferfaceSides vector
        for(moris::size_t iEl =0 ; iEl<tNumElem; iEl++)
        {
            if(mElementInterfaceSides(iEl,tLocGeomIndex) != std::numeric_limits<moris::size_t>::max())
            {
                if((moris_index)this->get_element_phase_index(iEl) == aPhaseIndex0 || (moris_index)this->get_element_phase_index(iEl) == aPhaseIndex1)
                {
                    tInterfaceCMInfo(tCount,0) = iEl;

                    tCellsIndexInSet(iEl) = tCount;

                    tInterfaceCMInfo(tCount,1) = mElementInterfaceSides(iEl,tLocGeomIndex);
                    tCount++;
                }
            }
        }

        // Size out extra space
        tInterfaceCMInfo.resize(tCount,2);

        // get the interface pairs with respect to this geometry
        bool tNoPair = true;
        moris::Matrix<moris::IndexMat> tInterfaceElementPairs;
        moris::Matrix<moris::IndexMat> tInterfacePairSideOrds;
        unzip_child_mesh_interface_get_interface_element_pairs(tLocGeomIndex,tNoPair,tInterfaceElementPairs,tInterfacePairSideOrds);

        MORIS_ASSERT(!tNoPair,"No pair found");

        Matrix<IndexMat> tCellsToKeepInSet(this->get_num_entities(EntityRank::ELEMENT),1,MORIS_INDEX_MAX);

        // iterate through element pairs
        uint tKeepCount = 0;
        for(moris::uint i = 0; i < tInterfaceElementPairs.n_cols(); i ++)
        {
            moris_index tCell0 = tInterfaceElementPairs(0,i);
            moris_index tCell1 = tInterfaceElementPairs(1,i);

            if(tCellsIndexInSet(tCell0) != MORIS_INDEX_MAX && tCellsIndexInSet(tCell1) != MORIS_INDEX_MAX )
            {
                tCellsToKeepInSet(tCell0) = 1;
                tCellsToKeepInSet(tCell1) = 1;
                tKeepCount = tKeepCount+2;
            }
        }

        moris::Matrix< moris::IdMat > tInterfaceInfo(tKeepCount,2);
        tCount = 0;
        for(moris::size_t iEl =0 ; iEl<tNumElem; iEl++)
        {
            if(tCellsToKeepInSet(iEl) == 1 && (moris_index)get_element_phase_index(iEl) == aPhaseIndex0)
            {
                moris_index tIndexInSet = tCellsIndexInSet(iEl);
                tInterfaceInfo(tCount,1) = tInterfaceCMInfo(tIndexInSet,1);

                if(aIndexFlag == 1)
                {
                    tInterfaceInfo(tCount,0) = mChildElementInds(tInterfaceCMInfo(tIndexInSet,0));
                }
                else if (aIndexFlag == 0)
                {
                    tInterfaceInfo(tCount,0) = mChildElementIds(tInterfaceCMInfo(tIndexInSet,0));
                }
                else
                {
                    tInterfaceInfo(tCount,0) = tInterfaceCMInfo(tIndexInSet,0);
                }
                tCount++;
            }
        }

        // convert to procs indices or global ids
        tInterfaceInfo.resize(tCount,2);

        return tInterfaceInfo;
    }

    // ---------------------------------------------------------------------------------

    moris::Memory_Map
    Child_Mesh::get_memory_usage()
    {
        moris::Memory_Map tMemoryMap;

        tMemoryMap.mMemoryMapData["mParentElementIndex"] = sizeof(mParentElementIndex);
        tMemoryMap.mMemoryMapData["mElementTopology"] = sizeof(mElementTopology);
        tMemoryMap.mMemoryMapData["mConnectivity"] = sizeof(mConnectivity);
        tMemoryMap.mMemoryMapData["mElementToNode"] = mElementToNode.capacity();
        tMemoryMap.mMemoryMapData["mElementEdgeParentInds"] = mElementEdgeParentInds.capacity();
        tMemoryMap.mMemoryMapData["mElementEdgeParentRanks"] = mElementEdgeParentRanks.capacity();
        tMemoryMap.mMemoryMapData["mElementFaceParentInds"] = mElementFaceParentInds.capacity();
        tMemoryMap.mMemoryMapData["mElementFaceParentRanks"] = mElementFaceParentRanks.capacity();
        tMemoryMap.mMemoryMapData["mElementInterfaceSides"] = mElementInterfaceSides.capacity();
        tMemoryMap.mMemoryMapData["mElementEdgeParentInds"] = mElementEdgeParentInds.capacity();
        tMemoryMap.mMemoryMapData["mElementEdgeParentInds"] = mElementEdgeParentInds.capacity();
        tMemoryMap.mMemoryMapData["mSpatialDimension"] = sizeof(mSpatialDimension);
        tMemoryMap.mMemoryMapData["mGeometryIndex"] = mGeometryIndex.capacity();
        tMemoryMap.mMemoryMapData["mChildElementIds"] = mChildElementIds.capacity();
        tMemoryMap.mMemoryMapData["mChildElementInds"] = mChildElementInds.capacity();
        tMemoryMap.mMemoryMapData["mVertices"] = mVertices.capacity();
        tMemoryMap.mMemoryMapData["mNodeIds"] = mNodeIds.capacity();
        tMemoryMap.mMemoryMapData["mNodeInds"] = mNodeInds.capacity();
        tMemoryMap.mMemoryMapData["mNodeParentRank"] = mNodeParentRank.capacity();
        tMemoryMap.mMemoryMapData["mNodeParentInds"] = mNodeParentInds.capacity();
        tMemoryMap.mMemoryMapData["mNodeIndsToCMInd"] = mNodeIndsToCMInd.size() * 2 * sizeof(moris::size_t);
        tMemoryMap.mMemoryMapData["mNodeParametricCoord"] = mNodeParametricCoord.capacity();
        tMemoryMap.mMemoryMapData["mFaceToNode"] = mFaceToNode.capacity();
        tMemoryMap.mMemoryMapData["mNodeToFace"] = mNodeToFace.capacity();
        tMemoryMap.mMemoryMapData["mFaceToElement"] = mFaceToElement.capacity();
        tMemoryMap.mMemoryMapData["mElementToFace"] = mElementToFace.capacity();
        tMemoryMap.mMemoryMapData["mFaceParentInds"] = mFaceParentInds.capacity();
        tMemoryMap.mMemoryMapData["mFaceParentInds"] = mFaceParentInds.capacity();
        tMemoryMap.mMemoryMapData["mFaceParentRanks"] = mFaceParentRanks.capacity();
        tMemoryMap.mMemoryMapData["mEdgeToNode"] = mEdgeToNode.capacity();
        tMemoryMap.mMemoryMapData["mNodeToEdge"] = mNodeToEdge.capacity();
        tMemoryMap.mMemoryMapData["mEdgeToElement"] = mEdgeToElement.capacity();
        tMemoryMap.mMemoryMapData["mElementToEdge"] = mElementToEdge.capacity();
        tMemoryMap.mMemoryMapData["mEdgeParentInds"] = mEdgeParentInds.capacity();
        tMemoryMap.mMemoryMapData["mEdgeParentRanks"] = mEdgeParentRanks.capacity();
        tMemoryMap.mMemoryMapData["mElementToElement"] = mElementToElement.capacity();
        tMemoryMap.mMemoryMapData["mIntersectConnectivity"] = mIntersectConnectivity.capacity();
        tMemoryMap.mMemoryMapData["mIntersectedCMNodeIndex"] = mIntersectedCMNodeIndex.capacity();
        tMemoryMap.mMemoryMapData["mIntersectedEdges"] = mIntersectedEdges.capacity();
        tMemoryMap.mMemoryMapData["mEdgeOnInterface"] = mEdgeOnInterface.capacity();
        tMemoryMap.mMemoryMapData["mElementPhaseIndices"] = mElementPhaseIndices.capacity();
        tMemoryMap.mMemoryMapData["mElementBinIndex"] = mElementBinIndex.capacity();
        tMemoryMap.mMemoryMapData["mBinBulkPhase"] = mBinBulkPhase.capacity();
        tMemoryMap.mMemoryMapData["mElementBinIndex"] = mElementBinIndex.capacity();
        tMemoryMap.mMemoryMapData["mBinBulkPhase"] = mBinBulkPhase.capacity();
        tMemoryMap.mMemoryMapData["mSubPhaseBinId"] = mSubPhaseBinId.capacity();
        tMemoryMap.mMemoryMapData["mSubPhaseBinIndices"] = mSubPhaseBinIndices.capacity();
        tMemoryMap.mMemoryMapData["mSubPhaseBins"] = mSubPhaseBins.capacity();
        tMemoryMap.mMemoryMapData["mSubphaseBasisIndices"] = moris::internal_capacity(mSubphaseBasisIndices);
        tMemoryMap.mMemoryMapData["mSubphaseBasisEnrichmentLevel"] = mSubphaseBasisEnrichmentLevel.capacity();
        tMemoryMap.mMemoryMapData["mDoubleSideSetSubphaseInds"] = mDoubleSideSetSubphaseInds.capacity();
        tMemoryMap.mMemoryMapData["mDoubleSideSetCellPairs"] = mDoubleSideSetCellPairs.capacity();
        tMemoryMap.mMemoryMapData["mDoubleSideSetFacetPairs"] = mDoubleSideSetFacetPairs.capacity();

        return tMemoryMap;
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Child_Mesh::pack_interface_sides_loc_inds() const
    {
        // Loop bound and sizing
        moris::size_t tNumElem = get_num_entities(EntityRank::ELEMENT);

        moris::Matrix< moris::IdMat > tInterfaceSideSetInfo(tNumElem,2);

        // Keep track of the number of interface sides
        moris::size_t tCount = 0;

        // Iterate over each element and if the element has an interface side it will be in mElementInferfaceSides vector
        for(moris::size_t iEl =0 ; iEl<tNumElem; iEl++)
        {
            //TODO: NOTE THIS WILL NOT WORK WITH MULTI-MATERIAL YET (AT LEAST NOT SPLIT THEM UP)
            if(mElementInterfaceSides(iEl) != std::numeric_limits<moris::size_t>::max())
            {
                tInterfaceSideSetInfo(tCount,0) = mChildElementInds(iEl);
                tInterfaceSideSetInfo(tCount,1) = mElementInterfaceSides(iEl);
                tCount++;
            }
        }

        // Size out space
        tInterfaceSideSetInfo.resize(tCount,2);

        return tInterfaceSideSetInfo;
    }

    // ---------------------------------------------------------------------------------

    enum CellTopology
    Child_Mesh::template_to_cell_topology(enum TemplateType aTemplateType)
    {
        switch(aTemplateType)
        {
            case(TemplateType::REGULAR_SUBDIVISION_HEX8):
            case(TemplateType::TET_4):
            {
                return CellTopology::TET4;
                break;
            }
            case(TemplateType::HEX_8):
                          {
                return CellTopology::HEX8;
                break;
                          }
            case(TemplateType::REGULAR_SUBDIVISION_QUAD4):
            case(TemplateType::TRI_3):
            {
                return CellTopology::TRI3;
                break;
            }
            case(TemplateType::QUAD_4):
                          {
                return CellTopology::QUAD4;
                break;
                          }
            default:
            {
                MORIS_ASSERT(0, "Invalid topology for this function at present.");
            }
        }

        return CellTopology::INVALID;
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::generate_face_connectivity_and_ancestry(moris::Matrix< moris::IndexMat > const & aElementToNodeLocal)
    {

        xtk::create_faces_from_element_to_node(
                mConnectivity.get(),
                this->get_num_entities(EntityRank::NODE),
                aElementToNodeLocal,
                mElementToFace,
                mFaceToNode,
                mNodeToFace,
                mFaceToElement);

        // Convert face to node from cm indices to proc indices
        mFaceToNode = this->convert_to_proc_indices(mFaceToNode);

        mHasFaceConn = true;

        this->setup_face_ancestry();

        mNodeToFace.set_size(0,0);
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::generate_edge_connectivity_and_ancestry(moris::Matrix< moris::IndexMat > const & aElementToNodeLocal)
    {
        create_edges_from_element_to_node(
                mElementTopology,
                get_num_entities(EntityRank::NODE),
                aElementToNodeLocal,
                mElementToEdge,
                mEdgeToNode,
                mNodeToEdge,
                mEdgeToElement);

        // convert back from local to proc indices
        mEdgeToNode = convert_to_proc_indices(mEdgeToNode);

        mHasEdgeConn = true;

        mEdgeOnInterface.resize(1,mEdgeToNode.n_rows());

        setup_edge_ancestry();

        mNodeToEdge.set_size(0,0);
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::generate_element_to_element_connectivity()
    {
        if(mElementTopology == CellTopology::TET4)
        {
            MORIS_ASSERT(mHasFaceConn,
                    "Face connectivity needs to be consistent prior to generating element to element. It is not necessarily consistent because mHasFaceConn is set to false.");

            mHasElemToElem                = true;
            moris::size_t tNumFacePerElem = 4;
            Matrix<IndexMat> tElementToElementFacet;
            mElementToElement            = generate_element_to_element(mFaceToElement,
                    mNumElem,
                    tNumFacePerElem,
                    std::numeric_limits<moris::moris_index>::max(),
                    tElementToElementFacet);
        }
        else if(mElementTopology == CellTopology::TRI3)
        {
            MORIS_ASSERT(mHasEdgeConn,
                    "Edge connectivity needs to be consistent prior to generating element to element. It is not necessarily consistent because mHasEdgeConn is set to false.");

            mHasElemToElem                = true;
            moris::size_t tNumEdgePerElem = 3;
            mElementToElement                = generate_element_to_element_2D(mEdgeToElement,
                    mNumElem,
                    tNumEdgePerElem,
                    std::numeric_limits<moris::moris_index>::max());
        }
        else
        {
            MORIS_ASSERT(0, "Element to element connectivity generation requires either TET4 or TRI3 cell topology.");
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::set_up_proc_local_to_cm_local_node_map()
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
                MORIS_ERROR(tBreaker!=0," Attempted to add duplicate node in constructor, are your node indices correct?");
            }
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::add_nodes_to_map(moris::Matrix< moris::IndexMat > const & aNodesToAdd)
    {
        moris::size_t tNumNodes = get_num_entities(EntityRank::NODE);
        moris::size_t tNumNewNodes = aNodesToAdd.n_cols();

        for( moris::size_t i = 0; i<tNumNewNodes; i++)
        {
            MORIS_ASSERT(mNodeIndsToCMInd.find(aNodesToAdd(i)) == mNodeIndsToCMInd.end(),"node already in this child mesh");
            {
                mNodeIndsToCMInd[aNodesToAdd(i)] = i+tNumNodes;
            }
        }
    }

    // ---------------------------------------------------------------------------------
    /*
     * Takes the element face ordinal ancestry and creates a face ancestry
     */

    void
    Child_Mesh::setup_face_ancestry()
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

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::setup_edge_ancestry()
    {
        moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
        moris::size_t tNumEdgePerElement = mElementToEdge.n_cols();
        moris::size_t tNumEdges = get_num_entities(EntityRank::EDGE);

        mEdgeParentInds  = moris::Matrix< moris::IndexMat >(1,tNumEdges);
        mEdgeParentRanks = moris::Matrix< moris::DDSTMat >(1,tNumEdges);

        for( moris::size_t i = 0; i<tNumElements; i++)
        {
            for(moris::size_t j = 0; j<tNumEdgePerElement; j++)
            {
                moris::size_t tEdgeIndex = mElementToEdge(i,j);

                mEdgeParentInds(0,tEdgeIndex) = mElementEdgeParentInds(i,j);
                mEdgeParentRanks(0,tEdgeIndex) = mElementEdgeParentRanks(i,j);
            }
        }
    }

    // ---------------------------------------------------------------------------------
    moris::Matrix< moris::IndexMat >
    Child_Mesh::convert_to_cm_local_indices(moris::Matrix< moris::IndexMat > const & aEntityRankToNode) const
    {
        moris::size_t tNumRows = aEntityRankToNode.n_rows();
        moris::size_t tNumCols = aEntityRankToNode.n_cols();

        moris::Matrix< moris::IndexMat > tLocalEntityRankToNode(tNumRows,tNumCols);

        for(moris::size_t i = 0; i < tNumRows; i++)
        {
            for(moris::size_t j = 0; j < tNumCols; j++)
            {
                auto tIter = mNodeIndsToCMInd.find(aEntityRankToNode(i,j));

                MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

                tLocalEntityRankToNode(i,j) = tIter->second;
            }
        }

        return tLocalEntityRankToNode;
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Child_Mesh::convert_to_proc_indices(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
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

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Child_Mesh::convert_cm_loc_to_glob_ids(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
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

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Child_Mesh::convert_proc_to_glob_ids(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
    {
        moris::size_t tNumRows = aEntityRankToNodeLocal.n_rows();
        moris::size_t tNumCols = aEntityRankToNodeLocal.n_cols();

        moris::Matrix<  moris::IdMat  > tProcEntityRankToNode(tNumRows,tNumCols);

        for(moris::size_t i = 0; i < tNumRows; i++)
        {
            for(moris::size_t j = 0; j < tNumCols; j++)
            {
                moris::uint tNodeCMIndex = this->get_cm_local_node_index(aEntityRankToNodeLocal(i,j));

                tProcEntityRankToNode(i,j) = mNodeIds(0,tNodeCMIndex);
            }
        }

        return tProcEntityRankToNode;
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::modify_child_mesh_internal(enum TemplateType aTemplate)
    {
        if(aTemplate == TemplateType::HIERARCHY_TET4)
        {
            // construct intersection connectivity information
            create_intersect_connectivity();

            // Add new node inheritance information based on edge intersect information
            add_intersect_conn_node_inheritance();

            // Verify the topology before modifying
            MORIS_ASSERT(verify_tet4_topology(
                    this->get_element_to_node(),
                    this->get_element_to_edge(),
                    this->get_element_to_face(),
                    this->get_edge_to_node(),
                    this->get_face_to_node()),
                    "The generated mesh has an invalid topology");
            MORIS_ASSERT(mIntersectConnectivity.n_rows() == mNumElem, "There needs to be a row for each element in the child mesh in the intersect connectivity");

            // Iterate over number of elements (only ones that existed at the beginning of the modification)
            moris::size_t tNumExistElem = mNumElem;

            // Container for all the templates to add
            Cell<Mesh_Modification_Template> tTemplatesToAdd(tNumExistElem);

            // Keep track of the number of intersected elements
            moris::size_t tNumIntersected = 0;
            moris::size_t tNumNewElem     = 0;

            // Get child mesh local connectivities
            moris::Matrix< moris::IndexMat > tEdgeToNodeCMLoc = get_edge_to_node_local();
            moris::Matrix< moris::IndexMat > tElemToNodeCMLoc = get_element_to_node_local();
            moris::Matrix< moris::IndexMat > tElemToEdgeCMLoc = get_element_to_edge();

            for(moris::size_t iE = 0; iE<tNumExistElem; iE++)
            {
                if(mIntersectConnectivity(iE,0) == 3 || mIntersectConnectivity(iE,0) == 4)
                {
                    moris::size_t tPermutationId = std::numeric_limits<moris::size_t>::max();

                    moris::Matrix< moris::IndexMat > tIntersectConnRow = mIntersectConnectivity.get_row(iE);

                    // Sort nodes by edge intersection
                    moris::Matrix< moris::IndexMat > tSortedNodes
                    = sort_nodes(aTemplate, tIntersectConnRow, tEdgeToNodeCMLoc, iE, tPermutationId);

                    // Select more specific tet4 hierarchy template
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
                    moris::Matrix< moris::IndexMat > tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    moris::Matrix< moris::IndexMat > tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentFaceRanks = mElementFaceParentRanks.get_row(iE);

                    // Convert the child mesh local indices to process local indices
                    tSortedNodes = convert_to_proc_indices(tSortedNodes);

                    // Setup template with this information
                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template(
                            tElementsAncestry(0,0),
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

                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template(
                            tElementsAncestry(0,0),
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
                    moris::Matrix< moris::IndexMat > tElementsAncestry({{mParentElementIndex}}); 
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    moris::Matrix< moris::IndexMat > tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentFaceRanks = mElementFaceParentRanks.get_row(iE);

                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template(
                            tElementsAncestry(0,0),
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
                    std::cout << "Invalid connectivity for nodal hierarchy template, (should be 1,2,3 or 4 nodes)\n";
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
            cleanup_intersect_connectivity();

            // Verify the topology before modifying
            MORIS_ASSERT(verify_tet4_topology(
                    this->get_element_to_node(),
                    this->get_element_to_edge(),
                    this->get_element_to_face(),
                    this->get_edge_to_node(),
                    this->get_face_to_node()),
                    "The generated mesh has an invalid topology in modify_child_mesh()");
        }

        else if(aTemplate == TemplateType::REGULAR_SUBDIVISION_HEX8)
        {
            Mesh_Modification_Template tRegSubTemplate(
                    mParentElementIndex,
                    0,
                    mNodeInds,
                    mNodeParentInds,
                    mNodeParentRank,
                    mElementEdgeParentInds,
                    mElementEdgeParentRanks,
                    mElementFaceParentInds,
                    mElementFaceParentRanks,
                    TemplateType::REGULAR_SUBDIVISION_HEX8);

            // Since this template takes HEX8 elements to TET4 elements, we need to resize the
            // connectivity information
            mElementToNode.resize(1,4);
            mElementEdgeParentInds.resize(1,6);
            mElementEdgeParentRanks.resize(1,6);
            mElementFaceParentInds.resize(1,4);
            mElementFaceParentRanks.resize(1,4);

            // Allocate space for new elements (should add 23 elements)
            allocate_more_elements(tRegSubTemplate.mNumNewElem - tRegSubTemplate.mNumElemToReplace);

            // Insert the template in the mesh
            insert_child_mesh_template(tRegSubTemplate);

            mElementTopology = CellTopology::TET4;

            // generate face, edge and element to element connectivity
            generate_connectivities(true,true,true);

            // Verify the topology before modifying
            MORIS_ASSERT(verify_tet4_topology(this->get_element_to_node(),
                    this->get_element_to_edge(),
                    this->get_element_to_face(),
                    this->get_edge_to_node(),
                    this->get_face_to_node()),
                    "The generated mesh has an invalid topology in modify_child_mesh()");

        }
        else if(aTemplate == TemplateType::REGULAR_SUBDIVISION_QUAD4)
        {
            Mesh_Modification_Template tRegSubTemplate(
                    mParentElementIndex,
                    0,
                    mNodeInds,
                    mNodeParentInds,
                    mNodeParentRank,
                    mElementEdgeParentInds,
                    mElementEdgeParentRanks,
                    {{}},
                    {{}},
                    TemplateType::REGULAR_SUBDIVISION_QUAD4);

            // Since this template takes QAUD4 elements to TRI3 elements, we need to resize the
            // connectivity information
            mElementToNode.resize(1,3);
            mElementEdgeParentInds.resize(1,4);
            mElementEdgeParentRanks.resize(1,4);

            // Allocate space for new elements (should add 3 elements)
            allocate_more_elements(tRegSubTemplate.mNumNewElem - tRegSubTemplate.mNumElemToReplace);

            // Insert the template in the mesh
            insert_child_mesh_template(tRegSubTemplate);

            mElementTopology = CellTopology::TRI3;

            // generate face, edge and element to element connectivity
            generate_connectivities(false,true,true);
        }
        else if(aTemplate == TemplateType::CONFORMAL_TRI3)
        {
            // construct intersection connectivity information
            create_intersect_connectivity();

            // Add new node inheritance information based on edge intersect information
            add_intersect_conn_node_inheritance();

            // Iterate over number of elements (only ones that existed at the beginning of the modification)
            moris::size_t tNumExistElem = mNumElem;

            // Container for all the templates to add
            Cell<Mesh_Modification_Template> tTemplatesToAdd(tNumExistElem);

            // Keep track of the number of intersected elements
            moris::size_t tNumIntersected = 0;
            moris::size_t tNumNewElem     = 0;

            for(moris::size_t iE = 0; iE<tNumExistElem; iE++)
            {

                if(mIntersectConnectivity(iE,0) == 2)
                {
                    auto tIntersectConnRow =  mIntersectConnectivity.get_row(iE);

                    // edge 0
                    moris_index tEdge0 = tIntersectConnRow(4);

                    // edge 1
                    moris_index tEdge1 = tIntersectConnRow(5);

                    moris::Matrix< moris::IndexMat > tEdgeOrdinals = get_edge_ordinal_from_element_and_edge_indices(iE,{{tEdge0,tEdge1}});

                    if(tEdgeOrdinals(1) < tEdgeOrdinals(0))
                    {
                        moris_index tSwap = tEdge0;
                        tIntersectConnRow(4) = tIntersectConnRow(5);
                        tIntersectConnRow(5) = tSwap;

                        tSwap = tIntersectConnRow(1);
                        tIntersectConnRow(1) = tIntersectConnRow(2);
                        tIntersectConnRow(2) = tSwap;
                    }

                    moris::size_t tPermutationId = tEdgeOrdinals(0) + tEdgeOrdinals(1);

                    // Get parent element information
                    moris::Matrix< moris::IndexMat > tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);

                    // get nodes on this (proc indices)
                    moris::Matrix<moris::IndexMat> tNodesInTemplate(1,5);
                    tNodesInTemplate({0,0},{0,2}) = this->get_element_to_node().get_row(iE);
                    tNodesInTemplate(3) = mNodeInds(tIntersectConnRow(1));
                    tNodesInTemplate(4) = mNodeInds(tIntersectConnRow(2));

                    // Setup template with this information
                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template
                            (tElementsAncestry(0,0),
                                    iE,
                                    tNodesInTemplate,
                                    tParentEdgeInds,
                                    tParentEdgeRanks,
                                    {{}},
                                    {{}},
                                    TemplateType::CONFORMAL_TRI3,
                                    tPermutationId);

                    // Increment the count of number of intersected elements and number of new elements
                    tNumNewElem = tNumNewElem + tTemplatesToAdd(tNumIntersected).mNumNewElem - tTemplatesToAdd(tNumIntersected).mNumElemToReplace;
                    tNumIntersected++;

                }
                else if(mIntersectConnectivity(iE,0) == 1)
                {
                    auto tIntersectConnRow =  mIntersectConnectivity.get_row(iE);

                    // edge 0
                    moris_index tEdge0 = tIntersectConnRow(4);

                    moris::Matrix< moris::IndexMat > tEdgeOrdinals = this->get_edge_ordinal_from_element_and_edge_indices(iE,{{tEdge0}});

                    moris::size_t tPermutationId = tEdgeOrdinals(0)+10;

                    // Get parent element information
                    moris::Matrix< moris::IndexMat > tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);

                    // get nodes on this (proc indices)
                    moris::Matrix<moris::IndexMat> tNodesInTemplate(1,4);
                    tNodesInTemplate({0,0},{0,2}) = this->get_element_to_node().get_row(iE);
                    tNodesInTemplate(3) = mNodeInds(tIntersectConnRow(1));

                    // Setup template with this information
                    tTemplatesToAdd(tNumIntersected) = Mesh_Modification_Template
                            (tElementsAncestry(0,0),
                                    iE,
                                    tNodesInTemplate,
                                    tParentEdgeInds,
                                    tParentEdgeRanks,
                                    {{}},
                                    {{}},
                                    TemplateType::CONFORMAL_TRI3,
                                    tPermutationId);

                    // Increment the count of number of intersected elements and number of new elements
                    tNumNewElem = tNumNewElem + tTemplatesToAdd(tNumIntersected).mNumNewElem - tTemplatesToAdd(tNumIntersected).mNumElemToReplace;
                    tNumIntersected++;
                    
                }
                else if(mIntersectConnectivity(iE,0) == 0)
                {
                    continue;
                }
                else
                {
                    std::cout<<"mIntersectConnectivity(iE,0) = "<<mIntersectConnectivity(iE,0)<<std::endl;
                    std::cout<<"Parent Cell Index = "<<this->get_parent_element_index()<<std::endl;
                    
                    Matrix<IndexMat> tNodeIds = this->get_node_ids();
                    moris::print(tNodeIds,"Node Ids");
                    MORIS_ERROR(0,"Unsupported case in 2D conformal tet template");
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
            generate_connectivities(false,true,true);

            // Clear the intersection connectivity
            cleanup_intersect_connectivity();

        }
        else
        {
            MORIS_ASSERT(0, "Template type not currently supported in modify_child_mesh.");
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::create_intersect_connectivity()
    {
        // Get entities of dimension edge connected to elements
        uint tNumIntersectedEdges = mIntersectedEdges.size();

        MORIS_ASSERT(tNumIntersectedEdges == mIntersectedCMNodeIndex.size(),
                "Dimension mismatch between cell of intersected edges and cell of intersected node indices");

        // Allocate Intersect connectivity matrix
        moris::size_t tNumElements = get_num_entities(EntityRank::ELEMENT);
        moris::size_t tNumEdgesToElem = get_element_to_edge().n_cols();
        mIntersectConnectivity = moris::Matrix< moris::IndexMat >(tNumElements, tNumEdgesToElem * 2 + 1, 0); // Needs to be zero for use column

        moris::Matrix< moris::IndexMat > const tEdgeToElem = get_edge_to_element();

        for(moris::size_t iIE = 0; iIE <tNumIntersectedEdges; iIE++)
        {
            moris::moris_index tCMNodeInd = mIntersectedCMNodeIndex(iIE);
            moris::moris_index tCMEdgeInd = mIntersectedEdges(iIE);

            moris::size_t tNumElemsConnected = 0;
            for(moris::size_t i = 0; i<tEdgeToElem.n_cols(); i++)
            {
                if(tEdgeToElem(tCMEdgeInd, i)== std::numeric_limits<moris::moris_index>::max())
                {
                    break;
                }
                tNumElemsConnected++;
            }

            for(moris::size_t i = 0; i < tNumElemsConnected; i++)
            {
                moris::size_t tElemCMInd = tEdgeToElem(tCMEdgeInd, i);

                MORIS_ASSERT(mIntersectConnectivity(tElemCMInd, 0) < (moris::moris_index)this->get_element_to_edge().n_cols(),
                        "Entity corresponding to provided aDInd has exceeded allocated space");
                MORIS_ASSERT(tElemCMInd < mIntersectConnectivity.n_rows(),
                        "aDInd is outside of bounds. Has auxiliary connectivity been initialized?");

                mIntersectConnectivity(tElemCMInd, mIntersectConnectivity(tElemCMInd, 0) + 1) = tCMNodeInd;
                mIntersectConnectivity(tElemCMInd, mIntersectConnectivity(tElemCMInd, 0) + mElementToEdge.n_cols() + 1) = tCMEdgeInd;
                mIntersectConnectivity(tElemCMInd, 0)++;
            }
        }
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Child_Mesh::sort_nodes(
            enum TemplateType                        aTemplate,
            moris::Matrix< moris::IndexMat > const & aIntConnectivity,
            moris::Matrix< moris::IndexMat > const & aEdgeToNodeCMLoc,
            moris::size_t const &                      aElementIndex,
            moris::size_t &                            aPermutation)
    {
        //Locate highest node in intersection connectivity
        switch(aTemplate)
        {
            case TemplateType::HIERARCHY_TET4:
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
                        std::cout << "Duplicate node not found, invalid edge intersection configuration";

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
                        std::cout << "Node 2 not found, invalid edge intersection configuration";

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
                            // if the edge h and
                            if(tNodesMH(0, i) == tNodesH(0, 0))
                            {
                                tSortedNodes(0, 1) = tNodesH(0, 1);  // High nodes independent node
                                tSortedNodes(0, 0) = tNodesMH(0, i); // Shared Node
                                tSortedNodes(0, 2) = tNodesMH(0, j); // Mid highs ind node
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

                        MORIS_ERROR(tSuccess == 1, "Sorting to find nodes 1,2,3 unsuccessful");

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

                        MORIS_ERROR(tSuccess == 1, "Sorting to find node 4 unsuccessful");

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
                        MORIS_ERROR(tSuccess == 1, "Sorting to find nodes 1,2,3 unsuccessful");

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
                        MORIS_ASSERT(tSuccess = 1, "Sorting to find node 4 unsuccessful");

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
                        MORIS_ERROR(tSuccess == 1, "Sorting to find nodes 1,2,3 unsuccessful");

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
                        MORIS_ASSERT(tSuccess = 1, "Sorting to find node 4 unsuccessful");
                    }
                    else
                        std::cout << "Sorting Failed (invalid flagging). Did a node appear twice?";
                    return tSortedNodes;
                }

                else
                {
                    std::cout << "SORTING NOT COMPLETED! Check to see if this function is called for a non-intersected element";
                    moris::Matrix< moris::IndexMat > dummy(1, 1);
                    return dummy;
                }

                break;
            }
            default:
            {
                std::cout << "Sorting for specified template type not implemented";

                moris::Matrix< moris::IndexMat > dummy(1, 1);
                return dummy;

                break;
            }
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::get_intersection_permutation(
            moris::Matrix< moris::IndexMat > const & aOrderedEdgeOrdinals,
            moris::size_t                          & aPermutation)
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
            std::cout<<"Permutation rule not implemented";
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::reindex_template_parent_information(Mesh_Modification_Template & aMeshModTemplate)
    {
        // Reindex parent ordinals for edges
        moris::size_t tNumEdgePerElem = aMeshModTemplate.mNewParentEdgeOrdinals.n_cols();
        moris::size_t tNumFacePerElem = aMeshModTemplate.mNewParentFaceOrdinals.n_cols();
        moris::size_t tNumElems       = aMeshModTemplate.mNumNewElem;
        moris::size_t tNumNodes       = aMeshModTemplate.mNodeInds.numel();
        moris::Matrix< moris::IndexMat > tParentElem({{aMeshModTemplate.mParentElemInd}});
        moris::Matrix< moris::DDSTMat > tParentElemRank({{3}});

        // Place ptrs in cell to avoid if statement checking rank in for loop
        Cell<moris::Matrix< moris::IndexMat >*> tParentEntitiesInds({& aMeshModTemplate.mParentNodeInds,
            & aMeshModTemplate.mParentEdgeInds,
            & aMeshModTemplate.mParentFaceInds,
            & tParentElem});

        Cell<moris::Matrix< moris::DDSTMat >*> tParentEntitiesRanks({& aMeshModTemplate.mParentNodeRanks,
            & aMeshModTemplate.mParentEdgeRanks,
            & aMeshModTemplate.mParentFaceRanks,
            & tParentElemRank});

        for(moris::size_t i = 0; i<tNumElems; i++)
        {
            // Reindex edges
            for(moris::size_t j = 0; j<tNumEdgePerElem; j++)
            {
                moris::size_t tEdgeParentRank = aMeshModTemplate.mNewParentEdgeRanks(i,j);
                aMeshModTemplate.mNewParentEdgeRanks(i,j) = (*tParentEntitiesRanks(tEdgeParentRank))(0,aMeshModTemplate.mNewParentEdgeOrdinals(i,j));
                aMeshModTemplate.mNewParentEdgeOrdinals(i,j) = (*tParentEntitiesInds(tEdgeParentRank))(0,aMeshModTemplate.mNewParentEdgeOrdinals(i,j));
            }

            // Reindex faces
            for(moris::size_t j = 0; j<tNumFacePerElem; j++)
            {
                moris::size_t tFaceParentRank = aMeshModTemplate.mNewParentFaceRanks(i,j);
                aMeshModTemplate.mNewParentFaceRanks(i,j) = (*tParentEntitiesRanks(tFaceParentRank))(0,aMeshModTemplate.mNewParentFaceOrdinals(i,j));
                aMeshModTemplate.mNewParentFaceOrdinals(i,j) = (*tParentEntitiesInds(tFaceParentRank))(0,aMeshModTemplate.mNewParentFaceOrdinals(i,j));
            }
        }

        // reindex node ancestry if there is any
        if(aMeshModTemplate.mHasNodeInheritance)
        {
            for( moris::size_t i = 0; i<tNumNodes; i++)
            {
                moris::size_t tNodeParentRank = aMeshModTemplate.mNewNodeParentRanks(i);
                moris::size_t tNodeParentOrd  = aMeshModTemplate.mNewNodeParentOrdinals(i);
                aMeshModTemplate.mNewNodeParentRanks(i)    = (*tParentEntitiesRanks(tNodeParentRank))(0,tNodeParentOrd);
                aMeshModTemplate.mNewNodeParentOrdinals(i) = (*tParentEntitiesInds(tNodeParentRank))(0,tNodeParentOrd);
            }
        }

        aMeshModTemplate.mIsReindexed = true;
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::cleanup_intersect_connectivity()
    {
        mIntersectConnectivity.set_size(0,0);
        mIntersectedEdges.clear();
        mIntersectedCMNodeIndex.clear();
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::add_intersect_conn_node_inheritance()
    {
        // Get entities of dimension edge connected to elements
        uint tNumIntersectedEdges = mIntersectedEdges.size();

        MORIS_ASSERT(tNumIntersectedEdges == mIntersectedCMNodeIndex.size(),
                "Dimension mismatch between cell of intersected edges and cell of intersected node indices");

        for( moris::uint i = 0 ; i < tNumIntersectedEdges; i++)
        {
            moris::moris_index tEdgeParentInd  = get_entity_parent_entity_proc_ind(EntityRank::EDGE,mIntersectedEdges(i));
            moris::moris_index tEdgeParentRank = get_entity_parent_entity_rank(EntityRank::EDGE,mIntersectedEdges(i));

            set_node_inheritance(mIntersectedCMNodeIndex(i),tEdgeParentInd,tEdgeParentRank);
        }
    }

    // ---------------------------------------------------------------------------------

    void
    Child_Mesh::construct_subphase_bins(
            moris::moris_index                     & aSubPhaseIndex,
            moris::Matrix< moris::IndexMat > const & aElementSubPhase)
    {
        // Number of bins corresponds to the maxmimum value in the element sub phase vector
        moris::size_t tNumBins = aElementSubPhase.max() + 1;
        moris::size_t tNumElements = this->get_num_entities(EntityRank::ELEMENT);

        // Initialize member variables
        // Element Sub-phase bins
        mBinBulkPhase = Cell<moris::moris_index>(tNumBins);
        mElementBinIndex = aElementSubPhase.copy();

        moris::Matrix< moris::DDSTMat > tBinSizeCounter(1,tNumBins,0);
        mSubPhaseBinId                = moris::Matrix<moris::IndexMat>(1,tNumBins);
        mSubPhaseBinIndices           = Cell<moris::moris_index>(tNumBins);
        mSubPhaseBins                 = Cell<moris::Matrix< moris::IndexMat >>(tNumBins);
        mSubphaseBasisEnrichmentLevel = Cell<Cell< moris_index >>(tNumBins);
        mSubphaseBasisIndices         = Cell<Cell< moris_index >>(tNumBins);

        for(moris::size_t i = 0; i<tNumBins; i++)
        {
            mSubPhaseBins(i) = moris::Matrix< moris::IndexMat >(1,tNumElements);

            if(i == 0)
            {
                mSubPhaseBinIndices(0) = mParentElementIndex;
            }
            else
            {
                mSubPhaseBinIndices(i) = aSubPhaseIndex;
                aSubPhaseIndex++;
            }
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
    // ---------------------------------------------------------------------------------

    void Child_Mesh::identify_hanging_nodes( const moris::Cell< moris_index > & aTransitionFacetIndices )
    {
        moris::Cell< moris_index> tHangingNodeIndices;

        for( uint Ik = 0; Ik < aTransitionFacetIndices.size(); Ik ++ )
        {
            moris_index tFacetIndex = aTransitionFacetIndices( Ik );

            for( uint Ii = 0; Ii < mNodeInds.numel(); Ii++ )
            {
                if( (this->get_entity_parent_entity_rank( EntityRank::NODE, Ii ) == (moris_index) this->get_facet_rank() ) &&
                    ( this->get_entity_parent_entity_proc_ind( EntityRank::NODE, Ii ) == tFacetIndex ) )
                {
                    tHangingNodeIndices.push_back( mNodeInds( Ii ) );

                    break;
                }

            }
        }

        mHangingNodes.set_size( tHangingNodeIndices.size(), 1 );

        for( uint Ii = 0; Ii < tHangingNodeIndices.size(); Ii++ )
        {
            mHangingNodes( Ii ) = tHangingNodeIndices( Ii );
        }
    }

    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------
}

