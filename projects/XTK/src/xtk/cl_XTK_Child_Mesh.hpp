/*
 * cl_XTK_Child_Mesh_Test.hpp
 *
 *  Created on: Jun 21, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_CHILD_MESH_HPP_
#define SRC_XTK_CL_XTK_CHILD_MESH_HPP_
#include <unordered_map>


// Linear algebra includes
#include "cl_Matrix.hpp"
#include "fn_isvector.hpp"
#include "fn_iscol.hpp"
#include "fn_trans.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_Cell.hpp"

#include "cl_XTK_Enums.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_XTK_Output_Options.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_create_faces_from_element_to_node.hpp"
#include "fn_create_edges_from_element_to_node.hpp"
#include "cl_XTK_Child_Mesh_Modification_Template.hpp"


// MTK includes
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_verify_tet_topology.hpp"

#include "assert.hpp"

using namespace moris;

namespace xtk
{

class Child_Mesh
{
public:

    Child_Mesh();

    Child_Mesh(moris::moris_index                      aParentElementIndex,
                    moris::Matrix< moris::IndexMat > & aNodeInds,
                    moris::Matrix< moris::IndexMat > & aElementNodeParentInds,
                    moris::Matrix< moris::DDSTMat >  & aElementNodeParentRanks,
                    moris::Matrix< moris::IndexMat > & aElementToNode,
                    moris::Matrix< moris::IndexMat > & aElementEdgeParentInds,
                    moris::Matrix< moris::DDSTMat >  & aElementEdgeParentRanks,
                    moris::Matrix< moris::IndexMat > & aElementFaceParentInds,
                    moris::Matrix< moris::DDSTMat >  & aElementFaceParentRanks,
                    moris::Matrix< moris::DDSTMat >  & aElementInferfaceSides );


    Child_Mesh(Mesh_Modification_Template & aMeshModTemplate);

    ~Child_Mesh(){}


    // Declare iterator type
    typename std::unordered_map<moris::size_t,moris::size_t>::iterator Map_Iterator;


    // --------------------------------------------------------------
    // Functions to access connectivity
    // --------------------------------------------------------------

    /*
     * Get number of a entity of a given rank
     */
    moris::size_t
    get_num_entities(enum EntityRank aEntityRank) const;

    /*
     * return element to node connectivity matrix (processor local indices)
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_node() const;
    /*
     * return element to node connectivity matrix (processor local indices)
     */
    moris::Matrix< moris::IdMat >
    get_element_to_node_glob_ids(moris::moris_index aCMElemIndex) const;
    /*
     * Compute and return the node to element connectivity. WARNING:
     * this computes the connectivity so use it sparingly.
     */
    moris::Matrix< moris::IndexMat >
    get_node_to_element_local() const;

    moris::moris_index
    get_element_node_ordinal(moris::moris_index aCMLocElemIndex,
                             moris::moris_index aProcLocNodeIndex);


    /*
     * Converts the existing element to node connectivity (which contains proc local indices)
     * to child mesh local indices
     */
    moris::Matrix< moris::IndexMat >
    get_element_to_node_local() const;


    /*
     * Converts the existing element to node connectivity (which contains proc local indices)
     * to child mesh local indices
     */
    moris::Matrix< moris::IdMat >
    get_element_to_node_global() const;

    /*
     * Return edge to node
     */
    moris::Matrix< moris::IndexMat > &
    get_edge_to_node();

    moris::Matrix< moris::IndexMat > const &
    get_edge_to_node() const;

    /*
     * Return edge to node
     */
    moris::Matrix< moris::IndexMat >
    get_edge_to_node_local() const;

    /*
     * Return element to edge
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_edge() const;

    /*
     * Return edges connected to elements
     */
    moris::Matrix< moris::IndexMat > const &
    get_edge_to_element() const;

    /*
     * return faces to node
     */
    moris::Matrix< moris::IndexMat > const &
    get_face_to_node() const;

    moris::Matrix< moris::IndexMat > const &
    get_face_to_element() const;

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    moris::Matrix< moris::IndexMat >
    get_face_to_node_local() const;

    /*
     * return element to face connectivity matrix (cm local element indices)
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_face() const
    {
        MORIS_ASSERT(mHasFaceConn,"Face connectivity has not been generated with call to generate_face_connectivity");
        return mElementToFace;
    }

    /*
     * Return element to element connectivity (cm local element indices)
     */
    moris::Matrix< moris::IndexMat > const &
    get_element_to_element() const;

    /*
     * Return the parametric coordinates of the nodes (ordered by cm local index)
     */
    moris::Matrix< moris::DDRMat > const &
    get_parametric_coordinates() const;

    /*
     * Return the parametric coordinates of a node using a node index
     */
    moris::Matrix< moris::DDRMat >
    get_parametric_coordinates(moris::moris_index aNodeIndex) const;


    /*!
     *
     */
    moris::mtk::Geometry_Type
    get_child_geometry_type() const;


    // Functions to access ancestry

    /*!
     * Get and entities parent entity rank
     */
    moris::moris_index
    get_entity_parent_entity_rank(enum EntityRank aEntityRank,
                                  moris::moris_index aCMLocIndex);

    /*!
     * Get and entities parent entity rank
     */
    moris::moris_index
    get_entity_parent_entity_proc_ind(enum EntityRank  aEntityRank,
                                      moris::moris_index aCMLocIndex);


    moris::moris_index
    get_parent_element_index() const;

    moris::Matrix< moris::IndexMat > const &
    get_face_parent_inds() const;

    moris::Matrix< moris::DDSTMat > const &
    get_face_parent_ranks() const;

    moris::Matrix< moris::IndexMat > const &
    get_edge_parent_inds() const;

    moris::Matrix< moris::DDSTMat > const &
    get_edge_parent_ranks() const;

    moris::Matrix< moris::IndexMat > const &
    get_node_indices() const;


    moris::Matrix< moris::IdMat > const &
    get_node_ids() const;

    moris::Matrix< moris::IdMat > const &
    get_element_ids() const;


    moris::Matrix< moris::IndexMat > const &
    get_element_inds() const;


    moris::Matrix<moris::IndexMat>
    convert_to_proc_local_elem_inds(moris::Matrix< moris::IndexMat > aConnOfElementCMIndices);

    // Function to access ordinals
    moris::Matrix< moris::IndexMat >
    get_edge_ordinal_from_element_and_edge_indices(moris::moris_index const & aElementIndex,
                                                   moris::Matrix< moris::IndexMat > const & aEdgeIndices) const;

    // Function to access ordinals
    moris::moris_index
    get_edge_ordinal_from_element_and_edge_indices(moris::moris_index const & aElementIndex,
                                                   moris::moris_index const & aEdgeIndex) const;


    // Function to access ordinals
    moris::moris_index
    get_face_ordinal_from_element_and_face_index(moris::moris_index const & aElementIndex,
                                                 moris::moris_index         aFaceIndex) const;
    /*
     * Returns the child element and face ordinal connected to a provided parent face
     */
    void
    get_child_elements_connected_to_parent_face(moris::moris_index         const & aParentFaceIndex,
                                                moris::Matrix< moris::IdMat >    & aChildElemsIdsOnFace,
                                                moris::Matrix< moris::IndexMat > & aChildElemsCMIndOnFace,
                                                moris::Matrix< moris::IndexMat > & aChildElemOnFaceOrdinal) const;

    // --------------------------------------------------------------
    // Functions to modify the mesh
    // --------------------------------------------------------------

    /*
     * Modify child mesh by selecting template using Intersection connectivity
     * to determine correct template and insert the template
     */
    void
    modify_child_mesh(enum TemplateType aTemplate);


    void
    initialize_unzipping();

    void
    finalize_unzipping();


    /*
     * Returns the child mesh local element index not the processor local index.
     */
    void
    unzip_child_mesh_interface_get_interface_element_pairs(moris::uint aGeometryIndex,
                                                           bool & aNoPairFoundFlag,
                                                           moris::Matrix<moris::IndexMat> & aInterfaceElementPairs,
                                                           moris::Matrix<moris::IndexMat> & aInterfacePairSideOrds);

    void
    unzip_child_mesh_interface(moris::moris_index                       aGeometryIndex,
                               moris::Matrix< moris::IndexMat > const & aInterfaceElementPairsCMIndex,
                               moris::Matrix< moris::IndexMat > const & aElementUsingZippedNodes,
                               moris::Matrix< moris::IndexMat > const & aInterfaceNodeIndices,
                               moris::Matrix< moris::IndexMat > const & aUnzippedInterfaceNodeIndices,
                               moris::Matrix< moris::IndexMat > const & aUnzippedInterfaceNodeIds);


    void
    convert_tet4_to_tet10_child();

    /*
     * Add more nodes to the existing node indices list
     */
    void
    add_node_indices(moris::Matrix< moris::IndexMat > const & aNewNodeInds);

    /*
     * Add more nodes to the existing node indices list
     */
    void
    add_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds);


    /*
     * sets the node ids. note this overwrites existing mNodeIds data
     */

    void
    set_node_ids(moris::Matrix< moris::IdMat > const & aNewNodeIds);

    /*
     * Sets the globally unique element Ids for the child mesh. This is important for mesh from data creation
     * @param[in] aElementId - First element Id (note: this is incremented as the id is assigned)
     */
    void set_child_element_ids(moris::moris_id & aElementId);

    /*
     * Sets the processor unique element indices for the child mesh.
     * @param[in] aElementInd - First element Ind (note: this is incremented as the ind is assigned)
     */
    void set_child_element_inds(moris::moris_index & aElementInd);

    /*!
     * Add node parametric coordinate for a single node
     */
    void
    add_node_parametric_coordinate( moris::size_t aNodeIndex,
                                    moris::Matrix< moris::DDRMat > const & aParamCoord );

    /*!
     * Add node parametric coordinate for multiple nodes
     * @param[in] aNodeIndices - List of node indices with parametric coordinates
     * @param[in] aParamCoord  - Parametric coordinates for nodes (note: row 1 of aParamCoord is the Parametric coordinate of aNodeIndices(1))
     */
    void
    add_node_parametric_coordinate( moris::Matrix< moris::IndexMat> const & aNodeIndices,
                                    moris::Matrix< moris::DDRMat >  const & aParamCoord );

    void
    allocate_parametric_coordinates( moris::size_t aNumNewNodes,
                                     moris::size_t aDimOfParmCoord);

    /*
     * Resizes element related matrices
     */
    void
    allocate_more_elements(moris::size_t const & aNumMoreElements);

    void
    insert_child_mesh_template(Mesh_Modification_Template & aMeshModTemplate);

    /*
     * Convert the element to node, parent information to  the local indexing scheme rather than ordinal
     */
    void reindex_template(Mesh_Modification_Template & aMeshModTemplate);

    /*
     * Replace the information associated with an element
     */
    void
    replace_element(moris::size_t                    const & aElementIndexToReplace,
                    moris::size_t                    const & aRowIndex,
                    moris::Matrix< moris::IndexMat > const & aElementToNode,
                    moris::Matrix< moris::IndexMat > const & aElementEdgeParentInds,
                    moris::Matrix< moris::DDSTMat >  const & aElementEdgeParentRanks,
                    moris::Matrix< moris::IndexMat > const & aElementFaceParentInds,
                    moris::Matrix< moris::DDSTMat >  const & aElementFaceParentRanks,
                    moris::Matrix< moris::DDSTMat >  const & aElementInterfaceFaces);

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
                moris::Matrix< moris::DDSTMat> const & aElementInterfaceFaces);


    /*!
     * Add node inheritance
     * @param[in] aCMLocNodeIndex - Child mesh local node index
     * @param[in] aProcLocParentEntityInd - Processor local parent entity index
     * @param[in] aParentEntityRank - Parent entity rank
     */
    void
    set_node_inheritance(moris::moris_index aCMLocNodeIndex,
                         moris::moris_index aProcLocParentEntityInd,
                         moris::moris_index aParentEntityRank);


    /*
     * Generates face connectivities, edge connectivities, element to element graph
     */

    void
    generate_connectivities(bool aGenerateFaceConn,
                            bool aGenerateEdgeConn,
                            bool aGenerateElemToElem);




    /**
      * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
      *       - 1 means the provided aDPrime1Ind is an XTK index
      *
      * aDPrime2Ind must be XTK local index
      */
     void
     add_entity_to_intersect_connectivity(moris::moris_index aCMNodeInd,
                                          moris::moris_index aCMEdgeInd,
                                          moris::size_t aFlag);

     void
     mark_edge_as_on_interface(moris::moris_index aEdgeIndex);


      /**
       * Tells the child mesh where a node index will be placed once it has been communicated
       */
      void set_pending_node_index_pointers(Cell<moris::moris_index*>            aNodeIndPtr,
                                           moris::Matrix< moris::DDRMat > const & aNodeParamCoordinate);

      /**
       *  retrieve_pending_node_inds
       * XTK mesh retrieves pending node indices. Prior to this call the unique node assignments
       * must be complete via the request structure in XTK model
       */
      void retrieve_pending_node_inds();


      /*
       * Take the information of edges on the interface and
       * figure out which faces are on the interface, then
       * mark element edges as on interface.
       */
      void
      mark_interface_faces_from_interface_coincident_faces();


      // --------------------------------------------------------------
      // Functions for elemental phase information and subphase bins
      // --------------------------------------------------------------

      void
      initialize_element_phase_mat();


      void
      set_element_phase_index(moris::size_t aEntityIndex,
                              moris::size_t aEntityPhaseIndex);

      moris::size_t
      get_element_phase_index( moris::size_t const & aEntityIndex) const;


      moris::Matrix< moris::IndexMat > const &
      get_element_phase_indices() const
      {
          MORIS_ASSERT(mHasPhaseInfo,"Elemental phase information not set");
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
      pack_interface_sides( bool aIndexFlag = false,
                            moris::moris_index aPhaseIndex = MORIS_INDEX_MAX) const
      {
          // Loop bound and sizing
          moris::size_t tNumElem = get_num_entities(EntityRank::ELEMENT);

          moris::Matrix< moris::IdMat > tInterfaceSideSetInfo(tNumElem,2);

          // Keep track of the number of interface sides
          moris::size_t tCount = 0;

         // Iterate over each element and if the element has an interface side it will be in mElementInferfaceSides vector
         for(moris::size_t iEl =0 ; iEl<tNumElem; iEl++)
         {

             if(mElementInterfaceSides(iEl) != std::numeric_limits<moris::size_t>::max())
             {
                 if((moris_index)this->get_element_phase_index(iEl) == aPhaseIndex || aPhaseIndex == MORIS_INDEX_MAX)
                 {
                     if( aIndexFlag )
                     {
                         tInterfaceSideSetInfo(tCount,0) = mChildElementInds(iEl);
                     }
                     else
                     {
                         tInterfaceSideSetInfo(tCount,0) = mChildElementIds(iEl);
                     }
                     tInterfaceSideSetInfo(tCount,1) = mElementInterfaceSides(iEl);
                     tCount++;
                 }
             }
         }

         // Size out space
         tInterfaceSideSetInfo.resize(tCount,2);

         return tInterfaceSideSetInfo;
      }

      moris::Matrix< moris::IdMat >
      pack_interface_sides_loc_inds() const
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

private:
    // Parent element index
    moris::moris_index              mParentElementIndex;

    // Element To Node and Ancestry Information (This is the only data that is set with templates.
    // all other is generated with an algorithm)
    // All node connectivity is indexed by proc local indexs
    enum CellTopology              mElementTopology;
    moris::size_t                    mNumElem;
    moris::Matrix< moris::IndexMat > mElementToNode; /* proc inds*/
    moris::Matrix< moris::IndexMat > mElementEdgeParentInds;
    moris::Matrix< moris::DDSTMat  > mElementEdgeParentRanks;
    moris::Matrix< moris::IndexMat > mElementFaceParentInds;
    moris::Matrix< moris::DDSTMat  > mElementFaceParentRanks;
    moris::Matrix< moris::DDSTMat  > mElementInterfaceSides;

    // Child element information ---------------------------
    moris::Matrix< moris::IdMat >    mChildElementIds;
    moris::Matrix< moris::IndexMat > mChildElementInds;
    // No child element parent information is stored because
    // it is assumed that all are the parent element of this
    // child mesh.

    // Node information ------------------------------------
    moris::Matrix< moris::IdMat >   mNodeIds;
    moris::Matrix< moris::IndexMat> mNodeInds;
    moris::Matrix< moris::DDSTMat > mNodeParentRank;
    moris::Matrix< moris::IndexMat> mNodeParentInds;

    // Map where  key - proc local ind, val - local child mesh index
    std::unordered_map<moris::size_t, moris::size_t> mNodeIndsToCMInd;

    // Parametric coordinate relative to parent element
    moris::Matrix< moris::DDRMat >   mNodeParametricCoord;

    // Face Connectivity -----------------------------------
    bool mHasFaceConn;
    moris::Matrix< moris::IndexMat > mFaceToNode;/* proc inds*/
    moris::Matrix< moris::IndexMat > mNodeToFace;
    moris::Matrix< moris::IndexMat > mFaceToElement;
    moris::Matrix< moris::IndexMat > mElementToFace;
    moris::Matrix< moris::IndexMat > mFaceParentInds;
    moris::Matrix< moris::DDSTMat >  mFaceParentRanks;

    // Edge connectivity -----------------------------------
    bool mHasEdgeConn;
    moris::Matrix< moris::IndexMat > mEdgeToNode;/* proc inds*/
    moris::Matrix< moris::IndexMat > mNodeToEdge;
    moris::Matrix< moris::IndexMat > mEdgeToElement;
    moris::Matrix< moris::IndexMat > mElementToEdge;
    moris::Matrix< moris::IndexMat > mEdgeParentInds;
    moris::Matrix< moris::DDSTMat >  mEdgeParentRanks;

    // Element to Element graph ----------------------------
    bool mHasElemToElem;
    moris::Matrix< moris::IndexMat > mElementToElement;

    // Auxiliary connectivity data and pending nodes (mesh modification data)
    moris::Matrix<moris::IndexMat>  mIntersectConnectivity;
    moris::Cell< moris_index >      mIntersectedCMNodeIndex;
    moris::Cell< moris_index >      mIntersectedEdges;

    bool                             mHasCoincidentEdges;
    moris::Matrix< moris::IndexMat > mEdgeOnInterface;
    Cell<moris::moris_index*>        mPtrPendingNodeIndex;
    moris::Matrix< moris::DDRMat >   mPendingParamCoordinates;

    // Phase member variables -----------------------------
    bool                                   mHasPhaseInfo;
    moris::Matrix< moris::IndexMat >       mElementPhaseIndices;
    moris::Matrix< moris::IndexMat >       mElementBinIndex;
    Cell<moris::moris_index>               mBinBulkPhase;
    Cell<moris::Matrix< moris::IndexMat >> mSubPhaseBins;

    // Unzipping information
    bool mUnzippingFlag = false;

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

        mEdgeOnInterface.resize(1,mEdgeToNode.n_rows());


        setup_edge_ancestry();

    }

    void
    generate_element_to_element_connectivity()
    {
        MORIS_ASSERT(mHasFaceConn,
                     "Face connectivity needs to be consistent prior to generating element to element. It is not necessarily consistent because mHasFaceConn is set to false");
        mHasElemToElem          = true;
        moris::size_t tNumFacePerElem = 4;
        mElementToElement      = generate_element_to_element(mFaceToElement,
                                                              mNumElem,
                                                              tNumFacePerElem,
                                                              std::numeric_limits<moris::moris_index>::max());
    }

    void
    set_up_proc_local_to_cm_local_node_map()
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

            MORIS_ASSERT(mNodeIndsToCMInd.find(aNodesToAdd(i)) == mNodeIndsToCMInd.end(),"node already in this child mesh");
            {
                mNodeIndsToCMInd[aNodesToAdd(i)] = i+tNumNodes;
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
    convert_to_cm_local_indices(moris::Matrix< moris::IndexMat > const & aEntityRankToNode) const
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

    moris::moris_index
    get_cm_local_node_index(moris::moris_index aNodeProcIndex) const
    {
        auto tIter = mNodeIndsToCMInd.find(aNodeProcIndex);

        MORIS_ASSERT(tIter != mNodeIndsToCMInd.end(),"Node not in map, conversion to local indices failed");

        return tIter->second;
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
    convert_cm_loc_to_glob_ids(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
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

    moris::Matrix< moris::IdMat >
    convert_proc_to_glob_ids(moris::Matrix< moris::IndexMat > const & aEntityRankToNodeLocal) const
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


    void
    modify_child_mesh_internal(enum TemplateType aTemplate)
    {
        if(aTemplate == TemplateType::HIERARCHY_TET4)
        {
            // construct intersection connectivity information
            create_intersect_connectivity();

            // Add new node inheritance information based on edge intersect information
            add_intersect_conn_node_inheritance();

            // Verify the topology before modifying
            MORIS_ASSERT(verify_tet4_topology(this->get_element_to_node(),
                                              this->get_element_to_edge(),
                                              this->get_element_to_face(),
                                              this->get_edge_to_node(),
                                              this->get_face_to_node()),
                                              "The generated mesh has an invalid topology");
            MORIS_ASSERT(mIntersectConnectivity.n_rows() == mNumElem, "There needs to be a row for each element in the child mesh in the intersect connectivity");
            MORIS_ASSERT(aTemplate == TemplateType::HIERARCHY_TET4,"This function needs to be abstracted and has only been tested with HIERARCHY_TET4.");

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
                    moris::Matrix< moris::IndexMat > tElementsAncestry({{mParentElementIndex}}); // Not used
                    moris::Matrix< moris::IndexMat > tParentEdgeInds  = mElementEdgeParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks = mElementEdgeParentRanks.get_row(iE);
                    moris::Matrix< moris::IndexMat > tParentFaceInds  = mElementFaceParentInds.get_row(iE);
                    moris::Matrix< moris::DDSTMat >  tParentFaceRanks = mElementFaceParentRanks.get_row(iE);


                    // Convert the child mesh local indices to process local indices
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
                    std::cout << "Invalid connectivity for nodal hierarchy template, (should be 3 or 4 nodes)\n";
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

        }

        else if(aTemplate == TemplateType::REGULAR_SUBDIVISION_HEX8)
        {

            Mesh_Modification_Template tRegSubTemplate(mParentElementIndex,
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

            // generate face, edge and element to element connectivity
            generate_connectivities(true,true,true);

        }
    }


    void
    create_intersect_connectivity()
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
            std::cout<<"Permutation rule not implemented";
        }


    }


    void
    reindex_template_parent_information(Mesh_Modification_Template & aMeshModTemplate)
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


    /*
     * Clear member data associated with intersection connectivity
     */
    void
    cleanup_intersect_connectivity()
    {
        mIntersectConnectivity.resize(0,0);
        mIntersectedEdges.clear();
        mIntersectedCMNodeIndex.clear();
    }

    /*
     * Add nodes created during modification inheritance information
     */
    void
    add_intersect_conn_node_inheritance()
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
