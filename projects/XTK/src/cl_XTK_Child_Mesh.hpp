/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Child_Mesh.hpp
 *
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

#include "cl_Vector.hpp"

#include "cl_XTK_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_XTK_Output_Options.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_generate_element_to_element_2D.hpp"
#include "fn_create_faces_from_element_to_node.hpp"
#include "fn_create_edges_from_element_to_node.hpp"
#include "cl_XTK_Child_Mesh_Modification_Template.hpp"

#include "cl_MTK_Cell_Info_Factory.hpp"

// MTK includes
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_verify_tet_topology.hpp"

#include "cl_XTK_Side_Cluster.hpp"

#include "assert.hpp"

using namespace moris;

namespace moris
{
    class Memory_Map;
    namespace mtk
    {
        class Vertex;
        class Mesh;
    }    // namespace mtk
}    // namespace moris

namespace moris::xtk
{

    class Child_Mesh
    {
      public:
        // ----------------------------------------------------------------------------------

        Child_Mesh();

        // ----------------------------------------------------------------------------------

        Child_Mesh( moris::uint           aSpatialDimension,
                moris::moris_index        aParentElementIndex,
                Matrix< IndexMat >       &aNodeInds,
                Matrix< IndexMat >       &aElementNodeParentInds,
                Matrix< moris::DDSTMat > &aElementNodeParentRanks,
                Matrix< IndexMat >       &aElementToNode,
                Matrix< IndexMat >       &aElementEdgeParentInds,
                Matrix< moris::DDSTMat > &aElementEdgeParentRanks,
                Matrix< IndexMat >       &aElementFaceParentInds,
                Matrix< moris::DDSTMat > &aElementFaceParentRanks,
                Matrix< moris::DDSTMat > &aElementInferfaceSides );

        Child_Mesh( moris_index              aChildMeshIndex,
                moris::uint                  aSpatialDimension,
                moris::mtk::Cell            *aParentCell,
                moris::mtk::Mesh            *aBackgroundMesh,
                Vector< moris_index > const &aGeometryIndices );

        // ----------------------------------------------------------------------------------

        Child_Mesh( Mesh_Modification_Template &aMeshModTemplate );

        // ----------------------------------------------------------------------------------

        ~Child_Mesh();

        // ----------------------------------------------------------------------------------

        // Declare iterator type
        typename std::unordered_map< moris::size_t, moris::size_t >::iterator Map_Iterator;

        // --------------------------------------------------------------
        // Functions to access connectivity
        // --------------------------------------------------------------

        /*
         * Get number of a entity of a given rank
         */
        moris::size_t
        get_num_entities( mtk::EntityRank aEntityRank ) const;

        // ----------------------------------------------------------------------------------

        mtk::EntityRank
        get_facet_rank() const;

        // ----------------------------------------------------------------------------------

        mtk::EntityRank
        get_facet_rank_internal() const;

        // ----------------------------------------------------------------------------------

        /*!
         * Returns the connectivity pointer (i.e. a TET4 connectivity)
         */
        moris::mtk::Cell_Info const *
        get_cell_info() const;

        /*!
         * Returns the connectivity pointer (i.e. a TET4 connectivity)
         */
        std::shared_ptr< moris::mtk::Cell_Info >
        get_cell_info_sp() const;

        // ----------------------------------------------------------------------------------

        /*
         * return element to node connectivity matrix (processor local indices)
         */
        Matrix< IndexMat > const &
        get_element_to_node() const;

        // ----------------------------------------------------------------------------------

        /*
         * return element to node connectivity matrix (processor local indices)
         */
        Matrix< IdMat >
        get_element_to_node_glob_ids( moris::moris_index aCMElemIndex ) const;
        /*
         * Compute and return the node to element connectivity. WARNING:
         * this computes the connectivity so use it sparingly.
         */
        Matrix< IndexMat >
        get_node_to_element_local() const;

        // ----------------------------------------------------------------------------------

        moris::moris_index
        get_element_node_ordinal( moris::moris_index aCMLocElemIndex,
                moris::moris_index                   aProcLocNodeIndex );

        // ----------------------------------------------------------------------------------

        /*
         * Converts the existing element to node connectivity (which contains proc local indices)
         * to child mesh local indices
         */
        Matrix< IndexMat >
        get_element_to_node_local() const;

        // ----------------------------------------------------------------------------------

        /*
         * Converts the existing element to node connectivity (which contains proc local indices)
         * to child mesh local indices
         */
        Matrix< IdMat >
        get_element_to_node_global() const;

        // ----------------------------------------------------------------------------------

        /*
         * Return edge to node
         */
        Matrix< IndexMat > &
        get_edge_to_node();

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_edge_to_node() const;

        // ----------------------------------------------------------------------------------

        /*
         * Return edge to node
         */
        Matrix< IndexMat >
        get_edge_to_node_local() const;

        // ----------------------------------------------------------------------------------

        /*
         * Return element to edge
         */
        Matrix< IndexMat > const &
        get_element_to_edge() const;

        // ----------------------------------------------------------------------------------

        /*
         * Return edges connected to elements
         */
        Matrix< IndexMat > const &
        get_edge_to_element() const;

        // ----------------------------------------------------------------------------------

        /*
         * return faces to node
         */
        Matrix< IndexMat > const &
        get_face_to_node() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_facet_to_node() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_face_to_element() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_facet_to_element() const;

        // ----------------------------------------------------------------------------------

        /*
         * return element to face connectivity matrix (cm local element indices)
         */
        Matrix< IndexMat >
        get_face_to_node_local() const;

        // ----------------------------------------------------------------------------------

        /*
         * return element to face connectivity matrix (cm local element indices)
         */
        Matrix< IndexMat > const &
        get_element_to_face() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_element_to_facet() const;

        // ----------------------------------------------------------------------------------

        /*
         * get the facet ordinal given a cm cell index and cm facet index
         */
        moris_index
        get_cell_facet_ordinal( moris_index aCellIndex, moris_index aFacetIndex ) const;

        // ----------------------------------------------------------------------------------

        /*
         * Return element to element connectivity (cm local element indices)
         */
        Matrix< IndexMat > const &
        get_element_to_element() const;

        // ----------------------------------------------------------------------------------

        /*
         * Return the parametric coordinates of the nodes (ordered by cm local index)
         */
        Matrix< DDRMat > const &
        get_parametric_coordinates() const;

        // ----------------------------------------------------------------------------------

        /*
         * Return the parametric coordinates of a node using a node index
         */
        Matrix< DDRMat >
        get_parametric_coordinates( moris::moris_index aNodeIndex ) const;

        // ----------------------------------------------------------------------------------

        /*!
         *
         */
        moris::mtk::Geometry_Type
        get_child_geometry_type() const;

        // ----------------------------------------------------------------------------------

        moris::mtk::Interpolation_Order
        get_child_interpolation_order() const;

        // ----------------------------------------------------------------------------------
        // Functions to access ancestry
        // ----------------------------------------------------------------------------------

        /*!
         * Get and entities parent entity rank
         */
        moris::moris_index
        get_entity_parent_entity_rank( mtk::EntityRank aEntityRank,
                moris::moris_index                     aCMLocIndex );

        // ----------------------------------------------------------------------------------

        /*!
         * Get and entities parent entity rank
         */
        moris::moris_index
        get_entity_parent_entity_proc_ind( mtk::EntityRank aEntityRank,
                moris::moris_index                         aCMLocIndex );

        // ----------------------------------------------------------------------------------

        moris::moris_index
        get_parent_element_index() const;

        // ----------------------------------------------------------------------------------

        moris::moris_index
        get_child_mesh_index() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_face_parent_inds() const;

        // ----------------------------------------------------------------------------------

        Matrix< moris::DDSTMat > const &
        get_face_parent_ranks() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_facet_parent_inds() const;

        // ----------------------------------------------------------------------------------

        Matrix< moris::DDSTMat > const &
        get_facet_parent_ranks() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_edge_parent_inds() const;

        // ----------------------------------------------------------------------------------

        Matrix< moris::DDSTMat > const &
        get_edge_parent_ranks() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_node_indices() const;

        // ----------------------------------------------------------------------------------

        Matrix< IdMat > const &
        get_node_ids() const;

        // ----------------------------------------------------------------------------------

        Vector< moris::mtk::Vertex const * > const &
        get_vertices() const;

        // ----------------------------------------------------------------------------------

        Matrix< IdMat > const &
        get_element_ids() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_element_inds() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat >
        convert_to_proc_local_elem_inds( Matrix< IndexMat > aConnOfElementCMIndices );

        // ----------------------------------------------------------------------------------

        // Function to access ordinals
        Matrix< IndexMat >
        get_edge_ordinal_from_element_and_edge_indices( moris::moris_index const &aElementIndex,
                Matrix< IndexMat > const                                         &aEdgeIndices ) const;

        // ----------------------------------------------------------------------------------

        // Function to access ordinals
        moris::moris_index
        get_edge_ordinal_from_element_and_edge_indices( moris::moris_index const &aElementIndex,
                moris::moris_index const                                         &aEdgeIndex ) const;

        // ----------------------------------------------------------------------------------

        // Function to access ordinals
        moris::moris_index
        get_face_ordinal_from_element_and_face_index( moris::moris_index const &aElementIndex,
                moris::moris_index                                              aFacetIndex ) const;

        // ----------------------------------------------------------------------------------

        /*
         * Returns the child element and face ordinal connected to a provided parent face
         */
        void
        get_child_elements_connected_to_parent_facet(
                moris::moris_index const &aParentFacetIndex,
                moris::uint              &aNumElemsOnFace,
                Matrix< IdMat >          &aChildElemsIdsOnFacet,
                Matrix< IndexMat >       &aChildElemsCMIndOnFacet,
                Matrix< IndexMat >       &aChildElemOnFacetOrdinal ) const;

        // ----------------------------------------------------------------------------------

        moris::moris_index
        get_cm_local_node_index( moris::moris_index aNodeProcIndex ) const;

        // ----------------------------------------------------------------------------------

        bool
        has_interface_along_geometry( moris_index aGeomIndex ) const;

        // ----------------------------------------------------------------------------------

        moris_index
        get_local_geom_index( moris_index aGeomIndex ) const;

        // --------------------------------------------------------------
        // Functions to modify the mesh
        // --------------------------------------------------------------

        /*
         * Modify child mesh by selecting template using Intersection connectivity
         * to determine correct template and insert the template
         */
        void
        modify_child_mesh( enum TemplateType aTemplate );

        // ----------------------------------------------------------------------------------

        void
        initialize_unzipping();

        // ----------------------------------------------------------------------------------

        void
        finalize_unzipping();

        // ----------------------------------------------------------------------------------

        /*
         * Returns the child mesh local element index not the processor local index.
         */
        void
        unzip_child_mesh_interface_get_interface_element_pairs( moris::uint aGeometryIndex,
                bool                                                       &aNoPairFoundFlag,
                Matrix< IndexMat >                                         &aInterfaceElementPairs,
                Matrix< IndexMat >                                         &aInterfacePairSideOrds ) const;

        // ----------------------------------------------------------------------------------

        void
        unzip_child_mesh_interface( moris::moris_index aGeometryIndex,
                Matrix< IndexMat > const              &aInterfaceElementPairsCMIndex,
                Matrix< IndexMat > const              &aElementUsingZippedNodes,
                Matrix< IndexMat > const              &aInterfaceNodeIndices,
                Matrix< IndexMat > const              &aUnzippedInterfaceNodeIndices,
                Matrix< IndexMat > const              &aUnzippedInterfaceNodeIds );

        // ----------------------------------------------------------------------------------

        void
        convert_tet4_to_tet10_child();

        // ----------------------------------------------------------------------------------

        /*
         * Add more nodes to the existing node indices list
         */
        void
        add_node_indices( Matrix< IndexMat > const &aNewNodeInds );

        // ----------------------------------------------------------------------------------

        /*
         * Add more nodes to the existing node indices list
         */
        void
        add_node_ids( Matrix< IdMat > const &aNewNodeIds );

        // ----------------------------------------------------------------------------------

        void
        add_vertices( Vector< moris::mtk::Vertex const * > const &aVertices );

        // ----------------------------------------------------------------------------------

        /*
         * sets the node ids. note this overwrites existing mNodeIds data
         */

        void
        set_node_ids( Matrix< IdMat > const &aNewNodeIds );

        // ----------------------------------------------------------------------------------

        /*
         * Sets the globally unique element Ids for the child mesh. This is important for mesh from data creation
         * @param[in] aElementId - First element Id (note: this is incremented as the id is assigned)
         */
        void set_child_element_ids( moris::moris_id &aElementId );

        // ----------------------------------------------------------------------------------

        /*
         * Sets the processor unique element indices for the child mesh.
         * @param[in] aElementInd - First element Ind (note: this is incremented as the ind is assigned)
         */
        void set_child_element_inds( moris::moris_index &aElementInd );

        // ----------------------------------------------------------------------------------

        /*!
         * Add node parametric coordinate for a single node
         */
        void
        add_node_parametric_coordinate( moris::size_t aNodeIndex,
                Matrix< DDRMat > const               &aParamCoord );

        // ----------------------------------------------------------------------------------

        /*!
         * Add node parametric coordinate for multiple nodes
         * @param[in] aNodeIndices - List of node indices with parametric coordinates
         * @param[in] aParamCoord  - Parametric coordinates for nodes (note: row 1 of aParamCoord is the Parametric coordinate of aNodeIndices(1))
         */
        void
        add_node_parametric_coordinate( Matrix< IndexMat > const &aNodeIndices,
                Matrix< DDRMat > const                           &aParamCoord );

        // ----------------------------------------------------------------------------------

        void
        allocate_parametric_coordinates( moris::size_t aNumNewNodes,
                moris::size_t                          aDimOfParmCoord );

        // ----------------------------------------------------------------------------------

        /*
         * Resizes element related matrices
         */
        void
        allocate_more_elements( moris::size_t const &aNumMoreElements );

        // ----------------------------------------------------------------------------------

        void
        insert_child_mesh_template( Mesh_Modification_Template &aMeshModTemplate );

        // ----------------------------------------------------------------------------------

        void
        setup_template_interface_facets( Mesh_Modification_Template &aMeshModTemplate,
                Matrix< DDSTMat >                                   &aInterfaceSideOrdinals );

        // ----------------------------------------------------------------------------------

        /*
         * Convert the element to node, parent information to  the local indexing scheme rather than ordinal
         */
        void reindex_template( Mesh_Modification_Template &aMeshModTemplate );

        // ----------------------------------------------------------------------------------

        /*
         * Replace the information associated with an element
         */
        void
        replace_element( moris::size_t const   &aElementIndexToReplace,
                moris::size_t const            &aRowIndex,
                moris::size_t const            &aGeomIndex,
                Matrix< IndexMat > const       &aElementToNode,
                Matrix< IndexMat > const       &aElementEdgeParentInds,
                Matrix< moris::DDSTMat > const &aElementEdgeParentRanks,
                Matrix< IndexMat > const       &aElementFaceParentInds,
                Matrix< moris::DDSTMat > const &aElementFaceParentRanks,
                Matrix< moris::DDSTMat > const &aElementInterfaceFaces );

        // ----------------------------------------------------------------------------------

        /*
         * Replace the information associated with an element
         */
        void
        add_element( moris::size_t const       &aRowIndex,
                moris::size_t const            &aGeomIndex,
                Matrix< IndexMat > const       &aElementToNode,
                Matrix< IndexMat > const       &aElementEdgeParentInds,
                Matrix< moris::DDSTMat > const &aElementEdgeParentRanks,
                Matrix< IndexMat > const       &aElementFaceParentInds,
                Matrix< moris::DDSTMat > const &aElementFaceParentRanks,
                Matrix< moris::DDSTMat > const &aElementInterfaceFaces );

        // ----------------------------------------------------------------------------------

        /*!
         * Add node inheritance
         * @param[in] aCMLocNodeIndex - Child mesh local node index
         * @param[in] aProcLocParentEntityInd - Processor local parent entity index
         * @param[in] aParentEntityRank - Parent entity rank
         */
        void
        set_node_inheritance( moris::moris_index aCMLocNodeIndex,
                moris::moris_index               aProcLocParentEntityInd,
                moris::moris_index               aParentEntityRank );

        // ----------------------------------------------------------------------------------

        /*
         * Generates face connectivities, edge connectivities, element to element graph
         */

        void
        generate_connectivities( bool aGenerateFaceConn,
                bool                  aGenerateEdgeConn,
                bool                  aGenerateElemToElem );

        // ----------------------------------------------------------------------------------

        void
        add_new_geometry_interface( moris_index aChildIndex );

        // ----------------------------------------------------------------------------------

        /**
         * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
         *       - 1 means the provided aDPrime1Ind is an XTK index
         *
         * aDPrime2Ind must be XTK local index
         */
        void
        add_entity_to_intersect_connectivity( moris::moris_index aCMNodeInd,
                moris::moris_index                               aCMEdgeInd,
                moris::size_t                                    aFlag );

        // ----------------------------------------------------------------------------------

        void
        mark_edge_as_on_interface( moris::moris_index aEdgeIndex );

        void
        mark_facet_as_on_interface( moris_index aFacetCMIndex,
                moris_index                     aGeometryIndex );

        // ----------------------------------------------------------------------------------

        /*
         * Take the information of edges on the interface and
         * figure out which faces are on the interface, then
         * mark element edges as on interface.
         */
        void
        mark_interface_faces_from_interface_coincident_faces();

        // ----------------------------------------------------------------------------------

        // --------------------------------------------------------------
        // Functions for elemental phase information and subphase bins
        // --------------------------------------------------------------

        void
        initialize_element_phase_mat();

        // ----------------------------------------------------------------------------------

        void
        set_element_phase_index( moris::size_t aEntityIndex,
                moris::size_t                  aEntityPhaseIndex );

        // ----------------------------------------------------------------------------------

        moris::size_t
        get_element_phase_index( moris::size_t const &aEntityIndex ) const;

        // ----------------------------------------------------------------------------------

        moris_index
        get_element_subphase_index( moris::size_t const &aEntityIndex ) const;

        // ----------------------------------------------------------------------------------

        /*!
         * Returns the active bulk phases in this child mesh
         */
        Vector< moris_index >
        get_active_bulk_phases( Vector< moris_index > &aBulkPhasesIndex ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_element_phase_indices() const;

        // ----------------------------------------------------------------------------------

        moris_index
        get_element_subphase_id( moris::size_t const &aEntityIndex ) const;

        // ----------------------------------------------------------------------------------

        moris::size_t
        get_num_subphase_bins() const;

        // ----------------------------------------------------------------------------------

        Vector< moris::moris_index > const &
        get_subphase_bin_bulk_phase() const;

        // ----------------------------------------------------------------------------------

        Vector< Matrix< IndexMat > > const &
        get_subphase_groups() const;

        // ----------------------------------------------------------------------------------

        Vector< moris_index > const &
        get_subphase_indices() const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > const &
        get_subphase_ids() const;

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_loc_index( moris_index aSubPhaseIndex ) const;

        // ----------------------------------------------------------------------------------

        void
        get_subphases_attached_to_facet( moris_index aFacetIndex,
                Vector< moris_index >               &aSubPhaseCMIndex,
                Vector< moris_index >               &aRepresentativeChildCellInd,
                Vector< moris_index >               &aRepresentativeChildCellSideOrdinal ) const;

        // ----------------------------------------------------------------------------------

        /*
         * Constructs the internal double side sets for this child mesh
         * @return Whether there are interfaces between another child mesh
         */
        bool
        construct_internal_double_sides_between_subphases();

        // ----------------------------------------------------------------------------------

        bool
        has_inter_child_mesh_interfaces();

        // ----------------------------------------------------------------------------------

        uint
        get_num_double_side_interfaces() const;

        // ----------------------------------------------------------------------------------

        Vector< moris_index > const &
        get_double_side_interface_subphase_indices( moris_index aDblSideCMIndex ) const;

        // ----------------------------------------------------------------------------------

        Vector<Vector< moris_index > > const &
        get_double_side_interface_cell_pairs( moris_index aDblSideCMIndex ) const;

        // ----------------------------------------------------------------------------------

        Vector<Vector< moris_index > > const &
        get_double_side_interface_cell_pairs_facet_ords( moris_index aDblSideCMIndex ) const;

        // ----------------------------------------------------------------------------------

        void
        delete_double_sides_interface_sets();

        // ----------------------------------------------------------------------------------

        void
        print_double_sides_between_subphases( moris_index aVerboseLevel = 0 );

        // ----------------------------------------------------------------------------------

        /*
         * Sets the elemental subphase value in this child mesh
         */
        void
        set_elemental_subphase( moris::moris_index &aSubPhaseIndex,
                Matrix< IndexMat > const           &aElementSubPhase );

        // ----------------------------------------------------------------------------------

        void
        set_subphase_id( moris_id const &aSubphaseIndex,
                moris_id                &aSubphaseId );

        // ----------------------------------------------------------------------------------

        /*
         * Returns the elemental subphase values
         */
        Matrix< IndexMat > const &
        get_elemental_subphase_bin_membership() const;

        // ----------------------------------------------------------------------------------

        /*
         * Add a basis to subphase bin relationship
         */
        void
        add_basis_and_enrichment_to_subphase_group( moris_index aSubphaseBinIndex,
                moris_index                                     aBasisIndex,
                moris_index                                     aBasisEnrLev );

        // ----------------------------------------------------------------------------------

        Vector< moris_index > const &
        get_subphase_basis_indices( moris_index aSubphaseBin ) const;

        // ----------------------------------------------------------------------------------

        Vector< moris_index > const &
        get_subphase_basis_enrichment_levels( moris_index aSubphaseBin ) const;

        // --------------------------------------------------------------

        void
        reindex_cells( Vector< moris_index > &aOldIndexToNewCellIndex );

        // --------------------------------------------------------------
        // Functions IO
        // --------------------------------------------------------------

        void pack_child_mesh_by_phase( moris::size_t const &aNumPhases,
                Vector< Matrix< IdMat > >                  &aElementIds,
                Vector< Matrix< IdMat > >                  &aElementCMInds ) const;

        // ----------------------------------------------------------------------------------

        Matrix< moris::DDSTMat > const &
        get_cell_interface_side_ords();
        /*
         * Fora a given interface on geometry, pack the sides
         */
        Matrix< IdMat >
        pack_interface_sides( moris_index aGeometryIndex,
                moris_index               aPhaseIndex0,
                moris_index               aPhaseIndex1,
                moris_index               aIndexFlag = 0 ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IdMat >
        pack_interface_sides_loc_inds() const;

        // ----------------------------------------------------------------------------------

        mtk::CellTopology
        template_to_cell_topology( enum TemplateType aTemplateType );

        // ----------------------------------------------------------------------------------

        void
        mark_as_hmr_child_mesh() { mHMR = true; };

        // ----------------------------------------------------------------------------------

        void identify_hanging_nodes( const Vector< moris_index > &aTransitionFacetIndices );

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > get_hanging_nodes() const { return mHangingNodes; };

        // ----------------------------------------------------------------------------------
        // Memory Printing / Computing
        // ----------------------------------------------------------------------------------
        moris::Memory_Map
        get_memory_usage();

      private:
        // child mesh index
        moris::moris_index mChildMeshIndex;

        // Parent element index
        moris::moris_index mParentElementIndex;

        // Element To Node and Ancestry Information (This is the only data that is set with templates.
        // all other is generated with an algorithm)
        // All node connectivity is indexed by proc local indexs
        mtk::CellTopology                        mElementTopology;
        std::shared_ptr< moris::mtk::Cell_Info > mConnectivity;

        moris::size_t            mNumElem;
        Matrix< IndexMat >       mElementToNode; /* proc inds*/
        Matrix< IndexMat >       mElementEdgeParentInds;
        Matrix< moris::DDSTMat > mElementEdgeParentRanks;
        Matrix< IndexMat >       mElementFaceParentInds;
        Matrix< moris::DDSTMat > mElementFaceParentRanks;
        Matrix< moris::DDSTMat > mElementInterfaceSides;

        // spatial dimension
        moris::uint mSpatialDimension;

        // Geometries which intersect this child mesh
        Vector< moris_index > mGeometryIndex;

        // Child element information ---------------------------
        Matrix< IdMat >    mChildElementIds;
        Matrix< IndexMat > mChildElementInds;
        // No child element parent information is stored because
        // it is assumed that all are the parent element of this
        // child mesh.

        // Node information ------------------------------------
        Vector< moris::mtk::Vertex const * > mVertices;
        Matrix< IdMat >                      mNodeIds;
        Matrix< IndexMat >                   mNodeInds;
        Matrix< moris::DDSTMat >             mNodeParentRank;
        Matrix< IndexMat >                   mNodeParentInds;

        // Map where  key - proc local ind, val - local child mesh index
        std::unordered_map< moris::size_t, moris::size_t > mNodeIndsToCMInd;

        // Parametric coordinate relative to parent element
        Matrix< DDRMat > mNodeParametricCoord;

        // Face Connectivity -----------------------------------
        bool                     mHasFaceConn;
        Matrix< IndexMat >       mFaceToNode; /* proc inds*/
        Matrix< IndexMat >       mNodeToFace;
        Matrix< IndexMat >       mFaceToElement;
        Matrix< IndexMat >       mElementToFace;
        Matrix< IndexMat >       mFaceParentInds;
        Matrix< moris::DDSTMat > mFaceParentRanks;

        // Edge connectivity -----------------------------------
        bool                     mHasEdgeConn;
        Matrix< IndexMat >       mEdgeToNode; /* proc inds*/
        Matrix< IndexMat >       mNodeToEdge;
        Matrix< IndexMat >       mEdgeToElement;
        Matrix< IndexMat >       mElementToEdge;
        Matrix< IndexMat >       mEdgeParentInds;
        Matrix< moris::DDSTMat > mEdgeParentRanks;

        // Element to Element graph ----------------------------
        bool               mHasElemToElem;
        Matrix< IndexMat > mElementToElement;

        // Auxiliary connectivity data and pending nodes (mesh modification data)
        Matrix< IndexMat >    mIntersectConnectivity;
        Vector< moris_index > mIntersectedCMNodeIndex;
        Vector< moris_index > mIntersectedEdges;

        bool               mHasCoincidentEdges;
        Matrix< IndexMat > mEdgeOnInterface;

        bool mHasInterChildMeshInterface;

        // Phase member variables -----------------------------
        bool                         mHasPhaseInfo;
        Matrix< IndexMat >           mElementPhaseIndices;
        Matrix< IndexMat >           mElementBinIndex;
        Vector< moris::moris_index > mBinBulkPhase;

        Matrix< IndexMat >              mSubPhaseBinId;      /*glob id of subphase bin*/
        Vector< moris_index >           mSubPhaseBinIndices; /*proc index of subphase bin*/
        Vector< Matrix< IndexMat > >    mSubPhaseBins;
        Vector< Vector< moris_index > > mSubphaseBasisIndices;
        Vector< Vector< moris_index > > mSubphaseBasisEnrichmentLevel;

        // Double side set between subphases
        Vector<Vector< moris_index > >           mDoubleSideSetSubphaseInds;
        Vector<Vector<Vector< moris_index > > > mDoubleSideSetCellPairs;
        Vector<Vector<Vector< moris_index > > > mDoubleSideSetFacetPairs;

        Matrix< IndexMat > mHangingNodes;

        // Unzipping information
        bool mUnzippingFlag = false;

        // HMR Flag
        bool mHMR = false;

      private:
        // ----------------------------------------------------------------------------------
        void
        generate_face_connectivity_and_ancestry( Matrix< IndexMat > const &aElementToNodeLocal );

        // ----------------------------------------------------------------------------------

        void
        generate_edge_connectivity_and_ancestry( Matrix< IndexMat > const &aElementToNodeLocal );

        // ----------------------------------------------------------------------------------

        void
        generate_element_to_element_connectivity();

        // ----------------------------------------------------------------------------------

        void
        set_up_proc_local_to_cm_local_node_map();

        // ----------------------------------------------------------------------------------

        /*
         * Note: this function only adds nodes which do not already exist to the map
         * and ignores nodes which are already in the map.
         */
        void
        add_nodes_to_map( Matrix< IndexMat > const &aNodesToAdd );

        // ----------------------------------------------------------------------------------

        /*
         * Takes the element face ordinal ancestry and creates a face ancestry
         */
        void
        setup_face_ancestry();

        // ----------------------------------------------------------------------------------

        void
        setup_edge_ancestry();

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat >
        convert_to_cm_local_indices( Matrix< IndexMat > const &aEntityRankToNode ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat >
        convert_to_proc_indices( Matrix< IndexMat > const &aEntityRankToNodeLocal ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IdMat >
        convert_cm_loc_to_glob_ids( Matrix< IndexMat > const &aEntityRankToNodeLocal ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IdMat >
        convert_proc_to_glob_ids( Matrix< IndexMat > const &aEntityRankToNodeLocal ) const;

        // ----------------------------------------------------------------------------------

        void
        modify_child_mesh_internal( enum TemplateType aTemplate );

        // ----------------------------------------------------------------------------------

        void
        create_intersect_connectivity();

        // ----------------------------------------------------------------------------------

        /*
         * Edges needed to get permutation
         */
        Matrix< IndexMat >
        sort_nodes( enum TemplateType     aTemplate,
                Matrix< IndexMat > const &aIntConnectivity,
                Matrix< IndexMat > const &aEdgeToNodeCMLoc,
                moris::size_t const      &aElementIndex,
                moris::size_t            &aPermutation );

        // ----------------------------------------------------------------------------------

        /*
         * This assumes the edge ordinals provided are ordered in the same order as is found in the auxiliary connectivity
         */
        void
        get_intersection_permutation( Matrix< IndexMat > const &aOrderedEdgeOrdinals,
                moris::size_t                                  &aPermutation );

        // ----------------------------------------------------------------------------------

        void
        reindex_template_parent_information( Mesh_Modification_Template &aMeshModTemplate );

        // ----------------------------------------------------------------------------------

        /*
         * Clear member data associated with intersection connectivity
         */
        void
        cleanup_intersect_connectivity();

        // ----------------------------------------------------------------------------------

        /*
         * Add nodes created during modification inheritance information
         */
        void
        add_intersect_conn_node_inheritance();

        // ----------------------------------------------------------------------------------

        void
        construct_subphase_bins( moris::moris_index &aSubPhaseIndex,
                Matrix< IndexMat > const            &aElementSubPhase );

        // ----------------------------------------------------------------------------------

        template< typename Mat_Type >
        void
        row_vector_connectivity_check( Matrix< Mat_Type > &aConnectivity )
        {
            MORIS_ASSERT( moris::isvector( aConnectivity ), "Provided connectivity is not a vector" );

            if ( moris::iscol( aConnectivity ) )
            {
                aConnectivity = moris::trans( aConnectivity );
            }
        }
    };

}    // namespace moris::xtk
#endif /* SRC_XTK_CL_XTK_CHILD_MESH_HPP_ */
